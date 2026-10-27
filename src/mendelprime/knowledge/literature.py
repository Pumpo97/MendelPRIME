"""Bounded PubMed and OMIM retrieval as plain state and document records."""

from __future__ import annotations

import hashlib
import json
import math
import re
import ssl
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path

from mendelprime.domain import MendelPrimeError, Source, case_tuple, validate_gene, variant_tuple
from mendelprime.runtime.http import urlopen
from mendelprime.runtime.storage import (
    atomic_replace_bytes,
    open_regular,
    private_directory,
    read_regular_bounded,
)

SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
FETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNAVAILABLE_TTL_SECONDS = 24 * 60 * 60
_MAX_CACHE_BYTES = 8 * 1024 * 1024
_MAX_SEARCH_BYTES = 2 * 1024 * 1024
_MAX_FETCH_BYTES = 64 * 1024 * 1024
_NCBI_INTERVAL_SECONDS = 0.34
_MAX_RETRY_DELAY_SECONDS = 30
_RETRYABLE_HTTP = {408, 429, 500, 502, 503, 504}
_SAFE_FETCH_ERRORS = {
    "pubmed_record_missing_or_unsupported",
    "article_or_abstract_missing",
    "pubmed_cache_write_failed",
    "pubmed_fetch_failed",
    "pubmed_http_error",
    "pubmed_invalid_response",
    "pubmed_response_too_large",
    "pubmed_tls_error",
    "pubmed_transport_error",
    "source_citation_identifier_unresolvable",
}


def new_pubmed(cache_root, timeout_seconds=45):
    """Create provider state; tests may replace network, clock, and sleep callables."""
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    root = Path(cache_root)
    private_directory(root)
    return {
        "cache_root": root,
        "timeout_seconds": float(timeout_seconds),
        "request_lock": threading.Lock(),
        "next_request_at": 0.0,
        "monotonic": time.monotonic,
        "rate_sleep": time.sleep,
    }


def search(state, query, limit):
    if not isinstance(query, str) or not query.strip():
        raise ValueError("PubMed query must be nonempty")
    if type(limit) is not int or not 1 <= limit <= 10_000:
        raise ValueError("PubMed search limit must be between 1 and 10000")
    relative = Path("search") / f"{_hash(query)}_{limit}.json"
    cached = _cached_pmids(_read_cache(state, relative), query)
    if cached is not None:
        return cached
    params = urllib.parse.urlencode(
        {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": limit,
            "sort": "relevance",
        }
    )
    def decode(content):
        payload = json.loads(content.decode("utf-8"))
        if not isinstance(payload, dict) or not isinstance(payload.get("esearchresult"), dict):
            raise ValueError("Expected a PubMed search response")
        if any(name in block for block in (payload, payload["esearchresult"]) for name in ("error", "ERROR")):
            raise ValueError("PubMed returned an error document")
        raw = payload["esearchresult"].get("idlist")
        if not isinstance(raw, list):
            raise ValueError("Expected PubMed identifiers")
        return tuple(str(value) for value in raw if str(value).isdigit())

    pmids = _network_bytes(state, f"{SEARCH_URL}?{params}", _MAX_SEARCH_BYTES, decode=decode)
    _write_cache(state, relative, {"query": query, "pmids": pmids})
    return pmids


def fetch_many(state, pmids):
    unique = tuple(dict.fromkeys(_require_pmid(pmid) for pmid in pmids))
    publications, errors, missing = {}, {}, []
    now = state.get("clock", time.time)()
    for pmid in unique:
        cached = _read_cache(state, Path("articles") / f"PMID_{pmid}.json")
        publication = _publication_from_cache(cached)
        if publication is not None and publication["pmid"] == pmid:
            publications[pmid] = publication
        elif _cached_unavailable(cached, pmid, now):
            errors[pmid] = "article_or_abstract_missing"
        else:
            missing.append(pmid)

    for start in range(0, len(missing), 100):
        batch = missing[start : start + 100]
        try:
            fetched = _fetch(state, batch)
        except Exception as error:
            resource = isinstance(error, MendelPrimeError) and error.category == "resource"
            errors.update({pmid: error.code if resource else "pubmed_fetch_failed" for pmid in batch})
            continue
        for pmid in batch:
            publication = fetched.get(pmid)
            if publication is None:
                errors[pmid] = "pubmed_record_missing_or_unsupported"
                continue
            relative = Path("articles") / f"PMID_{pmid}.json"
            if not publication["abstract"].strip():
                errors[pmid] = "article_or_abstract_missing"
                try:
                    _write_cache(state, relative, _unavailable(pmid, state.get("clock", time.time)()))
                except MendelPrimeError as error:
                    errors[pmid] = error.code
                continue
            publications[pmid] = publication
            try:
                _write_cache(state, relative, publication)
            except MendelPrimeError as error:
                errors[pmid] = error.code
    return publications, errors


def parse_pubmed_xml(text):
    publications = {}
    root = ET.fromstring(text)
    if root.tag != "PubmedArticleSet" or root.find(".//ERROR") is not None:
        raise ET.ParseError("Expected PubMed records, not a provider error document")
    for record in root:
        book = record.tag == "PubmedBookArticle"
        if record.tag not in {"PubmedArticle", "PubmedBookArticle"}:
            continue
        container = record.find("BookDocument" if book else "MedlineCitation")
        if container is None:
            continue
        article = container if book else container.find("Article")
        if article is None:
            continue
        pmid = _xml_text(container.find("PMID"))
        if not pmid.isdigit():
            continue
        abstract = (_xml_text(node) for node in article.findall("Abstract/AbstractText"))
        kinds = (_xml_text(node) for node in record.findall(".//PublicationTypeList/PublicationType"))
        kinds = tuple(dict.fromkeys(filter(None, kinds)))
        if book and "Book Chapter" not in kinds:
            kinds += ("Book Chapter",)
        publications[pmid] = {
            "pmid": pmid,
            "title": _xml_text(article.find("ArticleTitle")),
            "abstract": " ".join(filter(None, abstract)),
            "journal": _xml_text(article.find("Book/BookTitle" if book else "Journal/Title")),
            "year": _xml_text(
                article.find("Book/PubDate/Year" if book else "Journal/JournalIssue/PubDate/Year")
            ),
            "publication_types": kinds,
        }
    return publications


def load_omim(records_path, mim2gene_path, text_limit=6000):
    records = {}
    with open_regular(Path(records_path)) as (stream, _status):
        for raw in stream:
            try:
                value = json.loads(raw.decode("utf-8", errors="ignore"))
            except (json.JSONDecodeError, RecursionError):
                continue
            if isinstance(value, dict) and (omim_id := str(value.get("omim_id") or "")):
                records[omim_id] = value
    entries = defaultdict(list)
    with open_regular(Path(mim2gene_path)) as (stream, _status):
        for raw in stream:
            line = raw.decode("utf-8", errors="ignore")
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(parts) < 4 or parts[1] != "gene":
                continue
            for symbol in re.split(r",\s*", parts[3]):
                gene = symbol.strip().upper()
                if gene and parts[0] not in entries[gene]:
                    entries[gene].append(parts[0])
    return {
        "records": records,
        "entries": {key: tuple(value) for key, value in entries.items()},
        "text_limit": max(300, text_limit),
    }


def molecular_genetics(omim, gene, limit=1):
    documents = []
    for omim_id in omim["entries"].get(gene, ()):
        record = omim["records"].get(omim_id, {})
        serialized = json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        version = "sha256:" + _hash(serialized)
        for index, paragraph in enumerate((record.get("sections") or {}).get("molecular_genetics", ()), 1):
            if not isinstance(paragraph, dict):
                continue
            raw_text = str(paragraph.get("text") or "")
            normalized = " ".join(raw_text.split())
            if len(normalized) < 40:
                continue
            prompt_text = normalized[: omim["text_limit"]]
            pmids = tuple(str(value) for value in paragraph.get("pmids", ()) if str(value).isdigit())
            source = _source(
                "omim",
                "gene_disease",
                f"https://omim.org/entry/{omim_id}#molecular-genetics-{index}",
                version,
                raw_text,
                gene,
            )
            documents.append(
                evidence_document(
                    source,
                    f"OMIM:{omim_id}:molecular_genetics:{_hash(prompt_text)[:12]}",
                    _omim_title(record, omim_id),
                    prompt_text,
                    provenance=("omim_gene_entry",),
                    cited_pmids=pmids,
                    section="molecular_genetics",
                )
            )
            if len(documents) == limit:
                return tuple(documents)
    return tuple(documents)


def cited_pmids(omim, gene, limit=10):
    pmids = []
    for document in molecular_genetics(omim, gene):
        _append_pmids(pmids, document["cited_pmids"], limit)
    for omim_id in omim["entries"].get(gene, ()):
        for variant in omim["records"].get(omim_id, {}).get("allelic_variants", ()):
            if not isinstance(variant, dict):
                continue
            for part in variant.get("sub_parts", ()):
                if isinstance(part, dict):
                    _append_pmids(pmids, part.get("pmids", ()), limit)
        if len(pmids) >= limit:
            break
    return tuple(pmids)


def select_articles(pubmed, genes, paper_scores, excluded_pmids=()):
    excluded = {_pmid(value) for value in excluded_pmids} - {""}
    selected = {}
    for gene in genes:
        scored = _scored_papers(paper_scores.get(gene, ()), excluded)
        articles, _errors = fetch_many(pubmed, tuple(pmid for pmid, _score in scored))
        article = next(
            (
                articles[pmid]
                for pmid, _score in scored
                if pmid in articles and articles[pmid]["title"] and articles[pmid]["abstract"]
            ),
            None,
        )
        route = "amelie"
        if article is None:
            query = (
                f'("{gene}"[Title/Abstract]) AND '
                "(clinical OR patient OR patients OR syndrome OR phenotype OR case)"
            )
            hits = search(pubmed, query, 5 + len(excluded))
            pmids = tuple(pmid for pmid in hits if _pmid(pmid) not in excluded)[:5]
            articles, _errors = fetch_many(pubmed, pmids)
            article = next(
                (
                    articles[pmid]
                    for pmid in pmids
                    if pmid in articles and articles[pmid]["title"] and articles[pmid]["abstract"]
                ),
                None,
            )
            route = "pubmed_fallback"
        if article is None:
            selected[gene] = {"gene": gene, "pmid": "", "title": "", "source": None, "route": "unavailable"}
        else:
            selected[gene] = {
                "gene": gene,
                "pmid": article["pmid"],
                "title": article["title"],
                "source": _publication_source(article, gene),
                "route": route,
            }
    return selected


def retrieve_evidence(
    pubmed,
    omim,
    profile,
    genes,
    seed_pmids,
    must_review,
    excluded_pmids=(),
    citation_seeds=(),
    *,
    variant_aware=False,
    variant_search_terms_by_gene=None,
):
    if any(
        seed.get("gene") not in genes
        or seed.get("route")
        not in {
            "clinvar_citation",
            "clingen_citation",
        }
        for seed in citation_seeds
    ):
        raise ValueError("variant citation seeds must belong to admitted genes and known routes")
    excluded = {_pmid(value) for value in excluded_pmids} - {""}
    prepared, all_pmids = {}, []
    labels = tuple(term["label"] for term in profile["terms"] if term["present"])
    for gene in genes:
        validate_gene(gene)
        provenance = defaultdict(set)
        for value in must_review.get(gene, ()):
            if pmid := _pmid(value):
                provenance[pmid].add("must_review")
        variant_pmids = set()
        for seed in citation_seeds:
            if seed["gene"] == gene and (pmid := _pmid(seed.get("pmid"))):
                if variant_aware and pmid not in variant_pmids and len(variant_pmids) >= 5:
                    continue
                variant_pmids.add(pmid)
                provenance[pmid].add(seed["route"])
                if variant_aware:
                    provenance[pmid].update(("variant", "clinvar_variant_exact"))
        papers = sorted(
            enumerate(seed_pmids.get(gene, ())),
            key=lambda item: (-_paper(item[1])[1], item[0], _paper(item[1])[0]),
        )
        for _index, paper in papers[:5]:
            if pmid := _pmid(_paper(paper)[0]):
                provenance[pmid].add("amelie")
        omim_context = molecular_genetics(omim, gene)
        if variant_aware:
            for document in omim_context:
                for pmid in document.get("cited_pmids", ()):
                    if pmid := _pmid(pmid):
                        provenance[pmid].add("omim_molecular_genetics")
        for pmid in cited_pmids(omim, gene, 10):
            provenance[pmid].add("omim")
        errors, searched = [], []
        queries = evidence_queries(gene, labels)
        terms = ()
        if variant_aware:
            terms = tuple(dict.fromkeys((variant_search_terms_by_gene or {}).get(gene, ())))[:4]
            variant_queries = tuple(
                f'("{gene}"[Title/Abstract]) AND ("{term.replace(chr(34), "")}"[Title/Abstract])'
                for term in terms
                if isinstance(term, str) and term.strip()
            )
            queries = (*variant_queries, *queries)
        for query in queries:
            try:
                hits = search(pubmed, query, 30)
            except Exception as error:
                errors.append(_safe_error(error, "pubmed_search_failed"))
                continue
            for value in hits:
                if (pmid := _pmid(value)) and pmid not in searched:
                    searched.append(pmid)
        for pmid in searched[:20]:
            provenance[pmid].add("pubmed_search")
        for pmid in excluded:
            provenance.pop(pmid, None)
        prepared[gene] = provenance, omim_context, errors, terms
        all_pmids.extend(provenance)

    publications, fetch_errors = fetch_many(pubmed, tuple(all_pmids))
    output = {}
    for gene, (provenance, omim_context, errors, terms) in prepared.items():
        documents = [
            _article_document(gene, publications[pmid], tuple(sorted(routes)), 3500)
            for pmid, routes in provenance.items()
            if pmid in publications and publications[pmid]["abstract"].strip()
        ]
        fetch_issues = (
            f"PMID:{pmid}:{_safe_fetch_error(error)}"
            for pmid, error in fetch_errors.items()
            if pmid in provenance
        )
        output[gene] = {
            "gene": gene,
            "documents": _select_variant_documents(documents, 10)
            if variant_aware
            else _select_documents(documents, 10),
            "omim_context": omim_context,
            "retrieval_errors": tuple((*errors, *fetch_issues)),
        }
        if variant_aware:
            output[gene]["variant_expansion"] = {
                "retrieval_strategy": "september_vcf_v1",
                "variant_search_terms": terms,
                "coverage": (),
                "issues": (),
            }
            if profile.get("case_key") is not None:
                output[gene]["variant_expansion"]["case_key"] = case_tuple(profile["case_key"])
    return output


def _select_variant_documents(documents, limit):
    """Select VCF documents by route priority and per-route limits."""
    selected, seen = [], set()
    for route, maximum in (
        ("clinvar_variant_exact", 3),
        ("omim_matched_allelic_variant", 2),
        ("omim_molecular_genetics", 2),
        ("verification_search", 5),
        ("must_review", 1),
        ("variant", 2),
        ("amelie", 3),
        ("omim", 3),
        ("pubmed_search", None),
        (None, None),
    ):
        added = 0
        for document in documents:
            if len(selected) >= limit or maximum is not None and added >= maximum:
                break
            if (route is None or route in document["provenance"]) and document["prompt_id"] not in seen:
                selected.append(document)
                seen.add(document["prompt_id"])
                added += 1
    priority = (
        "clinvar_variant_exact",
        "omim_matched_allelic_variant",
        "omim_molecular_genetics",
        "amelie",
        "pubmed_search",
    )
    return tuple(
        sorted(
            selected,
            key=lambda document: next(
                (index for index, route in enumerate(priority) if route in document["provenance"]),
                len(priority),
            ),
        )
    )


def _variant_search_terms(target, *, report=False):
    primary = target["primary"]
    terms = tuple(
        dict.fromkeys(
            value.split(":")[-1]
            for row in (primary, *target.get("aliases", ()))
            for name in ("hgvs_c", "hgvs_p")
            if isinstance((value := row.get(name)), str) and value.strip()
        )
    )
    if not report:
        return terms
    # Parentheses and the p. prefix are HGVS syntax, not how most papers spell a variant.
    protein_codes = dict(
        zip(
            "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter".split(),
            "ARNDCQEGHILKMFPSTWYV*",
            strict=True,
        )
    )
    normalized = [term.removeprefix("p.").strip("()") if term.startswith("p.") else term for term in terms]
    normalized.extend(
        re.sub("|".join(protein_codes), lambda m: protein_codes[m[0]], term)
        for term in tuple(normalized)
        if not term.startswith("c.")
    )
    return tuple(dict.fromkeys((*terms, *normalized)))


def _variant_article_key(document):
    parent = document.get("metadata", {}).get("parent_pmid", "")
    if parent.startswith("PMID:"):
        return parent
    source = document["source"]
    if source.provider == "pmc" and source.kind in {"derived", "annotation"}:
        try:
            lineage = json.loads(source.content)
        except (ValueError, RecursionError):
            lineage = {}
        if isinstance(lineage, dict) and lineage.get("record_type") == "pmc_article_lineage_v1":
            return lineage.get("pmid", "")
    identifier = document["prompt_id"]
    if identifier.startswith("PMID:"):
        return identifier.split("#", 1)[0]
    locator = document["source"].locator
    match = re.search(r"(?:pubmed\.ncbi\.nlm\.nih\.gov/|PMID:)(\d+)", locator)
    return f"PMID:{match[1]}" if match else ""


def _pmc_resource(pubmed, pmid):
    relative = Path("ranking_pmc") / f"PMID_{pmid}.json"
    cached = _read_cache(pubmed, relative)
    if isinstance(cached, dict) and cached.get("pmid") == pmid:
        if cached.get("status") == "not_available":
            return cached
        if cached.get("status") == "found" and isinstance(cached.get("xml"), str):
            try:
                _check_pmc_identity(cached["xml"], pmid, cached.get("pmcid", ""))
            except (MendelPrimeError, ET.ParseError):
                pass
            else:
                return cached
    params = urllib.parse.urlencode({"ids": pmid, "format": "json"})
    raw = _network_bytes(
        pubmed, "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?" + params, _MAX_SEARCH_BYTES
    )
    try:
        payload = json.loads(raw)
        records = payload["records"]
        pmcid = next((row.get("pmcid", "") for row in records if str(row.get("pmid")) == pmid), "")
    except (KeyError, ValueError, TypeError, AttributeError) as error:
        raise _resource_error("pmc_invalid_response") from error
    if not pmcid:
        result = {"pmid": pmid, "status": "not_available"}
    else:
        if not re.fullmatch(r"PMC\d+", pmcid):
            raise _resource_error("pmc_invalid_identity")
        query = urllib.parse.urlencode({"db": "pmc", "id": pmcid, "retmode": "xml"})
        xml = _network_bytes(pubmed, FETCH_URL + "?" + query, _MAX_CACHE_BYTES).decode("utf-8")
        _check_pmc_identity(xml, pmid, pmcid)
        result = {"pmid": pmid, "pmcid": pmcid, "status": "found", "xml": xml}
    _write_cache(pubmed, relative, result)
    return result


def _main_pmc_article(xml):
    root = ET.fromstring(xml)
    articles = (
        (root,)
        if root.tag.rsplit("}", 1)[-1] == "article"
        else tuple(node for node in root if node.tag.rsplit("}", 1)[-1] == "article")
    )
    if len(articles) != 1:
        raise _resource_error("pmc_article_identity_mismatch")
    return articles[0]


def _check_pmc_identity(xml, pmid, pmcid):
    if not isinstance(pmcid, str) or not re.fullmatch(r"PMC\d+", pmcid):
        raise _resource_error("pmc_invalid_identity")
    article = _main_pmc_article(xml)
    front = next((node for node in article if node.tag.rsplit("}", 1)[-1] == "front"), ())
    metadata = next((node for node in front if node.tag.rsplit("}", 1)[-1] == "article-meta"), ())
    identifiers = tuple(node for node in metadata if node.tag.rsplit("}", 1)[-1] == "article-id")
    actual = {
        "".join(node.itertext()).strip()
        for node in identifiers
        if node.get("pub-id-type") in {"pmc", "pmcid"}
    }
    if not actual or not any(value.removeprefix("PMC") == pmcid.removeprefix("PMC") for value in actual):
        raise _resource_error("pmc_article_identity_mismatch")
    article_pmids = {
        "".join(node.itertext()).strip() for node in identifiers if node.get("pub-id-type") == "pmid"
    }
    if article_pmids and article_pmids != {pmid}:
        raise _resource_error("pmc_parent_pmid_mismatch")


def _pmc_passages(xml, gene, terms, limit, chars, *, report=False):
    root = _main_pmc_article(xml)
    body = next((node for node in root if node.tag.rsplit("}", 1)[-1] == "body"), None)
    license_text = " ".join(
        _xml_text(node) for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "license"
    )
    if body is None:
        return (), license_text
    parents = {child: node for node in root.iter() for child in node}
    gene_pattern = re.compile(r"(?<![A-Za-z0-9])" + re.escape(gene) + r"(?![A-Za-z0-9])", re.I)
    raw = []
    for index, node in enumerate(body.iter()):
        match_content = None
        tag = node.tag.rsplit("}", 1)[-1]
        ancestors, parent = [], parents.get(node)
        while parent is not None:
            ancestors.append(parent)
            parent = parents.get(parent)
        tags = {item.tag.rsplit("}", 1)[-1] for item in ancestors}
        if tags.intersection({"supplementary-material", "ref-list"}):
            continue
        section = next((item for item in ancestors if item.tag.rsplit("}", 1)[-1] == "sec"), None)
        title = _xml_text(section.find("title")) if section is not None else ""
        if re.search(r"references|acknowledg|author contributions|funding|competing interests", title, re.I):
            continue
        if tag == "p" and not tags.intersection({"table-wrap", "fig", "caption"}):
            content, kind = _xml_text(node), "paragraph"
            siblings = list(parents[node]) if node in parents else [node]
            position = siblings.index(node)
            labels = set(
                re.findall(r"\b(?:patient|proband|individual|family)\s+([A-Z0-9][A-Z0-9.-]*)", content, re.I)
            )
            window = []
            for neighbor in siblings[max(0, position - 1) : position + 2]:
                if neighbor.tag.rsplit("}", 1)[-1] != "p":
                    continue
                neighbor_text = _xml_text(neighbor)
                neighbor_labels = set(
                    re.findall(
                        r"\b(?:patient|proband|individual|family)\s+([A-Z0-9][A-Z0-9.-]*)",
                        neighbor_text,
                        re.I,
                    )
                )
                if (
                    neighbor is node
                    or not neighbor_labels
                    or (labels and labels.intersection(neighbor_labels))
                ):
                    window.append(neighbor_text)
            if len(window) > 1:
                content, kind = " ".join(window), "paragraph_window"
        elif tag == "tr":
            content, kind = " | ".join(_xml_text(cell) for cell in node), "table_row"
            match_content = content
            if report:
                table = next((a for a in ancestors if a.tag.rsplit("}", 1)[-1] == "table-wrap"), None)
                if table is not None:
                    context = [
                        _xml_text(a)
                        for a in table.iter()
                        if a.tag.rsplit("}", 1)[-1] in {"label", "caption", "thead", "table-wrap-foot"}
                    ]
                    for row in table.iter("tr"):
                        if row is not node and (
                            row.find("th") is not None
                            or re.match(
                                r"(?:WT(?:\b|\d)|wild[ -]?type)", _xml_text(next(iter(row), None)), re.I
                            )
                        ):
                            context.append(" | ".join(_xml_text(cell) for cell in row))
                    content = "\n".join((content, *dict.fromkeys(context)))
        elif tag == "caption" and "fig" in tags:
            content, kind = _xml_text(node), "figure_caption"
        else:
            continue
        if len(content) < (40 if kind == "table_row" else 80):
            continue
        match_text = match_content if report and match_content is not None else content
        matches = tuple(
            term
            for term in terms
            if (
                re.search(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", match_text, re.I)
                if report
                else term.casefold() in content.casefold()
            )
        )
        gene_match = bool(gene_pattern.search(match_text))
        if not matches and not gene_match:
            continue
        score = (
            6 * len(matches)
            + 4 * gene_match
            + min(3, len(re.findall(r"patient|proband|family|segregat|functional|RNA|trans", content, re.I)))
        )
        score += int(bool(re.search(r"result|clinical|genetic", title, re.I))) * 2
        score += 1.5 if kind == "table_row" else 0.5 if kind == "figure_caption" else 0
        raw.append((score, index, title, content, kind, matches))
    selected, token_sets = [], []
    for score, index, title, content, kind, matches in sorted(raw, key=lambda row: (-row[0], row[1])):
        tokens = set(re.findall(r"\w+", content.casefold()))
        if any(len(tokens & other) / max(1, len(tokens | other)) >= 0.8 for other in token_sets):
            continue
        selected.append(
            {
                "locator": f"passage-{index}",
                "section": title,
                "text": content if report else content[:chars],
                "content_type": kind,
                "query_hits": matches,
            }
        )
        token_sets.append(tokens)
        if len(selected) == limit:
            break
    return tuple(selected), license_text


def context_pmids(source):
    """Read explicit original citation identifiers, never arbitrary numbers."""
    groups = re.findall(r"\bPMIDs?\b\s*:?\s*([0-9][\d,\s]*)", source.content, re.I)
    identifiers = [value for group in groups for value in re.findall(r"\d+", group)]
    identifiers.extend(re.findall(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", source.content))
    try:
        row = json.loads(source.content).get("row", {})
    except (ValueError, AttributeError):
        row = {}
    if isinstance(row, dict):
        identifiers.extend(
            str(value)
            for key, value in row.items()
            if re.fullmatch(r"(?:Haploinsufficiency|Triplosensitivity) PMID\d+", key) and str(value).isdigit()
        )
    return tuple(dict.fromkeys(identifiers))


def report_citations(gene, targets, seeds, contexts):
    keys = {variant_tuple(t["primary"]["key"], MendelPrimeError("variant_key_required")) for t in targets}
    plan = {
        str(s["pmid"]): ("clinvar_citation", "post_allele_citation")
        for s in seeds
        if s.get("gene") == gene
        and variant_tuple(s["variant_key"], MendelPrimeError("variant_key_required")) in keys
        and str(s["pmid"]).isdigit()
    }
    for pmid in _context_citation_origins(gene, contexts):
        plan[pmid] = tuple(dict.fromkeys((*plan.get(pmid, ()), "post_gene_context_citation")))
    return plan


def _context_citation_origins(gene, contexts):
    origins = defaultdict(list)
    for document in contexts:
        source = document["source"]
        if (
            source.provider not in {"clingen", "genereviews"}
            or source.kind != "gene_disease"
            or source.genes != (gene,)
            or source.variant_keys
        ):
            continue
        identifiers = context_pmids(source)
        if document.get("cited_pmids") and tuple(document["cited_pmids"]) != identifiers:
            raise MendelPrimeError("report_citation_binding_invalid", field=gene)
        for identifier in identifiers:
            origins[identifier].append(source)
    return {identifier: tuple(rows) for identifier, rows in origins.items()}


def complete_report_documents(pubmed, gene, evidence, targets, seeds, case_key, *, fetched=None, contexts=()):
    """Complete explicit citations; preserve existing sources and ranking evidence."""
    current = evidence or {"gene": gene, "documents": (), "omim_context": (), "retrieval_errors": ()}
    plan = report_citations(gene, targets, seeds, contexts)
    origins = _context_citation_origins(gene, contexts)
    documents = {doc["source"].source_id: doc for doc in current["documents"]}
    existing = {
        _variant_article_key(doc).removeprefix("PMID:"): doc
        for doc in documents.values()
        if doc["source"].provider == "pubmed"
    }
    missing = tuple(pmid for pmid in plan if pmid not in existing)
    publications, errors = fetch_many(pubmed, missing) if fetched is None and missing else fetched or ({}, {})
    errors = {pmid: error for pmid, error in errors.items() if pmid in missing}
    resolution = {source.source_id: source for source in current.get("citation_resolution", ())}
    for pmid, error in tuple(errors.items()):
        if error != "pubmed_record_missing_or_unsupported" or pmid not in origins:
            continue
        errors[pmid] = "source_citation_identifier_unresolvable"
        for origin in origins[pmid]:
            receipt = _citation_resolution_source(gene, pmid, origin, case_key)
            resolution[receipt.source_id] = receipt
    for pmid, routes in plan.items():
        if pmid in existing:
            doc = {
                **existing[pmid],
                "provenance": tuple(dict.fromkeys((*existing[pmid]["provenance"], *routes))),
            }
        elif pmid in publications:
            doc = _article_document(gene, publications[pmid], routes, None)
        else:
            continue
        documents[doc["source"].source_id] = doc
    expanded = expand_variant_evidence(
        pubmed,
        gene,
        {**current, "documents": tuple(documents.values())},
        targets,
        max_articles=3,
        purpose="report",
        case_key=case_key,
        preferred_articles=frozenset(
            f"PMID:{pmid}" for pmid, routes in plan.items() if "post_gene_context_citation" in routes
        ),
    )
    return {
        **expanded,
        "citation_resolution": tuple(resolution.values()),
        "retrieval_errors": tuple(
            dict.fromkeys(
                (
                    *expanded.get("retrieval_errors", ()),
                    *(f"PMID:{pmid}:{_safe_fetch_error(error)}" for pmid, error in errors.items()),
                )
            )
        ),
    }


def _citation_resolution_source(gene, identifier, origin, case_key):
    content = json.dumps(
        {
            "record_type": "source_citation_resolution_v1",
            "provider": origin.provider,
            "gene": gene,
            "source_id": origin.source_id,
            "identifier_type": "PMID",
            "raw_identifier": identifier,
            "status": "unresolvable",
            "error_code": "source_citation_identifier_unresolvable",
            "not_medical_evidence": True,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return Source(
        provider="mendelprime",
        kind="derived",
        locator=f"method:source-citation-resolution/{gene}/PMID-{identifier}",
        version="source-citation-resolution-v1",
        content=content,
        genes=(gene,),
        case_key=case_tuple(case_key),
        parent_source_ids=(origin.source_id,),
    )


def expand_variant_evidence(
    pubmed,
    gene,
    evidence,
    targets,
    *,
    max_articles=2,
    max_passages=4,
    passage_chars=1800,
    previous=None,
    excluded_pmids=(),
    selected_variant_keys=None,
    case_key=None,
    purpose="ranking",
    preferred_articles=frozenset(),
):
    """Bounded variant queries and PMC passages; never fetch supplementary files.

    Targets are prepared records for one gene. Coverage travels with the returned
    evidence row, so a second pass reuses successful and unavailable outcomes.
    """
    validate_gene(gene)
    if purpose not in {"ranking", "report"}:
        raise ValueError("unknown literature purpose")
    if preferred_articles and purpose != "report":
        raise ValueError("curated article priority is report-only")
    for value, maximum in (
        (max_articles, 3 if purpose == "report" else 2),
        (max_passages, 4),
        (passage_chars, 1800),
    ):
        if type(value) is not int or not 1 <= value <= maximum:
            raise ValueError("ranking literature budget exceeds restored bounds")
    if evidence["gene"] != gene or any(target["primary"]["gene"] != gene for target in targets):
        raise MendelPrimeError("gene_scope_mismatch", field="variant_expansion")
    all_targets = tuple(targets)
    if selected_variant_keys is not None:
        by_key = {
            variant_tuple(
                target["primary"]["key"], MendelPrimeError("variant_key_required", field="variant_expansion")
            ): target
            for target in targets
        }
        if any(key not in by_key for key in selected_variant_keys):
            raise MendelPrimeError("selected_variant_missing", field="variant_expansion")
        targets = tuple(by_key[key] for key in selected_variant_keys)
    case_keys = {target["primary"]["source"].case_key for target in all_targets}
    if case_key is None:
        case_key = evidence.get("variant_expansion", {}).get("case_key")
    if case_key is not None:
        case_keys.add(case_tuple(case_key))
    if len(case_keys) > 1 or None in case_keys:
        raise MendelPrimeError("target_case_scope_required", field="variant_expansion")
    case_key = next(iter(case_keys), None)
    excluded = {_pmid(value) for value in excluded_pmids} - {""}
    documents = [
        document
        for document in evidence["documents"]
        if _variant_article_key(document).removeprefix("PMID:") not in excluded
    ]
    if previous is not None:
        if previous["gene"] != gene:
            raise MendelPrimeError("gene_scope_mismatch", field="variant_expansion.previous")
        for document in previous["documents"]:
            if _variant_article_key(document).removeprefix("PMID:") in excluded:
                continue
            if all(item["source"].source_id != document["source"].source_id for item in documents):
                documents.append(document)
    previous_coverage = tuple((previous or evidence).get("variant_expansion", {}).get("coverage", ()))
    coverage, issues = list(previous_coverage), []
    terms = [gene, "functional assay", "segregation", "RNA", "in trans"]
    terms.extend(
        term for target in all_targets for term in _variant_search_terms(target, report=purpose == "report")
    )
    existing = {_variant_article_key(document) for document in documents}
    pending_pmids = []
    for target in tuple(targets)[: 8 if purpose == "report" else 2]:
        key = variant_tuple(
            target["primary"]["key"], MendelPrimeError("variant_key_required", field="variant_expansion")
        )
        identifiers = _variant_search_terms(target, report=purpose == "report")
        quoted = " OR ".join('"' + term.replace('"', "") + '"[Title/Abstract]' for term in identifiers)
        query = f'"{gene}"[Title/Abstract] AND ({quoted})' if quoted else ""
        row = {
            "resource": "variant_pubmed",
            "identity": {"gene": gene, "variant_key": key},
            "query": query,
            "article_keys": (),
            "status": "no_identifiers",
        }
        prior = next(
            (
                item
                for item in previous_coverage
                if item.get("resource") == "variant_pubmed"
                and item.get("identity") == row["identity"]
                and item.get("query") == query
                and item["status"] != "failed"
            ),
            None,
        )
        if prior is not None:
            continue
        initial = evidence.get("variant_expansion", {})
        if (
            purpose == "ranking"
            and initial.get("retrieval_strategy") == "september_vcf_v1"
            and set(identifiers).intersection(initial.get("variant_search_terms", ()))
            and not any(
                not str(error).endswith("article_or_abstract_missing")
                for error in evidence.get("retrieval_errors", ())
            )
        ):
            row["status"] = "handled_by_initial_retrieval"
            coverage.append(row)
            continue
        if query:
            try:
                hits = tuple(pmid for pmid in search(pubmed, query, max_articles) if pmid not in excluded)
                row.update(
                    status="ok" if hits else "no_results", article_keys=tuple(f"PMID:{pmid}" for pmid in hits)
                )
                pending_pmids.extend(pmid for pmid in hits if f"PMID:{pmid}" not in existing)
            except Exception as error:
                row["status"] = "failed"
                issues.append(_safe_error(error, "variant_search_failed"))
        coverage = [
            item
            for item in coverage
            if not (
                item.get("resource") == "variant_pubmed"
                and item.get("identity") == row["identity"]
                and item.get("query") == query
            )
        ]
        coverage.append(row)
    pmids = tuple(dict.fromkeys(pending_pmids))[:max_articles]
    publications, errors = fetch_many(pubmed, pmids) if pmids else ({}, {})
    for pmid, publication in publications.items():
        documents.append(_article_document(gene, publication, ("variant_targeted",), 3500))
    issues.extend(f"PMID:{pmid}:{_safe_fetch_error(error)}" for pmid, error in errors.items())
    for row in coverage:
        if row.get("resource") == "variant_pubmed" and any(
            key.removeprefix("PMID:") in errors
            and errors[key.removeprefix("PMID:")] != "article_or_abstract_missing"
            for key in row.get("article_keys", ())
        ):
            row["status"] = "failed"
    abstract_documents = [document for document in documents if document["source"].provider == "pubmed"]
    ranked = sorted(
        enumerate(abstract_documents),
        key=lambda row: (
            purpose == "report"
            and not bool(
                {"variant_targeted", "post_allele_citation", "clinvar_variant_exact"}.intersection(
                    row[1]["provenance"]
                )
            ),
            _variant_article_key(row[1]) not in preferred_articles,
            -sum(term.casefold() in row[1]["source"].content.casefold() for term in terms),
            -len(re.findall(r"patients?|famil(?:y|ies)|cohort|case series|phenotyp", row[1]["title"], re.I))
            if _variant_article_key(row[1]) in preferred_articles
            else 0,
            row[0],
        ),
    )
    selected = []
    for _, document in ranked:
        key = _variant_article_key(document)
        if key and key not in {item[0] for item in selected}:
            selected.append((key, document))
        if len(selected) == max_articles:
            break
    for article_key, abstract in selected:
        query_terms = tuple(sorted(set(terms)))
        prior = next(
            (
                row
                for row in previous_coverage
                if row.get("resource") == "pmc"
                and row.get("article_key") == article_key
                and row["status"] != "failed"
                and (row["status"] == "not_available" or row.get("query_terms") == query_terms)
                and (purpose != "report" or row.get("extraction_version") == "report-table-context-v2")
            ),
            None,
        )
        if prior is not None:
            continue
        row = {
            "resource": "pmc",
            "article_key": article_key,
            "query_terms": query_terms,
            "status": "not_available",
            "source_ids": (),
            **({"extraction_version": "report-table-context-v2"} if purpose == "report" else {}),
        }
        try:
            article = _pmc_resource(pubmed, article_key.removeprefix("PMID:"))
            if article["status"] == "found":
                passages, license_text = _pmc_passages(
                    article["xml"],
                    gene,
                    tuple(terms),
                    max_passages,
                    passage_chars,
                    **({"report": True} if purpose == "report" else {}),
                )
                if purpose == "report" and passages:
                    documents = [
                        d
                        for d in documents
                        if not (d["source"].provider == "pmc" and _variant_article_key(d) == article_key)
                    ]
                ids = []
                for passage in passages:
                    source = Source(
                        provider="pmc",
                        kind="literature",
                        locator=f"https://pmc.ncbi.nlm.nih.gov/articles/{article['pmcid']}/#{passage['locator']}",
                        version="sha256:"
                        + hashlib.sha256(article["xml"].encode()).hexdigest()
                        + (";report-table-context-v2" if purpose == "report" else ""),
                        content=passage["text"],
                        genes=(gene,),
                    )
                    fulltext = evidence_document(
                        source,
                        f"PMC:{article['pmcid'].removeprefix('PMC')}#{passage['locator']}",
                        abstract["title"],
                        passage["text"],
                        provenance=("verification_search", "pmc_fulltext"),
                        query_hits=passage["query_hits"],
                        section=passage["section"],
                    )
                    fulltext["metadata"] = {
                        "pmcid": article["pmcid"],
                        "parent_pmid": article_key,
                        "parent_source_id": abstract["source"].source_id,
                        "content_type": passage["content_type"],
                        "license": license_text,
                    }
                    if all(item["source"].source_id != source.source_id for item in documents):
                        documents.append(fulltext)
                    ids.append(source.source_id)
                    link_content = json.dumps(
                        {
                            "record_type": "pmc_article_lineage_v1",
                            "pmid": article_key,
                            "pmcid": article["pmcid"],
                            "abstract_source_id": abstract["source"].source_id,
                            "passage_source_id": source.source_id,
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    link = Source(
                        provider="pmc",
                        kind="derived" if case_key is not None else "annotation",
                        locator=f"method:vcf-ranking/pmc-lineage/{article['pmcid']}/{passage['locator']}",
                        version="pmc-article-lineage-v1"
                        + (";report-table-context-v2" if purpose == "report" else ""),
                        content=link_content,
                        genes=(gene,),
                        case_key=case_key,
                        parent_source_ids=(abstract["source"].source_id, source.source_id)
                        if case_key is not None
                        else (),
                    )
                    if all(item["source"].source_id != link.source_id for item in documents):
                        documents.append(
                            evidence_document(
                                link,
                                f"PMC_LINK:{article['pmcid']}#{passage['locator']}",
                                "Source-bound PMC article lineage",
                                link_content,
                                provenance=("pmc_article_lineage",),
                                section="article_lineage",
                            )
                        )
                    ids.append(link.source_id)
                row.update(status="expanded" if passages else "no_relevant_passages", source_ids=tuple(ids))
        except Exception as error:
            row["status"] = "failed"
            issues.append(_safe_error(error, "pmc_expansion_failed"))
        coverage = [
            item
            for item in coverage
            if not (item.get("resource") == "pmc" and item.get("article_key") == article_key)
        ]
        coverage.append(row)
    registry = {}
    for document in (*documents, *evidence.get("omim_context", ())):
        source = document["source"]
        if source.genes and gene not in source.genes:
            raise MendelPrimeError("document_gene_scope_mismatch", field="variant_expansion")
        previous_source = registry.setdefault(source.source_id, source)
        if previous_source != source:
            raise MendelPrimeError("source_id_collision", field="variant_expansion")
    if any(parent not in registry for source in registry.values() for parent in source.parent_source_ids):
        raise MendelPrimeError("parent_source_missing", field="variant_expansion")
    merged = {
        **evidence,
        "documents": tuple(documents),
        "retrieval_errors": tuple(dict.fromkeys((*evidence.get("retrieval_errors", ()), *issues))),
    }
    return {
        **merged,
        "variant_expansion": {
            **evidence.get("variant_expansion", {}),
            "coverage": tuple(coverage),
            "issues": tuple(dict.fromkeys(issues)),
        },
    }


def retrieve_reflection(pubmed, genes, queries, excluded_pmids):
    excluded = {_pmid(value) for value in excluded_pmids} - {""}
    plans, all_pmids = {}, []
    for gene in genes:
        validate_gene(gene)
        counts, matched, errors = defaultdict(int), defaultdict(list), []
        for phrase in queries.get(gene, ()):
            try:
                pmids = search(pubmed, reflection_query(gene, phrase), 10)
            except Exception as error:
                errors.append(f"reflection_search:{_safe_error(error, 'pubmed_reflection_search_failed')}")
                pmids = ()
            for value in pmids:
                pmid = _pmid(value)
                if not pmid or pmid in excluded:
                    continue
                counts[pmid] += 1
                if phrase not in matched[pmid]:
                    matched[pmid].append(phrase)
                if pmid not in all_pmids:
                    all_pmids.append(pmid)
            pubmed.get("sleep", time.sleep)(0.34)
        plans[gene] = counts, matched, errors

    publications, fetch_errors = fetch_many(pubmed, tuple(all_pmids))
    output = {}
    for gene in genes:
        counts, matched, retrieval_errors = plans[gene]
        retrieval_errors = list(retrieval_errors)
        for pmid in counts:
            if pmid not in fetch_errors:
                continue
            code = _safe_fetch_error(fetch_errors[pmid])
            if pmid not in publications and code != "article_or_abstract_missing":
                retrieval_errors.append(f"PMID:{pmid}:{code}")
        candidates = [
            publications[pmid]
            for pmid in counts
            if pmid in publications and publications[pmid]["abstract"].strip()
        ]
        candidates.sort(key=lambda article: _reflection_priority(article, gene, counts))
        documents = tuple(
            _article_document(
                gene,
                article,
                ("verification_search", "top2_reflection"),
                None,
                tuple(matched[article["pmid"]]),
            )
            for article in candidates[:6]
        )
        output[gene] = {
            "gene": gene,
            "documents": documents,
            "omim_context": (),
            "retrieval_errors": tuple(retrieval_errors),
        }
    return output


def evidence_queries(gene, labels):
    specific = " OR ".join(f'"{label}"' for label in labels[:4] if label)
    broader = " OR ".join(f'"{label}"' for label in labels[:10] if label)
    clinical = "(clinical OR patient OR patients OR syndrome OR phenotype OR case)"
    return tuple(
        query
        for query in (
            f'("{gene}"[Title/Abstract]) AND ({specific}) AND {clinical}' if specific else "",
            f'("{gene}"[Title/Abstract]) AND ({broader}) AND {clinical}' if broader else "",
            f'("{gene}"[Title/Abstract]) AND {clinical}',
        )
        if query
    )


def reflection_query(gene, phrase):
    safe = phrase.replace('"', "")[:180]
    parts = tuple(
        dict.fromkeys(
            value.strip()
            for value in re.split(r"\b(?:with|and|plus|versus)\b|[,;/]", safe, flags=re.IGNORECASE)
            if value.strip()
        )
    )[:3]
    if len(parts) > 1:
        clinical = " AND ".join(f'"{value}"[Title/Abstract]' for value in parts)
    else:
        words = tuple(
            dict.fromkeys(
                word
                for word in re.findall(r"[A-Za-z0-9-]+", safe)
                if len(word) > 2 and word.lower() not in {"the", "for", "from", "without"}
            )
        )[:6]
        clinical = " OR ".join(f'"{word}"[Title/Abstract]' for word in words)
    return (
        f'("{gene}"[Title/Abstract]) AND ({clinical}) AND '
        "(patient OR patients OR clinical OR phenotype OR syndrome OR case OR family)"
    )


def _fetch(state, pmids):
    if not pmids:
        return {}
    body = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}).encode("ascii")
    return _network_bytes(
        state, urllib.request.Request(FETCH_URL, data=body, method="POST"), _MAX_FETCH_BYTES,
        decode=lambda content: parse_pubmed_xml(content.decode("utf-8")),
    )


def _network_bytes(state, request, limit, *, decode=None):
    opener = state.get("urlopen", urlopen)
    sleeper = state.get("sleep", time.sleep)
    invalid_retried = False
    for attempt in range(6):
        try:
            _wait_for_request_turn(state)
            with opener(request, timeout=state["timeout_seconds"]) as response:
                content = response.read(limit + 1)
            if not isinstance(content, bytes):
                raise _resource_error("pubmed_invalid_response")
            if len(content) > limit:
                raise _resource_error("pubmed_response_too_large")
            return decode(content) if decode is not None else content
        except (ValueError, ET.ParseError, RecursionError) as error:
            if isinstance(error, MendelPrimeError) and error.code != "pubmed_invalid_response":
                raise
            if invalid_retried or attempt == 5:
                raise _resource_error("pubmed_invalid_response", attempts=attempt + 1) from error
            invalid_retried = True
            sleeper(2)
        except urllib.error.HTTPError as error:
            retryable = error.code in _RETRYABLE_HTTP
            if retryable and attempt < 5:
                delay = _http_retry_delay(error, attempt)
                error.close()
                if error.code == 429:
                    _defer_requests(state, delay)
                else:
                    sleeper(delay)
                continue
            error.close()
            raise _resource_error(
                "pubmed_http_error", attempts=attempt + 1, http_status=error.code, retryable=retryable
            ) from error
        except (OSError, TimeoutError, urllib.error.URLError) as error:
            reason = error.reason if isinstance(error, urllib.error.URLError) else error
            if isinstance(reason, (ssl.SSLCertVerificationError, ssl.SSLError)):
                raise _resource_error("pubmed_tls_error", attempts=attempt + 1, retryable=False) from error
            if attempt == 0:
                sleeper(2)
                continue
            raise _resource_error("pubmed_transport_error", attempts=attempt + 1, retryable=True) from error
        except Exception as error:
            raise _resource_error("pubmed_transport_error", attempts=attempt + 1) from error
    raise _resource_error("pubmed_transport_error", attempts=6, retryable=True)


def _wait_for_request_turn(state):
    lock = state.get("request_lock")
    if lock is None:
        return
    clock = state.get("monotonic", time.monotonic)
    with lock:
        now = clock()
        delay = max(0.0, state.get("next_request_at", 0.0) - now)
        if delay:
            state.get("rate_sleep", time.sleep)(delay)
            now = clock()
        state["next_request_at"] = max(now, state.get("next_request_at", 0.0)) + _NCBI_INTERVAL_SECONDS


def _defer_requests(state, seconds):
    lock = state.get("request_lock")
    if lock is None:
        state.get("sleep", time.sleep)(seconds)
        return
    clock = state.get("monotonic", time.monotonic)
    with lock:
        state["next_request_at"] = max(state.get("next_request_at", 0.0), clock() + seconds)


def _http_retry_delay(error, attempt):
    value = (getattr(error, "headers", None) or {}).get("Retry-After", "")
    return (
        min(int(value), _MAX_RETRY_DELAY_SECONDS)
        if value.isdigit()
        else min(2 ** (attempt + 1), _MAX_RETRY_DELAY_SECONDS)
    )


def _read_cache(state, relative):
    path = state["cache_root"] / relative
    try:
        private_directory(state["cache_root"], create=False)
        private_directory(path.parent, create=False)
        raw, _before, _after = read_regular_bounded(path, _MAX_CACHE_BYTES, follow_symlinks=False)
        return json.loads(raw.decode("utf-8"))
    except (OSError, OverflowError, UnicodeError, json.JSONDecodeError):
        return None


def _write_cache(state, relative, value):
    target = state["cache_root"] / relative
    try:
        private_directory(state["cache_root"])
        private_directory(target.parent)
        content = json.dumps(value, ensure_ascii=False).encode("utf-8")
        atomic_replace_bytes(target, content, sync=False)
    except Exception as error:
        raise _resource_error("pubmed_cache_write_failed") from error


def _publication_from_cache(value):
    if not isinstance(value, dict) or not isinstance(value.get("publication_types", ()), (list, tuple)):
        return None
    publication = {
        "pmid": str(value.get("pmid") or ""),
        "title": str(value.get("title") or ""),
        "abstract": str(value.get("abstract") or ""),
        "journal": str(value.get("journal") or ""),
        "year": str(value.get("year") or ""),
        "publication_types": tuple(str(item) for item in value.get("publication_types", ())),
    }
    return publication if publication["pmid"].isdigit() and publication["abstract"].strip() else None


def _publication_source(publication, gene):
    snapshot = json.dumps(publication, ensure_ascii=True, sort_keys=True).encode()
    return _source(
        "pubmed",
        "literature",
        f"https://pubmed.ncbi.nlm.nih.gov/{publication['pmid']}/",
        "sha256:" + hashlib.sha256(snapshot).hexdigest(),
        publication["abstract"],
        gene,
    )


def _article_document(gene, publication, provenance, prompt_limit, query_hits=()):
    prompt_text = " ".join(publication["abstract"].split())
    if prompt_limit is not None:
        prompt_text = prompt_text[:prompt_limit]
    return evidence_document(
        _publication_source(publication, gene),
        f"PMID:{publication['pmid']}",
        publication["title"],
        prompt_text,
        provenance=provenance,
        publication_types=publication["publication_types"],
        query_hits=query_hits,
    )


def evidence_document(
    source,
    prompt_id,
    title,
    prompt_text,
    *,
    provenance=(),
    publication_types=(),
    query_hits=(),
    cited_pmids=(),
    section="",
):
    return {
        "source": source,
        "prompt_id": prompt_id,
        "title": title,
        "prompt_text": prompt_text,
        "provenance": tuple(provenance),
        "publication_types": tuple(publication_types),
        "query_hits": tuple(query_hits),
        "cited_pmids": tuple(cited_pmids),
        "section": section,
    }


def _select_documents(documents, limit):
    selected, seen = [], set()
    for route, maximum in (
        ("must_review", 1),
        ("amelie", 3),
        ("omim", 3),
        ("clinvar_citation", None),
        ("clingen_citation", None),
        ("pubmed_search", None),
    ):
        added = 0
        for document in documents:
            if len(selected) == limit or maximum is not None and added == maximum:
                break
            if route in document["provenance"] and document["prompt_id"] not in seen:
                selected.append(document)
                seen.add(document["prompt_id"])
                added += 1
    for document in documents:
        if len(selected) == limit:
            break
        if document["prompt_id"] not in seen:
            selected.append(document)
            seen.add(document["prompt_id"])
    return tuple(sorted(selected, key=lambda item: item["prompt_id"].removeprefix("PMID:")))


def _reflection_priority(article, gene, hits):
    text = f"{article['title']} {article['abstract']}"
    mentioned = bool(re.search(rf"(?<![A-Z0-9]){re.escape(gene)}(?![A-Z0-9])", text, re.IGNORECASE))
    patient = bool(
        {value.lower() for value in article["publication_types"]}
        & {"case reports", "clinical study", "observational study"}
    )
    return -hits.get(article["pmid"], 0), -int(mentioned), -int(patient), article["pmid"]


def _scored_papers(values, excluded):
    candidates = values if excluded else values[:5]
    selected = (
        (_pmid(_paper(item)[0]) if excluded else str(_paper(item)[0]).strip(), _paper(item)[1])
        for item in candidates
    )
    return tuple(
        sorted(
            ((pmid, score) for pmid, score in selected if not excluded or pmid and pmid not in excluded),
            key=lambda item: item[1],
            reverse=True,
        )[:5]
    )


def _paper(value):
    return (value.get("pmid", ""), value.get("score")) if isinstance(value, Mapping) else value[:2]


def _append_pmids(target, values, limit):
    for value in values if isinstance(values, (list, tuple)) else ():
        pmid = str(value).strip()
        if pmid.isdigit() and pmid not in target:
            target.append(pmid)
        if len(target) == limit:
            return


def _cached_pmids(value, query):
    if not isinstance(value, dict) or value.get("query") != query:
        return None
    raw = value.get("pmids")
    if not isinstance(raw, list) or any(not str(item).isdigit() for item in raw):
        return None
    return tuple(str(item) for item in raw)


def _unavailable(pmid, timestamp):
    return {
        "status": "unavailable",
        "pmid": pmid,
        "reason": "article_or_abstract_missing",
        "timestamp": timestamp,
    }


def _cached_unavailable(value, pmid, now):
    if not isinstance(value, dict):
        return False
    timestamp = value.get("timestamp")
    return bool(
        value.get("status") == "unavailable"
        and value.get("pmid") == pmid
        and pmid.isdigit()
        and value.get("reason") == "article_or_abstract_missing"
        and type(timestamp) in {int, float}
        and math.isfinite(timestamp)
        and 0 <= timestamp <= now
        and now - timestamp < UNAVAILABLE_TTL_SECONDS
    )


def _safe_error(error, fallback):
    return error.code if isinstance(error, MendelPrimeError) and error.category == "resource" else fallback


def _safe_fetch_error(error):
    return error if isinstance(error, str) and error in _SAFE_FETCH_ERRORS else "pubmed_fetch_failed"


def _resource_error(code, **metadata):
    return MendelPrimeError(code, category="resource", source="pubmed", **metadata)


def _source(provider, kind, locator, version, content, gene):
    return Source(
        provider=provider, kind=kind, locator=locator, version=version, content=content, genes=(gene,)
    )


def _omim_title(record, omim_id):
    value = str(record.get("page_title") or "").replace("Entry - ", "")
    return value.replace(" - OMIM - (OMIM.ORG)", "").strip() or f"OMIM {omim_id}"


def _xml_text(node):
    return " ".join("".join(node.itertext()).split()) if node is not None else ""


def _hash(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _require_pmid(value):
    pmid = str(value).strip()
    if not pmid.isdigit():
        raise ValueError("numeric PMID required")
    return pmid


def _pmid(value):
    text = str(value).strip()
    text = text[5:].strip() if text.upper().startswith("PMID:") else text
    return text if text.isdigit() else ""
