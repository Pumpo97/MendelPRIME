"""Phenotype profiles, resources, and native candidate acquisition."""

import errno
import hashlib
import hmac
import http.client
import io
import json
import logging
import math
import os
import re
import secrets
import socket
import ssl
import time
from pathlib import Path
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

from mendelprime import domain
from mendelprime.knowledge.hpo import ancestors, consulted_terms, depth
from mendelprime.runtime.http import tls_context
from mendelprime.runtime.storage import (
    atomic_replace_bytes,
    open_regular,
    private_directory,
    read_regular_bounded,
    stat_fingerprint,
)

MAX_RESPONSE_BYTES, _MAX_CACHE_BYTES = 8 * 1024 * 1024, 64 * 1024 * 1024
_P2G_URL, _AMELIE_URL = ("https://phen2gene.wglab.org/api", "https://amelie.stanford.edu/api/gene_list_api/")
_PSEUDONYM_DOMAIN = b"mendelprime.amelie.patient.v1\0"
_GENE = re.compile(r"[A-Z0-9][A-Z0-9._-]*", re.ASCII)
_P2D_FILES = (
    "data/matrix/lin_similarity_matrix.json",
    "data/association/Disease-Hpo/disease2hpo20210413.json",
    "data/association/Gene-Disease/disease2genecard2021.json",
    "data/obo/hp20210413.obo",
)
_RETRY_HTTP = {408, 409, 429, 500, 502, 503, 504}
_RETRY_ERRNO = {
    getattr(errno, name)
    for name in (
        "ECONNABORTED",
        "ECONNREFUSED",
        "ECONNRESET",
        "EHOSTDOWN",
        "EHOSTUNREACH",
        "ENETDOWN",
        "ENETRESET",
        "ENETUNREACH",
        "EPIPE",
        "ETIMEDOUT",
    )
    if hasattr(errno, name)
}


def require_positive_hpo(patient):
    domain.require(
        any(term["present"] for term in patient.phenotypes), "phenotypes", "positive_hpo_required_for_ranking"
    )


def build_profile(patient, hpo):
    terms = []
    for observation in patient.phenotypes:
        identifier = observation["hpo_id"]
        record = hpo["terms"].get(identifier)
        ic, level = hpo["gene_ic"].get(identifier, 0.0), depth(hpo, identifier)
        parents = record["parent_ids"] if record else ()
        terms.append(
            {
                "hpo_id": identifier,
                "label": observation["label"]
                or (record["name"] if record and record["name"] else identifier),
                "present": observation["present"],
                "information_content": ic,
                "depth": level,
                "priority_score": 0.55 * ic + 0.45 * min(level / 10.0, 1.0),
                "role": "supportive",
                "definition": record["definition"] if record else "",
                "parent_ids": parents,
                "parent_labels": tuple(hpo["terms"].get(parent, {}).get("name", "") for parent in parents),
                "ontology_known": record is not None,
            }
        )
    present = [term for term in terms if term["present"]]
    eligible = [
        term
        for term in present
        if term["hpo_id"] not in {"HP:0000001", "HP:0000005", "HP:0000118"}
        and not any(
            other["hpo_id"] != term["hpo_id"] and term["hpo_id"] in ancestors(hpo, other["hpo_id"])
            for other in present
        )
    ]

    def ordering(term):
        return -term["priority_score"], term["label"]

    selected = [term for term in eligible if term["priority_score"] >= 0.42]
    if not selected and eligible:
        selected = sorted(eligible, key=ordering)[:1]
    prioritized = {term["hpo_id"] for term in sorted(selected, key=ordering)[:6]}
    for term in terms:
        if term["hpo_id"] in prioritized:
            term["role"] = "specificity_prioritized"
    return domain.freeze_record(
        {
            "case_key": patient.key,
            "terms": terms,
            "clinical_description": patient.clinical_description,
            "ontology_manifest": hpo["ontology_manifest"],
            "statistics_manifest": hpo["statistics_manifest"],
            "consulted_hpo_ids": consulted_terms(hpo, tuple(term["hpo_id"] for term in terms)),
        }
    )


def _resource_error(code, **values):
    return domain.MendelPrimeError(code, category="resource", **values)


def _check(condition, error):
    if not condition:
        raise error


def _resource(path, *, as_json, errors="strict"):
    value = digest = failure = None
    try:
        with open_regular(path) as (stream, before):
            identity = stat_fingerprint(before)
            checksum = hashlib.sha256()
            while chunk := stream.read(1024 * 1024):
                checksum.update(chunk)
            _check(stat_fingerprint(os.fstat(stream.fileno())) == identity, OSError())
            digest = checksum.hexdigest()
            stream.seek(0)
            reader = io.TextIOWrapper(stream, encoding="utf-8", errors=errors)
            try:
                value = json.load(reader) if as_json else reader.read()
            finally:
                reader.detach()
            _check(stat_fingerprint(os.fstat(stream.fileno())) == identity, OSError())
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError, ValueError):
        failure = _resource_error("phenotype_resource_load_failed", source=str(path))
    if failure:
        raise failure
    return value, digest


def _parents(text):
    output, current = {}, ""
    for line in text.splitlines():
        if line == "[Term]":
            current = ""
        elif line.startswith("id:"):
            current = line.split(":", 1)[1].strip()
            output.setdefault(current, [])
        elif current and line.startswith("is_a:"):
            output[current].append(line.split(":", 1)[1].split("!", 1)[0].strip())
    return {term: tuple(values) for term, values in output.items()}


def _annotations(text):
    output = {}
    lines = iter(text.splitlines())
    header = tuple(item.strip().lower() for item in next(lines, "").split("\t"))
    gene_at = next(
        (header.index(name) for name in ("gene_symbol", "gene symbol", "symbol", "gene") if name in header), 1
    )
    hpo_at = next((header.index(name) for name in ("hpo_id", "hpo-id", "hpo id", "hpo") if name in header), 2)
    for line in lines:
        parts = line.split("\t")
        if len(parts) > max(gene_at, hpo_at):
            gene, term = parts[gene_at].strip().upper(), parts[hpo_at].strip()
            if gene and gene != "-" and term.startswith("HP:"):
                output.setdefault(gene, set()).add(term)
    return {gene: frozenset(terms) for gene, terms in output.items()}


def load_phenotype_resources(config):
    resources = config.get("resources", config)
    root, overlap_path = resources.get("phen2disease_root"), resources.get("genes_to_phenotype")
    _check(
        isinstance(root, Path) and isinstance(overlap_path, Path),
        _resource_error("phenotype_resources_required", field="resources"),
    )
    payloads, versions = [], []
    for relative in _P2D_FILES[:3]:
        payload, digest = _resource(root / relative, as_json=True)
        payloads.append(payload)
        versions.append((relative, digest))
    parent_text, digest = _resource(root / _P2D_FILES[3], as_json=False)
    versions.append((_P2D_FILES[3], digest))
    overlap_text, overlap_digest = _resource(overlap_path, as_json=False, errors="replace")
    similarity, raw_disease_hpo, raw_disease_genes = payloads
    disease_hpo = {
        str(disease): tuple(term for term in terms if str(term).startswith("HP:"))
        for disease, terms in raw_disease_hpo.items()
    }
    disease_genes = {
        str(disease): tuple(
            str(gene).upper() for gene in genes if str(gene).strip() and str(gene).strip() != "-"
        )
        for disease, genes in raw_disease_genes.items()
    }
    parents = _parents(parent_text)
    children = {}
    for child, values in parents.items():
        for parent in values:
            children.setdefault(parent, set()).add(child)
    descendants = {}
    for term in set(parents) | set(children):
        found, pending = set(), list(children.get(term, ()))
        while pending:
            child = pending.pop()
            if child not in found:
                found.add(child)
                pending.extend(children.get(child, ()))
        descendants[term] = found
    counts, total = {}, max(len(disease_hpo), 1)
    for terms in disease_hpo.values():
        for term in set(terms):
            counts[term] = counts.get(term, 0) + 1
    ic = {term: -math.log2(count / total) for term, count in counts.items()}
    return {
        "similarity": similarity,
        "disease_hpo": disease_hpo,
        "disease_genes": disease_genes,
        "descendants": descendants,
        "disease_ic": ic,
        "epsilon": max(ic.values(), default=1.0),
        "phen2disease_versions": tuple(versions),
        "annotations": _annotations(overlap_text),
        "overlap_versions": (("hpo_annotations", overlap_digest),),
    }


def _positive(profile):
    values = tuple(term["hpo_id"] for term in profile["terms"] if term["present"])
    _check(values, _resource_error("positive_hpo_required"))
    return values


def _limit(value):
    _check(type(value) is int and value >= 1, ValueError("limit must be a positive integer"))
    return value


def _version(method, resources):
    encoded = json.dumps(resources, separators=(",", ":")).encode("ascii")
    return f"{method}+resources.sha256:{hashlib.sha256(encoded).hexdigest()}"


def _selected(rows, limit):
    output = {}
    for raw_gene, score in rows:
        gene = raw_gene.strip().upper()
        if gene:
            output.setdefault(gene, score)
    return list(output.items())[:limit]


def _tool_record(tool_id, profile, entries, content, version):
    genes = tuple(item["gene"] for item in entries)
    domain.Ranking(genes=genes)
    _check(entries, _resource_error("empty_tool_ranking"))
    source = domain.Source(
        provider=tool_id,
        kind="tool_observation",
        locator=f"method:{tool_id}/ranking",
        version=version,
        content=content,
        genes=genes,
        case_key=profile["case_key"],
    )
    return domain.freeze_record({"tool_id": tool_id, "entries": entries, "sources": (source,), "issues": ()})


def rank_overlap(state, profile, limit=100):
    query = tuple(dict.fromkeys(_positive(profile)))
    query_set, rows, details = set(query), [], {}
    for gene, terms in state["annotations"].items():
        matched = tuple(sorted(query_set & terms))
        score = len(matched) + len(matched) / max(len(query_set), 1)
        rows.append((gene, score))
        details[gene] = {
            "hpo_ids": matched,
            "query_term_count": len(query_set),
            "gene_term_count": len(terms),
        }
    rows.sort(key=lambda item: (-item[1], item[0]))
    selected = _selected(rows, _limit(limit))
    entries = tuple({"gene": gene, "score": score, "context": details.get(gene)} for gene, score in selected)
    content = json.dumps(
        {"queryHPO": list(query), "resource_versions": state["overlap_versions"], "rows": entries},
        ensure_ascii=True,
        separators=(",", ":"),
    )
    return _tool_record(
        "local_overlap", profile, entries, content, _version("hpo-overlap-v1", state["overlap_versions"])
    )


def _specific(state, terms):
    values = set(terms)
    return tuple(term for term in terms if not state["descendants"].get(term, set()) & (values - {term}))


def _weighted(state, pairs):
    values = [
        (
            max(
                (
                    1.0 if term == other else float(state["similarity"].get(term, {}).get(other, 0.0))
                    for other in candidates
                ),
                default=0.0,
            ),
            state["disease_ic"].get(term, state["epsilon"]),
        )
        for term, candidates in pairs
    ]
    return sum(score * weight for score, weight in values) / (
        sum(weight for _, weight in values) or state["epsilon"]
    )


def rank_phen2disease(state, profile, limit=100):
    query = _positive(profile)
    valid = tuple(term for term in query if term in state["similarity"])
    patient = _specific(state, valid) or valid
    disease_scores = {
        disease: _weighted(state, ((term, terms) for term in patient))
        + _weighted(state, (*((term, terms) for term in patient), *((term, patient) for term in terms)))
        for disease, raw in state["disease_hpo"].items()
        if (terms := _specific(state, tuple(term for term in raw if term in state["similarity"])))
    }
    scores, diseases = {}, {}
    for disease, genes in state["disease_genes"].items():
        score = disease_scores.get(disease, 0.0)
        for gene in genes:
            if score > scores.get(gene, float("-inf")):
                scores[gene], diseases[gene] = score, disease
    rows = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    details = {
        gene: {"disease_id": diseases.get(gene, ""), "integrated_score": round(score, 6)}
        for gene, score in rows[: _limit(limit)]
    }
    selected = _selected(rows, limit)
    entries = tuple({"gene": gene, "score": score, "context": details.get(gene)} for gene, score in selected)
    content = json.dumps(
        {"queryHPO": list(query), "resource_versions": state["phen2disease_versions"], "rows": entries},
        ensure_ascii=True,
        separators=(",", ":"),
    )
    return _tool_record(
        "phen2disease", profile, entries, content, _version("phen2disease-v1", state["phen2disease_versions"])
    )


def _canonical(value):
    return json.dumps(value, ensure_ascii=True, allow_nan=False, sort_keys=True, separators=(",", ":"))


def _ca_bundle(external):
    value = external.get("ca_bundle")
    _check(
        value is None or (isinstance(value, Path) and value.is_file()),
        _resource_error("ca_bundle_unavailable", field="external.ca_bundle"),
    )
    return value


def _urlopen(request, timeout, ca_bundle, verify_tls):
    context = tls_context() or ssl.create_default_context(cafile=str(ca_bundle) if ca_bundle else None)
    if not verify_tls:
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE
    opener = urlrequest.OpenerDirector()
    for handler in (
        urlrequest.ProxyHandler(),
        urlrequest.UnknownHandler(),
        urlrequest.HTTPHandler(),
        urlrequest.HTTPDefaultErrorHandler(),
        urlrequest.HTTPSHandler(context=context),
        urlrequest.HTTPErrorProcessor(),
    ):
        opener.add_handler(handler)
    return opener.open(request, timeout=timeout)


def _failure(error):
    if isinstance(error, urlerror.HTTPError):
        return "tool_http_error", error.code, error.code in _RETRY_HTTP
    reason = error.reason if isinstance(error, urlerror.URLError) else error
    if isinstance(reason, (ssl.CertificateError, ssl.SSLError)):
        return "tool_tls_error", None, False
    if isinstance(reason, (TimeoutError, socket.timeout)) or (
        isinstance(reason, OSError) and reason.errno == getattr(errno, "ETIMEDOUT", None)
    ):
        return "tool_timeout", None, True
    retryable = isinstance(
        reason, (ConnectionError, socket.gaierror, http.client.IncompleteRead, http.client.RemoteDisconnected)
    ) or (isinstance(reason, OSError) and reason.errno in _RETRY_ERRNO)
    return ("tool_connection_error" if retryable else "tool_transport_error"), None, retryable


def _send_json(request, timeout, retries, ca_bundle, verify_tls):
    for attempt in range(retries + 1):
        error = payload = code = None
        try:
            with _urlopen(request, timeout, ca_bundle, verify_tls) as response:
                raw = response.read(MAX_RESPONSE_BYTES + 1)
            code = (
                "invalid_tool_response"
                if not isinstance(raw, bytes)
                else "tool_response_too_large"
                if len(raw) > MAX_RESPONSE_BYTES
                else "empty_tool_response"
                if not raw
                else None
            )
            if code is None:
                try:
                    payload = json.loads(
                        raw.decode("utf-8"), parse_constant=lambda _value: _check(False, ValueError())
                    )
                except (UnicodeError, json.JSONDecodeError, RecursionError, ValueError):
                    code = "invalid_tool_response_json"
        except Exception as caught:
            error = caught
        if error is not None:
            code, status, retryable = _failure(error)
            if isinstance(error, urlerror.HTTPError):
                error.close()
            if retryable and attempt < retries:
                time.sleep(min(2**attempt, 8))
                continue
            raise _resource_error(code, attempts=attempt + 1, http_status=status, retryable=retryable)
        _check(not code, _resource_error(code or "invalid_tool_response", attempts=attempt + 1))
        return payload
    raise _resource_error("tool_transport_error", attempts=retries + 1, retryable=True)


def _score(value):
    _check(not isinstance(value, bool), _resource_error("invalid_tool_response"))
    try:
        score = float(value)
    except (TypeError, ValueError, OverflowError):
        score = math.nan
    _check(math.isfinite(score), _resource_error("invalid_tool_response"))
    return score


def _response_gene(value):
    gene = str(value or "").strip().upper()
    _check(not gene or _GENE.fullmatch(gene), _resource_error("invalid_tool_response"))
    return gene or None


def _external_record(tool, profile, entries, endpoint, hpo_ids, genes, payload, **request_values):
    request = {"endpoint": endpoint, "hpo_ids": sorted(set(hpo_ids)), **request_values}
    if genes:
        request["genes"] = list(genes)
    content = _canonical({"semantic_request": request, "response": payload})
    return _tool_record(
        tool,
        profile,
        entries,
        content,
        "phen2gene-response-v1" if tool == "phen2gene" else "amelie-response-v1",
    )


def rank_phen2gene(profile, limit=100, *, external):
    hpo_ids, maximum = _positive(profile), _limit(limit)
    acknowledged = external.get("phen2gene_acknowledged", False)
    _check(type(acknowledged) is bool, TypeError("phen2gene_acknowledged must be bool"))
    _check(acknowledged, _resource_error("external_send_not_acknowledged"))
    ca = _ca_bundle(external)
    query = urlparse.urlencode({"HPO_list": ";".join(hpo_ids), "weight_model": "sk"})
    request = urlrequest.Request(f"{_P2G_URL}?{query}", headers={"Accept": "application/json"}, method="GET")
    payload = _send_json(request, 120.0, 0, ca, True)
    rows = payload.get("results") if isinstance(payload, dict) else None
    _check(isinstance(rows, list), _resource_error("invalid_tool_response"))
    entries, seen = [], set()
    for row in rows:
        _check(isinstance(row, dict), _resource_error("invalid_tool_response"))
        gene = _response_gene(row.get("Gene") or row.get("gene"))
        if gene and gene not in seen:
            seen.add(gene)
            entries.append(
                {"gene": gene, "score": _score(row.get("Score") or row.get("score") or 0.0), "context": None}
            )
    return _external_record(
        "phen2gene", profile, tuple(entries[:maximum]), _P2G_URL, hpo_ids, (), payload, weight_model="sk"
    )


def _amelie_entries(payload, limit):
    _check(isinstance(payload, list), _resource_error("invalid_tool_response"))
    order, first_scores, latest = [], {}, {}
    for item in payload:
        _check(
            isinstance(item, list) and len(item) >= 2 and isinstance(item[1], list),
            _resource_error("invalid_tool_response"),
        )
        gene = _response_gene(item[0])
        if gene is None:
            continue
        papers = tuple(
            {"pmid": str(paper[0]).strip(), "score": _score(paper[1])}
            for paper in item[1]
            if isinstance(paper, list) and len(paper) >= 2 and str(paper[0]).strip()
        )
        if gene not in first_scores:
            order.append(gene)
            first_scores[gene] = sum(paper["score"] for paper in papers) / len(papers) if papers else 0.0
        latest[gene] = papers
    return tuple(
        {"gene": gene, "score": first_scores[gene], "context": latest[gene]} for gene in order[:limit]
    )


def _private_file(status):
    return os.name != "posix" or (status.st_uid == os.getuid() and status.st_mode & 0o077 == 0)


def _read_cache(path, binding):
    try:
        private_directory(path.parent.parent, create=False)
        private_directory(path.parent, create=False)
        raw, status, _after = read_regular_bounded(path, _MAX_CACHE_BYTES, follow_symlinks=False)
        _check(_private_file(status), OSError())
        record = json.loads(raw.decode("utf-8"))
        _check(
            isinstance(record, dict)
            and set(record) == {"binding", "response"}
            and record["binding"] == binding,
            ValueError(),
        )
        _check(_amelie_entries(record["response"], max(len(record["response"]), 1)), ValueError())
        return record["response"]
    except (
        OSError,
        OverflowError,
        UnicodeError,
        json.JSONDecodeError,
        TypeError,
        ValueError,
        domain.MendelPrimeError,
    ):
        return None


def _write_cache(path, binding, payload):
    private_directory(path.parent.parent)
    private_directory(path.parent)
    raw = _canonical({"binding": binding, "response": payload}).encode("utf-8")
    atomic_replace_bytes(path, raw, mode=0o600)


def _secret(root):
    private_directory(root)
    path = root / ".pseudonym_key"
    for _attempt in range(100):
        try:
            with open_regular(path, follow_symlinks=False) as (stream, status):
                value = stream.read(65)
            _check(_private_file(status) and 32 <= len(value) <= 64, OSError())
            return value
        except FileNotFoundError:
            value = secrets.token_bytes(32)
            try:
                descriptor = os.open(
                    path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0), 0o600
                )
            except FileExistsError:
                time.sleep(0.01)
                continue
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(value)
                stream.flush()
                os.fsync(stream.fileno())
            os.chmod(path, 0o600)
            return value
    raise OSError


def rank_amelie(profile, seeds, cache_root, limit=100, *, external, secret):
    hpo_ids, maximum = _positive(profile), _limit(limit)
    _check(
        type(seeds) is tuple and all(isinstance(value, str) for value in seeds),
        TypeError("seeds must be a tuple of strings"),
    )
    genes = tuple(dict.fromkeys(value.strip().upper() for value in seeds if value.strip()))[:maximum]
    _check(
        genes and all(_GENE.fullmatch(gene) for gene in genes),
        ValueError("AMELIE requires canonical seed genes"),
    )
    acknowledged, verify_tls = (
        external.get("amelie_acknowledged", False),
        external.get("amelie_verify_tls", True),
    )
    _check(
        type(acknowledged) is bool and type(verify_tls) is bool,
        TypeError("AMELIE consent and TLS settings must be bool"),
    )
    binding = {
        "endpoint": _AMELIE_URL,
        "genes": list(genes),
        "hpo_ids": sorted(set(hpo_ids)),
        "schema": "amelie-response-v1",
        "verify_tls": verify_tls,
    }
    root = Path(cache_root)
    path = root / "by_request" / f"{hashlib.sha256(_canonical(binding).encode('ascii')).hexdigest()}.json"
    payload = _read_cache(path, binding)
    entries = _amelie_entries(payload, maximum) if payload is not None else None
    if payload is None:
        _check(acknowledged, _resource_error("external_send_not_acknowledged"))
        ca = _ca_bundle(external)
        if not verify_tls:
            logging.getLogger(__name__).warning(
                "AMELIE TLS certificate and hostname verification are disabled for %s", _AMELIE_URL
            )
        failure = None
        try:
            private = secret.encode("utf-8") if isinstance(secret, str) else secret
            private = private or _secret(root)
        except (OSError, TypeError, ValueError):
            failure = _resource_error("pseudonym_key_unavailable")
        if failure:
            raise failure
        patient = (
            "MG-"
            + hmac.new(
                private,
                _PSEUDONYM_DOMAIN + _canonical(list(profile["case_key"])).encode("ascii"),
                hashlib.sha256,
            ).hexdigest()[:24]
        )
        form = {"patientName": patient, "phenotypes": ",".join(hpo_ids), "genes": ",".join(genes)}
        request = urlrequest.Request(
            _AMELIE_URL,
            data=urlparse.urlencode(form).encode("ascii"),
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            method="POST",
        )
        payload = _send_json(request, 120.0, 1, ca, verify_tls)
        entries = _amelie_entries(payload, maximum)
        _check(entries, _resource_error("empty_tool_ranking"))
        try:
            _write_cache(path, binding, payload)
        except OSError:
            pass
    return _external_record("amelie", profile, entries, _AMELIE_URL, hpo_ids, genes, payload)


def _valid_record(record, expected, profile):
    if not hasattr(record, "get"):
        return "tool_identity_mismatch"
    if record.get("tool_id") != expected or not record.get("entries"):
        return "tool_identity_mismatch" if record.get("tool_id") != expected else "empty_tool_ranking"
    genes = tuple(item["gene"] for item in record["entries"])
    try:
        domain.Ranking(genes=genes)
    except domain.MendelPrimeError:
        return "tool_identity_mismatch"
    sources = record.get("sources", ())
    if not sources or any(
        not isinstance(source, domain.Source)
        or source.case_key != profile["case_key"]
        or source.provider != expected
        or source.kind != "tool_observation"
        or set(source.genes) != set(genes)
        for source in sources
    ):
        return "tool_provenance_mismatch"
    return ""


def collect_candidates(state, profile, cache_root, external, secret):
    rankings, failures = [], []
    functions = (
        ("phen2disease", lambda: rank_phen2disease(state, profile, 150)),
        ("phen2gene", lambda: rank_phen2gene(profile, 150, external=external)),
        ("local_overlap", lambda: rank_overlap(state, profile, 150)),
    )

    def invoke(name, call):
        try:
            record = call()
            code = _valid_record(record, name, profile)
            _check(not code, _resource_error(code or "tool_execution_failed"))
            rankings.append(record)
        except domain.MendelPrimeError as error:
            code = error.code if re.fullmatch(r"[a-z][a-z0-9_]*", error.code) else "tool_execution_failed"
            failures.append(domain.freeze_record({"tool_id": name, "code": code}))
        except Exception:
            failures.append(domain.freeze_record({"tool_id": name, "code": "tool_execution_failed"}))

    for name, call in functions:
        invoke(name, call)
    genes = {record["tool_id"]: tuple(item["gene"] for item in record["entries"]) for record in rankings}
    seeds = domain.source_diverse_rrf(genes, 100, 60, 10, ("phen2disease", "phen2gene", "local_overlap"))
    invoke(
        "amelie",
        lambda: rank_amelie(
            profile, seeds, Path(cache_root) / "amelie", 100, external=external, secret=secret
        ),
    )
    genes = {record["tool_id"]: tuple(item["gene"] for item in record["entries"]) for record in rankings}
    consensus = domain.source_diverse_rrf(
        genes, 100, 30, 10, ("phen2disease", "phen2gene", "local_overlap", "amelie")
    )
    return domain.freeze_record(
        {
            "rankings": tuple(rankings),
            "amelie_seeds": seeds,
            "consensus": consensus,
            "failures": tuple(failures),
        }
    )
