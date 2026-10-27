"""Deterministic semantic, dense, and hybrid similar-case retrieval."""

from __future__ import annotations

import hashlib
import json
import re
import threading
from collections import defaultdict
from importlib import import_module
from pathlib import Path

from mendelprime.domain import Source, validate_gene
from mendelprime.knowledge.hpo import ancestors

_HPO_ID = re.compile(r"HP:\d{7}", re.ASCII)
_INSTRUCTION = (
    "Instruct: Given a patient's phenotypes and clinical description, "
    "retrieve patient cases likely caused by the same Mendelian disease or causal gene.\nQuery:"
)


def load_index(root):
    """Load the RareA corpus and its aligned dense index."""
    numpy = import_module("numpy")
    root = Path(root)
    records = _load_source(root)
    directory = root / "rarea"
    try:
        identifiers = json.loads((directory / "record_ids.json").read_text(encoding="utf-8"))
        matrix = numpy.load(directory / "embeddings.npy", allow_pickle=False, mmap_mode="r")
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("invalid similar-case index for rarea") from error
    if (
        not isinstance(identifiers, list)
        or any(not isinstance(value, str) for value in identifiers)
        or identifiers != [record["record_id"] for record in records]
        or getattr(matrix, "ndim", None) != 2
        or matrix.shape[0] != len(records)
    ):
        raise ValueError("similar-case index mismatch for rarea")
    return {
        "root": root,
        "sources": ("rarea",),
        "records": records,
        "records_by_source": {"rarea": records},
        "embeddings": {"rarea": matrix},
        "semantic_hpo": None,
        "semantic_terms": (),
        "lin_cache": {},
        "model_key": None,
        "models": None,
        "lock": threading.RLock(),
        "inference_lock": threading.Lock(),
    }


def _load_source(root):
    rows = []
    digest = hashlib.sha256()
    try:
        with (root / "rarea" / "records.jsonl").open("rb") as stream:
            for raw_line in stream:
                digest.update(raw_line)
                if not raw_line.strip():
                    continue
                value = json.loads(raw_line.decode("utf-8"))
                if not isinstance(value, dict):
                    raise ValueError("record must be an object")
                rows.append(value)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("invalid similar-case corpus for rarea") from error
    if not rows:
        raise ValueError("empty similar-case corpus for rarea")
    version = "sha256:" + digest.hexdigest()
    records = tuple(_case_record(row, "rarea", version) for row in rows)
    identifiers = [record["record_id"] for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate similar-case record for rarea")
    return records


def _case_record(row, source, version):
    record_id = str(row.get("record_id") or "")
    values = row.get("candidate_genes") or (row.get("gene"),)
    values = (values,) if isinstance(values, str) else values
    if not isinstance(values, (list, tuple)):
        raise ValueError("invalid similar-case genes")
    genes = tuple(
        dict.fromkeys(
            str(value).strip().upper() for value in values if value is not None and str(value).strip()
        )
    )
    for gene in genes:
        validate_gene(gene)
    raw_phenotypes = row.get("phenotype_evidence") or ()
    metadata = row.get("metadata") or {}
    if not isinstance(metadata, dict):
        raise ValueError("invalid similar-case metadata")
    if not isinstance(raw_phenotypes, (list, tuple)) or any(
        not isinstance(value, dict) for value in raw_phenotypes
    ):
        raise ValueError("invalid similar-case phenotype evidence")
    phenotypes = tuple(
        {
            "label": str(value.get("hpo_label") or value.get("matched_concept") or ""),
            "assertion": str(value.get("assertion") or "present"),
            "subject": str(value.get("subject") or "patient"),
        }
        for value in raw_phenotypes
    )
    retrieval_text = str(row.get("retrieval_text") or row.get("case_text") or "")
    case_text = str(row.get("case_text") or retrieval_text)
    provider = str(row.get("source") or record_id.split(":", 1)[0] or source).strip().lower()
    group_id = str(row.get("group_id") or record_id)
    hpo_ids = _text_values(row.get("hpo_ids") or ())
    pmid_values = (
        row.get("pmids")
        or row.get("source_pmids")
        or metadata.get("pmids")
        or metadata.get("source_pmids")
        or ()
    )
    pmids = _text_values(pmid_values)
    if not record_id or not group_id or not genes or not case_text or provider != source:
        raise ValueError("invalid similar-case record")
    if any(not _HPO_ID.fullmatch(value) for value in hpo_ids):
        raise ValueError("invalid similar-case HPO identifier")
    if any(not item["assertion"] or not item["subject"] for item in phenotypes):
        raise ValueError("invalid similar-case phenotype evidence")
    return {
        "record_id": record_id,
        "provider": provider,
        "group_id": group_id,
        "candidate_genes": genes,
        "case_text": case_text,
        "retrieval_text": retrieval_text,
        "gene_link_type": str(row.get("gene_link_type") or ""),
        "hpo_ids": hpo_ids,
        "pmids": pmids,
        "phenotypes": phenotypes,
        "corpus_version": version,
        "gene_case_rank": _gene_case_rank(metadata.get("gene_case_rank")),
    }


def _text_values(value):
    if not isinstance(value, (list, tuple)):
        raise ValueError("similar-case text sequence required")
    return tuple(dict.fromkeys(str(item).strip() for item in value if str(item).strip()))


def _gene_case_rank(value):
    try:
        rank = int(value or 1)
    except (TypeError, ValueError) as error:
        raise ValueError("invalid similar-case gene case rank") from error
    if rank < 1:
        raise ValueError("invalid similar-case gene case rank")
    return rank


def semantic_hits(state, hpo, profile, gene_ranked=False, gene_limit=10, record_pool=50):
    """Rank cases with bidirectional IC-weighted Lin similarity."""
    records, case_terms, lin_cache = _semantic_context(state, hpo)
    patient = _prune(hpo, {term["hpo_id"] for term in profile["terms"] if term["present"]})
    information_content = hpo["disease_ic"]
    fallback = max(information_content.values(), default=1.0)

    def lin(left, right):
        key = (left, right) if left <= right else (right, left)
        if key not in lin_cache:
            common = (ancestors(hpo, left) | {left}) & (ancestors(hpo, right) | {right})
            mica = max((information_content.get(term, 0.0) for term in common), default=0.0)
            denominator = information_content.get(left, 0.0) + information_content.get(right, 0.0)
            lin_cache[key] = 2.0 * mica / denominator if denominator else float(left == right)
        return lin_cache[key]

    def score(case):
        if not patient or not case:
            return 0.0

        def weighted(source, target):
            weights = [information_content.get(term, fallback) for term in source]
            total = sum(
                max((lin(term, other) for other in target), default=0.0) * weight
                for term, weight in zip(source, weights, strict=True)
            )
            return total, sum(weights)

        patient_total, patient_weight = weighted(patient, case)
        case_total, case_weight = weighted(case, patient)
        return patient_total / (patient_weight or fallback) + (patient_total + case_total) / (
            patient_weight + case_weight or fallback
        )

    scored = [
        (score(terms), record["record_id"], index)
        for index, (record, terms) in enumerate(zip(records, case_terms, strict=True))
    ]
    grouped = {}
    for item in scored:
        group = records[item[2]]["group_id"]
        if group not in grouped or (-item[0], item[1]) < (-grouped[group][0], grouped[group][1]):
            grouped[group] = item
    ranked = sorted(grouped.values(), key=lambda item: (-round(item[0], 8), item[1]))
    if not gene_ranked:
        return tuple(
            _semantic_hit(
                records[index],
                rank,
                value,
                records[index]["candidate_genes"][0],
                records[index]["gene_case_rank"],
            )
            for rank, (value, _, index) in enumerate(ranked[:10], 1)
        )
    best = {}
    for value, record_id, index in ranked:
        for gene in records[index]["candidate_genes"]:
            if gene not in best or (-value, record_id) < (-best[gene][0], best[gene][1]):
                best[gene] = (value, record_id, index)
    selected = sorted(best, key=lambda gene: (-round(best[gene][0], 8), gene))[:gene_limit]
    output = []
    for gene_rank, gene in enumerate(selected, 1):
        case_rank = 0
        for record_rank, (value, _, index) in enumerate(ranked, 1):
            record = records[index]
            if gene not in record["candidate_genes"] or (
                record_rank > record_pool and index != best[gene][2]
            ):
                continue
            case_rank += 1
            hit = _semantic_hit(record, record_rank, value, gene, case_rank)
            hit.update(
                rerank_rank=gene_rank,
                rerank_score=best[gene][0],
                semantic_gene_rank=gene_rank,
                semantic_record_rank=record_rank,
            )
            output.append(hit)
    return tuple(output)


def _semantic_context(state, hpo):
    with state["lock"]:
        if state["semantic_hpo"] is not hpo:
            state["semantic_hpo"] = hpo
            state["semantic_terms"] = tuple(
                _prune(hpo, set(record["hpo_ids"])) for record in state["records"]
            )
            state["lin_cache"] = {}
        return state["records"], state["semantic_terms"], state["lin_cache"]


def _prune(hpo, terms):
    return tuple(
        sorted(
            term
            for term in terms
            if not any(term != other and term in ancestors(hpo, other) for other in terms)
        )
    )


def _semantic_hit(record, rank, score, gene, gene_case_rank):
    return _hit(
        record,
        gene,
        retrieval_method="phen2disease_bidirectional_ic_lin",
        retrieval_rank=rank,
        retrieval_score=score,
        rerank_rank=rank,
        rerank_score=score,
        gene_case_rank=gene_case_rank,
        semantic_record_rank=rank,
    )


def dense_hits(state, profile, embedding_model, cross_encoder, device="cuda"):
    """Run local-only embedding retrieval, cross reranking, and weighted gene RRF."""
    models = _models(state, embedding_model, cross_encoder, device)
    with state["inference_lock"]:
        query = _embed(models, _query_text(profile), device)
        candidates = []
        for provider in state["sources"]:
            scores = state["embeddings"][provider] @ query
            candidates.extend(
                {
                    "provider": provider,
                    "index": index,
                    "score": float(score),
                    "record": state["records_by_source"][provider][index],
                }
                for index, score in enumerate(scores)
            )
        ranked = sorted(candidates, key=lambda value: (-value["score"], value["provider"], value["index"]))
        selected = []
        groups = set()
        for candidate in ranked:
            if candidate["record"]["group_id"] not in groups:
                groups.add(candidate["record"]["group_id"])
                selected.append(candidate)
            if len(selected) == 50:
                break
        labels = (term["label"] for term in profile["terms"] if term["present"] and term["label"])
        patient_text = "; ".join(dict.fromkeys(labels))
        pairs = [(patient_text, _positive_labels(candidate["record"])) for candidate in selected]
        scores = _cross_scores(models, pairs, device)
    with_ranks = [(candidate, rank) for rank, candidate in enumerate(selected, 1)]
    cross_ranked = sorted(
        zip(with_ranks, scores, strict=True), key=lambda value: (-round(value[1], 5), value[0][1])
    )
    dense_ranks = _gene_ranks(candidate["record"] for candidate in selected)
    cross_ranks = _gene_ranks(candidate["record"] for (candidate, _), _ in cross_ranked)
    records = defaultdict(list)
    seen = defaultdict(set)
    for record_rank, ((candidate, _), cross_score) in enumerate(cross_ranked, 1):
        for gene in candidate["record"]["candidate_genes"]:
            if candidate["record"]["group_id"] not in seen[gene]:
                seen[gene].add(candidate["record"]["group_id"])
                records[gene].append((candidate, record_rank, cross_score))
    fusion = {gene: 1 / (60 + rank) + 0.5 / (60 + cross_ranks[gene]) for gene, rank in dense_ranks.items()}
    genes = sorted(fusion, key=lambda gene: (-fusion[gene], dense_ranks[gene], gene))[:10]
    dense_order = {id(candidate): rank for candidate, rank in with_ranks}
    output = []
    for gene_rank, gene in enumerate(genes, 1):
        for case_rank, (candidate, cross_rank, cross_score) in enumerate(records[gene], 1):
            output.append(
                _hit(
                    candidate["record"],
                    gene,
                    retrieval_method="dense_cross_encoder_rrf",
                    retrieval_rank=dense_order[id(candidate)],
                    retrieval_score=candidate["score"],
                    rerank_rank=gene_rank,
                    rerank_score=fusion[gene],
                    gene_case_rank=case_rank,
                    dense_gene_rank=dense_ranks[gene],
                    cross_encoder_gene_rank=cross_ranks[gene],
                    cross_encoder_record_rank=cross_rank,
                    cross_encoder_score=cross_score,
                )
            )
    return tuple(output)


def _models(state, embedding_model, cross_encoder, device):
    key = (str(embedding_model), str(cross_encoder), device)
    with state["lock"]:
        if state["models"] is None:
            embedding_path = _local_model_directory(embedding_model, "embedding_model")
            cross_path = _local_model_directory(cross_encoder, "cross_encoder")
            torch = import_module("torch")
            transformers = import_module("transformers")
            tokenizer = transformers.AutoTokenizer.from_pretrained(
                embedding_path, padding_side="left", local_files_only=True
            )
            dtype = torch.bfloat16 if device.startswith("cuda") else torch.float32
            encoder = (
                transformers.AutoModel.from_pretrained(embedding_path, dtype=dtype, local_files_only=True)
                .to(device)
                .eval()
            )
            cross_tokenizer = transformers.AutoTokenizer.from_pretrained(cross_path, local_files_only=True)
            cross = (
                transformers.AutoModelForSequenceClassification.from_pretrained(
                    cross_path, local_files_only=True
                )
                .to(device)
                .eval()
            )
            state["models"] = {
                "torch": torch,
                "tokenizer": tokenizer,
                "encoder": encoder,
                "cross_tokenizer": cross_tokenizer,
                "cross": cross,
            }
            state["model_key"] = key
        if state["model_key"] != key:
            raise ValueError("dense retrieval state is already bound to different models")
        return state["models"]


def _local_model_directory(value, field):
    path = Path(value)
    if not path.is_absolute() or not path.is_dir():
        raise ValueError(f"{field} must be an existing local directory")
    return str(path)


def _embed(models, text, device):
    torch = models["torch"]
    encoded = models["tokenizer"](
        [_INSTRUCTION + text], padding=True, truncation=True, max_length=1024, return_tensors="pt"
    ).to(device)
    with torch.inference_mode():
        hidden = models["encoder"](**encoded).last_hidden_state
        if bool((encoded["attention_mask"][:, -1].sum() == 1).item()):
            value = hidden[:, -1]
        else:
            value = hidden[
                torch.arange(hidden.shape[0], device=hidden.device), encoded["attention_mask"].sum(dim=1) - 1
            ]
        return torch.nn.functional.normalize(value, p=2, dim=1)[0].cpu().float().numpy()


def _cross_scores(models, pairs, device):
    if not pairs:
        return []
    torch = models["torch"]
    encoded = models["cross_tokenizer"](
        pairs, padding=True, truncation=True, max_length=128, return_tensors="pt"
    ).to(device)
    with torch.inference_mode():
        values = models["cross"](**encoded).logits.reshape(-1).cpu().float().tolist()
        return [float(value) for value in values]


def hybrid_hits(dense_hits, semantic_gene_hits):
    """Fuse dense and semantic gene ranks with equal-route reciprocal rank fusion."""
    routes = {"deeprare_global": dense_hits, "phen2disease": semantic_gene_hits}
    ranks = {route: _hit_gene_ranks(hits) for route, hits in routes.items()}
    scores = defaultdict(float)
    source_count = defaultdict(int)
    best_rank = {}
    for values in ranks.values():
        for gene, rank in values.items():
            scores[gene] += 1 / (60 + rank)
            source_count[gene] += 1
            best_rank[gene] = min(best_rank.get(gene, 10**9), rank)
    genes = tuple(
        sorted(scores, key=lambda gene: (-scores[gene], -source_count[gene], best_rank[gene], gene))[:10]
    )
    selected = set(genes)
    merged = {}
    for route, hits in routes.items():
        for hit in hits:
            route_rank = ranks[route].get(hit["gene"])
            if hit["gene"] not in selected or route_rank is None:
                continue
            key = hit["gene"], hit["group_id"]
            if key not in merged:
                merged[key] = {
                    "primary": hit,
                    "primary_route_rank": route_rank,
                    "routes": set(),
                    "ranks": {},
                    "records": defaultdict(list),
                }
            value = merged[key]
            if _hit_priority(hit, route_rank) < _hit_priority(value["primary"], value["primary_route_rank"]):
                value["primary"] = hit
                value["primary_route_rank"] = route_rank
            value["routes"].add(route)
            value["ranks"][route] = route_rank
            value["records"][route].append(hit["record_id"])
    output = []
    for gene_rank, gene in enumerate(genes, 1):
        values = (value for (item_gene, _), value in merged.items() if item_gene == gene)
        values = sorted(
            values,
            key=lambda value: (
                -len(value["routes"]),
                min(value["ranks"].values()),
                value["primary"]["gene_case_rank"],
                value["primary"]["retrieval_rank"],
                value["primary"]["source"].provider,
                value["primary"]["record_id"],
            ),
        )
        for case_rank, value in enumerate(values, 1):
            hit = dict(value["primary"])
            hit.update(
                retrieval_method="deeprare_phen2disease_gene_rrf",
                rerank_rank=gene_rank,
                rerank_score=scores[gene],
                gene_case_rank=case_rank,
                hybrid_routes=tuple(sorted(value["routes"])),
                hybrid_route_gene_ranks=tuple(sorted(value["ranks"].items())),
                hybrid_route_record_ids=tuple(
                    (route, tuple(dict.fromkeys(identifiers)))
                    for route, identifiers in sorted(value["records"].items())
                ),
            )
            output.append(hit)
    return tuple(output)


def _hit(record, gene, **values):
    digest = hashlib.sha256(record["record_id"].encode("utf-8")).hexdigest()
    source = Source(
        provider=record["provider"],
        kind="similar_case",
        locator=f"resource:similar-cases/{record['provider']}/{digest}",
        version=record["corpus_version"],
        content=record["case_text"],
        genes=(gene,),
    )
    hit = {
        "record_id": record["record_id"],
        "group_id": record["group_id"],
        "gene": gene,
        "gene_link_type": record["gene_link_type"],
        "hpo_ids": record["hpo_ids"],
        "pmids": record["pmids"],
        "source": source,
        "dense_gene_rank": None,
        "cross_encoder_gene_rank": None,
        "cross_encoder_record_rank": None,
        "cross_encoder_score": None,
        "semantic_gene_rank": None,
        "semantic_record_rank": None,
        "hybrid_routes": (),
        "hybrid_route_gene_ranks": (),
        "hybrid_route_record_ids": (),
    }
    hit.update(values)
    return hit


def _gene_ranks(records):
    output = {}
    for record in records:
        for gene in record["candidate_genes"]:
            output.setdefault(gene, len(output) + 1)
    return output


def _hit_gene_ranks(hits):
    ranks = {}
    for hit in hits:
        if hit["rerank_rank"] <= 10:
            ranks[hit["gene"]] = min(ranks.get(hit["gene"], 10**9), hit["rerank_rank"])
    return ranks


def _hit_priority(hit, route_rank):
    return (
        route_rank,
        hit["gene_case_rank"],
        hit["retrieval_rank"],
        hit["source"].provider,
        hit["record_id"],
    )


def _query_text(profile):
    present = tuple(term for term in profile["terms"] if term["present"])
    ordered = sorted(present, key=lambda term: (-term["information_content"], term["hpo_id"]))
    specific = {term["hpo_id"] for term in ordered[:8]}
    terms = [
        ("Specific phenotype: " if term["hpo_id"] in specific else "")
        + f"{term['hpo_id']} {term['label']}"
        + (f": {term['definition']}" if term["definition"] else "")
        for term in present
    ]
    sections = ["Patient phenotypes: " + "; ".join(terms)]
    absent = tuple(term for term in profile["terms"] if not term["present"])
    if absent:
        labels = (f"{term['hpo_id']} {term['label']}" for term in absent)
        sections.append("Explicitly absent phenotypes: " + "; ".join(labels))
    if profile["clinical_description"].strip():
        sections.append("Clinical description: " + profile["clinical_description"].strip())
    return "\n".join(sections)


def _positive_labels(record):
    labels = [
        value["label"]
        for value in record["phenotypes"]
        if value["assertion"].lower() == "present" and value["subject"].lower() == "patient"
    ]
    if not any(labels):
        text = re.split(
            r"\b(?:Absent|Explicitly absent) phenotypes\s*:",
            record["retrieval_text"],
            maxsplit=1,
            flags=re.IGNORECASE,
        )[0]
        text = re.sub(r"^\s*(?:Patient|Present) phenotypes\s*:\s*", "", text, flags=re.IGNORECASE)
        without_ids = re.sub(r"\s*\[HP:\d+\]|\bHP:\d+\b", "", text, flags=re.IGNORECASE)
        labels = re.split(r"[;,\n]+", without_ids)
    normalized = (" ".join(str(value or "").split()) for value in labels)
    return "; ".join(dict.fromkeys(value for value in normalized if value))


def unique_upper_gene_groups(hits):
    groups = defaultdict(set)
    for hit in hits:
        gene = hit["gene"].upper()
        if hit["group_id"] not in groups[gene]:
            groups[gene].add(hit["group_id"])
            yield hit


def unique_case_hits(hits):
    records = set()
    groups = defaultdict(set)
    for hit in hits:
        record = hit["record_id"], hit["gene"]
        if record not in records and hit["group_id"] not in groups[hit["gene"]]:
            records.add(record)
            groups[hit["gene"]].add(hit["group_id"])
            yield hit


__all__ = [
    "dense_hits",
    "hybrid_hits",
    "load_index",
    "semantic_hits",
    "unique_case_hits",
    "unique_upper_gene_groups",
]
