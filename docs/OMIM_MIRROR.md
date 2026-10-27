# Local OMIM mirror format

MendelPRIME uses a **locally mirrored, normalized OMIM JSONL snapshot** for fast
lookup. This is the application's storage format, not the official OMIM API
response schema. Mirroring and preprocessing preserve entry identity, section
text and linked PMIDs while avoiding repeated OMIM web/API requests at runtime.
The snapshot is loaded once per KnowledgeBase instance and reused across its
cases. Other evidence services, such as PubMed, can still require network access.

The supported deployment route is local files. No OMIM API key or API adapter is required.
Users supply the locally prepared snapshot appropriate to their authorized use.

## One-record demo

The source includes:

```text
examples/resources/omim/omim_demo.jsonl
examples/resources/omim/mim2gene_demo.txt
```

The JSONL contains exactly one record for ASPM (MIM 605481). It uses one short
sentence excerpt from the existing local mirror, with the PMID explicitly
present in that sentence. Its companion mapping contains the matching gene/MIM
row. The excerpt is not the complete entry, and the empty allelic_variants list
does not assert that the actual OMIM entry has no allelic variants. The added
demo_metadata field labels this restricted scope and is not used for ranking.

Do not use this small demo as the OMIM knowledge base for a cohort or interpret
missing genes as a biological result. It is a runnable data-format example,
not a distribution of the full mirrored database or an accuracy benchmark.

## Reader contract

Each nonempty line is one JSON object. Keep text and citation bindings together:

| Field | Meaning and current reader use |
|---|---|
| omim_id | String entry identifier; used to join the gene mapping |
| page_title | Human-readable title |
| sections.molecular_genetics | Ordered paragraph objects containing original text and associated PMIDs |
| sections.* | Other mirrored sections may be retained; current molecular retrieval does not automatically consume every section |
| allelic_variants | Optional variant records; the citation-seed reader consumes sub_parts[].pmids |
| scraped_at / other provenance | May retain the actual snapshot acquisition information; do not invent a date when it is unknown |

The mim2gene table is tab-delimited. Its gene rows supply the MIM identifier in
column 1, entry type `gene` in column 2 and approved gene symbol in column 4.
Preserve the matching map when preparing or transferring a snapshot. Never
attach all references from an entry to each paragraph: citations must remain
bound to the text from which they were obtained.

## Configure a complete local mirror

In the selected HPO/VCF resources.toml, keep:

```toml
[resources]
omim_jsonl = "core/omim/omim_clean_full.jsonl"
mim2gene = "core/omim/mim2gene.txt"
```

These paths resolve beneath resource_root. HPO plus clinical description reuses
the HPO resource manifest. Supply the full snapshot separately and retain its
provenance/version; changing the database content can change retrieved evidence.
The application records content-derived source versions for traceability.
The small demonstration does not imply permission to redistribute a full OMIM
mirror; source-data terms remain separate from the application source release.
