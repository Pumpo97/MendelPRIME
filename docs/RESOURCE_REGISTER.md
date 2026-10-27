# Resource register

This catalog covers all 46 resource path/command fields; these are not 46 independent downloads.

A = prepared project asset; B = upstream ready-made file; C = raw-source preparation needed; D = tool installation/configuration.
MANE, the AI-CURA protocol/workbook and ClinGen SVI are report-only resources.
VCF ranking needs only the AI-CURA derived knowledge and pinned manifest, supplied
in `mendelprime-vcf-ranking-guidance-v1.tar.gz`. Other molecular ranking resources
remain required. The JSON catalog's `required_when` records these distinctions.
Install the required files on the destination host. Upstream raw downloads may need conversion to the runtime format.
Sizes exclude sidecars unless stated; nested directories/aliases overlap and must not be summed. Detailed bytes, indexes and closure files are in [RESOURCE_REGISTER.json](../resources/RESOURCE_REGISTER.json).

| Config field | Type | Target beneath resource_root | Current payload | Preparation path | Declaration |
|---|---|---|---:|---|---|
| `resources.hpo_obo` | A | `core/hpo/hp.obo` | 10.06 MB | [P00](RESOURCE_PREPARATION.md) | HPO project resource attachment |
| `resources.hpo_stats` | A | `core/hpo/hpo_stats.json` | 5.88 MB | [P00](RESOURCE_PREPARATION.md) | Preserve all existing values |
| `resources.genes_to_phenotype` | A | `core/hpo/genes_to_phenotype.txt` | 20.53 MB | [P00](RESOURCE_PREPARATION.md) | Same frozen HPO snapshot |
| `resources.phen2disease_root` | B | `phenotype/phen2disease` | 3.86 GB | [P01](RESOURCE_PREPARATION.md) | Four upstream precomputed files; extraction/layout, not new inference |
| `resources.omim_jsonl` | C | `core/omim/omim_clean_full.jsonl` | 528.16 MB | [P02](RESOURCE_PREPARATION.md) | Current mirror is usable; raw pages are not this schema |
| `resources.mim2gene` | B | `core/omim/mim2gene.txt` | 983.16 KB | [P02](RESOURCE_PREPARATION.md) | Version-matched official table; no model inference |
| `resources.similar_case_root` | A | `phenotype/similar-cases/index` | 92.41 MB | [P00](RESOURCE_PREPARATION.md) | Records, IDs and vectors are inseparable |
| `resources.embedding_model` | B | `phenotype/models/qwen3-embedding-4b` | 8.06 GB | [P03](RESOURCE_PREPARATION.md) | Pinned weight shards plus configuration/tokenizer |
| `resources.cross_encoder` | A | `phenotype/models/medcpt-cross-encoder` | 438.89 MB | [P03](RESOURCE_PREPARATION.md) | Published support asset is converted safetensors |
| `genotype.java` | D | `tools/java/bin/java` | 332.90 MB | [P04](RESOURCE_PREPARATION.md) | Measured size is the complete installed Java root |
| `genotype.exomiser_jar` | D | `genotype/exomiser/exomiser-cli-14.0.0.jar` | 45.17 MB | [P04](RESOURCE_PREPARATION.md) | Measured size includes launcher and adjacent libraries |
| `genotype.exomiser_properties` | D | `genotype/exomiser/application-2406.properties` | 4.80 KB | [P04](RESOURCE_PREPARATION.md) | Configure paths; do not change analysis policy |
| `genotype.lirical_jar` | D | `genotype/lirical/lirical-cli-2.4.1.jar` | 30.93 MB | [P04](RESOURCE_PREPARATION.md) | Built ZIP is available; data is separate |
| `genotype.lirical_data` | B | `genotype/lirical/data` | 363.07 MB | [P04](RESOURCE_PREPARATION.md) | Complete download directory with transcript serializations |
| `genotype.gene_map` | B | `genotype/lirical/data/hgnc_complete_set.txt` | 16.90 MB | [P04](RESOURCE_PREPARATION.md) | Already part of lirical_data |
| `genotype.assemblies[0].exomiser_genome_data` | B | `genotype/exomiser/data/2406_hg19` | 28.39 GB | [P04](RESOURCE_PREPARATION.md) | Official 2406 hg19 data; shared by LIRICAL |
| `genotype.assemblies[0].exomiser_phenotype_data` | B | `genotype/exomiser/data/2406_phenotype` | 12.79 GB | [P04](RESOURCE_PREPARATION.md) | Not interchangeable with the project HPO package |
| `genotype.assemblies[0].lirical_exomiser_data` | B | `genotype/exomiser/data/2406_hg19` | 28.39 GB | [P04](RESOURCE_PREPARATION.md) | Alias of Exomiser genome data, not another download |
| `genotype.assemblies[0].lirical_variants_db` | B | `genotype/exomiser/data/2406_hg19/2406_hg19_variants.mv.db` | 27.69 GB | [P04](RESOURCE_PREPARATION.md) | Contained in Exomiser genome data |
| `genotype.assemblies[0].lirical_clinvar_db` | B | `genotype/exomiser/data/2406_hg19/2406_hg19_clinvar.mv.db` | 113.10 MB | [P04](RESOURCE_PREPARATION.md) | Native mv.db, not project SQLite |
| `variant_knowledge.clinvar_db` | A | `variant-knowledge/clinvar/clinvar_20260117.sqlite` | 4.33 GB | [P05](RESOURCE_PREPARATION.md) | Already processed; raw TSV/VCF needs conversion |
| `variant_knowledge.clingen_validity` | B | `variant-knowledge/clingen/gene_disease_validity.csv` | 1.12 MB | [P06](RESOURCE_PREPARATION.md) | Reader accepts official preamble; manifest still required |
| `variant_knowledge.clingen_vcep` | B | `variant-knowledge/clingen/variant_pathogenicity.tsv` | 31.54 MB | [P06](RESOURCE_PREPARATION.md) | Current export is TSV; validate actual delimiter and headers |
| `variant_knowledge.clingen_manifest` | C | `variant-knowledge/clingen/resource_manifest.json` | 725 B | [P06](RESOURCE_PREPARATION.md) | Generate hashes binding both exports and release date |
| `explanation_resources.reference_fasta` | C | `assessment/grch37/reference/hg19.fa` | 3.20 GB | [P07](RESOURCE_PREPARATION.md) | Extract matching FASTA and create matching fai; no coordinate conversion |
| `explanation_resources.vep_command` | D | `tools/vep/vep` | runtime-dependent | [P08](RESOURCE_PREPARATION.md) | Executable plus Perl/container/runtime/mount closure |
| `explanation_resources.vep_cache` | B | `assessment/grch37/vep/cache` | 15.62 GB | [P08](RESOURCE_PREPARATION.md) | Indexed 104_GRCh37 cache; no new variant scoring at preparation |
| `explanation_resources.vep_fasta` | C | `assessment/grch37/vep/final_hg19.fa` | 3.15 GB | [P07](RESOURCE_PREPARATION.md) | Separate existing file; do not assume byte identity to reference_fasta |
| `explanation_resources.clinvar_protein_disease_db` | A | `assessment/grch37/clinvar-protein/clinvar_protein_acmg_disease_20260117.sqlite` | 97.93 MB | [P05](RESOURCE_PREPARATION.md) | Preserve disease fields and comparator selection |
| `explanation_resources.gnomad_sites` | B | `assessment/grch37/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz` | 63.15 GB | [P09](RESOURCE_PREPARATION.md) | Native v2.1.1 VCF plus matching index |
| `explanation_resources.gnomad_coverage` | C | `assessment/grch37/gnomad/gnomad.exomes.coverage.summary.tsv.bgz` | 1.69 GB | [P09](RESOURCE_PREPARATION.md) | Native v2.1 table plus coordinate index; no missing-value imputation |
| `explanation_resources.gnomad_constraint` | B | `assessment/grch37/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` | 4.61 MB | [P09](RESOURCE_PREPARATION.md) | Official v2.1.1 gene table; parser reads compressed text |
| `explanation_resources.spliceai_snv` | B | `assessment/grch37/predictions/spliceai_scores.masked.snv.hg19.vcf.gz` | 28.63 GB | [P10](RESOURCE_PREPARATION.md) | Masked v1.3 hg19 precomputed VCF plus index |
| `explanation_resources.spliceai_indel` | B | `assessment/grch37/predictions/spliceai_scores.masked.indel.hg19.vcf.gz` | 68.86 GB | [P10](RESOURCE_PREPARATION.md) | Same score model/version as SNV file |
| `explanation_resources.revel` | C | `assessment/grch37/predictions/new_tabbed_revel.tsv.gz` | 652.22 MB | [P10](RESOURCE_PREPARATION.md) | Raw CSV -> exact nine-column TSV -> BGZF -> tabix |
| `explanation_resources.mane_summary` | B | `assessment/transcripts/MANE.GRCh38.v1.5.summary.txt.gz` | 1.12 MB | [P11](RESOURCE_PREPARATION.md) | GRCh38 columns are transcript context, not hg19 coordinates |
| `explanation_resources.autopvs1_root` | D | `assessment/grch37/autopvs1` | 3.29 GB | [P11](RESOURCE_PREPARATION.md) | Pinned source plus disclosed local FASTA compatibility patch |
| `explanation_resources.ai_cura_protocol` | C | `guidelines/ai-cura/AI-CURA_Protocol.md` | 16.27 KB | [P12](RESOURCE_PREPARATION.md) | Derived protocol text, not a raw PDF |
| `explanation_resources.ai_cura_workbook` | B | `guidelines/ai-cura/adz4172_tables_s1_to_s13.xlsx` | 127.29 KB | [P12](RESOURCE_PREPARATION.md) | Preserve the official XLSX bytes |
| `explanation_resources.ai_cura_knowledgebase` | C | `guidelines/ai-cura/ai_cura_knowledgebase.json` | 75.86 KB | [P12](RESOURCE_PREPARATION.md) | Workbook tables -> schema-bound JSON |
| `explanation_resources.ai_cura_manifest` | C | `guidelines/ai-cura/resource_manifest.json` | 2.09 KB | [P12](RESOURCE_PREPARATION.md) | Exact public and original pins supported; do not arbitrarily rewrite |
| `explanation_resources.clingen_guideline_manifest` | C | `guidelines/clingen-svi-20261021/manifest.json` | 1.40 KB (closure 359.57 KB) | [P12](RESOURCE_PREPARATION.md) | Original guidance -> text extraction -> relative manifest |
| `explanation_resources.clingen_dosage_table` | B | `report/clingen-dosage/clingen_dosage_detailed_GRCh37_20260924.tsv` | 248.80 KB | [P13](RESOURCE_PREPARATION.md) | Optional GRCh37 report enrichment |
| `explanation_resources.genereviews_index` | C | `report/genereviews/genereviews_index.tsv` | 104 B | [P13](RESOURCE_PREPARATION.md) | Current configured subset contains three genes; chapters are fetched online |
| `explanation_resources.post_evidence_catalog` | C | `report/curated-literature/lmx1b-family-series-20261024/catalog.json` | 25.12 KB (closure 992.56 KB) | [P13](RESOURCE_PREPARATION.md) | Manual allele/quote/family verification; all referenced files required |
| `explanation_resources.kinship_executable` | D | `tools/plink2/plink2` | 42.49 MB | [P13](RESOURCE_PREPARATION.md) | Optional kinship QC executable, not variant evidence by itself |

## Scope declarations

- HPO and processed RareA vectors are designated project data attachments; actual upload is a separate step.
- OMIM remains the local mirror interface. Its production snapshot is user-supplied; the source contains a one-record format demo.
- ClinVar SQLite is already prepared and can be installed directly from the provided asset. Raw upstream ClinVar files are not drop-in replacements.
- Online LLM, phenotype and evidence services are listed separately in RESOURCE_PREPARATION.md; they have no fixed local payload size.
- VEP command size excludes its image/Perl/runtime closure. AutoPVS1 size includes its currently configured reference dependency and overlaps the shared reference.
- The GeneReviews index is a three-gene subset, and the curated literature catalog is a scoped evidence collection, not comprehensive databases.
- Keep resource versions, scientific values and source bindings consistent with their manifests.

## Shipped preparation scripts

Conversion, indexing, pinned-resource preparation and compatibility installation are documented in [preparation/README.md](../preparation/README.md). The JSON register includes the available commands. Prepared assets can be reused directly. These utilities do not run an LLM or rank patients; manually reviewed evidence still requires verified source bindings.

## HPO statistics

`hpo_stats.json` stores annotation frequencies and information content (IC).
It is separate from `hp.obo`, which defines terms and ontology relationships.
The snapshot has 19,077 entries in each IC map, 5,251 genes, 9,120 diseases and
329,279 annotation records. These are resource counts, not patient prevalences.

| Field | Meaning |
|---|---|
| `hpo_gene_counts`, `hpo_disease_counts` | Annotated gene/disease counts including ancestor propagation |
| `direct_hpo_gene_counts`, `direct_hpo_disease_counts` | Counts before propagation |
| `ic_gene`, `ic_disease` | Gene/disease annotation information content consumed by the reader |

```text
IC_gene(h)    = -ln((gene_count(h) + 1) / (total_genes + 1))
IC_disease(h) = -ln((disease_count(h) + 1) / (total_diseases + 1))
```

The +1 terms smooth zero counts. Broad terms have lower IC; annotation coverage
also affects specificity. For example, HP:0001250 has 2,111 associated genes and
`ic_gene=0.9109735905965348`; HP:0000118 has 5,242 and `ic_gene=0.001715102850154058`.

`tools/phenotype.py::build_profile` combines gene IC and ontology depth to mark
up to six specificity-prioritized hints. The tools still receive the complete
positive HPO input. RareA query construction uses these hints, while its ontology
matching uses disease IC and shared ancestors. Phen2Disease calculates IC from
its own dataset; HPOOverlap uses overlap counts. The stored `ic_case`,
`hpo_case_counts` and `case_statistics` fields are not consumed by this reader.

The HPO package installs `hp.obo`, `genes_to_phenotype.txt` and `hpo_stats.json`
under `core/hpo/`. Its ontology is `hp/releases/2025-03-03`; the package date
identifies the project snapshot, not the upstream ontology release. Of 329,339
annotation rows, 60 reference 39 terms outside this ontology; statistics use the
329,279 known-term rows. Preserve the supplied values and version match.
Configure `hpo_stats` even with RareA disabled because profile construction uses
gene IC. Source basenames are provenance labels, not destination filesystem paths.
