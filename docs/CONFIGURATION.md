# Configuration by input and output

This guide follows the actual configuration loader, Runner, Medical Tools and
assessment-resource readers. It describes the minimum profiles and their dependencies.

The [installer](INSTALL.md) generates a private `run.toml`/`resources.toml` pair
under WORK/configs/phenotype-rank, phenotype-report, vcf-rank or vcf-report.
The `minimal/*.toml` paths below identify source templates, not the generated
paths used by the README. Existing manual configurations remain supported.

## 1. Choose the branch and output

| Input | Ranking profile | Report | Local resources |
|---|---|---|---|
| HPO | `minimal/phenotype.toml` | Add `--explanation` | Common HPO/OMIM plus gene annotations and Phen2Disease |
| HPO + clinical text | Same profile | Add `--explanation` | Exactly the same as HPO |
| HPO + VCF, trio or singleton | `minimal/vcf_rank.toml` | Select `minimal/vcf_report.toml` instead | Common HPO/OMIM, native tools and molecular resources; reports add VEP, MANE and full guidance |

A Case and a Cohort use the same profile and resource installation. The input
record selects the scientific branch; there is no `--mode hpo`, `--trio` or
`--singleton` CLI flag. A VCF cohort may mix trio and proband-only cases. For a
cohort mixing VCF and non-VCF records, configure the union of both resource sets;
none of the minimal profiles is a shortcut for that mixed-resource deployment.

`ranking.similar_cases = false` disables RareA only in the phenotype branch.
`ranking.reflection = false` disables optional Reflection. Both are explicit in
all minimum profiles. The regular `configs/hpo`, `clinical` and `vcf` templates
describe expanded installations; they are not the minimum-deployment examples.

### What the report switch actually changes

| Work | HPO/clinical, report off | HPO/clinical, report on | VCF, report off | VCF, report on |
|---|---|---|---|---|
| Ranking and its literature evidence | Yes | Yes | Yes | Yes |
| Pre-ranking molecular preparation and ACMG evidence | No | No | Yes | Yes |
| Local VEP with `ranking_vep_enabled=false` | No | No | No | Yes, post-ranking |
| Post-ranking ACMG for selected Top5-gene variants | No | No | No | Yes |
| Report selection and writing | No | Yes | No | Yes |
| Formal files per case | Two | Two | Two | Two |

The report switch is not a request to rerank genes. With reporting off, the
Markdown is empty and provenance is retained for a later `explain` invocation.
VCF report annotation does not enable ranking VEP: that separate switch remains
false in both minimum profiles. Conversely, explicitly enabling ranking VEP
requires the VEP resources even without a report.

## 2. Where each setting belongs

```text
Case JSON/YAML or Cohort JSONL -> patient HPO/text, VCF, build and family roles
minimal/*.toml                -> endpoint, execution and scientific switches
minimal/*.resources.toml      -> installed paths and actual resource versions
environment                  -> API key and explicit TLS compatibility policy
CLI                          -> input, config, output and occasional overrides
```

`resource_file` is relative to its run TOML. `resource_root` is relative to its
resource TOML unless absolute. Resource entries resolve under that root; absolute
entries remain absolute. The first `vep_command` item is similarly resolved;
other arguments are literal except `{workdir}`. TOML never expands `$RES` or `~`.

The split files are the recommended interface. A legacy single TOML remains
readable: omit `resource_file`, put supported sections together, and use absolute
resource paths or paths relative to that TOML. Do not copy `resource_root` into
the legacy single-file form. Unknown settings are rejected.

## 3. Run parameters and their scope

Values below are the **minimum templates' effective values**, not a promise that
omitting every field uses those values. Keep the provided profile intact except
for the documented choices.

| Parameter in run TOML | Minimum value | Applies to / meaning |
|---|---|---|
| `resource_file` | Matching `.resources.toml` | One selected resource manifest |
| `llm.base_url` | Replace the `.invalid` placeholder | Your provider's base URL, not the final `/chat/completions` URL |
| `llm.model` | `<MODEL_ID>` in source templates | Supply a real ID; the installer requires it and writes it explicitly |
| `llm.wire_api` | `chat_completions` | Supported alternatives: `chat_completions`, `responses`; match provider |
| `llm.transport_retries` | `1` | Transport retry setting, subject to stage-specific attempt budgets; not structural-output repair |
| `llm.api_key_env` | Implicit `MENDELPRIME_API_KEY` | Name of the environment variable, never the key value |
| `execution.lanes` | `1` | Concurrent cases within one process |
| `execution.max_llm_concurrency` | `4` | Shared ceiling on active LLM calls across those cases; not a host-wide ceiling |
| `execution.queue_timeout_seconds` | `1800` | Waiting for an LLM budget slot, **not HTTP request timeout** |
| `execution.cache_choices` | `false` | Reuse of eligible LLM responses; does not disable database/literature/tool caches |
| `execution.device` | `cpu` | Retrieval device when retrieval is enabled; does not control Java or disable retrieval |
| `workers.sources` | HPO `3`; VCF `2` | Concurrent source-hypothesis LLM jobs, not raw medical-tool processes |
| `workers.briefs` | `4` | Per-case gene Brief jobs |
| `workers.pairwise` | `4` | Per-case Pairwise jobs within total LLM budget and retry policy |
| `workers.literature` | VCF `4` | VCF literature interpretation and report evidence worker pool; not HPO PubMed retrieval concurrency |
| `workers.post_acmg` | VCF report `4` | Concurrent selected-variant assessments; not four extra global LLM slots |
| `external.phen2gene_acknowledged` | HPO `false` | Set true only after authorizing phenotype upload |
| `external.amelie_acknowledged` | HPO `false` | Set true only after authorizing phenotype/candidate-gene upload |
| `external.amelie_verify_tls` | HPO `true` | Service-specific TLS setting |
| `ranking.similar_cases` | `false` | Optional RareA branch, HPO/clinical only |
| `ranking.reflection` | `false` | Optional additional ranking verification |
| `genotype.max_heap` | VCF `12g` | Maximum Java heap per native-tool process; not total host RAM |
| `genotype.timeout_seconds` | VCF `1800` | Timeout per native-tool invocation, not whole case or LLM |
| `genotype.lirical_transcript_db` | Implicit `REFSEQ` | Native LIRICAL transcript selection; leave unchanged for this profile |
| `explanation_resources.ranking_vep_enabled` | VCF `false` | Local VEP before ranking; independent of post-ranking report annotation |
| `explanation_resources.vep_timeout_seconds` | Implicit `300` | Local VEP invocation timeout when used |
| `explanation_resources.pvs1_timeout_seconds` | Implicit `300` | AutoPVS1 invocation timeout |
| `output.explanation` | `false`; VCF report profile `true` | Report and post-ranking ACMG, not pre-ranking VCF evidence |

HPO templates omit VCF-only worker settings. VCF templates omit phenotype-upload
settings because the native VCF branch does not invoke Phen2Gene or AMELIE.
The rank-only template omits the unused post-ranking worker setting. These
omissions use the configuration defaults.

Important schema defaults differ from minimum profiles: similarity=true,
device=cuda, cache_choices=true, total LLM concurrency=8, source workers=4,
Brief/Pairwise workers=8. Do not remove those explicit minimum overrides.
Per-stage HTTP timeouts, token limits, Prompt text, candidate counts and ranking
formulas are implementation settings, not arbitrary TOML knobs. Do not invent
`timeout`, `top_k` or extra model-routing fields in these files.

Optional service settings are `external.ca_bundle` (custom CA file) and
`external.amelie_pseudonym_secret_env` (defaults to
`MENDELPRIME_AMELIE_PSEUDONYM_SECRET`). If no pseudonym secret is supplied, the
phenotype client creates a private cache secret. The current phenotype clients
do not use a separate configured Phen2Gene/AMELIE API key.

CLI `--lanes`, `--max-llm-concurrency` and `--[no-]explanation` override TOML;
omitting them keeps TOML values. `--input`, `--config`, `--output-root` are the
only required `run` options. The README additionally sets `--cache-root` to
make private work placement explicit; otherwise it defaults to `.cache` beside
the output root. `--audit-root` is optional and can contain patient evidence.

## 4. Resource parameters by profile

All names below are literal supported keys. Their paths, sizes and acquisition
recipes are in the [46-field register](RESOURCE_REGISTER.md). The count is not
the number of downloads: several keys refer to the same installed dataset.
`genotype.assemblies` is an array of tables, written as `[[genotype.assemblies]]`.

| Resource section and fields | HPO / clinical | VCF rank | VCF report |
|---|---|---|---|
| `resources`: `hpo_obo`, `hpo_stats`, `omim_jsonl`, `mim2gene` | Required | Required | Required |
| `resources`: `genes_to_phenotype`, `phen2disease_root` | Required | Not used | Not used |
| `resources`: `similar_case_root`, `embedding_model`, `cross_encoder` | Only if RareA enabled | Not used | Not used |
| `genotype`: `java`, `exomiser_jar`, `exomiser_properties`, `lirical_jar`, `lirical_data`, `gene_map` | No | Required | Required |
| `genotype`: `exomiser_version`, `lirical_version` | No | 14.0.0 (default), 2.4.1 (explicit) | Same |
| `genotype.assemblies[]`: `assembly`, `exomiser_release`, `phenotype_release` | No | GRCh37, 2406, 2406 | Same |
| `genotype.assemblies[]`: `exomiser_genome_data`, `exomiser_phenotype_data`, `lirical_exomiser_data`, `lirical_variants_db`, `lirical_clinvar_db` | No | Required; LIRICAL reuses genome files | Same |
| `variant_knowledge`: `clinvar_db`, `clingen_validity`, `clingen_vcep`, `clingen_manifest` | No | Required | Required |
| `explanation_resources`: `assembly`, `reference_fasta`, `clinvar_protein_disease_db`, `clinvar_protein_release` | No | GRCh37 + required resources/release | Same |
| `explanation_resources`: `gnomad_sites`, `gnomad_coverage`, `population_release`, `gnomad_constraint`, `constraint_release` | No | Required | Required |
| `explanation_resources`: `spliceai_snv`, `spliceai_indel`, `spliceai_release`, `revel`, `revel_release` | No | Required | Required |
| `explanation_resources`: `autopvs1_root` | No | Required | Required |
| `explanation_resources`: `ai_cura_knowledgebase`, `ai_cura_manifest` | No | Required, prepared two-file package | Required |
| `explanation_resources`: `mane_summary`, `mane_release`, `ai_cura_protocol`, `ai_cura_workbook`, `clingen_guideline_manifest` | No | Not read or checked | Required |
| `explanation_resources`: `vep_command`, `vep_cache`, `vep_fasta`, `vep_cache_version`, `vep_tool_version` | No | Only if ranking VEP explicitly enabled | Required: 104 cache / 104.3 tool |
| `explanation_resources`: `clingen_dosage_table`, `genereviews_index`, `post_evidence_catalog`, `kinship_executable` | No | Not needed | Optional enrichment, absent from minimum |

A six-path HPO manifest means nine actual files: the common four, gene
annotations and Phen2Disease's four-file tree. A VCF manifest has many more
companions: Java libraries, native `lib/` directories, FASTA `.fai`, five tabix
indexes, VEP cache contents, AutoPVS1 configuration/data and every file referenced
by active guideline manifests. Do not count only the visible TOML filenames.

The rank-only AI-CURA reader verifies the pinned manifest and derived JSON bytes,
including the declared workbook identity. It does not open the original workbook
or protocol, or load unused S9/ClinGen SVI text into ranking. Reports retain the
full capsule verification. This changes resource loading, not ACMG thresholds.
Inactive report resources do not enter ranking cache fingerprints. A saved
rank-only result can still be passed to `explain` after configuring report resources.

Every explicit release label must describe the installed file. MANE's GRCh38
summary is a transcript registry, not VCF liftover. Only GRCh37 is qualified here.
The local OMIM interface uses no OMIM API key. GeneReviews/curated evidence are
optional scoped resources, not full database coverage merely because configured.

## 5. Network, tools and cold start

| Component | HPO / clinical | VCF | Operation |
|---|---|---|---|
| Application LLM | Yes | Yes | Configured endpoint/key; patient evidence can leave the host |
| Phen2Gene and AMELIE | Yes | No | Remote phenotype tools; explicit upload acknowledgements |
| PubMed / PMC | Yes | Yes | Search, article retrieval and source binding; cached when available |
| ClinGen detail pages | No | Evidence-dependent | Additional curation/mechanism text |
| Ensembl GRCh37 REST | No | Report-dependent fallback | Exact transcript/HGVS checks when local annotation needs reconciliation |
| GeneReviews / NCBI Bookshelf | No | Optional report enrichment | Index maps genes to chapters; chapter text is fetched online |
| Java / Exomiser / LIRICAL | No | Local processes | Configured executable, JAR and data; no Docker or web server required |
| AutoPVS1 / local VEP | No | Local processes | AutoPVS1 for molecular evidence; VEP according to switches above |

There is no single `offline=true` flag that preserves this full pipeline.
Service addresses other than the configurable LLM endpoint are currently
implementation-defined, not general TOML endpoint options. Doctor validates
dependency/layout/startup contracts, not service availability or API capacity.

Cold execution means a new process and empty application cache, not empty
installed databases or cleared OS/provider caches. Report-on can take longer
because of selected-variant annotation, post-ranking evidence retrieval, two-step
ACMG judgments and report calls. A cold HPO run also loads the large Phen2Disease
JSON. Cohorts share resource loading. Hardware/time requirements are not inferred
from archive sizes; Java heap settings are not a total-memory guarantee.
The Phen2Disease similarity matrix is about 3.6 GB of JSON and is loaded before
the first LLM request. The `Loading shared HPO and evidence resources` message
includes this work. Use one Cohort invocation to share resource loading across cases.

Keep HTTPS verification on unless explicitly accepting a documented exception.
`MENDELPRIME_VERIFY_TLS=0` disables authentication checks for application HTTPS
transports, not scientific validation or separate installers/native programs.
It is an environment setting, not another TOML parameter. Prefer a valid CA
configuration or the narrower AMELIE-only setting when sufficient.

If AMELIE fails certificate validation, this is a service/trust-chain issue,
not a missing scientific database. Prefer fixing the certificate/CA chain.
For an explicitly authorized compatibility exception, set
`external.amelie_verify_tls = false` in the existing `[external]` run-TOML
section; leave other HTTPS verification enabled. The application logs this
exception. Neither installation nor a successful doctor check proves that the
remote service's certificate or API is currently usable.

## 6. Resource layout and preflight

- FASTA files require adjacent `.fai` indexes. gnomAD sites/coverage, SpliceAI
  SNV/indel and REVEL require matching `.tbi` or `.csi` indexes.
- VEP needs its complete runtime and `homo_sapiens/104_GRCh37` cache. The
  [VEP guide](VEP_RUNTIME.md) documents executable prefixes, container mounts and FASTA.
- Exomiser and LIRICAL launcher JARs require the adjacent `lib/` dependencies
  declared in their `META-INF/MANIFEST.MF`. A launcher JAR alone is insufficient.
- Exomiser data directories/releases from TOML are passed as Spring command-line
  overrides. Relocate any other custom paths retained in the properties file.
- AutoPVS1 needs its configured source/data, MaxEnt matrices and reference/index.
  Keep required HG19/HG38 annotation tables even for a GRCh37 run; HG38 FASTA is
  optional. Its `DEFAULT.vep_cache` is not used by the selected PVS1 path.
- Keep the AI-CURA pinned manifest unchanged. ClinGen SVI text paths and reviewed
  evidence files resolve relative to their owning manifests; copy complete directories.
- RareA requires aligned `records.jsonl`, `record_ids.json` and `embeddings.npy`,
  plus both model directories. See [package installation](RESOURCE_ASSETS.md).

`doctor` needs no patient or API key. It checks dependencies, layouts and indexes;
`--check-tools` additionally starts native version/help checks and an AutoPVS1
import. For containerized VEP, use `--work-root` with an existing mounted directory.
For a self-contained resource export, `--portable-root` rejects paths and symlinks
escaping that bundle. Shared system installations can omit that export check.
Doctor does not verify network availability, model capacity or scientific performance.
