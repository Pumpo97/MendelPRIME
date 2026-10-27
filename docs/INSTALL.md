# Installing MendelPRIME

Use the current source checkout, Linux x86-64 and Python 3.12. Installation is
separate from inference: `install.py` never submits a patient or calls an LLM.
The resource archives remain pinned to Release 1.0.0; its original source
archive predates the installer. Clone the current branch for these commands.

## Prerequisites

- Python 3.12 with venv/pip. With an existing Conda/Miniforge installation, run
  `conda create -n mendelprime-bootstrap -c conda-forge python=3.12 pip`, then
  `conda activate mendelprime-bootstrap` before using the installer.
- VCF requires Java 17: a JRE for prepared programs, a JDK only for source builds.
  Pass `--java` with the executable, or put a complete Java installation on PATH.
- VCF also needs libpq for the psycopg2 dependency of HGVS. A source build needs
  a C compiler, Python/PostgreSQL development headers and pg_config on PATH.
  Alternatively add libpq to the user-owned Conda environment and keep its bin
  on PATH during installation. Python extras do not install OS shared libraries.
- No root access or replacement of system Java is required.

## Capability selection

Set WORK, MENDELPRIME_BASE_URL and MENDELPRIME_MODEL as in the
[README](../README.md#quick-start). Use an actual model ID supported by your
endpoint; examples intentionally do not select a specific model.

| Installer selection | Python group | Configuration under WORK/configs |
|---|---|---|
| `--input-type hpo` or `clinical` | Base | `phenotype-rank/run.toml` and `resources.toml` |
| Either above plus `--report` | Same base | `phenotype-report/run.toml` and `resources.toml` |
| `--input-type vcf` | VCF | `vcf-rank/run.toml` and `resources.toml` |
| `--input-type vcf --report` | Same VCF | `vcf-report/run.toml` and `resources.toml` |
| HPO/clinical plus `--with-rarea` | Retrieval | Same directory name with `-rarea` appended |

HPO and clinical share resources. Case/Cohort and trio/proband-only do not
require different installations. Mixed VCF/non-VCF cohorts need a manually
configured union; the installer selects one resource family.

```bash
# Read-only plan; no endpoint, key or installed databases needed:
python3.12 install.py --input-type vcf --report --work "$WORK" --plan

# Install with already downloaded Release archives:
python3.12 install.py --input-type vcf --work "$WORK" --archives /absolute/downloads --java /absolute/jdk17/bin/java

# Reuse installed resources without downloading or extracting archives:
python3.12 install.py --input-type vcf --work "$WORK" --resources /absolute/resources

# Generate the TOML pair only; no environment installation or doctor:
python3.12 install.py --input-type hpo --work "$WORK" --configure-only
```

`--download-resources` downloads only the selected prepared archives, verifies
their SHA256 against resources/ASSETS.json, and extracts them. `--archives`
imports the same files locally. With both, missing archives are downloaded there.
With neither, resources must already be installed. VCF uses the Release LIRICAL
ZIP by default, not a source build.

The installer creates/reuses WORK/.venv, installs the selected Python group with
the existing constraints from a temporary build copy, then runs pip check and
doctor with native-tool checks. Build artifacts stay outside the source checkout.
Missing tools, libraries or data produce a nonzero exit; installing Python
packages alone is not a successful full deployment. Missing resources are not
silently disabled.

## Configuration and repeated installation

Paths are relative to resource_root unless absolute. TOML does not expand shell
variables. The installer writes an absolute root and discovers Java on PATH when
available. Scientific defaults, Prompts and tool analysis policies are unchanged.

| User setting | Generated setting |
|---|---|
| `--base-url` / `MENDELPRIME_BASE_URL` | run.toml: `llm.base_url` |
| `--model` / `MENDELPRIME_MODEL` | run.toml: `llm.model` |
| `--wire-api` | `llm.wire_api`: chat_completions by default, or responses |
| `--resources` | resources.toml: top-level `resource_root` |
| `--java` | resources.toml: `genotype.java` |
| `--report` | run.toml: `output.explanation`; selects VCF report resources |
| `--with-rarea` | `ranking.similar_cases`; adds three retrieval paths |
| `--ack-phenotype-upload` | `external.phen2gene_acknowledged` and `amelie_acknowledged`; use only with authorization |
| API key | Runtime environment `MENDELPRIME_API_KEY`; never saved by the installer |
| Concurrency | run.toml: `execution.lanes`, `max_llm_concurrency`, `workers`; defaults are 1 lane / 4 global slots |

Existing configuration pairs are reused, never overwritten. Edit them directly;
conflicting explicit installer arguments fail instead of being ignored. Ranking
and report profiles have separate directories. Identical resource files are
retained; different existing files cause an error. Run one installer at a time.
Keep WORK, archives and installed resources outside CODE.
Archive imports use temporary extraction space and verify existing payloads;
allow space for both the archive and expanded data. They do not delete an older
installation or replace a different database snapshot.

For HPO/clinical, adding `--explanation` at runtime with the rank profile is also
sufficient once resources are ready. VCF needs the report resource profile, not
just a runtime flag.

## Resource-to-parameter map

The last column names fields in the generated resources.toml. Exact versions,
companion indexes and all 46 fields remain in the [register](RESOURCE_REGISTER.md).

| Needed by | Resource and acquisition | Preparation/install route | TOML fields |
|---|---|---|---|
| All | HPO Release package | Installer extracts unchanged | `resources.hpo_obo`, `hpo_stats`; phenotype also `genes_to_phenotype` |
| All | Authorized complete OMIM mirror + matching mapping | [Schema](OMIM_MIRROR.md), [OMIM preparation](../preparation/README.md#omim-and-genereviews) | `resources.omim_jsonl`, `mim2gene` |
| HPO/clinical | Phen2Disease upstream dataset and similarity matrix | [Four-file layout](RESOURCE_PREPARATION.md#p01-phen2disease-upstream-precomputed-data-not-another-inference-run); no upstream ranking run | `resources.phen2disease_root` |
| VCF | Java 17 and complete Exomiser 14.0.0 distribution | [Native installation](RESOURCE_PREPARATION.md#p04-java-exomiser-and-lirical-install-full-programs-and-separate-data); retain lib and relocate properties | `genotype.java`, `exomiser_jar`, `exomiser_properties` |
| VCF | LIRICAL 2.4.1 Release ZIP | Installer preserves JAR + lib; [manual route](LIRICAL.md) | `genotype.lirical_jar`, `lirical_version` |
| VCF | Exomiser 2406 genome/phenotype archives; LIRICAL data/HGNC | [Data acquisition](LIRICAL.md#install-the-data-separately); genome .mv.db files are shared | `genotype.assemblies[0].exomiser_genome_data`, `exomiser_phenotype_data`, `lirical_exomiser_data`, `lirical_variants_db`, `lirical_clinvar_db`; `genotype.lirical_data`, `gene_map` |
| VCF | Both ClinVar SQLite databases in Release | Installer extracts; [raw conversion](../preparation/README.md#clinvar-databases) only for separately acquired snapshots | `variant_knowledge.clinvar_db`; `explanation_resources.clinvar_protein_disease_db`, `clinvar_protein_release` |
| VCF | ClinGen validity and VCEP exports | [Acquire and bind exports](../preparation/README.md#clingen-tables-and-guidance) | `variant_knowledge.clingen_validity`, `clingen_vcep`, `clingen_manifest` |
| VCF | Verified GRCh37 reference | [FASTA/index preparation](../preparation/README.md#fasta-references) | `explanation_resources.reference_fasta` plus adjacent .fai |
| VCF | gnomAD sites/coverage/constraint, SpliceAI SNV/indel, REVEL | [Acquisition and BGZF/tabix preparation](../preparation/README.md#revel-and-populationprediction-indexes) | `explanation_resources.gnomad_sites`, `gnomad_coverage`, `gnomad_constraint`, `spliceai_snv`, `spliceai_indel`, `revel` and their release fields |
| VCF | AutoPVS1 pinned source/data and reference | [Compatibility installation](../preparation/README.md#autopvs1-compatibility-installation) | `explanation_resources.autopvs1_root` |
| VCF | AI-CURA ranking knowledge + manifest in Release | Installer extracts; original workbook not needed for ranking | `explanation_resources.ai_cura_knowledgebase`, `ai_cura_manifest` |
| VCF reports only | VEP 104.3/cache 104/FASTA; MANE | [VEP](VEP_RUNTIME.md), [MANE](RESOURCE_PREPARATION.md#p11-mane-and-autopvs1) | `explanation_resources.vep_command`, `vep_cache`, `vep_fasta`, `vep_cache_version`, `vep_tool_version`, `mane_summary`, `mane_release` |
| VCF reports only | AI-CURA workbook/protocol; ClinGen SVI originals/text | [AI-CURA capsule](../preparation/README.md#ai-cura-capsule), [SVI](../preparation/README.md#clingen-tables-and-guidance) | `explanation_resources.ai_cura_protocol`, `ai_cura_workbook`, `clingen_guideline_manifest` |
| Optional HPO/clinical retrieval | RareA and retrieval-support packages; matching Qwen weights | Installer extracts packages; [complete query model](RESOURCE_ASSETS.md#complete-the-query-model) separately | `resources.similar_case_root`, `embedding_model`, `cross_encoder` |
| Optional VCF report enrichment | GeneReviews, dosage, reviewed literature, PLINK2 | [Scoped resources](RESOURCE_PREPARATION.md#p13-optional-report-enrichment-and-kinship-qc) | `explanation_resources.genereviews_index`, `clingen_dosage_table`, `post_evidence_catalog`, `kinship_executable` |

**VCF ranking still requires pre-ranking molecular evidence/ACMG.** Not all
explanation_resources fields are report-only; use the table and selected profile.
HPO reports require no variant databases. The one-record OMIM demo is not a
substitute for a complete database.

## Offline packages and optional retrieval

First obtain a complete Python 3.12 / Linux wheelhouse including build
dependencies, and the selected Release archives:

```bash
python3.12 install.py --input-type vcf --work "$WORK" --offline --wheelhouse /absolute/wheels --archives /absolute/downloads
```

This disables pip indexes and installer downloads, not network use during
inference. The application still needs LLM and evidence services. Doctor checks
local resources/startup, not service availability or model quality.

For containerized VEP, mount the selected WORK and resource directories at their
configured paths. The installer puts doctor work under WORK/cache/doctor; it
does not assume that an arbitrary host temporary directory is visible inside VEP.

`--with-rarea` installs the retrieval extra. Select the appropriate PyTorch wheel
source for CPU/CUDA when necessary and configure execution.device in the
generated TOML. CPU retrieval still needs model weights. The source extras hpo
and full mean Python dependency groups, not installed native tools/databases.

## Manual installation and source builds

The existing manual interface remains: create a Python 3.12 environment, then
install `.` for phenotype without RareA, `.[hpo]` for retrieval, `.[vcf]` for
VCF, or `.[full]` for both, using `-c constraints/runtime.txt`. Select one group.
Dependency-only files are under requirements/. Copy the desired templates and
replace endpoint/model placeholders and resource paths before running.

To build the source distribution/wheel, install the tools pinned by
constraints/build.txt and run `python -m build --no-isolation`. Keep deployment
files when using a wheel: it installs runtime and Prompts, not install.py,
configs, demos, preparation scripts, native tools or resources. The single
inference entry remains `python main.py` / `mendelprime`.
