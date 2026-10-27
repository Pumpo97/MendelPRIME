# Resource origin, preparation and declarations

This guide separates **prepared resources** from **raw downloads that need
processing**. It covers all 46
resource path/command fields in the HPO and VCF profiles; clinical input reuses
the HPO resources. The measured item-by-item table is [RESOURCE_REGISTER.md](RESOURCE_REGISTER.md)
and the machine-readable record is [RESOURCE_REGISTER.json](../resources/RESOURCE_REGISTER.json).
Aliases and nested files overlap: 46 fields are not 46 independent downloads.

## Categories and publication decision

- **A, prepared project asset:** extract the supplied archive, verify hashes and
  configure its paths. Do not redo the original transformation during inference.
- **B, upstream ready-made/precomputed:** no scientific recalculation is needed,
  but obtain the correct release, all companion files, and the required layout.
- **C, preparation required from raw inputs:** explicit conversion, indexing,
  extraction or derived metadata is necessary. A renamed raw file is not enough.
- **D, install/configure:** deploy a tool and its dependencies; this is distinct
  from transforming scientific data.

HPO, RareA, ClinVar, AI-CURA ranking guidance and the LIRICAL distribution are
provided as separate resource packages. Retrieval support contains the complete
MedCPT reranker and Qwen configuration/tokenizer files, not the Qwen weight
shards. These packages are separate from the source repository.

Executable preparation utilities are now shipped under preparation/. Start with
[their README](../preparation/README.md) for original acquisition links, commands,
outputs and failure behavior. They run outside inference and need no old project
checkout. Existing reviewed/pinned resources are preserved; manual biological
curation is not replaced with an automatic generator.
The shipped scripts use `pysam` for BGZF/tabix; standalone HTSlib executables are
needed only if choosing the optional manual commands below. Preparation can
also require GNU sort, Poppler pdftotext, Git/JDK and the upstream Maven Wrapper.
These utilities are not all installed by the Python extras and are unnecessary
at inference time when their verified prepared outputs are already supplied.

## P00. HPO and RareA: distribute the prepared assets

**Input:** the verified project snapshots. **Action:** install the archives in
[RESOURCE_ASSETS.md](RESOURCE_ASSETS.md), without editing their contents.

- HPO: `core/hpo/{hp.obo,hpo_stats.json,genes_to_phenotype.txt}`. The ontology is
  `hp/releases/2025-03-03`; statistic/count consistency and byte identity were
  checked. Preserve the documented 60 annotation rows outside this ontology.
- RareA: `phenotype/similar-cases/index/rarea/{records.jsonl,record_ids.json,embeddings.npy}`.
  Keep all 5,867 records in their original vector-row order. Do not publish only
  the vectors and omit the records needed for source-bound explanations.

**Acceptance:** archive/payload hashes, row/ID alignment, correct model binding
and real local retrieval. The small HPO pack does not contain Phen2Disease's
independent data. Embeddings do not by themselves anonymize underlying text;
preserve source attribution and the applicable corpus publication conditions.

## P01. Phen2Disease: upstream precomputed data, not another inference run

**Origin:** [official repository](https://github.com/ZhuLab-Fudan/Phen2Disease),
including its linked precomputed similarity matrix. The application consumes
only these four paths under `phenotype/phen2disease/`:

```text
data/matrix/lin_similarity_matrix.json
data/association/Disease-Hpo/disease2hpo20210413.json
data/association/Gene-Disease/disease2genecard2021.json
data/obo/hp20210413.obo
```

**Path:** obtain the matching repository data and matrix archive -> decompress
the matrix JSON -> arrange these four paths -> record hashes -> open them with
the existing phenotype reader. Current selected payload is about 3.86 GB.
No model training or patient analysis from the upstream pipeline is needed for
this deployment. Do not replace its 2021 ontology/mappings with the common HPO
snapshot merely because both are HPO resources. This four-file dataset is not
included in the project's prepared resource archives; obtain it from the upstream source.

Example installation after obtaining the matrix from the link in the upstream
README. Set RAW and RES to your private download/resource directories; the target
must be unused. This copies only the consumed resources, not upstream test cases:

```bash
set -e
git clone https://github.com/ZhuLab-Fudan/Phen2Disease.git "$RAW/Phen2Disease"
git -C "$RAW/Phen2Disease" rev-parse HEAD
# Download lin_similarity_matrix.json.gz using the upstream README's link
# and place it at "$RAW/lin_similarity_matrix.json.gz".
test -s "$RAW/lin_similarity_matrix.json.gz"
test ! -e "$RES/phenotype/phen2disease"
mkdir -p "$RES/phenotype/phen2disease/data/matrix"
for item in association/Disease-Hpo/disease2hpo20210413.json association/Gene-Disease/disease2genecard2021.json obo/hp20210413.obo; do
  test -s "$RAW/Phen2Disease/data/$item" || exit 1
  mkdir -p "$RES/phenotype/phen2disease/data/$(dirname "$item")"
  cp "$RAW/Phen2Disease/data/$item" "$RES/phenotype/phen2disease/data/$item"
done
gzip -dc "$RAW/lin_similarity_matrix.json.gz" > "$RES/phenotype/phen2disease/data/matrix/lin_similarity_matrix.json"
```

The block uses `set -e` so a failed check/decompression stops installation.
Retain the printed upstream commit and checksums of all four
files. A new main-branch checkout is not a declaration of historical snapshot
identity. No `pip install Phen2Disease` or upstream patient-ranking command is
needed. Doctor checks the layout; an ordinary HPO run loads the actual tables.

## P02. OMIM: local mirror conversion, plus an ordinary gene map

**Input:** an authorized local OMIM mirror and a corresponding mim2gene table.
**Path:** retain entry identifiers and source text -> group paragraphs by their
actual section -> bind PMIDs to the actual paragraph/allelic-variant reference
-> write one object per JSONL line -> join by MIM ID using mim2gene.

**Output:** `core/omim/omim_clean_full.jsonl` and `core/omim/mim2gene.txt`.
The latter is a directly readable tab-delimited mapping; the former is a project
normalization format, not raw HTML or an official API response. Current mirror
payload is about 528 MB. See [OMIM_MIRROR.md](OMIM_MIRROR.md) for the exact schema,
one-record demo and runnable reader test. The `omim` preparation command validates
and installs that normalized mirror. It does not scrape arbitrary HTML or infer
reference bindings; users provide their existing complete local mirror. No OMIM API
implementation is required. Never treat the demo as full cohort coverage.

## P03. Retrieval models: download exact weights or reuse the support asset

**Qwen:** download the two weight shards at the fixed revision documented in
RESOURCE_ASSETS.md; combine with the matching configuration/tokenizer sidecars.
This is a model installation, not regeneration of document vectors. Verify the
model checksum list before loading the RareA index.

**MedCPT:** use the supplied support archive for exact reproduction. Its weights
are saved as safetensors; the upstream repository also offers a PyTorch binary.
If separately preparing a conversion, load the trusted upstream snapshot with
the supported Transformers version, save via `save_pretrained(...,
safe_serialization=True)`, retain the tokenizer/configuration and compare every
state-dict tensor before/after. That conversion is optional for the generic
reader, but an independently downloaded model is not declared equivalent to the
verified support asset without checking it. Do not substitute another model,
query instruction or embedding dimensionality during deployment.

## P04. Java, Exomiser and LIRICAL: install full programs and separate data

- Java: install a complete Java 17 runtime; use a JDK for source builds. Configure
  the absolute executable. A copied `bin/java` without its runtime libraries is
  insufficient; no system Java replacement or root execution is required.
- Exomiser: unpack the complete 14.0.0 distribution, preserving its declared
  `lib/` class path. Configure the properties file and the 2406 genome and
  phenotype directories. The adapter supplies resource-path overrides; do not
  alter the scientific analysis settings while relocating other custom paths.
- Exomiser data: unpack the official **2406_hg19** and **2406_phenotype** releases
  from the [upstream data release](https://github.com/exomiser/Exomiser/discussions/562).
  These are ready-made native databases, not files to rebuild from the project
  HPO pack or ClinVar SQLite. The selected payloads are about 28.39 GB and 12.79 GB.
- LIRICAL: install the supplied 2.4.1 ZIP or follow the tested [official-source
  build](LIRICAL.md#optional-source-build). The default [installer](INSTALL.md)
  uses the prepared Release ZIP. Its data directory is separate. The upstream `download`
  command obtains its HPO/HGNC/Orphanet/MedGen/transcript data, but a fresh download
  is not automatically the historical snapshot; record its versions/hashes.

LIRICAL reuses the Exomiser genome directory and its variants/ClinVar `.mv.db`
files. `lirical_exomiser_data`, `lirical_variants_db` and `lirical_clinvar_db`
are not three additional datasets to download. `gene_map` is already inside the
LIRICAL data directory. **Acceptance:** class-path checks, actual version/startup
checks and a real VCF run with the existing native parser, not just a Java banner.

### Exomiser installation and TOML binding

Get the full `exomiser-cli-14.0.0-distribution.zip` from the
[official release](https://github.com/exomiser/Exomiser/releases/tag/14.0.0), not
the web/REST JAR. Download both 2406 databases from the data-release page above.
Verify checksums from the selected release where provided, inspect the archive
layout, and unpack into a new versioned resource directory. The final tree must
retain the CLI's libraries and separate database directories:

```text
genotype/exomiser/
  exomiser-cli-14.0.0.jar
  lib/...
  application.properties
  data/2406_hg19/2406_hg19_variants.mv.db
  data/2406_hg19/2406_hg19_clinvar.mv.db
  data/2406_hg19/...transcript and other native files...
  data/2406_phenotype/...
```

The extracted directory may instead keep its upstream versioned name; absolute
TOML paths support this without copying large files. With the downloaded ZIP in
RAW and a new installation directory, inspect and unpack it:

```bash
set -e
unzip -l "$RAW/exomiser-cli-14.0.0-distribution.zip"
mkdir -p "$RES/genotype"
test ! -e "$RES/genotype/exomiser-cli-14.0.0"
unzip "$RAW/exomiser-cli-14.0.0-distribution.zip" -d "$RES/genotype"
export EXOMISER_HOME=/absolute/path/to/the/extracted/exomiser-cli-14.0.0
export JAVA_HOME=/absolute/path/to/complete/jdk17
"$JAVA_HOME/bin/java" -version
"$JAVA_HOME/bin/java" -jar "$EXOMISER_HOME/exomiser-cli-14.0.0.jar" --version
```

The Java directory must be a complete Linux x86-64 installation, not a copied
executable. A Java 17 JRE is enough to run the tools; a JDK is needed only when
building LIRICAL. Use an administrator-provided installation or the official
[Temurin 17 distributions](https://adoptium.net/temurin/releases/?version=17),
retaining the runtime libraries and matching host architecture.

Set `genotype.java` and `genotype.exomiser_jar` to those actual absolute paths in
the resource TOML. Point `genotype.exomiser_properties` at the complete supplied
`application.properties`, or at a private relocated copy. The template filename
`application-2406.properties` is a destination convention, not a separately
downloaded database. Do not create an empty properties file merely to pass an
existence check. Keep other upstream settings; relocate custom optional paths.

Set `genotype.assemblies[0]` to GRCh37, releases `2406`, and the actual genome and
phenotype directories after unpacking their archives. MendelPRIME passes these
locations as Spring command-line overrides to `java -jar`; its adapter writes
the per-case analysis YAML and selects the declared samples/PED internally.
LIRICAL is also invoked by `java -jar`, using the same genome databases and its
own data directory. No Java server, Python-to-Java bridge or independent Offline
framework needs installation. Run the minimal VCF doctor with `--check-tools`
after both tool installations and all mandatory molecular resources are ready.

## P05. ClinVar: supplied SQLite is ready; raw upstream files are not

The recommended user path is to extract the provided two-database archive.
For reconstruction, [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
provides the source files, but they must be mapped into the runtime schemas:

```text
variant_summary.txt.gz + var_citations.txt[.gz], same recorded snapshot
  -> variants / identifiers / citations / metadata tables and indexes
  -> variant-knowledge/clinvar/clinvar_20260117.sqlite

that SQLite + the matching variant_summary.txt.gz
  -> GRCh37 pathogenic/likely-pathogenic comparator rows
  -> transcript, coding HGVS, protein HGVS, codon and disease fields
  -> assessment/grch37/clinvar-protein/clinvar_protein_acmg_disease_20260117.sqlite
```

The historical main builder uses **PositionVCF, ReferenceAlleleVCF and
AlternateAlleleVCF**, not generic HGVS-position columns. Preserve VariationID,
AlleleID, assembly, gene, classification, review status and citations. The
protein builder uses the matching summary file to retain PhenotypeIDS and
PhenotypeList; omitting it can lose disease context. Its fixed comparator filter
is part of the existing resource recipe, not permission to invent pathogenicity.

**Implementation status:** the preserved builders now live in preparation/ and
are invoked by `prepare_resources.py clinvar` and `clinvar-protein`. They have
no old-project import/path dependency. The protein command requires the matching
summary and checks the release; existing outputs are not overwritten. Use the
supplied SQLite when exact existing file hashes are required. **Acceptance:** SQLite
integrity/schema/index checks plus exact-allele and same-transcript/codon reader
queries with source/citation checks. Do not rename a VCF to `.sqlite` or replace
the Exomiser `.mv.db` with this database.

## P06. ClinGen exports: use original rows and add a binding manifest

**Input:** [gene-validity CSV](https://search.clinicalgenome.org/kb/gene-validity/download)
and [VCEP classification export](https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download?type=csv).
The current validity file includes an official preamble; the reader already
locates its header, so removing the preamble is unnecessary. The current VCEP
file is tab-separated. Inspect the response, not only a URL's `type=csv` label;
if an export is actually comma-separated, use a CSV parser/writer to change the
delimiter while preserving quoted text, headers and every field, never string
replacement over clinical summaries.

**Output:** `variant-knowledge/clingen/{gene_disease_validity.csv,
variant_pathogenicity.tsv,resource_manifest.json}`. The metadata manifest must
have schema `mgagent_clingen_resource_v1`, an ISO `release_date`, and
`files.<filename>.sha256` for both exports. Generate hashes after any necessary
delimiter conversion and retain original URLs. **Acceptance:** required columns
and field lengths in `knowledge/variants.py`, manifest hashes and lookup tests.
The 2026-10-27 live export is not the same snapshot as the configured 2026-08-31
files; a data update requires a new resource version, not a silent replacement.

## P07. References: preserve sequence identity, generate matching indexes

**Input:** the validated GRCh37/hg19 FASTA, not an arbitrary file with the same
basename. **Path:** decompress into the configured `.fa` -> verify contigs/build
and representative REF alleles -> generate `.fai` with `pysam.faidx` or samtools
-> record FASTA/index hashes together. Current paths are
`assessment/grch37/reference/hg19.fa` and
`assessment/grch37/vep/final_hg19.fa`, each with an adjacent `.fai`.

These existing files have different byte sizes (about 3.20 GB and 3.15 GB).
Do not assume they are interchangeable or deduplicate them by name alone; check
sequence content, wrapping, contigs and the consuming tool's requirements.
BGZF FASTA also needs its appropriate compression index if that route is used.
Never fix an unexplained stale-index warning with `touch`; rebuild in a new
location using the `fasta` preparation command when the original data/index
binding cannot be established.

For mitochondrial records, also follow the [reference-identity boundary](INPUTS.md#vcf-inputs):
the older hg19 chrM and rCRS are not interchangeable. Replacing a FASTA or changing
its filename does not convert the VCF coordinates.

## P08. VEP: install a runtime, not just a cache directory

**Path:** obtain the pinned 104.3 executable/image and dependencies -> unpack a
matching **homo_sapiens/104_GRCh37** cache -> configure FASTA/index -> expose the
actual executable/container prefix as `vep_command` -> test mounts and JSON
annotation output. Preserve the current Ensembl cache type; an unrelated RefSeq
cache from another tool's instructions is not automatically equivalent.

The argument-vector prefix must provide its own Perl/Python/container/image
store and mounts. Only its first executable is resource-root-resolved; other
arguments are literal except `{workdir}`. Mount both data and work directories.
Do not duplicate scientific VEP options already supplied by the adapter.
Install the VEP runtime separately from the resource cache. Verify it with
`doctor --check-tools --explanation` followed by real annotation with the normal
adapter. VEP offline mode does not make report literature retrieval offline.

## P09. gnomAD: sites are native VCF; coverage needs an index

- Sites: use native **2.1.1 exomes GRCh37** VCF plus its matching `.tbi`/`.csi`,
  retaining REF/ALT, FILTER, AC, AN, AF and population fields. Do not replace
  them with an ANNOVAR-style summary or fill missing AF/AN values with zero.
- Coverage: obtain the **2.1 exome coverage** summary used alongside sites 2.1.1.
  The current header starts with `chrom`, `pos` and includes `over_20`. The
  reader accepts the original non-comment header. For a sorted BGZF file without
  an index, a prep-only tabix command is:

```bash
tabix -S 1 -s 1 -b 2 -e 2 gnomad.exomes.coverage.summary.tsv.bgz
```

This command assumes the first line is the single non-comment header and data
positions are 1-based; do not use it blindly with another layout. Preserve a
valid supplied index instead of rebuilding unnecessarily. Plain gzip is not
BGZF. Coverage output stays a TSV, not a VCF or fabricated AF table.
- Constraint: use the official **2.1.1** compressed gene table with the header
  expected by `assessment_resources.py`; no tabix is required by this reader.

**Outputs:** the three configured files under `assessment/grch37/gnomad/`, with
site/coverage indexes. **Acceptance:** target-allele query, coverage positions
and `over_20` values, missing/AN=0 behavior, and exact header/version checks.
See [gnomAD's GRCh37 resource definitions](https://broadinstitute.github.io/gnomad_methods/api_reference/resources/grch37/gnomad.html).

## P10. SpliceAI and REVEL: distinguish native predictions from converted tables

**SpliceAI:** obtain both masked **v1.3 hg19** SNV and indel VCFs plus matching
indexes from the [upstream resource](https://github.com/Illumina/SpliceAI).
Preserve VCF columns and the allele/gene-specific SpliceAI INFO values. No
TensorFlow installation or model rerun is necessary when using these tables.
Do not replace them with SpliceAI2, GRCh38 or unmasked scores as a format fix.

**REVEL path:** [official score CSV](https://sites.google.com/site/revelgenomics/downloads)
-> parse CSV -> retain/reorder exactly these nine columns -> sort by contig and
numeric hg19 position -> write tab-delimited text with the `#chr` header -> BGZF
compress -> tabix index on column 1/2/2 (1-based coordinates).

```text
chr  hg19_pos  grch38_pos  ref  alt  aaref  aaalt  REVEL  Ensembl_transcriptid
```

The separators in the actual file must be tabs. Preserve scores, missing values,
alleles and semicolon-separated transcript lists; do not use GRCh38 positions
as hg19 coordinates or silently discard unparseable rows. Write/audit failures
before accepting a new prepared file. **Output:**
`assessment/grch37/predictions/new_tabbed_revel.tsv.gz` plus `.tbi`.
**Acceptance:** exact header and allele queries through `annotation.py`.
The current file is already converted. The shipped `revel` command implements
the stable, bounded-memory conversion and index construction; `coverage` and
`vcf-index` handle the corresponding indexed resource formats. See preparation/
README.md for commands and required input formats.

## P11. MANE and AutoPVS1

**MANE:** download the official [1.5 summary](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/)
and preserve the compressed TSV header. Its GRCh38 coordinates identify
transcript context; they do not convert the input VCF to GRCh38. Configure
`assessment/transcripts/MANE.GRCh38.v1.5.summary.txt.gz` unchanged.
This resource is required for VCF reports, not rank-only execution. AutoPVS1 below
remains a ranking dependency.

**AutoPVS1 path:** obtain the [upstream source](https://github.com/JiguangPeng/autopvs1)
at the verified commit recorded in RESOURCE_REGISTER.json -> include its data,
MaxEnt matrices and Python support -> configure local paths in config.ini ->
supply HG19 FASTA/index and both HG19/HG38 annotation tables required by eager
imports -> preserve the compatibility change described below -> test the
isolated import and real PVS1 callback. The runtime's Python dependencies are
covered by the VCF package group; a source clone alone does not supply FASTA.

**Disclosed local modification:** the verified checkout changes the two eager
FASTA loads in read_data.py to a helper returning `Fasta(path)` only when the
file exists. MendelPRIME still requires the HG19 reference before running the
GRCh37 path. This permits an unconfigured HG38 FASTA to remain absent, without
changing PVS1 decision rules. The base commit is
`7fb1be97667e5ef576f81bf2fabbddcf9a4c7594`; the modified file hash is recorded in
the register. An unmodified clone needs its own complete eager-load dependencies.
The compatibility patch is part of the external tool snapshot, not silently
present in upstream source. Preserve licensing/attribution and record the patch.
Changing config.ini locations must not alter decision tables or thresholds.

At that pinned commit, the existing imports already provide `os` and `Fasta`.
The disclosed replacement for the two eager assignments is:

```python
def optional_fasta(path):
    return Fasta(path) if os.path.isfile(path) else None

genome_hg19 = optional_fasta(config['HG19']['genome'])
genome_hg38 = optional_fasta(config['HG38']['genome'])
```

This is the previously validated external compatibility change, not a new rule.
The shipped `autopvs1` installer verifies the pinned tracked files, applies only
this replacement, checks the resulting hash and performs an isolated import.

## P12. AI-CURA and ClinGen guidance: derived, pinned resource capsules

For rank-only deployment, extract `mendelprime-vcf-ranking-guidance-v1.tar.gz`.
It supplies the existing derived JSON and pinned manifest; no original XLSX,
protocol or SVI text is needed at ranking runtime. The full recipes below are
for rebuilding that knowledge or enabling VCF reports.

**AI-CURA input:** the [paper/supplement](https://doi.org/10.1126/scitranslmed.adz4172)
and original S1-S13 workbook. **Path:** retain workbook bytes -> extract the
S8/S9/S10 tables into the current `ai_cura_knowledgebase_v1` JSON -> retain the
reviewed protocol Markdown -> bind all source/derived artifacts in the current
manifest. **Outputs:** four configured files under `guidelines/ai-cura/`.
The workbook itself is ready-made; the JSON, protocol and manifest are derived.

**Preparation:** the table extractor is included in
preparation/. The `aicura` command reconstructs the knowledge base from the exact
workbook and installs the unchanged reviewed protocol with the pinned public
manifest supplied as small preparation assets. The public manifest replaces
machine paths with basenames without changing artifact hashes, rule boundaries
or versions. The reader accepts exactly its hash or the original legacy hash;
the report loader checks the protocol, workbook and derived knowledge. Ranking
checks the pinned manifest and derived knowledge against the declared workbook
identity without opening the report-only originals. Arbitrary
equivalent-looking manifests are not accepted. The actual manifest hash is
retained in provenance.

**ClinGen SVI path:** obtain the original PM2, PM3 and PS2/PM6 guidance versions
identified in the manifest -> extract UTF-8 text with the recorded
`pdftotext -layout` procedure -> preserve original documents -> create relative
manifest entries `name`, `url`, original `sha256`, `text_file`, `text_sha256`.
The active reader requires the manifest plus referenced text files; originals
are copied by the shipped `svi` command for traceability. Different extraction output is a new
text snapshot, not a reason to bypass hashes. **Output:** complete
`guidelines/clingen-svi-20261021/` directory. Text-only copies are not documents
with verified provenance unless those bindings are retained.

## P13. Optional report enrichment and kinship QC

- ClinGen dosage: obtain the GRCh37 tabular snapshot from the official
  [ClinGen dosage resource](https://search.clinicalgenome.org/kb/gene-dosage), preserving comments,
  gene identifiers, scores and dates. Place it at the configured dosage path.
- GeneReviews: preserve a gene/chapter mapping with exact columns
  `NBK_id`, `GR_shortname`, `genesymbol`. The current configured file is a
  **three-gene subset**, not the complete official index. To reproduce it,
  select the same rows from the matching official mapping; expanding coverage
  is a resource change. Chapters are still fetched from NCBI on first use and
  cached; an index is not an offline full-text collection.
- Curated original evidence: original ClinVar/XML/PDF/JATS -> stable UTF-8 text
  -> manual genomic/allele/transcript identity and quote verification -> optional
  independently checked family-series ledger -> `verified_original_evidence_v1`
  catalog. Include every referenced original and text file with hashes and
  relative paths. Do not insert Gold grades or precomputed ACMG conclusions.
  This is a scoped reviewed collection, not an automatic universal conversion.
  The current catalog has three records covering two genes; its full referenced
  file closure is about 0.99 MB, not merely the 25 KB catalog JSON.
- PLINK2: download a Linux build from the [official installation page](https://www.cog-genomics.org/plink/2.0/),
  unpack it under a user-owned tools directory, run `plink2 --version`, and set
  `explanation_resources.kinship_executable` to its absolute executable. Existing
  kinship QC uses actual samples; no cohort truth/resource table substitutes
  for family observations. It is optional report enrichment when configured.

**Acceptance:** validate the actual index/catalog schema, complete file closure,
quote substrings and allele identity. For GeneReviews, also check the fetched
chapter's expected mechanism/inheritance/table sections; page layout or network
failures are not evidence of absent biology. Keep optional scope explicit.
The `genereviews` command performs explicit subset selection. `extract-text`
prepares PDF/XML text views; `curated` requires --ack-reviewed and reuses the
unchanged runtime source/allele/family/quote validators before publishing a full
catalog closure. These scripts do not fabricate the required human review.

## Online dependencies (not downloadable resource archives)

Python 3.12 and the selected `base`/`hpo`/`vcf`/`full` package group are runtime
prerequisites listed in INSTALL.md and requirements/. Native libpq/build headers
may also be needed. A site-specific CA bundle, if configured, belongs to the
private deployment rather than a scientific data archive. Patient VCF/PED and
phenotype inputs are supplied by the user; Gold labels are evaluation-only.

| Service | Role | Required declaration |
|---|---|---|
| Configured LLM | Ranking and optional report calls | User-provided key/model/endpoint; concurrency is per process |
| Phen2Gene | HPO-derived gene candidates | External upload acknowledgement and network access |
| AMELIE | HPO/literature-supported candidates | Explicit acknowledgement and correct TLS trust; no silent insecure fallback |
| PubMed/PMC | Article search, abstracts and supported full texts | Network/rate limits; unavailable papers are not negative evidence |
| ClinGen public pages | Source-bound curation detail | Online fetch/cache may still occur despite local export tables |
| GeneReviews/NCBI Books | Mechanism chapter text | Index-only deployment is not offline text coverage |
| GRCh37 Ensembl REST | Supported independent transcript/sequence checks | Assembly-specific endpoint and exact identity validation |

## Maintainer and user declarations

1. **No blanket raw-download claim:** A/B files still require the correct release
   and sidecars. C/D items need the documented preparation or complete installed
   snapshot. File existence and doctor alone do not prove scientific coverage.
2. **No silent preprocessing:** do not change coordinates, score values, candidate
   scope, missing-value semantics or clinical classifications to make a reader
   accept a file. New data versions are recorded and validated separately.
3. **Resource locations:** configure paths relative to resource_root or use
   explicit absolute paths valid on the destination host. Do not recreate paths
   from another installation.
4. **No double download accounting:** reuse nested/shared Exomiser/LIRICAL data
   and verified shared references; sizes in the register cannot be summed into
   one unique requirement. The current two reference FASTAs are not declared
   interchangeable without validation.
5. **No implied completeness or offline guarantee:** three-gene indexes, scoped
   curation and online evidence fetches are disclosed. The complete variant
   resource closure is not supplied by the small project attachments alone.
6. **Attribution:** retain upstream notices, source references, versions and
   checksums. Follow the data and tool distribution conditions when sharing resources.
