# MendelPRIME

MendelPRIME is an agent system for causal-gene prioritization in Mendelian
disorders. It integrates heterogeneous patient information with specialized
medical tools and knowledge retrieval to compare competing genetic hypotheses.
The system accepts HPO terms alone, HPO terms with clinical descriptions, or
HPO terms with genomic variants.

In the cohorts evaluated in our study, MendelPRIME achieved state-of-the-art
causal-gene prioritization performance among the compared methods in both
phenotype-driven and joint phenotype-variant evaluations. It produces stage-level
and final gene rankings. Optional traceable
reports summarize the case, explain supporting and conflicting evidence, and
compare leading candidates with references to databases and literature. For VCF
inputs, reports also include variant-level ACMG interpretation. Case and Cohort
workflows share the same interface, with resource profiles for ranking alone or
ranking with reports.

[Quick start](#quick-start) | [Demos](examples/README.md) |
[Resources](#resources-and-configuration) | [Documentation](#documentation)

![Figure 1. MendelPRIME architecture, supported inputs and outputs.](docs/assets/figure1.png)

**Figure 1. Overview of MendelPRIME.** Heterogeneous patient inputs (A) are
processed through the Harness, OTEM and CEMR with medical tools and knowledge
resources (B), producing gene priorities, evidence profiles and traceable
sources (C).

## Supported inputs

All three inputs support **Case** and **Cohort** processing through `python main.py`.

| Input | Patient information | Included demo |
|---|---|---|
| HPO | Positive HPO identifiers; explicitly absent terms are optional | [Case](examples/hpo/case.json) / [Cohort](examples/hpo/cohort.jsonl) |
| HPO + clinical description | HPO identifiers and clinical text | [Case](examples/clinical/case.json) / [Cohort](examples/clinical/cohort.jsonl) |
| HPO + VCF | HPO identifiers, VCF, assembly and sample roles | [Case](examples/vcf/case.json) / [Cohort](examples/vcf/cohort.jsonl), with the simulated VCF included |

Clinical text supplements HPO input. VCF supports **trio** and **proband-only**
analysis with **GRCh37/hg19** resources. Sample roles must be specified; column
order is not a pedigree. See [input schemas](docs/INPUTS.md).

The HPO and VCF demos share one published phenotype profile[^phenopackets].
The clinical demo is a separate published case with a short narrative containing
biochemical measurements and ultrasound findings[^clinical-demo]. The VCF uses a public GIAB
background[^giab] with a simulated variant, not clinical sequencing from that
patient. It was selected from a completed evaluation with the causal gene at
Top1. [Demo provenance and reference results](examples/README.md#data-and-reference-result)
are separate from inference inputs.

## Tools and knowledge

| Mode | Tools and resources |
|---|---|
| Shared | HPO[^hpo], local OMIM[^omim], PubMed[^pubmed]/PMC[^pmc], and a configured LLM |
| HPO / HPO+clinical | Built-in HPOOverlap and Phen2Disease[^phen2disease]; online Phen2Gene[^phen2gene] and AMELIE[^amelie] |
| VCF ranking | Local Exomiser[^exomiser], LIRICAL[^lirical] and AutoPVS1[^autopvs1]; ClinVar[^clinvar], ClinGen[^clingen], reference FASTA, gnomAD[^gnomad], SpliceAI[^spliceai], REVEL[^revel] and AI-CURA[^aicura] ranking knowledge |
| VCF report | Ranking resources plus VEP[^vep], MANE[^mane], full AI-CURA guidance and ClinGen SVI; optional GeneReviews[^genereviews] and PLINK2[^plink2] enrichment |
| Optional RareA retrieval | RareA[^rarea] records and vectors, Qwen3-Embedding[^qwen] and MedCPT[^medcpt] |

Exomiser and LIRICAL run through Java. HPO/clinical input needs no Java or variant
databases. Reports are optional; RareA and Reflection are disabled in the minimal
profiles. LLM, literature and configured evidence services require network access.

## Quick start

### Install

Use **Linux x86-64 and Python 3.12**. Set `CODE` to the extracted source directory;
keep private configuration, inputs and results outside it.

```bash
export CODE=/absolute/path/to/MendelPRIME
export WORK="$HOME/mendelprime-work"
test ! -e "$WORK" || exit 1
umask 077
mkdir -p "$WORK"
python3.12 -m venv "$WORK/.venv"
source "$WORK/.venv/bin/activate"
cd "$CODE"
python -m pip install -c constraints/runtime.txt .
python -m pip check
```

For VCF, replace the install command with `python -m pip install -c
constraints/runtime.txt '.[vcf]'`. Native tools and data are separate. See
[installation](docs/INSTALL.md) for Conda, compiler and native-library requirements.
Keep the source directory when installing a wheel: it contains configurations,
demos and preparation scripts that the wheel does not install.

### Prepare resources and configuration

The smallest HPO installation needs the **HPO resource package**, a complete
local **OMIM mirror + mim2gene**, and the **Phen2Disease dataset**. The one-record
OMIM example only illustrates its schema; it is not a replacement database.
Follow the [resource setup table](#resources-and-configuration) before running.

```bash
export CFG="$WORK/configs"
export INPUT="$CODE/examples"
export OUT="$WORK/outputs"
export RES="$WORK/resources"
cp -R "$CODE/configs" "$CFG"
mkdir -p "$OUT" "$RES" "$WORK/cache"
```

Edit the copied files, not the repository templates:

| Setting | Where to edit |
|---|---|
| Resource directory | `minimal/phenotype.resources.toml`: set `resource_root` to the actual absolute directory |
| LLM endpoint and model | `minimal/phenotype.toml`: edit the existing `[llm]` section |
| Phenotype-service consent | Set `[external].phen2gene_acknowledged` and `amelie_acknowledged` to `true` only after authorizing phenotype uploads |

TOML paths do not expand `$RES` or `~`. For VCF, use the matching resource file
and supply actual Java/JAR/data paths; report mode also requires VEP mounts.
Existing compatible resources can be reused in place.

```toml
[llm]
base_url = "https://your-provider.example/v1"
model = "gpt-6.1-sol"
wire_api = "chat_completions"
api_key_env = "MENDELPRIME_API_KEY"
transport_retries = 1
```

Import the key through the environment, never through public files:

```bash
read -r -s -p "LLM API key: " MENDELPRIME_API_KEY
export MENDELPRIME_API_KEY
printf '\n'
```

`.env` files are not loaded automatically. HTTPS verification is enabled by
default; see [network configuration](docs/CONFIGURATION.md#5-network-tools-and-cold-start).

### Run the HPO demo

```bash
python main.py doctor --config "$CFG/minimal/phenotype.toml" --profile hpo
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/case.json" --output-root "$OUT/hpo_case" --cache-root "$WORK/cache/hpo_case"
```

The minimal profiles already set **one case lane and four global LLM slots**.
No extra concurrency flags are required. `doctor` checks setup, not model quality
or API availability. A demo still needs its selected resources and a working API.

## Run a Case or Cohort

`run` is a subcommand of `main.py`, not a second script. A Case is one JSON/YAML
record; a Cohort JSONL has one record per line. Each included Cohort contains the
same single demo as its Case file; append your own unique cases for a real batch.

<details>
<summary>HPO Cohort, clinical Case/Cohort, and VCF Case/Cohort commands</summary>

```bash
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/cohort.jsonl" --output-root "$OUT/hpo_cohort" --cache-root "$WORK/cache/hpo_cohort"

python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/case.json" --output-root "$OUT/clinical_case" --cache-root "$WORK/cache/clinical_case"
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/cohort.jsonl" --output-root "$OUT/clinical_cohort" --cache-root "$WORK/cache/clinical_cohort"

python main.py doctor --config "$CFG/minimal/vcf_rank.toml" --profile vcf --check-tools
python main.py run --config "$CFG/minimal/vcf_rank.toml" --input "$INPUT/vcf/case.json" --output-root "$OUT/vcf_case" --cache-root "$WORK/cache/vcf_case"
python main.py run --config "$CFG/minimal/vcf_rank.toml" --input "$INPUT/vcf/cohort.jsonl" --output-root "$OUT/vcf_cohort" --cache-root "$WORK/cache/vcf_cohort"
```

</details>

| Requested output | Profile and option |
|---|---|
| HPO / clinical ranking | `configs/minimal/phenotype.toml` |
| HPO / clinical with report | Same profile; add `--explanation`, with no extra variant resources |
| VCF ranking | `configs/minimal/vcf_rank.toml` |
| VCF with report | `configs/minimal/vcf_report.toml`; install the additional report resources |

[Complete demo commands](examples/README.md#commands) cover all Case/Cohort and
report combinations. Use a new output directory for each run. VCF paths resolve
relative to the input manifest, so the supplied demo needs no path rewriting.

## Outputs

Each case produces exactly two files:

```text
outputs/hpo_case/demo_hpo/demo_001/
  result.json
  traceable_explanation.md
```

`result.json` contains stage rankings, the final ranking, status and evidence
provenance. `traceable_explanation.md` is empty when reports are off. Reports can
also be generated later from saved results without reranking:

```bash
python main.py explain --config "$CFG/minimal/vcf_report.toml" --input "$OUT/vcf_case" --output-root "$OUT/vcf_explained" --cache-root "$WORK/cache/vcf_explained"
python main.py evaluate --input "$OUT/vcf_case" --truth "$INPUT/truth.csv"
```

Use the phenotype profile for saved HPO/clinical results. Gold labels are used
only by `evaluate`. Existing case outputs are not overwritten. Exit codes are
`0` for success, `2` for degradation and `1` for failure or invalid setup.
[Recovery and report controls](docs/EXECUTION.md) describe partial failures.

## Resources and configuration

Source code, demo data and large resource packages are distinct. Download the six
prepared archives from the [1.0.0 release](https://github.com/Pumpo97/MendelPRIME/releases/tag/1.0.0),
along with [ASSETS.json](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/ASSETS.json)
and [SHA256SUMS](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/SHA256SUMS).
External resources remain subject to their own access and licensing terms.
The six prepared archives are **not** a complete installation. Start with the
[resource-set checklist](resources/README.md#choose-the-smallest-resource-set)
and obtain only the additional files needed for your selected profile.

| Capability | Ready-made packages | Obtain or prepare separately |
|---|---|---|
| HPO / clinical, with or without reports | HPO ontology, statistics and gene annotations | Full OMIM JSONL + `mim2gene`; Phen2Disease four-file dataset; LLM and phenotype-service access |
| VCF ranking | HPO, ClinVar SQLite, AI-CURA ranking guidance, LIRICAL CLI | OMIM; Java 17; Exomiser 14.0.0 and 2406 genome/phenotype data; LIRICAL data including HGNC[^hgnc]; ClinGen validity/VCEP; GRCh37 FASTA/indexes; gnomAD, SpliceAI, REVEL; AutoPVS1 source/data |
| Add VCF reports | Reuse the ranking resources | VEP 104.3/cache/FASTA; MANE; full AI-CURA protocol/workbook; ClinGen SVI guidance |
| Add RareA retrieval | RareA records/IDs/vectors and retrieval support | Matching Qwen weight shards; `.[hpo]` or `.[full]` Python dependencies |

**VCF ranking still uses pre-ranking ACMG when reports are off.** The report
switch controls post-ranking ACMG and narrative generation. Ranking-only does
not load VEP, MANE or the original AI-CURA protocol/workbook unless ranking VEP
is explicitly enabled. It still needs gnomAD, SpliceAI, REVEL and AutoPVS1;
Exomiser summary annotations do not replace these resources.

- [Prepared archives and extraction](docs/RESOURCE_ASSETS.md): HPO, ClinVar, guidance, RareA, retrieval support and LIRICAL.
- [Acquisition and version requirements](docs/RESOURCE_PREPARATION.md) and [preprocessing commands](preparation/README.md): references, predictors, ClinGen and AutoPVS1.
- [OMIM mirror schema](docs/OMIM_MIRROR.md): normalized local JSONL; no OMIM API key is used.
- [LIRICAL installation/build](docs/LIRICAL.md) and [VEP installation/mounts](docs/VEP_RUNTIME.md): full programs, not isolated JARs or executables.
- [Parameter reference](docs/CONFIGURATION.md) and [resource inventory](docs/RESOURCE_REGISTER.md): effective defaults, paths and companion indexes.

## Documentation

| Task | Guide |
|---|---|
| Try a demo or prepare patient inputs | [Demos](examples/README.md), [input schema](docs/INPUTS.md) |
| Install or configure a deployment | [Installation](docs/INSTALL.md), [configuration](docs/CONFIGURATION.md) |
| Obtain tools and databases | [Resource setup](#resources-and-configuration) |
| Generate reports, evaluate or recover failures | [Execution](docs/EXECUTION.md) |

## Resources and references

We thank the contributors of the tools and databases below, including RareA and
the original case-report authors. Cite the resources used in your analysis and
follow their version-specific licenses. MendelPRIME also uses Python[^python],
Java[^java], optional udocker[^udocker] and the libraries listed in `pyproject.toml`.

[^hpo]: Human Phenotype Ontology. [HPO project](https://hpo.jax.org/).
[^omim]: Online Mendelian Inheritance in Man. [OMIM](https://www.omim.org/).
[^pubmed]: NCBI PubMed. [Literature database](https://pubmed.ncbi.nlm.nih.gov/).
[^pmc]: PubMed Central. [Full-text archive](https://pmc.ncbi.nlm.nih.gov/).
[^phen2disease]: Phen2Disease. [Official repository](https://github.com/ZhuLab-Fudan/Phen2Disease).
[^phen2gene]: Phen2Gene. [Official repository](https://github.com/WGLab/Phen2Gene).
[^amelie]: AMELIE. [Service](https://amelie.stanford.edu/) and [publication](https://pubmed.ncbi.nlm.nih.gov/32434849/).
[^exomiser]: Exomiser. [Official repository](https://github.com/exomiser/Exomiser).
[^lirical]: LIRICAL. [Official repository](https://github.com/TheJacksonLaboratory/LIRICAL).
[^autopvs1]: AutoPVS1. [Official repository](https://github.com/JiguangPeng/autopvs1).
[^clinvar]: ClinVar and its submitters. [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
[^clingen]: ClinGen, its expert panels and SVI working group. [ClinGen](https://clinicalgenome.org/).
[^gnomad]: Genome Aggregation Database and contributing studies. [gnomAD](https://gnomad.broadinstitute.org/).
[^spliceai]: SpliceAI. [Official repository](https://github.com/Illumina/SpliceAI).
[^revel]: REVEL. [Official resource](https://sites.google.com/site/revelgenomics/).
[^aicura]: AI-CURA. [Paper and supplementary material](https://doi.org/10.1126/scitranslmed.adz4172).
[^vep]: Ensembl Variant Effect Predictor. [Documentation](https://www.ensembl.org/info/docs/tools/vep/index.html).
[^mane]: MANE, NCBI RefSeq and EMBL-EBI. [Matched transcript resource](https://www.ncbi.nlm.nih.gov/refseq/MANE/).
[^qwen]: Qwen3-Embedding-4B. [Model repository](https://huggingface.co/Qwen/Qwen3-Embedding-4B).
[^medcpt]: MedCPT Cross-Encoder. [NCBI model repository](https://huggingface.co/ncbi/MedCPT-Cross-Encoder).
[^rarea]: RareA similar-case corpus. [Prepared records and vector index](docs/RESOURCE_ASSETS.md#asset-contents); original references are retained in the records.
[^plink2]: PLINK 2. [Documentation](https://www.cog-genomics.org/plink/2.0/).
[^genereviews]: GeneReviews. [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK1116/).
[^hgnc]: HUGO Gene Nomenclature Committee. [HGNC](https://www.genenames.org/).
[^python]: Python. [Python 3.12 documentation](https://docs.python.org/3.12/).
[^java]: OpenJDK. [JDK 17 project](https://openjdk.org/projects/jdk/17/).
[^udocker]: udocker, an optional VEP container runtime. [Official repository](https://github.com/indigo-dc/udocker).
[^phenopackets]: Robinson and contributors. [384 Phenopackets, version 1.0](https://doi.org/10.5281/zenodo.3905420), CC-BY-4.0.
[^clinical-demo]: Vajro et al. (2018), Patient 3. [Published clinical and laboratory observations](https://doi.org/10.1186/s13023-017-0757-3), CC-BY-4.0.
[^giab]: Genome in a Bottle. [NIST project](https://www.nist.gov/programs-projects/genome-bottle); [public NIST VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf).
