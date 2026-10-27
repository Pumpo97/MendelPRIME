# Install the prepared resource assets

The resource packages include a **frozen HPO snapshot, prepared ClinVar
SQLite and a prebuilt RareA vector index**. OMIM remains local. No OMIM API key, new ClinVar
server, or user-side XML conversion is required by this release.

The archives are data-release attachments, separate from the source repository.
Download them from the [1.0.0 release](https://github.com/Pumpo97/MendelPRIME/releases/tag/1.0.0),
with its [resource manifest](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/ASSETS.json)
and [SHA256SUMS](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/SHA256SUMS).
The code archive alone is not a complete resource bundle.

## Asset contents

| Archive | Purpose | Paths after extraction into the resource root |
|---|---|---|
| [mendelprime-hpo-20260728-v1-public.tar.gz](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/mendelprime-hpo-20260728-v1-public.tar.gz) | Preserved ontology/annotations/statistic values, with public metadata | `core/hpo/hp.obo`; `core/hpo/hpo_stats.json`; `core/hpo/genes_to_phenotype.txt` |
| [mendelprime-clinvar-20260117-public.tar.gz](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/mendelprime-clinvar-20260117-public.tar.gz) | Preserved clinical data with sanitized build metadata | `variant-knowledge/clinvar/clinvar_20260117.sqlite`; `assessment/grch37/clinvar-protein/clinvar_protein_acmg_disease_20260117.sqlite` |
| [mendelprime-vcf-ranking-guidance-v1.tar.gz](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/mendelprime-vcf-ranking-guidance-v1.tar.gz) | Prepared AI-CURA knowledge for ranking, about 17 KB compressed | `guidelines/ai-cura/{ai_cura_knowledgebase.json,resource_manifest.json}` |
| [mendelprime-similar-cases-rarea-20261027-public.tar.gz](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/mendelprime-similar-cases-rarea-20261027-public.tar.gz) | Unchanged RareA corpus, identifier order and document embeddings | `phenotype/similar-cases/index/rarea/{records.jsonl,record_ids.json,embeddings.npy}` |
| [mendelprime-retrieval-support-20260815-public.tar.gz](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/mendelprime-retrieval-support-20260815-public.tar.gz) | Unchanged embedding sidecars and deployed safetensors MedCPT reranker | `phenotype/models/qwen3-embedding-4b/`; `phenotype/models/medcpt-cross-encoder/` |
| [lirical-cli-2.4.1-distribution.zip](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/lirical-cli-2.4.1-distribution.zip) | Complete upstream tool distribution, including notices | Install as described in [LIRICAL.md](LIRICAL.md); its data is separate |

Each project data tarball contains a versioned manifest and a `.sha256` payload
list in `manifests/`; LIRICAL keeps its upstream ZIP layout. Paths are relative
to the resource root.
The two SQLite files are distinct and both are needed for full VCF assessment.
They are derived lookup databases, not substitutes for Exomiser/LIRICAL's native
`2406_hg19_clinvar.mv.db`. Use the matching archive and payload checksums; verify
files before configuring them for inference.

The HPO archive preserves the project's `formal_20260728_v1` snapshot. Its OBO
declares `hp/releases/2025-03-03`; the bundle date is not an upstream HPO release
date. The annotation header does not declare a separate release, so its exact
bytes are versioned by SHA256. See [HPO statistics](RESOURCE_REGISTER.md#hpo-statistics) for the
count/IC consistency checks and the preserved historical annotation boundary.

The vector index contains **5,867 RareA records**, with **2,560
dimensions**. Identifiers, record ordering and vector rows are inseparable.
Do not replace `records.jsonl` with a reordered spreadsheet or regenerate only
one file. The records remain necessary to interpret a retrieved vector; a lone
`embeddings.npy` does not provide the original evidence or gene associations.
The reader loads only `rarea/`; other corpus folders are not consumed, and a
record declaring a different provider is rejected. The extraction target for a
new installation should contain only the documented RareA payload, not an older
multi-corpus index. A change in retrieval corpus can change downstream rankings;
preserve the documented corpus when reproducing a configuration.

## Verify and extract

Place the archives and `SHA256SUMS` in one download directory. Install into a new
resource root, or back up an existing installation before deliberately replacing
any resource. The examples assume no conflicting files already exist.

```bash
cd /absolute/path/to/downloads
sha256sum -c SHA256SUMS
export RES=/absolute/path/to/MendelPRIME-release-resources
mkdir -p "$RES"
tar -xzf mendelprime-hpo-20260728-v1-public.tar.gz -C "$RES"
# Add ClinVar only for VCF:
tar -xzf mendelprime-clinvar-20260117-public.tar.gz -C "$RES"
tar -xzf mendelprime-vcf-ranking-guidance-v1.tar.gz -C "$RES"
# Add these two only when enabling RareA retrieval:
tar -xzf mendelprime-similar-cases-rarea-20261027-public.tar.gz -C "$RES"
tar -xzf mendelprime-retrieval-support-20260815-public.tar.gz -C "$RES"
cd "$RES"
sha256sum -c manifests/mendelprime-hpo-20260728-v1.sha256
sha256sum -c manifests/mendelprime-clinvar-20260117.sha256
sha256sum -c manifests/mendelprime-vcf-ranking-guidance-v1.sha256
sha256sum -c manifests/mendelprime-similar-cases-rarea-20261027.sha256
sha256sum -c manifests/mendelprime-retrieval-support-20260815.sha256
```

The archives have no enclosing resource-root directory. The `-C` target above
is exactly the `resource_root` to put in each selected `resources.toml`.
The extracted files are self-contained; no external symlinks are required.

Only run extraction/payload-check commands for the assets selected by your
profile. `sha256sum -c SHA256SUMS` assumes the whole download set is present.
For a single downloaded asset, without requiring newer sha256sum options:

```bash
grep '  mendelprime-hpo-20260728-v1-public.tar.gz$' SHA256SUMS | sha256sum -c -
```

VCF ranking uses the two-file guidance package without the original AI-CURA XLSX,
protocol text, MANE or ClinGen SVI files. Its JSON is unchanged from the full
capsule and is verified against the pinned manifest. For later report generation,
add the original workbook and protocol using the [full capsule recipe](../preparation/README.md#ai-cura-capsule),
plus MANE, ClinGen SVI and VEP. Use a new preparation root if a destination already
exists; preparation utilities do not silently overwrite installed resources.

## Complete the query model

Precomputed document vectors avoid rebuilding the corpus, but a new patient's
query still needs the compatible embedding model and the cross-encoder. The
support archive includes the complete deployed MedCPT cross-encoder and small
Qwen sidecars; it deliberately excludes the approximately 8 GB embedding weights.

After installing the `hpo` or `full` Python group, download the two official
Qwen weight shards at the revision associated with this index:

```bash
export RES=/absolute/path/to/MendelPRIME-release-resources
python - <<'PY'
import os
from pathlib import Path
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="Qwen/Qwen3-Embedding-4B",
    revision="5cf2132abc99cad020ac570b19d031efec650f2b",
    allow_patterns=["model-*.safetensors"],
    local_dir=str(Path(os.environ["RES"]) / "phenotype/models/qwen3-embedding-4b"),
)
PY
cd "$RES"
sha256sum -c manifests/mendelprime-retrieval-models.sha256
```

Alternatively copy matching model files from an existing installation and run
the same checksum command. `model_compatibility.json` in the index also records
the verified runtime-file hashes. Do not substitute another embedding model
just because its output has 2,560 dimensions. The deployed MedCPT copy uses
safetensors; its upstream repository also offers a PyTorch `.bin` distribution,
which is not byte-identical to this prepared asset.

Upstream models: [Qwen3-Embedding-4B](https://huggingface.co/Qwen/Qwen3-Embedding-4B)
and [MedCPT-Cross-Encoder](https://huggingface.co/ncbi/MedCPT-Cross-Encoder).
Their own terms apply separately to the source code and corpus data release.
An index does not make underlying clinical text anonymous or grant redistribution
rights. Before public upload, check source-corpus terms and privacy, and retain
record groups/PMIDs and the evaluation exclusion policy. The prepared archives
do not replace the source-data permissions required for redistribution.

## Configure and check

In both HPO and VCF resource manifests, set `resource_root` to the absolute path
used above. Existing relative paths already match the archives. No new schema
or ranking parameter is needed. For HPO/clinical retrieval, keep:

```toml
[ranking]
similar_cases = true
```

Set `execution.device` to the tested GPU environment or `cpu`; CPU still requires
the models and can be much slower. Setting `similar_cases = false` omits the
entire branch and its model/index requirements, rather than running only part of
the original retrieval method. HPO statistics and the ontology are still needed.

```bash
python main.py doctor --config "$CFG/hpo/run.toml" --profile hpo
python main.py doctor --config "$CFG/minimal/vcf_report.toml" --profile vcf --check-tools
```

The first command deliberately checks the expanded RareA profile described in
this section. Configure its private resource manifest and device first. If you
are not enabling retrieval, use `$CFG/minimal/phenotype.toml` instead; do not
install neural models only to satisfy the wrong profile's doctor checks.

Run these from the MendelPRIME source directory after installing the remaining
resources listed in [the configuration matrix](CONFIGURATION.md#4-resource-parameters-by-profile). The prepared archives do **not**
include a complete OMIM mirror, Phen2Disease, native-tool databases, gnomAD, SpliceAI,
REVEL, MANE, VEP or the complete report guidance bundle. OMIM keeps the existing
`resources.omim_jsonl` and `resources.mim2gene` paths. Install LIRICAL from its
[pinned official source](LIRICAL.md), separately from its scientific databases.
The ranking-guidance archive supplies only the derived JSON and pinned manifest,
not the original report-assessment guidance. The source package includes a one-record [OMIM mirror example](OMIM_MIRROR.md)
for reader/format testing; it is not a substitute for the full local mirror.

Keep versioned resource packages immutable during a cohort run. Verify hashes
and test the actual readers before accepting a relocated bundle. An older
filesystem timestamp alone does not prove an index mismatch; neither does a
successful file copy prove the original pair was scientifically correct. Never
use `touch` to conceal an uninvestigated index warning.
