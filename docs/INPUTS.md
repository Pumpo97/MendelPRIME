# Patient inputs

All three input types use `python main.py run`. A Case file contains one record;
a Cohort file contains multiple records. They enter the same Runner and Harness.
The [runnable demos](../examples/README.md) include a shared HPO/simulated-VCF
example and a separate clinical case with laboratory values and imaging findings.
The VCF is a complete simulated spike-in dataset, not the patient's sequencing data.

## Shared fields

| Field | Meaning |
|---|---|
| `schema_version` | `mendelprime_case_v1` |
| `cohort` | Output group; defaults to `cases` |
| `case_id` | Stable identifier, unique within its cohort |
| `phenotypes.present` | Positive HPO IDs; at least one positive term is required for inference |
| `phenotypes.absent` | Explicitly negative HPO IDs, not unreported features |
| `clinical_description` | Optional verbatim clinical text, kept separate from HPO polarity |
| `genomic` | Optional VCF metadata; its presence selects the genotype branch |

HPO identifiers are sufficient; labels are optional. An alternative phenotype
representation is a list of objects with `hpo_id`, optional `label`, and boolean
`present`. Do not provide conflicting labels or positive/negative states for the
same identifier. Missing clinical findings are not negative observations.

Do not place Gold genes, expected rankings, or reference variant classifications
in the patient record. Gold data belongs only to the independent `evaluate`
command. Unknown input fields are rejected rather than silently ignored.

## HPO and clinical inputs

`examples/hpo/case.json` contains HPO only. The corresponding clinical example
adds `clinical_description`; it still requires HPO input. Both use the phenotype
ranking branch. The clinical description is evidence, not a replacement schema
or a request to change the ranking instructions. Remove personal identifiers
from free text before sending it to a configured external LLM service.

## VCF inputs

The `genomic` object includes:

| Field | Meaning |
|---|---|
| `vcf_path` | A VCF or supported compressed VCF path |
| `assembly` | Declare `GRCh37` explicitly for the currently validated resource profile; `hg19` is its alias |
| `family_mode` | `trio`, `proband_only`, or the supported `auto` resolver |
| `sample_roles` | Exact VCF sample IDs under `proband`, `father`, and `mother` |
| `pedigree_path` | Optional PED input, when available |
| `sex` | `MALE`, `FEMALE`, or `UNKNOWN` |
| `assay_type` | Such as `WES`, `WGS`, or `PANEL` |

For a trio, explicitly map the three roles. For a singleton, supply only the
proband and select `proband_only`. The same mode can deliberately analyze only
the proband from a family VCF. A PED file is not required when the explicit role
mapping provides the needed identities; runtime preparation handles tool input
generation. Do not infer family relationships from an arbitrary VCF column order.

Use either `sample_roles` or the legacy `proband_sample_id`, not both. When a PED
is also supplied, its relationships must agree with the metadata. State the
genome build instead of relying on incomplete VCF header hints. GRCh38 parsing
support alone does not establish that a complete GRCh38 resource set has been
installed or that its full ranking/report chain has been validated.

Mitochondrial reference identity is checked from the complete FASTA sequence,
not inferred from `chrM` or its length. The rCRS sequence (`NC_012920.1`) supports
reference-bound annotation. The older UCSC hg19 `chrM` (`NC_001807.4`) is retained
without automatic coordinate conversion or an rCRS ClinVar join; unavailable
annotations remain explicit evidence limitations. Unknown reference sequences
and REF-allele mismatches remain failures. Do not replace the mitochondrial
FASTA sequence without correspondingly validating the VCF coordinates. See the
[UCSC reference notes](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/).

The VCF demo includes `demo.vcf.gz` and its tabix index beside the input manifests.
It contains one simulated proband, not a trio. Use actual sample roles for your
own family data; do not add fictitious parents to the demo. Provenance and
evaluation-only labels are documented separately in [examples](../examples/README.md).

## Case and Cohort layouts

A single `.json`, `.yaml`, or `.yml` can contain one record or an array. A `.jsonl`
contains one record per nonempty line and is recommended for cohort manifests:

```text
inputs/
  cohort.jsonl
  case_001/family.vcf
  case_002/proband.vcf
```

Paths in a manifest resolve relative to that manifest, not the shell's current
directory. Alternatively use a directory of per-case `case.yaml` files:

```text
inputs/
  case_001/case.yaml
  case_001/family.vcf
  case_002/case.yaml
  case_002/proband.vcf
```

Directory discovery searches for `case.yaml`. If none exist, it tries the root
`cohort.jsonl`. Do not mix both conventions and expect both sets to run. Input
manifests must be regular files, not symlinks. Duplicate `(cohort, case_id)` keys
and case-insensitive filename collisions are rejected before inference.

Malformed input is a preflight error. A failure during an individual case is
recorded for that case while the remaining cohort jobs continue. Existing case
output directories are never overwritten or implicitly resumed. To retry, use
an input manifest containing only the intended cases and a new output root.
