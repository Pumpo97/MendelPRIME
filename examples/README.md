# Runnable demos

Complete the [installation and resource setup](../README.md#quick-start) first.
These examples supply patient inputs, not the external tools, databases or API
credentials. Run from the source directory with the `CFG`, `INPUT`, `OUT` and
`WORK` variables defined in that setup.

## Files

| Path | Contents |
|---|---|
| `hpo/case.json` | Five positive HPO terms |
| `clinical/case.json` | A separate published case: six positive HPO terms and a concise clinical description with laboratory and ultrasound findings |
| `vcf/case.json` | The same terms plus a GRCh37 simulated singleton VCF |
| Each mode's `cohort.jsonl` | One-record Cohort equivalent of its Case; add unique records for a larger batch |
| `vcf/demo.vcf.gz` and `.tbi` | Full BGZF VCF and tabix index, approximately 14.7 MB together |
| `truth.csv` | Evaluation-only Gold labels for the three input modes |
| `provenance.json` | Data attribution, transformations, checksums and recorded reference result |

Keep each VCF manifest beside `demo.vcf.gz` and its index. Its relative path and
explicit proband ID are already correct. The VCF contains one sample, so use
`proband_only`; it cannot demonstrate family segregation or trio inference.
For your own trio, follow the [sample-role schema](../docs/INPUTS.md#vcf-inputs).

## Commands

Use one command for the desired mode. Every output/cache name below is distinct.

### HPO

```bash
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/case.json" --output-root "$OUT/hpo_case" --cache-root "$WORK/cache/hpo_case"
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/cohort.jsonl" --output-root "$OUT/hpo_cohort" --cache-root "$WORK/cache/hpo_cohort"
```

With reports, reuse the same resources:

```bash
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/case.json" --output-root "$OUT/hpo_case_report" --cache-root "$WORK/cache/hpo_case_report" --explanation
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/hpo/cohort.jsonl" --output-root "$OUT/hpo_cohort_report" --cache-root "$WORK/cache/hpo_cohort_report" --explanation
```

### HPO plus clinical description

```bash
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/case.json" --output-root "$OUT/clinical_case" --cache-root "$WORK/cache/clinical_case"
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/cohort.jsonl" --output-root "$OUT/clinical_cohort" --cache-root "$WORK/cache/clinical_cohort"
```

With reports:

```bash
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/case.json" --output-root "$OUT/clinical_case_report" --cache-root "$WORK/cache/clinical_case_report" --explanation
python main.py run --config "$CFG/minimal/phenotype.toml" --input "$INPUT/clinical/cohort.jsonl" --output-root "$OUT/clinical_cohort_report" --cache-root "$WORK/cache/clinical_cohort_report" --explanation
```

### HPO plus VCF

Install the VCF dependencies, native tools and ranking resources first. Ranking
uses the smaller resource profile and does not need VEP or report-only guidance:

```bash
python main.py doctor --config "$CFG/minimal/vcf_rank.toml" --profile vcf --check-tools
python main.py run --config "$CFG/minimal/vcf_rank.toml" --input "$INPUT/vcf/case.json" --output-root "$OUT/vcf_case" --cache-root "$WORK/cache/vcf_case"
python main.py run --config "$CFG/minimal/vcf_rank.toml" --input "$INPUT/vcf/cohort.jsonl" --output-root "$OUT/vcf_cohort" --cache-root "$WORK/cache/vcf_cohort"
```

For reports, install the additional report resources and configure their paths:

```bash
python main.py doctor --config "$CFG/minimal/vcf_report.toml" --profile vcf --check-tools
python main.py run --config "$CFG/minimal/vcf_report.toml" --input "$INPUT/vcf/case.json" --output-root "$OUT/vcf_case_report" --cache-root "$WORK/cache/vcf_case_report"
python main.py run --config "$CFG/minimal/vcf_report.toml" --input "$INPUT/vcf/cohort.jsonl" --output-root "$OUT/vcf_cohort_report" --cache-root "$WORK/cache/vcf_cohort_report"
```

The profiles set one case lane and four total LLM slots, with Reflection and
RareA disabled. Changing the manifest does not change these settings. Within
each input mode, Case and Cohort contain the same example, not independent
benchmark observations. The clinical example differs from the HPO/VCF example.

## Outputs and evaluation

For the VCF Case command:

```text
outputs/vcf_case/demo_vcf/demo_001/
  result.json
  traceable_explanation.md
```

The report file is empty for ranking-only runs. Evaluate a completed run using
the separate truth file:

```bash
python main.py evaluate --input "$OUT/vcf_case" --truth "$INPUT/truth.csv"
```

Replace the input directory with `hpo_case` or `clinical_case` to evaluate those
runs. Only cases present in that result directory are scored; extra truth rows
do not enlarge the denominator. Do not pass `truth.csv` or `provenance.json` as
patient input. To compare separate runs, evaluate each output root separately.

## Data and reference result

### HPO and simulated VCF

The five observations are macrocephaly (`HP:0000256`), cerebral atrophy
(`HP:0002059`), subdural hemorrhage (`HP:0100309`), glutaric aciduria
(`HP:0003150`) and dystonia (`HP:0001332`). They come from Patient 5 in
[Abdul Wahab et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27672653/), as encoded
in [384 Phenopackets, version 1.0](https://doi.org/10.5281/zenodo.3905420)
by Robinson and contributors (CC-BY-4.0).

The simulated VCF pairs those phenotypes with sample `NIST7035` from the public
[GIAB NIST exome VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf).
It retains all **416,690 records** from the evaluated simulation, including the
homozygous spike-in at `GRCh37:19:13007113:G:A`. It is not a VCF from the patient
in the paper, nor a small file containing only the causal allele. No variants
were removed. The spike-in has `GT=1/1` but no
DP/GQ measurements; missing quality and parental evidence have not been invented.

The demo only removes machine-specific command headers and replaces the
machine-specific reference path with `GRCh37`. Variant records and genotypes
are unchanged, and the tabix index is regenerated. Upstream data attribution
remains separate from the license of MendelPRIME's own code.

| Recorded simulated-VCF reference | Gold gene | Rank |
|---|---|---:|
| Exomiser | GCDH | 1 |
| LIRICAL | GCDH | 1 |
| Final gene ranking | GCDH | 1 |

This reference comes from the completed EXO0001 evaluation, which used
GPT-5.5 and also included an AMELIE VCF branch. It is **not** a fresh execution
of the current minimal profile. The current VCF profile uses Exomiser and
LIRICAL; models, resources and live retrieval can change the final ranking.
No Top1 result is asserted for the HPO-only demo. Gene labels and reference
results are kept outside all six inference manifests.

### Clinical description

The clinical example is Patient 3 from [Vajro et al. (2018)](https://doi.org/10.1186/s13023-017-0757-3),
PMID 29321044. Its 61-word clinical passage covers age and course, fluctuating
transaminases, elevated ALP and CK, ultrasound steatosis, antithrombin-III
activity of 68%, and ceruloplasmin of 8 mg/dL (reference range 20-46).
The paired six-term HPO profile is preserved. Only line wrapping and whitespace
were normalized; no measurements or clinical findings were added.

The archived HPO-plus-clinical result ranked **TMEM199 first**, with status
`ok`. This is a historical reference, not a new evaluation of this release.
The causal gene, article identifiers, molecular findings and final syndrome name
are excluded from inference inputs. This case uses `demo_clinical/demo_002`;
its evaluation label is in `truth.csv`, separate from the inputs.

The clinical source is available under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
We credit Vajro and colleagues; provenance records the source and formatting changes.
