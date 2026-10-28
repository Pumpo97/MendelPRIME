# Execution

Install the appropriate dependencies and configure the resources first. Commands
below run from the source checkout; the installed `mendelprime` command invokes
the same `mendelprime.runtime.cli.main` implementation as `python main.py`.

Standard Case/Cohort commands are in the [README](../README.md#run-a-case-or-cohort).
Use its CODE/CFG/INPUT/OUT/WORK directory setup. All execution controls and
mode-specific dependencies are defined in [CONFIGURATION.md](CONFIGURATION.md).
This page covers private work directories, saved reports, evaluation and recovery.

## Private work and outputs

Use a fresh `--cache-root` for an application-level cold start. Reference data and
model weights must already be installed; this does not flush OS/provider caches.
Let MendelPRIME create the cache directory. If creating it yourself on Linux,
use `mkdir -m 700 CACHE` under an existing parent: caches and audits must belong
to the current user and must not be accessible by other users. An existing cache
with broader permissions is rejected rather than silently changing permissions.
Logs show case, ranking-stage and model-call progress. Public cases contain only:

```text
OUTPUT/cohort/case_id/result.json
OUTPUT/cohort/case_id/traceable_explanation.md
```

The Markdown file is empty when explanation is disabled. `result.json` contains
stage rankings, the final ranking, status, and versioned evidence provenance.
VCF evidence limits are separate from execution/contract failures. Successful
report generation does not erase a pre-existing ranking-stage degradation.

`--no-explanation` leaves the report empty and skips post-ranking ACMG; it does
not mean that degraded cases are forbidden from having reports. With reporting
enabled, completed variant assessments can be retained alongside an explicit
unclassified-variant list and a partial report. A failed assessment is never
silently replaced with a successful classification. Independent `explain`
preserves the original rankings and does not clear their historical failures.

## Independent reports and evaluation

```bash
python main.py explain --config "$CFG/vcf-report/run.toml" --input "$OUT/vcf_case" --output-root "$OUT/explained" --cache-root "$WORK/cache/explained"
python main.py explain --config "$CFG/vcf-report/run.toml" --input "$OUT/vcf_case" --output-root "$OUT/report_rewrite" --cache-root "$WORK/cache/report_rewrite" --rerender-report
python main.py evaluate --input "$OUT/vcf_cohort" --truth "$WORK/truth.csv"
```

`explain` accepts a saved case or result collection. By default it preserves an
existing report; `--rerender-report` regenerates its narrative while retaining
completed ACMG decisions. `--recompute-assessment` explicitly regenerates
post-ranking evidence/ACMG and the report. Neither option reruns gene ranking.
Legacy results without sufficient provenance cannot fabricate missing evidence.

For a transient post-ranking evidence-service failure in a VCF case, preserve the
first output and use the original canonical input to enable fresh family QC:

```bash
python main.py explain --config "$CFG/vcf-report/run.toml" --input "$OUT/first/cohort/case/result.json" --case-inputs "$WORK/inputs/vcf/case.json" --output-root "$OUT/post_recovery" --cache-root "$WORK/cache/post_recovery" --recompute-assessment
```

This rebuilds post-ranking ACMG/report evidence; it does not rerun Direct, tools,
Fusion, Admission, Global, Pairwise or Tail. Verify the new execution status and
compare all ranking stages before selecting the recovered package. Retain the
original degraded result; do not remove failure labels by editing its JSON.

Truth CSV fields are `cohort,case_id,gold_genes`; multiple genes may be
semicolon-separated or a JSON array. Evaluation emits per-case ranks and
cohort/overall metrics to stdout, and never injects Gold into inference.

## Failures and reruns

Exit codes are `0` for all cases ok, `2` if one or more are degraded, and `1` for
failed cases or invalid input/configuration. The final stdout JSON lists each
case's status and output path. A missing or failed resource is not negative
biological evidence. Source/allele/quotation checks are not bypassed on a retry.

If explanation-context validation fails after ranking, the completed rankings
are retained in a degraded two-file result, with an empty report and explicit
`explanation_context_invalid` diagnostics. Unvalidated provenance is not exported:
`provenance` is null and standalone `explain` refuses it until the evidence
context has been recovered and validated. This does not label failed evidence
or an unavailable report as successful.

Native regional candidates and transcript-specific annotations have different
scopes. When the same allele, transcript accession and coding HGVS identify a
different gene in ClinVar, the LIRICAL candidate and raw record are retained;
the foreign coding annotation is separated into `gene_scope` provenance rather
than assigned to that candidate. Other-gene ClinVar records remain in
`out_of_scope_records` and cannot supply its assertions or VCEP joins. Ambiguous
identities are not guessed. Database/query failures and invalid source bindings
remain failures; successful installation alone does not verify those bindings.

`report_passage_source_mismatch` means a selected passage was not offered under
the cited source. The existing single repair receives the offending claim,
citation and actual passage owners. Missing allele identity bindings are checked
in the same pass and included in that one repair, with scoped identity-source
options. These options are not automatic citations or pathogenicity evidence;
gene/allele eligibility, medical support and literal quotations still have to pass.
If that repair attempt also fails, report composition
does not run, even when ranking and ACMG have completed. This is an evidence
binding failure, not an installation failure. Preserve the failed result and
audit; an independent `explain` run can reuse completed assessments without
reranking, but must still satisfy the same citation checks.

Use a new output root for a rerun. To preserve completed work, select only failed
or unfinished cases in a new manifest; there is no automatic overwrite or
checkpoint-resume flag. Report-only failures can use `explain` on existing
results. On interruption, dispatch stops; active work can take time to release
resources. Run one foreground process first before scheduling a larger cohort.
