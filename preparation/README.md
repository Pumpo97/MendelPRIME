# Portable resource preparation

These source-distribution utilities prepare files **before** MendelPRIME runs.
They are not an alternative ranking pipeline and are not installed as another
console entry. Run them from the source checkout with Python 3.12 and the VCF
dependency group installed:

```bash
python -m pip install -c constraints/runtime.txt '.[vcf]'
python preparation/prepare_resources.py --help
export RAW=/absolute/path/to/downloaded-originals
export RES=/absolute/path/to/new-resource-root
```

The supplied HPO/statistics and RareA packages should be extracted unchanged,
not rebuilt by these scripts. No LLM, OMIM API, automatic web scraping or runtime
scientific-rule modification is involved. Large input/output checksums and
sorting can take time; progress stages are printed to stderr. The final stdout
line is a JSON result with a receipt path.

## Safety and dependencies

- Outputs are staged under RES and validated before publication. Existing files
  are never overwritten; use a new root or remove only your own discarded output
  after inspection. Failed preparation does not replace a working resource.
- Input fields/values are preserved; missing scores are not replaced with zero.
  Indexing requires sorted coordinates and the correct assembly. A declared
  GRCh37 label does not prove an arbitrary FASTA is GRCh37.
- Receipts record input/output checksums and sizes. Reference paths inside a
  reviewed catalog must stay within its directory. Scripts are relocatable;
  none import code from an old server project.
- GNU `sort` is used for bounded-memory stable REVEL sorting. Poppler
  `pdftotext` is needed for PDF text extraction unless supplying matching
  extracted SVI text. Indexing uses the installed `pysam`; separate bgzip/tabix
  executables are not required by these scripts.
- Work space must accommodate the new resource and temporary sorted text. These
  utilities do not fetch or approve licensed source data. Obtain originals from
  the official locations below and preserve their distribution terms.

## Obtain the original inputs

| Resource | Original acquisition | Important version boundary |
|---|---|---|
| ClinVar | [NCBI tab-delimited releases](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/): variant_summary.txt.gz and var_citations.txt | Match the two snapshots and declare their actual release, not a convenient old date |
| REVEL | [Official downloads](https://sites.google.com/site/revelgenomics/downloads); unpack the full CSV or CSV.gz | Requires the nine-column version including GRCh38 positions and transcript IDs; coordinates used here remain hg19 |
| gnomAD sites/coverage | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) | Sites 2.1.1 exomes GRCh37; coverage 2.1 exomes, not a GRCh38 or ANNOVAR summary |
| SpliceAI | [Official precomputed resources](https://github.com/Illumina/SpliceAI) | v1.3 masked hg19 SNV and indel VCFs; preserve access/usage conditions |
| Reference FASTA | The original provider of your verified GRCh37/hg19 reference, such as [Ensembl GRCh37](https://ftp.ensembl.org/pub/grch37/) | Use the selected snapshot and verify contigs/REF alleles; do not interchange the two configured FASTAs by name |
| ClinGen exports | [Validity CSV](https://search.clinicalgenome.org/kb/gene-validity/download), [VCEP export](https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download?type=csv) | Live URLs can change; preserve original files and the actual release date |
| AI-CURA workbook | Supplementary materials of [DOI 10.1126/scitranslmed.adz4172](https://doi.org/10.1126/scitranslmed.adz4172) | Original S1-S13 XLSX, SHA256 verified against the supplied immutable recipe |
| ClinGen SVI PDFs | Exact URLs in assets/clingen-svi/manifest.json | Original PDF and expected extracted-text hashes are pinned |
| OMIM | Your authorized local normalized mirror and the corresponding mim2gene mapping | No OMIM scraping/API command is supplied; use the documented mirror schema |
| GeneReviews | Official gene/chapter mapping from [NCBI Bookshelf GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/) | Header NBK_id, GR_shortname, genesymbol; explicit gene selection controls coverage |
| AutoPVS1 | [Official source](https://github.com/JiguangPeng/autopvs1) at commit 7fb1be97667e5ef576f81bf2fabbddcf9a4c7594 | All tracked files are checked against the supplied recipe before the known compatibility patch |
| Reviewed literature | Authorized original PDF/JATS/ClinVar records plus a manually reviewed catalog | Scripts verify bindings; they do not decide biological correspondence or invent quotes/grades |

A live download of a newer snapshot is a resource update, not proof of parity
with a historical run. The preparation recipes do not change runtime thresholds.

## ClinVar databases

Choose the **actual** source release. The date below is the existing snapshot
example, not an instruction to label new downloads as old data:

```bash
export RELEASE=2026-01-17
export STAMP="${RELEASE//-/}"
python preparation/prepare_resources.py clinvar --summary "$RAW/variant_summary.txt.gz" --citations "$RAW/var_citations.txt" --release "$RELEASE" --output-root "$RES"
python preparation/prepare_resources.py clinvar-protein --database "$RES/variant-knowledge/clinvar/clinvar_${STAMP}.sqlite" --summary "$RAW/variant_summary.txt.gz" --release "$RELEASE" --output-root "$RES"
```

Outputs are the coordinate database under `variant-knowledge/clinvar/` and the
protein/condition database under `assessment/grch37/clinvar-protein/`. The
builders preserve the schema, coordinate columns, comparator
filter and values. `--summary` is mandatory for the protein database so disease
identifiers and labels are not lost. SQLite build timestamps/layout can differ;
logical table equivalence, not byte equality of separately rebuilt SQLite files,
is the reconstruction criterion. Use the provided immutable databases when exact
historical file hashes are required.

## REVEL and population/prediction indexes

```bash
python preparation/prepare_resources.py revel --source "$RAW/revel.csv.gz" --sort-memory 256M --output-root "$RES"
python preparation/prepare_resources.py coverage --source "$RAW/gnomad.exomes.coverage.summary.tsv.bgz" --output-root "$RES"
python preparation/prepare_resources.py vcf-index --kind gnomad-sites --source "$RAW/gnomad.exomes.r2.1.1.sites.vcf.bgz" --output-root "$RES"
python preparation/prepare_resources.py vcf-index --kind spliceai-snv --source "$RAW/spliceai_scores.masked.snv.hg19.vcf.gz" --output-root "$RES"
python preparation/prepare_resources.py vcf-index --kind spliceai-indel --source "$RAW/spliceai_scores.masked.indel.hg19.vcf.gz" --output-root "$RES"
```

REVEL becomes `assessment/grch37/predictions/new_tabbed_revel.tsv.gz` plus its
index. All nine values are retained, including missing values and transcript
lists; sorting is stable within equal chromosome/position keys. Coverage and VCF
data go to the filenames already present in configs/vcf/resources.toml. Existing
BGZF bytes are copied without recompression; plain/gzip inputs are BGZF-compressed
without changing decompressed records. The indexing commands reject unsorted
input rather than silently reordering VCF/coverage observations. VCF inputs must
be sites-only eight-column resources, not family/patient genotype files.

## FASTA references

```bash
python preparation/prepare_resources.py fasta --source "$RAW/hg19.fa.gz" --assembly GRCh37 --kind reference --output-root "$RES"
python preparation/prepare_resources.py fasta --source "$RAW/final_hg19.fa.gz" --assembly GRCh37 --kind vep --output-root "$RES"
```

Outputs are `assessment/grch37/reference/hg19.fa` and
`assessment/grch37/vep/final_hg19.fa`, each with a newly matching .fai. This does
not perform liftover, rename contigs or validate a whole assembly against a
reference authority. Supply the verified source FASTA and preserve its hash.

## ClinGen tables and guidance

```bash
python preparation/prepare_resources.py clingen --validity "$RAW/gene_disease_validity.csv" --vcep "$RAW/variant_pathogenicity.tsv" --release 2026-08-31 --output-root "$RES"
python preparation/prepare_resources.py svi --source-dir "$RAW/clingen-svi" --output-root "$RES"
```

ClinGen export tables retain their fields; comma-delimited VCEP input is parsed
and written as TSV without replacing commas in clinical text. The command writes
`variant-knowledge/clingen/resource_manifest.json` and checks the actual reader's
header/row contracts. Guidance uses the three filenames/URLs in the supplied SVI
manifest, extracts with pdftotext -layout, verifies PDF/text hashes, and publishes
the complete `guidelines/clingen-svi-20261021/` folder. If another Poppler version
produces different text, it fails rather than weakening the binding. Alternatively
add `--text-dir "$RAW/verified-svi-text"` to supply the matching prepared text.

## AI-CURA capsule

Rank-only users can install the prepared
[two-file ranking-guidance archive](../docs/RESOURCE_ASSETS.md) directly. The
following command creates the full four-file capsule for VCF reports or resource
reconstruction; it requires the original workbook at preparation time. Reports
also require that workbook as part of their complete guidance capsule.

```bash
python preparation/prepare_resources.py aicura --workbook "$RAW/adz4172_tables_s1_to_s13.xlsx" --output-root "$RES"
```

This command checks the original workbook, reconstructs the existing S8/S9/S10
knowledge base, and requires its exact expected hash. It installs the unchanged
reviewed protocol and pinned public manifest from preparation/assets/ai-cura/, then
validates the four files with the actual runtime reader. Output is
`guidelines/ai-cura/`. It does **not** synthesize a new protocol from a paper or
change criteria/thresholds. The public manifest contains only provenance
basenames. Existing resources with the original pinned manifest remain readable.

## AutoPVS1 compatibility installation

```bash
git clone https://github.com/JiguangPeng/autopvs1.git "$RAW/autopvs1"
git -C "$RAW/autopvs1" checkout 7fb1be97667e5ef576f81bf2fabbddcf9a4c7594
python preparation/prepare_resources.py autopvs1 --source-dir "$RAW/autopvs1" --reference "$RES/assessment/grch37/reference/hg19.fa" --output-root "$RES"
```

The installer verifies the tracked source/data snapshot, copies it into
`assessment/grch37/autopvs1/`, applies only the known optional-FASTA compatibility
change, copies the explicitly supplied HG19 FASTA and creates its index. It also
accepts the already-patched verified snapshot. Unexpected source modifications
fail. HG38 annotation tables remain; absent HG38 FASTA is allowed by the existing
compatibility behavior, while the HG19 reference remains required. The isolated
import is tested before publishing. Upstream LICENSE/NOTICE files are retained.

## OMIM and GeneReviews

```bash
python preparation/prepare_resources.py omim --source "$RAW/omim_clean_full.jsonl" --mim2gene "$RAW/mim2gene.txt" --output-root "$RES"
python preparation/prepare_resources.py genereviews --source "$RAW/genereviews_official_index.tsv" --genes LMX1B GRIN2B CDKL5 --output-root "$RES"
```

OMIM validates and installs the existing normalized mirror under `core/omim/`;
it does not infer citation bindings from arbitrary HTML. The sample in
examples/resources/omim can exercise this command but is not cohort coverage.
GeneReviews writes `report/genereviews/genereviews_index.tsv`. The explicit gene
list above reproduces the current subset. Changing it changes resource coverage,
not merely file formatting; chapter text is still retrieved online at inference.
For this three-gene mapping without downloading the full index, use the supplied
`preparation/assets/genereviews_snapshot.tsv` as --source. It contains only the
same three mappings, not chapter text or full GeneReviews coverage.

## Original text and manually reviewed evidence

```bash
python preparation/prepare_resources.py extract-text --kind pdf --source "$RAW/article.pdf" --output-root "$RES"
python preparation/prepare_resources.py extract-text --kind xml --source "$RAW/identity.xml" --output-root "$RES"
python preparation/prepare_resources.py curated --catalog "$RAW/reviewed/catalog.json" --reference "$RES/assessment/grch37/reference/hg19.fa" --ack-reviewed --output-root "$RES"
```

Text extraction writes `prepared-text/<input-stem>.txt`. XML output preserves
element attributes, including identity/HGVS attributes; it is a serialized XML
text view, not a generated natural-language interpretation. Scanned PDFs with no
text fail and require separately reviewed OCR. Select and review quotes yourself.

The curated command accepts `verified_original_evidence_v1` with reviewed gene,
variant_key, identity/article metadata and quotes, plus optional family_series.
All original_file/text_file paths must stay inside the catalog's directory.
Missing source hashes can be filled from those files; existing mismatched hashes
fail. The unchanged runtime validator checks reference/allele identity, quotes,
family-series contracts and source bindings. Human biological review is still
required; --ack-reviewed is not a clinical certification. It does not invent
sources, grades or family counts. Output is
`report/curated-literature/reviewed/catalog.json` and its complete file closure;
set post_evidence_catalog to that path in the private configuration.

## Final checks

Configure the generated paths and declared releases in resources.toml. Existing
defaults match the documented snapshot names, except the generic reviewed catalog
path noted above. Run doctor with the appropriate native/report settings, then a
reader/tool check before changing an active cohort's resource version.
Do not put preparation output or receipts inside a formal case-results folder.
