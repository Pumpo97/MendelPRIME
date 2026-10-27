# External resource bundle

The code release does not contain large scientific datasets, model weights,
native-tool databases, or the VEP container. Place a prepared resource tree in a
sibling directory named `MendelPRIME-release-resources`, or edit each profile's
`resources.toml` to point to another installation.

Download the prepared archives and checksums from the
[1.0.0 release](https://github.com/Pumpo97/MendelPRIME/releases/tag/1.0.0).
Choose only the assets needed by your profile using the
[download and extraction guide](../docs/RESOURCE_ASSETS.md).

Expected top-level layout:

```text
core/{hpo,omim}
phenotype/{phen2disease,similar-cases,models}
genotype/{exomiser,lirical}
variant-knowledge/{clinvar,clingen}
assessment/grch37/{reference,vep,gnomad,predictions,autopvs1}
assessment/transcripts
guidelines/{ai-cura,clingen-svi-20261021}
report/{clingen-dosage,genereviews,curated-literature}
tools/{java,plink2}
```

Only install the directories needed by the selected input/report profile; this
tree is the union of capabilities, not a mandatory download list for every run.
See the [mode-specific matrix](../docs/CONFIGURATION.md#4-resource-parameters-by-profile)
and [acquisition/preparation guide](../docs/RESOURCE_PREPARATION.md).

## Choose the smallest resource set

| Run mode | Prepared archives to use | Additional preparation |
|---|---|---|
| HPO or HPO+clinical, with or without reports | HPO | Complete authorized OMIM mirror and matching mim2gene; Phen2Disease four-file dataset; configure LLM and phenotype-service access |
| VCF ranking without reports | HPO, ClinVar, VCF ranking guidance, LIRICAL program | Java, Exomiser program and both 2406 datasets, LIRICAL data, OMIM, ClinGen tables, reference/index, gnomAD, SpliceAI, REVEL and AutoPVS1 |
| Add VCF reports | Reuse the VCF ranking archives | VEP runtime/cache/reference, MANE, full AI-CURA capsule and ClinGen SVI guidance |
| Add RareA retrieval | RareA and retrieval support | Matching Qwen weight shards and optional retrieval Python dependencies |

The archives are not a complete installed environment. A LIRICAL program archive
does not contain its databases, and a retrieval-support archive does not contain
the Qwen weights. The one-record OMIM demo documents the format; it cannot replace
the full mirror. ClinVar SQLite cannot replace Exomiser/LIRICAL's native databases.

For each resource, the [preparation guide](../preparation/README.md) gives its
upstream source, conversion command and output location. Ready-made resources
and native-tool installation are covered in [RESOURCE_PREPARATION.md](../docs/RESOURCE_PREPARATION.md).
Set the selected resource TOML to those outputs, then run the corresponding
`doctor` command from the [main README](../README.md#run-a-case-or-cohort).

Small source-bound report resources may be distributed separately from the large
installed libraries. Keep manifests beside their relative files and preserve
hash-bound catalog directories intact. See [configuration](../docs/CONFIGURATION.md) for the complete
parameter/resource matrix, GRCh37 restriction, indexes, and VEP closure.

The frozen HPO ontology/statistics/annotation triple is distributed in
`mendelprime-hpo-20260728-v1-public.tar.gz`, with payload files under `core/hpo/`.
ClinVar's two prepared SQLite files and the aligned RareA vector index
are separate resource release assets, not files to commit into the source
repository. See [RESOURCE_ASSETS.md](../docs/RESOURCE_ASSETS.md) for their checksums,
extraction layout and model requirements, and [LIRICAL.md](../docs/LIRICAL.md)
for the official 2.4.1 source build.
OMIM remains local; no OMIM API key is required by the current code.
