# Installing MendelPRIME

MendelPRIME uses the distribution name `mendelprime` and exposes one console command,
`mendelprime`. Use a fresh virtual environment: another distribution that owns the same
`mendelprime` import package or console command must not be installed beside it.

## Supported platform

Use Python 3.12 on Linux x86-64. Package metadata limits Python to `>=3.12,<3.13`.
Other platforms require compatible Python packages, native executables and libraries.

## Install a capability group

From the unpacked source release, select the smallest group needed:

```bash
python3.12 -m venv .venv
. .venv/bin/activate

python -m pip install -c constraints/runtime.txt .
python -m pip install -c constraints/runtime.txt ".[hpo]"
python -m pip install -c constraints/runtime.txt ".[vcf]"
python -m pip install -c constraints/runtime.txt ".[full]"
```

Run only one of the four installation commands. The base installation (`.`) supports
complete HPO/clinical ranking and reports with RareA disabled, as well as
configuration, result evaluation and CLI discovery. Its direct Python dependency
is PyYAML 6.0.3. `hpo` adds RareA retrieval. `vcf` adds VCF parsing,
HGVS normalization, AutoPVS1 support, and HTML evidence parsing. `full` contains both application
paths.

The corresponding dependency-only files are `requirements/base.txt`, `requirements/hpo.txt`,
`requirements/vcf.txt`, and `requirements/full.txt`. They use `constraints/runtime.txt`, which records
the pinned dependency versions for Python 3.12 on Linux.

When `[ranking].similar_cases = false`, the base dependencies suffice for HPO and
clinical input; install `.[vcf]` to add VCF without installing GPU/model packages.
Do not confuse this switch with `device = "cpu"`: CPU alone still runs the
configured similarity models when the branch is enabled.

Keep the source archive (or its separate `configs/`, `constraints/`, `examples/`
and `docs/` deployment files) when installing a wheel. The wheel contains the
runtime and Prompts, not those deployment templates or any external resources.
No configuration is discovered from an old editable install or a Conda environment.

## PyTorch builds

The package pins the PyTorch API version but does not choose a CPU or accelerator wheel source.
Select the appropriate official PyTorch wheel for the deployment platform before installing the HPO
or full group, then install MendelPRIME with the runtime constraints. Match the
PyTorch build to the installed accelerator driver; the package requires PyTorch 2.9.1.

## Native and external dependencies

VCF analysis additionally requires externally managed scientific tools and datasets. These are not
embedded in the wheel:

- Java 17 for the complete Exomiser 14.0.0 and LIRICAL 2.4.1 distributions,
  including each launcher's adjacent `lib/` dependencies.
- VEP 104.3 with its cache and reference FASTA when VCF reports or ranking VEP are
  enabled. It is unnecessary for minimal VCF ranking-only or non-VCF input.
- A complete AutoPVS1 source and data tree. Its Python imports are covered by the `vcf` extra.
- PLINK2 only when optional kinship QC is enabled.
- The PostgreSQL client library (`libpq`) used by the `psycopg2` dependency of
  `hgvs`. If no suitable binary wheel is available, building `psycopg2` also needs
  a C compiler, Python development headers and PostgreSQL development headers/
  `pg_config`. A Python virtual environment isolates Python packages, not native
  shared libraries; record and verify those libraries on the deployment host.

For a user-owned Conda environment, installing `python=3.12 libpq pip` and
activating it puts `pg_config` on PATH. When invoking the Python executable by
absolute path without activation, prepend the environment's `bin` to PATH for
the pip build. Setting `PG_CONFIG` alone is not sufficient for psycopg2's build
discovery. Keep the generated environment outside the source tree.

If Python 3.12 is absent, use an existing Conda/Miniforge installation, or first
install Miniforge into a user-owned directory using the upstream installer and
verify its published checksum. Then use this alternative to the README's venv
block; do not stack the two environment recipes:

```bash
export CODE=/absolute/path/to/github
export WORK="$HOME/mendelprime-work"
test ! -e "$WORK"
umask 077
mkdir -p "$WORK"
conda create -y -n mendelprime -c conda-forge python=3.12 libpq pip
conda activate mendelprime
cd "$CODE"
python -m pip install -c constraints/runtime.txt '.[vcf]'
python -m pip check
```

This example installs all three input types without retrieval models; set
`similar_cases = false` in both phenotype run profiles. For enabled retrieval,
select the appropriate PyTorch build and install `.[full]` instead. A C compiler
and headers may still be required when a dependency must build from source.
Use an unused Conda environment name if `mendelprime` already exists. Continue
with the private WORK/config/input setup from the README after activation;
there is no need to create its `.venv` directory in this alternative.

`pysam` provides the VCF, BGZF, and tabix APIs used by MendelPRIME; host `samtools`, `bcftools`,
`bgzip`, and `tabix` commands are not Python-package prerequisites.

For LIRICAL, a JRE is enough to execute a built distribution; building the pinned
official source requires a JDK with `javac`, Git and network access to Maven
repositories. The [build instructions](LIRICAL.md) use the upstream Maven Wrapper
and keep all generated files outside MendelPRIME. Install the [ClinVar and
RareA assets](RESOURCE_ASSETS.md) separately. OMIM still uses local files.

## Reproducible build

Build in a clean Python 3.12 environment using the pinned backend and frontend versions:

```bash
python -m pip install -c constraints/build.txt build setuptools wheel
python -m build --no-isolation
```

The wheel must contain the `mendelprime` package, all ranking and explanation Prompt assets, and the
single `mendelprime = mendelprime.runtime.cli:main` entry point. Configuration, credentials, model
weights, patient data, databases, native tools, and scientific resources are deployment inputs and
must not be included in the distribution.

After installation, verify dependency consistency and CLI discovery without making model calls:

```bash
python -m pip check
mendelprime --help
python main.py doctor --config "$CFG/minimal/vcf_report.toml" --profile vcf --check-tools
```
