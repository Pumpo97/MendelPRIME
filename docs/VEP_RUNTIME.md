# Configuring the VEP runtime

VEP is unnecessary for phenotype inputs and for the minimal VCF rank-only
profile. The VCF report profile requires VEP **104.3**, its Ensembl
`homo_sapiens/104_GRCh37` cache and the matching FASTA/index. A cache archive is
not an executable. Install the runtime separately from MendelPRIME.

## Reuse a complete exported runtime

The following configuration uses udocker 1.3.17 with the image
`ensemblorg/ensembl-vep:release_104.3`. If reusing that resource export, preserve
both the udocker source directory and the complete image/container store. Copying
only its launcher or a container alias is insufficient. The store and its
native execution helpers must be usable by the invoking non-root user.

Verify the image and container in the new location before configuring inference:

```bash
python /absolute/tools/udocker-1.3.17/udocker/udocker --repo=/absolute/tools/udocker-home images
python /absolute/tools/udocker-1.3.17/udocker/udocker --repo=/absolute/tools/udocker-home ps
```

Use the actual container name shown by `ps`. The following argument-vector
template specifies the runtime. Replace every `/absolute/...` path and
`YOUR_VEP_CONTAINER`; TOML does not expand shell variables:

```toml
[explanation_resources]
vep_command = [
  "/absolute/environment/bin/python",
  "/absolute/tools/udocker-1.3.17/udocker/udocker",
  "--repo=/absolute/tools/udocker-home",
  "run", "--nobanner",
  "--volume=/absolute/resources:/absolute/resources",
  "--volume=/absolute/work:/absolute/work",
  "--workdir={workdir}",
  "--env=PATH=/opt/vep/src/ensembl-vep:/opt/vep/src/var_c_code:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
  "--env=PERL5LIB=/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules",
  "--env=OPT=/opt/vep", "--env=OPT_SRC=/opt/vep/src", "--env=LANG=en_US.UTF-8",
  "YOUR_VEP_CONTAINER", "/opt/vep/src/ensembl-vep/vep"
]
vep_cache = "/absolute/resources/assessment/grch37/vep/cache"
vep_fasta = "/absolute/resources/assessment/grch37/vep/final_hg19.fa"
vep_cache_version = "104"
vep_tool_version = "104.3"
```

Edit the existing `[explanation_resources]` section; do not add a duplicate.
The Python executable belongs to the new environment, not another deployment.
Mount the full private work root containing caches and temporary VEP inputs;
mounting only the input VCF directory is not enough. Paths must have the same
absolute names inside and outside this example container.

```bash
mkdir -p /absolute/work/doctor
python main.py doctor --config "$CFG/vcf-report/run.toml" --profile vcf --check-tools --explanation --work-root /absolute/work/doctor
```

Doctor launches the actual command and checks its version. A subsequent real
VCF run verifies annotation output. Do not add the scientific VEP arguments to
`vep_command`; the adapter supplies them.

Do not wrap a proot-based udocker run in `strace -f`: both use ptrace and the
nested tracing can fail with `ptrace(TRACEME): Operation not permitted` even
when an ordinary VEP annotation succeeds. Use MendelPRIME's `--audit-root` for
application events, and reproduce native-tool errors without an outer tracer.

## Acquiring a runtime on a new host

The [official VEP installation guide](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
describes native Perl, Docker and Singularity deployment. For this release,
retain **104.3** rather than installing whichever VEP version is newest.
The configuration uses the image tag
`ensemblorg/ensembl-vep:release_104.3`; record its resolved image digest when
acquiring it. For a site with an authorized Docker installation, a preparation
check is:

```bash
docker pull ensemblorg/ensembl-vep:release_104.3
docker image inspect --format '{{json .RepoDigests}}' ensemblorg/ensembl-vep:release_104.3
docker run --rm --user "$(id -u):$(id -g)" ensemblorg/ensembl-vep:release_104.3 vep --help
```

This is an acquisition/version check, not the application's complete mounted
command. Do not run patient analysis as root to bypass host permissions. Use
the site's container policy, or native/rootless installation.

Obtain the **Ensembl** `homo_sapiens_vep_104_GRCh37.tar.gz` indexed cache from
the [release-104 cache directory](https://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/)
or the corresponding official archived distribution. Verify the provider's
checksums and unpack into a new `assessment/grch37/vep/cache` directory. Confirm
it contains `homo_sapiens/104_GRCh37`, not a double-nested cache or GRCh38 cache.
Prepare the selected matching FASTA/index using the
[FASTA command](../preparation/README.md#fasta-references), then set all five VEP
fields from the minimum report resource manifest. Runtime, cache and FASTA are
three separate installations; the image does not supply the large reference data.

For a different container engine, adapt the command above to its mount
and user syntax, preserving the same absolute resource/work paths and executable.
The VEP adapter appends its own scientific arguments. Validate the actual command
with `doctor --check-tools` and a real annotation before accepting a new runtime.

## Native installations

A native VEP installation may instead use its absolute executable as the sole
command prefix. Supply its Perl modules, matching cache and FASTA independently.
Other container engines require their own equivalent mounts and entry point.
The source release does not include the image, udocker helpers or databases.
HTTPS compatibility settings in MendelPRIME do not configure
separate package managers, native tools or container engines.
