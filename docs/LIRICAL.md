# Build and install LIRICAL

The publication VCF profile selects official LIRICAL **2.4.1**. Its release fixes
omitted runtime-scoped dependencies in the distribution. MendelPRIME still calls
the Java CLI through its existing Medical Tools adapter; no service, Python
reimplementation or alternate ranking pipeline is introduced.

Sources: [release](https://github.com/TheJacksonLaboratory/LIRICAL/releases/tag/v2.4.1)
and [upstream build instructions](https://github.com/TheJacksonLaboratory/LIRICAL/blob/v2.4.1/docs/setup.rst).

## Build the pinned source

Use a Linux user-owned directory outside the MendelPRIME source tree. Install Git,
unzip and a **JDK 17** first. Set `JAVA_HOME` to its actual installation directory;
the system Java does not have to be replaced. `java` without `javac` is not enough
for compiling. Commands below require HTTPS access to GitHub and the Maven
repositories declared by upstream, including snapshot repositories if requested
by dependency resolution. No OMIM or LLM key is involved.

```bash
export JAVA_HOME=/absolute/path/to/jdk17
export PATH="$JAVA_HOME/bin:$PATH"
java -version
javac -version

mkdir -p "$HOME/mendelprime-build"
cd "$HOME/mendelprime-build"
git clone --depth 1 --branch v2.4.1 https://github.com/TheJacksonLaboratory/LIRICAL.git
cd LIRICAL
test "$(git rev-parse HEAD)" = 6e6add2df4023493aa05ea7718932db9ce97485c
./mvnw -B -Prelease -Dprotobuf.version=4.36.2 -Dmaven.repo.local="$HOME/mendelprime-build/m2" package
```

The wrapper downloads Maven; a system Maven installation is not required. The
command builds the official source and runs its tests. It does not build or
download Exomiser's genome databases. It does not require a locally modified
LIRICAL source tree or a manually assembled `lib/` directory.

Upstream declares a Protobuf version range. The command above pins the Maven
property to stable 4.36.2 for generated code and runtime, without changing Java source.

The distribution is `lirical-cli/target/lirical-cli-2.4.1-distribution.zip`.
Record its checksum and retain the complete ZIP. A pinned source revision alone
does not pin every transitive Maven download; preserve the built artifact and
its dependency hashes when reproducing a particular deployment.

## Install without overwriting an existing deployment

Set `RES` to a new resource root, or use a new versioned tool directory when an
existing deployment must remain available. The example refuses an existing
`genotype/lirical` directory. Other resources, including LIRICAL's data, are
separate from this program distribution.

```bash
RES=/absolute/path/to/MendelPRIME-release-resources
ZIP="$HOME/mendelprime-build/LIRICAL/lirical-cli/target/lirical-cli-2.4.1-distribution.zip"
test -f "$ZIP"
test ! -e "$RES/genotype/lirical"
mkdir -p "$RES/genotype"
unzip "$ZIP" -d "$RES/genotype"
mv "$RES/genotype/lirical-cli-2.4.1" "$RES/genotype/lirical"
"$JAVA_HOME/bin/java" -jar "$RES/genotype/lirical/lirical-cli-2.4.1.jar" --version
```

Keep the adjacent `lib/` directory. In the VCF `resources.toml`, set the following
fields in the existing `[genotype]` section, using the actual Java path:

```toml
[genotype]
java = "/absolute/path/to/jdk17/bin/java"
lirical_jar = "genotype/lirical/lirical-cli-2.4.1.jar"
lirical_version = "2.4.1"
lirical_data = "genotype/lirical/data"
gene_map = "genotype/lirical/data/hgnc_complete_set.txt"
```

Retain the other Exomiser, assembly and analysis fields. Old configs without an
explicit LIRICAL version keep the legacy default; set `lirical_version` explicitly
when adopting 2.4.1. A bare path change is not sufficient.

## Install the data separately

Prefer the fixed, tested LIRICAL data snapshot and Exomiser **2406** GRCh37 files
for this MendelPRIME profile. If acquiring LIRICAL data from upstream instead:

```bash
"$JAVA_HOME/bin/java" -jar "$RES/genotype/lirical/lirical-cli-2.4.1.jar" download -d "$RES/genotype/lirical/data"
```

This downloads the then-current HPO, disease/gene mappings and transcript files,
not necessarily the historical snapshot. Record versions and validate the result.
It does not download Exomiser's `2406_hg19_variants.mv.db` and
`2406_hg19_clinvar.mv.db`; obtain those from the
[Exomiser data release](https://github.com/exomiser/Exomiser/discussions/562).
MendelPRIME's prepared ClinVar SQLite databases cannot replace these Java-tool
databases. Required files and indexes are listed in [the resource register](RESOURCE_REGISTER.md)
and [preflight requirements](CONFIGURATION.md#6-resource-layout-and-preflight).

## Validate the actual tool

From the MendelPRIME source directory, after updating the profile:

```bash
python main.py doctor --config "$CFG/minimal/vcf_rank.toml" --profile vcf --check-tools
```

This checks the declared library class path and starts the configured executable.
Then run one real VCF case using the ordinary `python main.py run` entry. A version
banner alone does not prove genotype analysis works. Record the tool version
with your results; different releases can produce different outputs.
