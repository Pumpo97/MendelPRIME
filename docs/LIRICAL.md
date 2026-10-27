# Install LIRICAL

Use the **prepared LIRICAL 2.4.1 Release ZIP** by default. It contains the main
JAR and 62 library JARs; a Java 17 runtime is sufficient. Compiling is optional
and does not change how MendelPRIME invokes the tool.

## Install the Release program

The capability installer selects this ZIP automatically for either VCF mode:

```bash
python3.12 install.py --input-type vcf --work "$WORK" --download-resources --java /absolute/jdk17/bin/java
```

It installs the program under RES/genotype/lirical and writes its paths/version
to the generated resource TOML. Add `--report` for VCF reports; LIRICAL itself
is identical in both modes. Missing data is reported by doctor, not downloaded
silently. Reuse an existing archive with `--archives /absolute/downloads`.

For manual installation instead, download
[lirical-cli-2.4.1-distribution.zip](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/lirical-cli-2.4.1-distribution.zip)
and [SHA256SUMS](https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/SHA256SUMS)
into the same download directory:

```bash
export RES=/absolute/path/to/resources
export JAVA_HOME=/absolute/path/to/jdk17
export DOWNLOADS=/absolute/path/to/downloads
(
  set -e
  cd "$DOWNLOADS"
  grep '  lirical-cli-2.4.1-distribution.zip$' SHA256SUMS | sha256sum -c -
  test ! -e "$RES/genotype/lirical"
  test ! -e "$RES/genotype/lirical-cli-2.4.1"
  mkdir -p "$RES/genotype"
  unzip lirical-cli-2.4.1-distribution.zip -d "$RES/genotype"
  mv "$RES/genotype/lirical-cli-2.4.1" "$RES/genotype/lirical"
  "$JAVA_HOME/bin/java" -jar "$RES/genotype/lirical/lirical-cli-2.4.1.jar" --version
)
```

Keep the entire adjacent lib directory. In the existing genotype section of
your resource TOML, retain the other Exomiser/assembly fields and set:

```toml
[genotype]
java = "/absolute/path/to/jdk17/bin/java"
lirical_jar = "genotype/lirical/lirical-cli-2.4.1.jar"
lirical_version = "2.4.1"
lirical_data = "genotype/lirical/data"
gene_map = "genotype/lirical/data/hgnc_complete_set.txt"
```

## Install the data separately

The program ZIP does **not** contain its HPO, disease/gene mappings, HGNC or
transcript data. Reuse a matching complete data snapshot or obtain upstream data:

```bash
"$JAVA_HOME/bin/java" -jar "$RES/genotype/lirical/lirical-cli-2.4.1.jar" download -d "$RES/genotype/lirical/data"
```

A fresh download may differ from a prior snapshot; record its version/hashes.
This does not download Exomiser's genome databases. Install the official
[2406 hg19 and phenotype data](https://github.com/exomiser/Exomiser/discussions/562).
The same genome files are shared, not downloaded a second time for LIRICAL:

| Resource TOML field | Installed path beneath resource_root |
|---|---|
| genotype.assemblies[0].lirical_exomiser_data | genotype/exomiser/data/2406_hg19 |
| genotype.assemblies[0].lirical_variants_db | genotype/exomiser/data/2406_hg19/2406_hg19_variants.mv.db |
| genotype.assemblies[0].lirical_clinvar_db | genotype/exomiser/data/2406_hg19/2406_hg19_clinvar.mv.db |

The prepared ClinVar SQLite databases are for MendelPRIME's molecular evidence,
not substitutes for these Java-tool databases.

## Validate the installation

```bash
python main.py doctor --config "$CFG/vcf-rank/run.toml" --profile vcf --check-tools
```

Use vcf-report/run.toml for the report profile. Doctor checks class-path closure,
data layout and executable startup; it does not rank a case or test the LLM.

## Optional source build

Only use this route when deliberately rebuilding the program. Install Git and
JDK 17, set JAVA_HOME/PATH, and use a fresh directory outside MendelPRIME:

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

The Maven Wrapper downloads Maven and dependencies and runs upstream tests.
Its output is lirical-cli/target/lirical-cli-2.4.1-distribution.zip. Preserve that
artifact and its dependency hashes; a source revision does not pin every
transitive Maven download. Use the installation layout above, but verify against
your own build checksum rather than the prepared Release checksum.

Sources: [official release](https://github.com/TheJacksonLaboratory/LIRICAL/releases/tag/v2.4.1)
and [build instructions](https://github.com/TheJacksonLaboratory/LIRICAL/blob/v2.4.1/docs/setup.rst).
