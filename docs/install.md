# Install

Microsags v0.4 supports Linux x86-64. The recommended installation uses Conda
or Mamba to create an isolated environment from conda-forge and Bioconda.

## Install with Conda or Mamba

Clone Microsags, create the environment and install the program:

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags

conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

`environment.yml` installs the compiler and runtime dependencies. `install.sh`
builds the C++17 programs and installs them into the active environment.

## Install the current source without Git

This route does not require Git:

```bash
wget https://github.com/fuyucheng514-tech/cellbit/archive/refs/heads/main.tar.gz -O Microsags-main.tar.gz
mkdir Microsags-main
tar -xzf Microsags-main.tar.gz -C Microsags-main --strip-components=1
cd Microsags-main
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

## Build from source

A source build requires a C++17 compiler, CMake, pkg-config, zlib and HTSlib.
Runtime execution additionally requires the tools used by the selected input
route.

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
ctest --test-dir build --output-on-failure
cmake --install build --prefix "$HOME/.local"
"$HOME/.local/bin/microsags" --help
```

## Download the database

The source archive and the scientific database are separate Release assets.
Download the frozen packed database once:

```bash
mkdir -p "$HOME/microsags-data"
cd "$HOME/microsags-data"

wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
sha256sum -c Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
tar -xzf Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
```

## Configure Stage 3B databases once

CheckM2 and GTDB-Tk are called automatically by Assembly mode. Record their
database locations once so normal Assembly commands stay short:

```bash
mkdir -p "$HOME/.config/microsags"
cp config/paths.env.example "$HOME/.config/microsags/paths.env"
```

Edit the copied file:

```text
CHECKM2DB=/data/CheckM2/uniref100.KO.1.dmnd
GTDBTK_DATA_PATH=/data/GTDBTK/release
```

After this one-time configuration, users do not pass either database path on
every run. `--checkm2-database` and `--gtdbtk-data` remain available only as
explicit overrides.

## Build a database for a new GTDB release

If GTDB is updated, users can build a compatible database directly from the
new reference genomes and taxonomy table:

```bash
microsags sketch references/ -x taxonomy.csv -o database -t 32
```

`references/` must contain one GTDB genome FASTA per accession. Each filename
must contain its versioned `GCA_...` or `GCF_...` accession. `taxonomy.csv`
must map every reference accession to its GTDB taxonomy.

For example:

```text
references/
├── GCF_000001405.40_genomic.fna.gz
├── GCA_000002285.5_genomic.fna.gz
└── ...
```

```text
GCF_000001405.40,d__Bacteria;p__...;c__...;o__...;f__...;g__...;s__...
GCA_000002285.5,d__Bacteria;p__...;c__...;o__...;f__...;g__...;s__...
```

The taxonomy file has no header: column 1 is the versioned accession and
column 2 is its semicolon-delimited GTDB taxonomy. The accession set must match
the FASTA set exactly.

The command performs four steps automatically:

1. discover and validate the reference FASTA files;
2. generate DNA2bit sketches in parallel;
3. pack the sketches and bind them to the taxonomy table;
4. validate the one-to-one reference/taxonomy closure and write
   `COMPLETE.json` plus SHA-256 receipts.

The output is write-once. An existing `database/` is never overwritten. If the
command fails, the incomplete working directory is retained for diagnosis and
is not accepted by annotation mode. A completed database can be used directly:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

## Package-manager status

Microsags is not yet published as a Bioconda package. Therefore
`conda install -c bioconda microsags` is not currently supported. Use the
checked-in `environment.yml` and `install.sh` instead.
