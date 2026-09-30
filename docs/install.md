# Install

Microsags v0.5.2 supports Linux x86-64. The prebuilt download is the simplest
route: no compiler, Conda installation, or source checkout is required.

## Prebuilt Linux download (recommended)

Download and verify the program archive, then extract it into its final
location:

```bash
mkdir -p "$HOME/.local/share/microsags/v0.5.2"
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.5.2-linux1/Microsags-v0.5.2-linux-x86_64.tar.gz
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.5.2-linux1/Microsags-v0.5.2-linux-x86_64.tar.gz.sha256
sha256sum -c Microsags-v0.5.2-linux-x86_64.tar.gz.sha256
tar -xzf Microsags-v0.5.2-linux-x86_64.tar.gz -C "$HOME/.local/share/microsags/v0.5.2"
export PATH="$HOME/.local/share/microsags/v0.5.2/bin:$PATH"
conda-unpack
microsags --help
```

Run `conda-unpack` once after extraction. Despite its name, this command is
included inside the archive: you do not need to install Conda. Keep the
extracted directory in place after running it. Add the `export PATH=...` line
to your shell startup file if you want `microsags` available in new sessions.
The GTDB database is a separate download below. Assembly additionally needs
the external CheckM2 and GTDB-Tk programs and their databases.

## Build from a Git checkout with Conda or Mamba

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
It also includes libdeflate for accelerated reference-file decoding.

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

Optionally install the libdeflate development library for faster gzip decoding;
CMake detects it automatically. Without it, the zlib streaming reader remains
available and produces the same sketches. To explicitly disable the accelerator,
add `-DDNA2BIT_ENABLE_LIBDEFLATE=OFF` when configuring CMake.

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

The source code and the scientific database are separate downloads. Download
the frozen packed database once:

```bash
mkdir -p "$HOME/microsags-data"
cd "$HOME/microsags-data"

wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
sha256sum -c Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
tar -xzf Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
```

If GitHub transfers time out on a compute server, download the source archive,
database archive, and checksum file on a connected machine and copy them to the
server. Then use the archive-based install above and verify the database
checksum before extraction.

## Install and configure Stage 3B dependencies once

Assembly mode calls CheckM2 and GTDB-Tk as external programs. They are not
installed by Microsags' `environment.yml`; install them in separate environments
so their Python dependencies do not conflict:

```bash
mamba create -n microsags-checkm2 -c conda-forge -c bioconda checkm2=1.0.1
mamba create -n microsags-gtdbtk -c conda-forge -c bioconda gtdbtk=2.7.2
```

Their databases are separate downloads. Once both programs and databases are
available, record the absolute paths so normal Assembly commands stay short:

```bash
mkdir -p "$HOME/.config/microsags"
cp config/paths.env.example "$HOME/.config/microsags/paths.env"
```

Edit the copied file:

```text
CHECKM2_BIN=/absolute/path/to/microsags-checkm2/bin/checkm2
CHECKM2DB=/data/CheckM2/uniref100.KO.1.dmnd
GTDBTK_BIN=/absolute/path/to/microsags-gtdbtk/bin/gtdbtk
GTDBTK_DATA_PATH=/data/GTDBTK/release
```

After this one-time configuration, users do not pass these four paths on every
run. The `--checkm2`, `--gtdbtk`, `--checkm2-database`, and `--gtdbtk-data`
options remain available as explicit overrides. Microsags runs each external
program using its own environment's executables rather than the Microsags
runtime libraries.

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

From v0.5.1, `-t` controls persistent C++ sketch threads in one process,
instead of launching a program separately for every reference. The sketch
algorithm, parameters and database format are unchanged. Packing remains a
separate stage and is included in total database-construction time.

Version 0.5.2 also optimizes gzip decoding and sketch-counter memory access,
without changing sketch values or the database format. Existing databases do
not need to be rebuilt.

The output is write-once. An existing `database/` is never overwritten. If the
command fails, the incomplete working directory is retained for diagnosis and
is not accepted by annotation mode. A completed database can be used directly:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

## Package-manager status

Microsags is not published as a Bioconda package. Use the prebuilt Linux
download or build from source.
