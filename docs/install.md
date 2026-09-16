# Install

Microsags v0.2 supports Linux x86-64. The recommended installation uses Conda
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

Run installed commands directly:

```bash
microsags --help
microsags annotate --help
microsags assemble --help
```

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

## Download the k=17 GTDB232 database

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

The download is 1,389,486,830 bytes and expands to approximately 1.4 GB. It
contains 199,923 GTDB232 reference sketches generated with `k=17`,
`bit_len=55296` and `hash_type=0`. The published archive SHA256 is:

```text
a14eff367239fdcfa539b98d36eb46d24cd4388f98c07b5e5c2689fd996d0a31
```

Point Microsags at the extracted directory:

```bash
export MICROSAGS_DB="$HOME/microsags-data/dna2bit_gtdb232_packed_v1"
microsags annotate --input-type contigs SAGs/ -d "$MICROSAGS_DB" -o annotation_output
```

The packed index receipt binds its taxonomy and reference-manifest checksums.
Microsags fails closed when the database is missing or mismatched.

Stage 3B is launched by `microsags assemble`. CheckM2 and GTDB-Tk databases are
external data resources and must be supplied with `--checkm2-database` and
`--gtdbtk-data` (or their documented environment variables).

## Package-manager status

Microsags is not yet published as a Bioconda package. Therefore
`conda install -c bioconda microsags` is not currently supported. Use the
checked-in `environment.yml` and `install.sh` instead.
