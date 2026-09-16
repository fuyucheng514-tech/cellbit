# Install

Microsags v0.2 supports Linux x86-64. The recommended installation uses
[Pixi](https://pixi.sh/) to create an isolated environment from conda-forge and
Bioconda. The checked-in `pixi.lock` freezes the resolved Linux package set.

## Install with Pixi

Install Pixi using its official instructions, then clone Microsags:

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags

pixi install --frozen
pixi run install
pixi run verify
```

`pixi install --frozen` installs the compiler and runtime dependencies without
changing `pixi.lock`. `pixi run install` builds the C++17 programs and installs
them inside the project environment. `pixi run verify` checks the public
`microsags` command.

Run any installed command through Pixi:

```bash
pixi run microsags --help
pixi run microsags annotate --help
pixi run microsags assemble --help
```

## Install the current source without Git

This route does not require Git:

```bash
wget https://github.com/fuyucheng514-tech/cellbit/archive/refs/heads/main.tar.gz -O Microsags-main.tar.gz
mkdir Microsags-main
tar -xzf Microsags-main.tar.gz -C Microsags-main --strip-components=1
cd Microsags-main
pixi install --frozen
pixi run install
pixi run verify
```

## Install with Conda

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags

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
pixi run microsags annotate SAGs/ -d "$MICROSAGS_DB" -o annotation_output
```

The packed index receipt binds its taxonomy and reference-manifest checksums.
Microsags fails closed when the database is missing or mismatched.

Stage 3B is launched by `microsags assemble`. CheckM2 and GTDB-Tk databases are
external data resources and must be supplied with `--checkm2-database` and
`--gtdbtk-data` (or their documented environment variables).

## Package-manager status

Microsags is not yet published as a Bioconda package. Therefore
`conda install -c bioconda microsags` and `pixi global install microsags` are
not currently supported. Pixi installs the locked dependencies and builds the
checked-out source repository.
