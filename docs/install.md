# Install

Microsags v0.1 supports Linux x86-64. The recommended installation uses
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
pixi run sag-stage3b-tractor --help
```

## Install from the GitHub release

This route does not require Git:

```bash
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-v0.1.0-source.tar.gz
echo "97f5b6893dfb21ac6b9a58525de4c802abeb9693eaa18dd38c18f5a631f7ba8e  Microsags-v0.1.0-source.tar.gz" | sha256sum -c -
mkdir Microsags-v0.1.0-source
tar -xzf Microsags-v0.1.0-source.tar.gz -C Microsags-v0.1.0-source
cd Microsags-v0.1.0-source
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

## Scientific databases

The source archive deliberately excludes large scientific databases. A run
must provide a mutually compatible DNA2bit taxonomy table and packed index:

```bash
export MICROSAGS_DNA_TAX=/data/GTDB232/genome_taxonomy_1.csv
export MICROSAGS_DNA_PACKED_DB=/data/GTDB232/dna2bit-packed-index
```

The packed index receipt binds its taxonomy and reference-manifest checksums.
Microsags fails closed when the database is missing or mismatched.

Stage 3B additionally requires caller-generated CheckM2 quality evidence and a
bac120 marker table. CheckM2 and GTDB-Tk databases are external data resources,
not part of the Microsags source installation.

## Package-manager status

Microsags is not yet published as a Bioconda package. Therefore
`conda install -c bioconda microsags` and `pixi global install microsags` are
not currently supported. Pixi installs the locked dependencies and builds the
checked-out source repository.
