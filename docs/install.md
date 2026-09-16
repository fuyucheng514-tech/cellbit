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

## Package-manager status

Microsags is not yet published as a Bioconda package. Therefore
`conda install -c bioconda microsags` is not currently supported. Use the
checked-in `environment.yml` and `install.sh` instead.
