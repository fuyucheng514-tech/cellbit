# Install Microsags

Microsags currently supports Linux x86-64. The recommended installation uses
Conda or Mamba to create an isolated environment from conda-forge and Bioconda.

## Install with Conda or Mamba

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

`environment.yml` installs the compiler and runtime dependencies. `install.sh`
builds the C++17 sources and installs Microsags into the active environment.
Subsequent commands are invoked directly, for example:

```bash
microsags annotate examples/SAGs/*.fna \
  --database /data/Microsags-GTDB232-DNA2bit-k17-packed-v1 \
  --output output --threads 8
```

## Install from the GitHub release

```bash
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-v0.1.0-source.tar.gz
echo "97f5b6893dfb21ac6b9a58525de4c802abeb9693eaa18dd38c18f5a631f7ba8e  Microsags-v0.1.0-source.tar.gz" | sha256sum -c -
mkdir Microsags-v0.1.0-source
tar -xzf Microsags-v0.1.0-source.tar.gz -C Microsags-v0.1.0-source
cd Microsags-v0.1.0-source
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

## Build from source

Install a C++17 compiler, CMake, pkg-config, zlib and HTSlib, then run:

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
cmake --install build --prefix "$HOME/.local"
microsags --help
```

## Scientific databases

Large scientific databases are deliberately not stored in GitHub or inside the
Conda environment. Before a scientific run, provide:

- the GTDB R232 taxonomy table with `--dna-tax`;
- the matching, pre-built embedded packed DNA2bit index with
  `--dna-packed-db`;
- CheckM2 and GTDB-Tk databases only when the optional Stage 3B preparation or
  final quality evaluation is used.

Database versions and checksums are part of a run's reproducibility record, not
of the source installation.

## Package-manager status

Microsags is not yet published as a Bioconda package, so
`conda install -c bioconda microsags` is not currently supported. Installation
uses the checked-in `environment.yml` followed by `install.sh`.
