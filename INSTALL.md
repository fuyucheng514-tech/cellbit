# Install Microsags

Microsags currently supports Linux x86-64. The recommended installation uses
[pixi](https://pixi.sh/) to create an isolated, reproducible environment from
conda-forge and Bioconda.

## Install with pixi (recommended)

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
pixi install
pixi run install
pixi run microsags --help
```

`pixi install` resolves the pinned compiler and runtime dependencies.
`pixi run install` builds the C++17 sources and installs the commands into the
project's pixi environment. Subsequent commands can be run with `pixi run`, for
example:

```bash
pixi run microsags \
  --manifest examples/SAGs.example.tsv \
  --out output \
  --dna-tax /absolute/path/to/genome_taxonomy_1.csv \
  --dna-packed-db /absolute/path/to/packed-index \
  --threads 8 --memory-gb 64
```

## Install with conda

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
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
pixi environment. Before a scientific run, provide:

- the GTDB R232 taxonomy table with `--dna-tax`;
- the matching, pre-built embedded packed DNA2bit index with
  `--dna-packed-db`;
- CheckM2 and GTDB-Tk databases only when the optional Stage 3B preparation or
  final quality evaluation is used.

Database versions and checksums are part of a run's reproducibility record, not
of the source installation.

## Package-manager status

Microsags is not yet published in Bioconda. Therefore commands such as
`pixi global install microsags` and `conda install -c bioconda microsags` are
not advertised yet. The repository already contains the reproducible pixi and
conda source-install paths required before a Bioconda submission.
