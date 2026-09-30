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
Stage 3B additionally requires CheckM2 1.0.1 and GTDB-Tk 2.7.2, installed
separately because their Python environments differ. See the
[Install guide](docs/install.md#install-and-configure-stage-3b-dependencies-once)
for the commands and one-time executable/database path configuration.
Subsequent commands are invoked directly, for example:

```bash
microsags annotate --input-type contigs examples/SAGs/ \
  -d database -o annotation_output --threads 8
```

## Install the current source without Git

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

Annotation mode passes original FASTQ reads directly to the embedded DNA2bit
engine, without fastp. From v0.5.1, database construction uses persistent native
C++ threads; the annotation algorithm and database format are unchanged.

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

The Conda environment includes libdeflate. For manual builds, optionally install
the `libdeflate` development library before configuring CMake
to accelerate decoding of small gzip reference files. CMake detects it
automatically. Without it, Microsags retains its zlib reader and produces the
same sketches. To explicitly disable this optional path, configure with
`-DDNA2BIT_ENABLE_LIBDEFLATE=OFF`. Large files and multi-member gzip streams
continue to use the streaming reader; scientific parameters are unchanged.

## Scientific databases

Large scientific databases are deliberately not stored in GitHub or inside the
Conda environment. Before a scientific run, provide:

- the matching Microsags DNA2bit database with `-d database` in Annotation;
- CheckM2 and GTDB-Tk programs plus their databases when Assembly has
  unclassified SAGs and enters Stage 3B.

Database versions and checksums are part of a run's reproducibility record, not
of the source installation.

## Package-manager status

Microsags is not yet published as a Bioconda package, so
`conda install -c bioconda microsags` is not currently supported. Installation
uses the checked-in `environment.yml` followed by `install.sh`.
