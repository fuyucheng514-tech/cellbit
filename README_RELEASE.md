# DNA2bit-SAG Tractor original integrated release v0.1.0

This release contains the source-only original integrated software. It does not contain lake samples, FASTA/FASTQ files, result directories, build caches, or large databases.

## Install

```bash
conda env create -f environment.yml
conda activate dna2bit-sag-original
JOBS=8 PREFIX="$HOME/.local/dna2bit-sag-tractor" bash install.sh
export PATH="$HOME/.local/dna2bit-sag-tractor/bin:$PATH"
```

Without Conda, install a C++ compiler and CMake first, then run `bash install.sh`.

## Configure databases

Copy `config/paths.env.example` to `config/paths.env` and set the absolute paths to the teacher DNA2bit binary and GTDB232 database. These large external artifacts are intentionally not bundled in GitHub.

## Scope

This is the original integrated source line. It is separate from the frozen speed-v5 benchmark tree and from all lake test/result directories.
