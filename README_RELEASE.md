# Microsags source release v0.5.0

This release contains the source-only Microsags software. Annotation passes
original FASTQ reads directly to the embedded DNA2bit engine; no read trimming,
filtering, correction, or pre-assembly stage is invoked. The release does not
contain lake samples, FASTA/FASTQ files, result directories, build caches, or
large databases.

## Install

```bash
conda env create -f environment.yml
conda activate microsags
JOBS=8 PREFIX="$CONDA_PREFIX" bash install.sh
microsags --version
```

Without Conda, install a C++ compiler and CMake first, then run `bash install.sh`.

## Configure databases

Provide the packed DNA2bit database to Annotation mode with `-d database`.
Large scientific databases are intentionally not bundled in the source tree.

## Scope

This source release contains the validated v0.5.0 direct-read Annotation path.
Historical benchmark and lake result directories remain outside the repository.

`REPOSITORY_SHA256SUMS.txt` is the authoritative, checkout-relative checksum
inventory for this v0.5.0 source release. The older `SHA256SUMS.txt` is retained
only as historical provenance for the original archived package.
