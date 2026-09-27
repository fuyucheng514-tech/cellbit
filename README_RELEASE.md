# Microsags source release v0.5.2

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

For manual builds, see `INSTALL.md` for compiler and library dependencies.

## Configure databases

Provide the packed DNA2bit database to Annotation mode with `-d database`.
Large scientific databases are intentionally not bundled in the source tree.

## Scope

This release retains direct-read Annotation and the persistent native sketch
workers introduced in v0.5.1. It adds bounded libdeflate gzip decoding and
cache-friendly, overflow-guarded sketch counters. The Conda environment includes
libdeflate; manual builds still work without it using the existing zlib reader.
Sketch bytes, default parameters, database format, annotation rules and assembly
logic are unchanged. Existing packed databases do not need to be rebuilt.
Historical benchmark and lake result directories remain outside the repository.

Compatibility evidence and regression commands are recorded in
`tests/SKETCH_OPTIMIZATION_VALIDATION.md`. No comparative speed claim is made
from differently selected benchmark runs.

`REPOSITORY_SHA256SUMS.txt` is the authoritative, checkout-relative checksum
inventory for this v0.5.2 source release. The older `SHA256SUMS.txt` is retained
only as historical provenance for the original archived package.
