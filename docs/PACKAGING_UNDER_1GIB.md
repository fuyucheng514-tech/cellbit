# Contigs-start Lake packaging under 1 GiB

## Decision

Ship a **thin code package** and keep every release-specific database, Lake
input/evidence tree, and third-party environment outside it.  The final archive,
not an estimate or `du -h` rounding, must be at most `1,073,741,824` bytes.

The repository itself is about 0.6 MiB before compiled binaries.  The compiled
project executables are small (the server's stripped optimized Stage 3B
orchestrator was 289,080 bytes), so the 1 GiB limit is not a code-size problem.
It is impossible to satisfy by embedding the scientific references: the
CheckM2 DIAMOND file alone is 3,082,500,605 bytes, the existing skani R232
payload is described by the project audits as about 75 GB, and the experimental
`gtdb-ani-af` `SKETCHES.bin` observed in the server work is about 76 GB (with an
index-build temporary peak around 276 GB).  These observations are inventory
facts, not new benchmark results.

Recommended archive layout:

```text
cellbit-sag-linux-x86_64/
  bin/
    dna2bit-sag-pipeline
    sag-stage3b-tractor
    gtdb-ani-af
    cpp-subass
  libexec/
    stage3b_signed_leiden.py
  requirements-stage3b.txt
  docs/
  LICENSES/                 # only licenses/notices actually required
  BUILD_RECEIPT.json        # compiler/flags/source hashes
  RUNTIME_PATHS.example     # paths only; no database symlinks
```

Do not put symlinks to external databases in this tree.  Besides making the
archive non-self-contained, a later `tar --dereference` could silently ingest
the target.  Supply absolute paths in the launch configuration and bind them in
the run receipt.

## Runtime dependency boundary

For a manifest containing only contig FASTA inputs, `src/main.cpp` skips fastp
and SPAdes.  The actual 1–3A runtime is:

- the four project executables above (only `dna2bit-sag-pipeline` and
  `cpp-subass` are reached when 3B is not launched);
- the teacher binary `/home/data/shared/software/Dna2bit/dna2bit`;
- the external R232 bit database
  `/home/data/shared/software/Dna2bit/GTDB232/GTDB` and taxonomy
  `/home/data/shared/software/Dna2bit/GTDB232/genome_taxonomy_1.csv`;
- for labelled groups, the Flye prefix
  `/home/data/fyc/biosoft/miniconda3/envs/assemble`, specifically
  `bin/flye-modules`, `bin/flye-minimap2`, `bin/flye-samtools`,
  `lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg`, and the two
  polishing matrices in that config tree.

The optimized Stage 3B runtime additionally needs:

- the packaged `gtdb-ani-af` executable and the quality-passing
  Dna2bit-negative SAG FASTA files; current Stage3B uses exact SAG–SAG triangle
  and needs no GTDB ANI index, reference FASTA collection, or taxonomy argument;
- BLAST+ `makeblastdb` and `blastn`;
- a Python environment satisfying `requirements-stage3b.txt` and the packaged
  `stage3b_signed_leiden.py`;
- the same cpp-subass/Flye runtime described above.

CheckM2 1.0.1 and GTDB-Tk 2.7.2 are **not invoked by
`sag-stage3b-tractor`**.  They are conditional, external upstream evidence
generators for `--quality-manifest` and `--marker-map`.  A raw-evidence rebuild
therefore also needs the CheckM2 database
`/home/data/fyc/past/root_archive_20260806/checkm2_db/uniref100.KO.1.dmnd`
(exactly 3,082,500,605 bytes in the replay contract), the GTDB-Tk environment,
  and `/home/data/temp/release232`.  Historical Lake replay, as opposed to the
  current exact SAG–SAG path, additionally uses
`/home/data/temp/release232/skani/database` and the frozen Lake scripts/data.
Keep all of these out of the distributable.

The project binaries should remain dynamically linked.  In particular:

- `dna2bit-sag-pipeline` and `gtdb-ani-af`: zlib plus the normal C/C++ runtime;
- `sag-stage3b-tractor`: the normal C/C++ runtime;
- `cpp-subass`: htslib, zlib, threads, and whatever shared libraries that exact
  htslib build resolves (often compression, crypto, and curl libraries).

Do not guess or copy an entire system library tree.  Run `ldd` through the audit
below on the deployment host; install only unresolved non-system libraries in a
separate runtime prefix, or use the host's pinned module/conda environment.
BLAST+, Flye, dna2bit, the Leiden Python environment, CheckM2, and GTDB-Tk are
also external runtime prefixes and are not charged to the code archive.

## Reproducible build and release gate

Build/install into a new staging directory, strip only the four project
executables, and create the archive without dereferencing links.  The exact
compiler and archive program/version should be recorded in `BUILD_RECEIPT.json`.
A suitable sequence is:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=MinSizeRel \
  -DCMAKE_INSTALL_PREFIX="$PWD/release-root"
cmake --build build -j 16
cmake --install build
strip --strip-unneeded release-root/bin/{dna2bit-sag-pipeline,sag-stage3b-tractor,gtdb-ani-af,cpp-subass}
install -Dm755 python/stage3b_signed_leiden.py release-root/libexec/stage3b_signed_leiden.py
install -Dm644 requirements-stage3b.txt release-root/requirements-stage3b.txt
tar --sort=name --mtime='UTC 1970-01-01' --owner=0 --group=0 --numeric-owner \
  -czf cellbit-sag-linux-x86_64.tar.gz -C release-root .
```

Run the footprint audit against both the staging tree and final archive.  This
example selects the full contigs-start path; replace the placeholder with the
exact deployed Stage3B Python prefix:

```bash
python3 tools/audit_package_footprint.py \
  --package-root release-root \
  --archive cellbit-sag-linux-x86_64.tar.gz --require-archive \
  --external-runtime dna2bit=/home/data/shared/software/Dna2bit/dna2bit \
  --external-runtime flye=/home/data/fyc/biosoft/miniconda3/envs/assemble \
  --external-runtime blastn=/home/data/fyc/biosoft/miniconda3/envs/assemble/bin/blastn \
  --external-runtime makeblastdb=/home/data/fyc/biosoft/miniconda3/envs/assemble/bin/makeblastdb \
  --external-runtime stage3b_python=/ABSOLUTE/PINNED/STAGE3B/PYTHON/PREFIX \
  --external-reference dna2bit_r232=/home/data/shared/software/Dna2bit/GTDB232/GTDB \
  --external-reference dna2bit_taxonomy=/home/data/shared/software/Dna2bit/GTDB232/genome_taxonomy_1.csv \
  --ldd-binary bin/dna2bit-sag-pipeline \
  --ldd-binary bin/sag-stage3b-tractor \
  --ldd-binary bin/gtdb-ani-af \
  --ldd-binary bin/cpp-subass \
  --json-out PACKAGE_FOOTPRINT.json
```

For the optional raw-evidence rebuild add:

```bash
  --external-runtime checkm2=/ABSOLUTE/CHECKM2_1.0.1/PREFIX \
  --external-runtime gtdbtk=/home/data/fyc/biosoft/miniconda3/envs/gtdbtk \
  --external-reference checkm2_db=/home/data/fyc/past/root_archive_20260806/checkm2_db/uniref100.KO.1.dmnd \
  --known-bytes checkm2_db=3082500605 \
  --external-reference gtdb_r232=/home/data/temp/release232
```

`PACKAGE_FOOTPRINT.json` is the release report.  It contains the exact byte
total and SHA-256 of every regular package file, the exact archive byte count
and SHA-256, non-dereferenced sizes of every declared external component, and
resolved `ldd` libraries.  It fails if either package/archive exceeds 1 GiB, an
external component is placed inside the package, a package symlink escapes the
root, a known reference payload name/suffix (`SKETCHES.bin`, `POSTINGS.bin`,
`REFS.tsv`, `.bit`, `.dmnd`, `.mmi`, or `.msh`) appears in the package, a
required external path is missing, a known byte count drifts, or a dynamic
library is unresolved.  Use `--allow-missing-external` only for a workstation
plan-only inventory; never for a release gate.
