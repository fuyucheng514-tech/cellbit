# Microsags

Microsags is a command-line workflow for DNA2bit annotation and
species-guided subassembly of single-amplified genomes (SAGs).

## Install

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
```

## Annotation mode

Annotation accepts paired FASTQ, singleton FASTQ, or per-SAG contigs.
Original FASTQ reads are passed directly to the embedded DNA2bit engine.
No trimming, filtering, read correction, or SPAdes step is run.

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

The only published result is `annotation_output/annotations.tsv`, with one row
per input SAG and exactly two columns: `sag_id` and `species`. Rejected or
unmatched SAGs are reported as `UNCLASSIFIED`.

## Assembly mode

Assembly is the next step and accepts contigs only. If annotation used reads,
assemble every SAG independently first and preserve the same SAG identifiers.

```bash
microsags assemble SAG_contigs/ \
  -a annotation_output/annotations.tsv \
  -o assembly_output
```

Labelled SAGs enter Stage 3A. Unclassified SAGs enter Stage 3B. Both retain the
existing `cpp-subass`/Flye subassembly workflow.

To override only the Stage 3B Leiden resolution:

```bash
microsags assemble SAG_contigs/ \
  -a annotation_output/annotations.tsv \
  -o assembly_output \
  --leiden-resolution 0.18
```

If `--leiden-resolution` is omitted, the frozen default sweep is used.

For a real downloadable 3A/3B HQ example, see
[`examples/pa_hq_3a_3b/`](examples/pa_hq_3a_3b/README.md). Stage 3B needs
separately installed CheckM2 and GTDB-Tk programs and databases; configure
their paths once using the [Install guide](docs/install.md).

Assembly publishes only `stage3a.tsv`, final Stage 3A FASTA files,
`stage3b_clusters.tsv`, and final Stage 3B FASTA files. Temporary scientific
work is kept outside the result directory on failure and removed after a
successful publication.

## Build a database

```bash
microsags sketch references/ -x taxonomy.csv -o database -t 32
```

This mode is intended for a new GTDB release. `references/` contains one GTDB
genome FASTA per accession and `taxonomy.csv` maps those accessions to GTDB
taxonomy. Microsags sketches the references, builds the packed database,
validates that every reference has taxonomy, and writes a `COMPLETE.json`
receipt. The output directory is write-once and can be passed directly to
`microsags annotate -d database`.

In v0.5.1, reference sketching uses one native C++ process with `-t` persistent
threads. It no longer starts a separate program for each reference. DNA2bit
sketch parameters, packed database format and annotation rules are unchanged.

Version 0.5.2 further optimizes gzip decoding and sketch-counter memory access.
The Conda installation includes the optional libdeflate accelerator. Existing
databases and command lines remain compatible; no database rebuild is required.

## Documentation

Full documentation is available at
[fuyucheng514-tech.github.io/cellbit](https://fuyucheng514-tech.github.io/cellbit/).

The [performance page](https://fuyucheng514-tech.github.io/cellbit/performance/)
shows four separately downloadable SVG panels for simulated SAGs and 24 real
samples. The current annotation path avoids temporary per-SAG read receipts,
processes compressed contigs directly, and uses the full requested thread
budget for independent sketch tasks. Annotation labels remain unchanged in
the verified regression runs.

## Citation

A Microsags manuscript and formal citation are in preparation. DNA2bit and
external dependencies should be cited according to their original publications.
