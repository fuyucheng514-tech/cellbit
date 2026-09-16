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
FASTQ reads are quality-controlled with fastp and passed directly to DNA2bit;
SPAdes is not run.

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

## Documentation

Full documentation is available at
[fuyucheng514-tech.github.io/cellbit](https://fuyucheng514-tech.github.io/cellbit/).

## Citation

A Microsags manuscript and formal citation are in preparation. DNA2bit and
external dependencies should be cited according to their original publications.
