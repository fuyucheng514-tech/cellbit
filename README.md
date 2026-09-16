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

The principal result is `annotation_output/02_dna2bit/labels.tsv`.

## Assembly mode

Assembly is the next step and accepts contigs only. If annotation used reads,
assemble every SAG independently first and preserve the same SAG identifiers.

```bash
microsags assemble --input-type contigs SAG_contigs/ \
  --annotations annotation_output \
  -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Labelled SAGs enter Stage 3A. Unclassified SAGs enter Stage 3B. Both retain the
existing `cpp-subass`/Flye subassembly workflow.

To override only the Stage 3B Leiden resolution:

```bash
microsags assemble --input-type contigs SAG_contigs/ \
  --annotations annotation_output \
  -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database \
  --leiden-resolution 0.18
```

If `--leiden-resolution` is omitted, the frozen default sweep is used.

## Build a database

```bash
microsags sketch references/ -x taxonomy.csv -o database
```

## Documentation

Full documentation is available at
[fuyucheng514-tech.github.io/cellbit](https://fuyucheng514-tech.github.io/cellbit/).

## Citation

A Microsags manuscript and formal citation are in preparation. DNA2bit and
external dependencies should be cited according to their original publications.
