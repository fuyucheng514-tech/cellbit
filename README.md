# Microsags

Microsags is a command-line tool for species annotation and assembly of
single-amplified genomes (SAGs). It accepts assembled FASTA files or paired
FASTQ reads, assigns species with DNA2bit, and can run the complete Stage 3A
and Stage 3B assembly workflow.

## Features

- Accepts a directory, individual files, or a file list.
- Automatically detects assembled contigs and paired reads.
- Runs species annotation alone with `annotate`.
- Runs annotation, Stage 3A and Stage 3B with `assemble`.
- Builds a compatible DNA2bit database with `sketch`.
- Writes intermediate files, timing records and a completion receipt.

## Install

Microsags supports Linux x86-64. Conda or Mamba is recommended.

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

The scientific database is distributed separately from the source code. See
the [installation guide](https://fuyucheng514-tech.github.io/cellbit/install/)
for the download command.

## Quick start

Annotate SAGs:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

Run the complete Stage 3A and Stage 3B workflow:

```bash
microsags assemble SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Build a DNA2bit database:

```bash
microsags sketch references/ -x taxonomy.csv -o database
```

Use `-t` to set the number of threads and `-i fastq`, `-i singleton` or
`-i contigs` to force an input type. Without `-i`, Microsags recognizes `.fq`, `.fastq`, `.fa`,
`.fasta` and `.fna`, including gzip-compressed files.

## Input

For reads, provide one R1/R2 pair per SAG:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
└── SAG_0001_R2.fastq.gz
```

For assembled data, place one FASTA file per SAG in a directory:

```text
SAGs/
└── SAG_0001.fna
```

## Output

- `02_dna2bit/labels.tsv`: accepted species annotations;
- `03A_subassemble/`: species-guided Stage 3A assemblies;
- `stage3b/`: DNA2bit-negative Stage 3B results;
- `TIMING.tsv`: stage wall-clock times;
- `COMPLETE.json`: final completion status.

## Documentation

Documentation and examples are available at
[fuyucheng514-tech.github.io/cellbit](https://fuyucheng514-tech.github.io/cellbit/).

## Citation

A Microsags manuscript and formal citation are in preparation. DNA2bit and
external dependencies should be cited according to their original publications.
