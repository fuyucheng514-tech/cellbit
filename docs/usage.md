# Usage

## Annotate SAGs

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

Key output: `annotation_output/02_dna2bit/labels.tsv`.

## Assemble SAGs

```bash
microsags assemble SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

`assemble` runs the complete workflow: input preparation, DNA2bit annotation,
Stage 3A and Stage 3B.

## Build a database

```bash
microsags sketch references/ -x taxonomy.csv -o database
```

## Input examples

Paired reads:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
└── SAG_0001_R2.fastq.gz
```

Assembled SAGs:

```text
SAGs/
└── SAG_0001.fna
```

## Common options

| Option | Meaning |
|---|---|
| `-d, --database` | DNA2bit database directory |
| `-o, --output` | Output directory |
| `-t, --threads` | Number of threads |
| `-i, --input-type` | `auto`, `reads` or `contigs` |
| `-l, --file-list` | Text file containing one input path per line |
| `-h, --help` | Show command help |

Run `microsags COMMAND --help` for every available option.
