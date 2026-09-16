# Usage examples

Microsags has two operating modes:

1. **Annotation mode** assigns species labels to SAGs.
2. **Assembly mode** runs annotation followed by Stage 3A and Stage 3B.

Both modes accept paired FASTQ reads or assembled FASTA contigs.

## Input

### Paired FASTQ reads

Place one R1/R2 pair per SAG in the input directory:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
└── SAG_0001_R2.fastq.gz
```

Use:

```bash
--input-type fastq
```

FASTQ input is processed with fastp and SPAdes before annotation.

### Assembled contigs

Place one FASTA file per SAG in the input directory:

```text
SAGs/
└── SAG_0001.fna
```

A SAG FASTA may contain multiple contigs:

```text
>contig_1
ATGCGTACGTTAGCTAGCTAGCTGACTG...
>contig_2
GCTTACGATCGATCGGATCGATGCA...
```

Use:

```bash
--input-type contigs
```

Contig input skips fastp and SPAdes.

If `--input-type` is omitted, Microsags automatically recognizes
`.fq`, `.fastq`, `.fa`, `.fasta` and `.fna`, including gzip-compressed
files.

## Modes

### Annotation mode

Annotation mode assigns a DNA2bit species label to each eligible SAG and then
stops. It does not run Stage 3A or Stage 3B.

Paired FASTQ reads:

```bash
microsags annotate --input-type fastq SAGs/ -d database -o annotation_output
```

Assembled contigs:

```bash
microsags annotate --input-type contigs SAGs/ -d database -o annotation_output
```

Automatic input detection:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

Main outputs:

- `02_dna2bit/search_result.csv`: raw DNA2bit search results;
- `02_dna2bit/labels.tsv`: accepted species annotations;
- `03B_unclassified_pending.tsv`: rejected or no-hit SAGs;
- `COMPLETE.json`: completion status.

### Assembly mode

Assembly mode runs the complete workflow:

```text
input preparation
→ DNA2bit annotation
→ Stage 3A assembly of labelled SAGs
→ Stage 3B assembly of unlabelled SAGs
```

Paired FASTQ reads:

```bash
microsags assemble --input-type fastq SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Assembled contigs:

```bash
microsags assemble --input-type contigs SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Main outputs:

- `02_dna2bit/labels.tsv`: accepted species annotations;
- `03A_subassemble/`: Stage 3A groups and assemblies;
- `stage3b/`: Stage 3B clustering and assemblies;
- `TIMING.tsv`: stage wall-clock times;
- `COMPLETE.json`: final completion status.

## Common options

| Option | Meaning |
|---|---|
| `-i, --input-type` | `auto`, `fastq` or `contigs` |
| `-d, --database` | DNA2bit database directory |
| `-o, --output` | New output directory |
| `-t, --threads` | Number of worker threads |
| `-l, --file-list` | File containing one input path per line |
| `-h, --help` | Show command help |

## Completion check

A completed run writes `COMPLETE.json`. Check it with:

```bash
python -m json.tool assembly_output/COMPLETE.json
```
