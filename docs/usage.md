# Usage examples

Microsags provides two operating modes: **Annotation mode** and
**Assembly mode**.

## Annotation mode

Annotation mode assigns a DNA2bit species label to each eligible SAG and then
stops. It does not run Stage 3A or Stage 3B.

### Paired FASTQ reads

```bash
microsags annotate --input-type fastq SAGs/ -d database -o annotation_output
```

Example input with one R1/R2 pair per SAG:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
├── SAG_0001_R2.fastq.gz
├── SAG_0002_R1.fastq.gz
├── SAG_0002_R2.fastq.gz
├── SAG_0003_R1.fastq.gz
└── SAG_0003_R2.fastq.gz
```

FASTQ input is processed with fastp and SPAdes before annotation.

### Singleton FASTQ reads

```bash
microsags annotate --input-type singleton SAGs/ -d database -o annotation_output
```

Example input with one single-end FASTQ file per SAG:

```text
SAGs/
├── SAG_0001.fastq.gz
├── SAG_0002.fastq.gz
└── SAG_0003.fastq.gz
```

Singleton input is processed with fastp and SPAdes single-end mode before
annotation.

### Assembled contigs

```bash
microsags annotate --input-type contigs SAGs/ -d database -o annotation_output
```

Example input with one FASTA file per SAG:

```text
SAGs/
├── SAG_0001.fna
├── SAG_0002.fna
└── SAG_0003.fna
```

A SAG FASTA may contain multiple contigs:

```text
>contig_1
ATGCGTACGTTAGCTAGCTAGCTGACTG...
>contig_2
GCTTACGATCGATCGGATCGATGCA...
```

Contig input skips fastp and SPAdes.

If `--input-type` is omitted, Microsags automatically recognizes standard
FASTQ and FASTA filename extensions:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

Main outputs:

- `02_dna2bit/search_result.csv`: raw DNA2bit search results;
- `02_dna2bit/labels.tsv`: accepted species annotations;
- `03B_unclassified_pending.tsv`: rejected or no-hit SAGs;
- `COMPLETE.json`: completion status.

## Assembly mode

Assembly mode runs DNA2bit annotation followed by Stage 3A assembly of
labelled SAGs and Stage 3B assembly of unlabelled SAGs.

### Paired FASTQ reads

```bash
microsags assemble --input-type fastq SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Example input with one R1/R2 pair per SAG:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
├── SAG_0001_R2.fastq.gz
├── SAG_0002_R1.fastq.gz
├── SAG_0002_R2.fastq.gz
├── SAG_0003_R1.fastq.gz
└── SAG_0003_R2.fastq.gz
```

FASTQ input is processed with fastp and SPAdes before annotation and assembly.

### Singleton FASTQ reads

```bash
microsags assemble --input-type singleton SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Example input with one single-end FASTQ file per SAG:

```text
SAGs/
├── SAG_0001.fastq.gz
├── SAG_0002.fastq.gz
└── SAG_0003.fastq.gz
```

Singleton input is processed with fastp and SPAdes single-end mode before
annotation and assembly.

### Assembled contigs

```bash
microsags assemble --input-type contigs SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

Example input with one FASTA file per SAG:

```text
SAGs/
├── SAG_0001.fna
├── SAG_0002.fna
└── SAG_0003.fna
```

Contig input skips fastp and SPAdes.

Main outputs:

- `02_dna2bit/labels.tsv`: accepted species annotations;
- `03A_subassemble/`: Stage 3A groups and assemblies;
- `stage3b/`: Stage 3B clustering and assemblies;
- `TIMING.tsv`: stage wall-clock times;
- `COMPLETE.json`: final completion status.

## Common options

| Option | Meaning |
|---|---|
| `-i, --input-type` | `auto`, `fastq`, `singleton` or `contigs` |
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
