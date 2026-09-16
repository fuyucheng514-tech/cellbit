# Microsags

Microsags is a command-line tool for species annotation and assembly of
single-amplified genomes (SAGs). It accepts assembled FASTA files or paired
FASTQ reads and automatically selects the correct route.

## Main commands

### Species annotation

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

This produces DNA2bit species annotations and stops before SAG aggregation.

### Complete assembly

```bash
microsags assemble SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

This runs annotation followed by Stage 3A and Stage 3B.

### Build a database

```bash
microsags sketch references/ -x taxonomy.csv -o database
```

## Input

For assembled SAGs, put one `.fa`, `.fasta` or `.fna` file per SAG in
`SAGs/`. For reads, put paired `_R1` and `_R2` FASTQ files in the directory.
Gzip-compressed files are supported.

Use `-i contigs` or `-i reads` to select the route explicitly. If omitted,
Microsags detects the input type from standard filenames and validates the
sequence contents.

## Output

- `02_dna2bit/labels.tsv`: species annotations;
- `03A_subassemble/`: Stage 3A assemblies;
- `stage3b/`: Stage 3B results;
- `TIMING.tsv`: stage timings;
- `COMPLETE.json`: completion status.

Continue with [Install](install.md) or [Usage](usage.md).
