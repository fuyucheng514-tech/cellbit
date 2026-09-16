# Usage examples

Microsags has two consecutive modes. Run **Annotation mode** first. If genome
assembly is required, prepare one contig FASTA per SAG and then run
**Assembly mode** with the completed annotation directory.

## Annotation mode

Annotation mode assigns a DNA2bit species label to each SAG. It accepts paired
FASTQ reads, singleton FASTQ reads, or assembled contigs. It does not run
Stage 3A, Stage 3B or subassembly.

For FASTQ input, Microsags runs fastp and sends all cleaned reads from each SAG
directly to DNA2bit. It does not run SPAdes.

### Paired FASTQ reads

```bash
microsags annotate --input-type fastq SAGs/ -d database -o annotation_output
```

```text
SAGs/
├── SAG_0001_R1.fastq.gz
├── SAG_0001_R2.fastq.gz
├── SAG_0002_R1.fastq.gz
├── SAG_0002_R2.fastq.gz
├── SAG_0003_R1.fastq.gz
├── SAG_0003_R2.fastq.gz
└── ...
```

### Singleton FASTQ reads

```bash
microsags annotate --input-type singleton SAGs/ -d database -o annotation_output
```

```text
SAGs/
├── SAG_0001.fastq.gz
├── SAG_0002.fastq.gz
├── SAG_0003.fastq.gz
└── ...
```

### Assembled contigs

```bash
microsags annotate --input-type contigs SAGs/ -d database -o annotation_output
```

```text
SAGs/
├── SAG_0001.fna
├── SAG_0002.fna
├── SAG_0003.fna
└── ...
```

The principal result is
`annotation_output/02_dna2bit/labels.tsv`.

Example result (illustrative species names):

```text
sag_id      species
SAG_0001    Example_species_A
SAG_0002    Example_species_B
SAG_0003    Example_species_C
...         ...
```

## Prepare contigs for Assembly mode

Assembly mode accepts contigs only. If Annotation mode used FASTQ reads, first
assemble each SAG independently with SPAdes or another short-read assembler.
Keep the same SAG identifiers:

```text
Annotation SAG ID          Assembly input
SAG_0001             →     SAG_0001.fna
SAG_0002             →     SAG_0002.fna
SAG_0003             →     SAG_0003.fna
```

Microsags requires the contig SAG set to match the union of labelled and
unclassified SAGs in `annotation_output`. A missing, extra or duplicate SAG
causes a hard stop.

## Assembly mode

Assembly mode reads prior species annotations and per-SAG contigs. It does not
run fastp, SPAdes or DNA2bit again.

```bash
microsags assemble --input-type contigs SAG_contigs/ \
  --annotations annotation_output \
  -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

```text
SAG_contigs/
├── SAG_0001.fna
├── SAG_0002.fna
├── SAG_0003.fna
└── ...
```

Labelled SAGs are grouped by their DNA2bit species label and processed by
Stage 3A. Unclassified SAGs are processed by the Stage 3B evidence and
graph-clustering workflow. Both routes retain the existing
`cpp-subass`/Flye subassembly logic.

The Stage 3B Leiden resolution is optional and adjustable:

```bash
microsags assemble --input-type contigs SAG_contigs/ \
  --annotations annotation_output \
  -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database \
  --leiden-resolution 0.18
```

When this option is omitted, Microsags uses the frozen default Leiden sweep.
All other clustering rules remain unchanged.

Example output:

```text
assembly_output/
├── 03A_subassemble/
│   ├── Example_species_A/run/assembly.fasta
│   ├── Example_species_B/run/assembly.fasta
│   └── ...
├── stage3b/
│   └── result/
│       ├── 06_subassemble/cluster_00001/run/assembly.fasta
│       └── ...
└── COMPLETE.json
```

## Common options

| Option | Mode | Meaning |
|---|---|---|
| `-i, --input-type` | Annotation | `auto`, `fastq`, `singleton` or `contigs` |
| `-i, --input-type` | Assembly | `contigs` only |
| `-d, --database` | Annotation | DNA2bit database directory |
| `-a, --annotations` | Assembly | Completed Annotation-mode output directory |
| `-o, --output` | Both | New write-once output directory |
| `-t, --threads` | Both | Number of worker threads |
| `--leiden-resolution` | Assembly | Optional Stage 3B Leiden resolution |
| `-h, --help` | Both | Show command help |
