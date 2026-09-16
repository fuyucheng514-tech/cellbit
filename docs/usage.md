# Usage examples

## Input routes

The public CLI accepts sequence paths, a directory, or a path list and creates
the internal tab-separated manifest automatically. Each internal row represents
one SAG and has one of the following forms.

| Input route | Manifest row | Automatic processing |
|---|---|---|
| Paired reads | `SAG_ID<TAB>R1<TAB>R2` | fastp → SPAdes `--sc --careful` → length gate |
| Existing contigs | `SAG_ID<TAB>assembly_fasta` | preserve/decompress FASTA → length gate |

### Paired FASTQ reads

`SAGs/` can contain one R1/R2 pair per SAG:

```text
SAGs/
├── SAG_0001_R1.fastq.gz
└── SAG_0001_R2.fastq.gz
```

```bash
microsags annotate --input-type fastq SAGs/ -d database -o annotation_output
```

### Assembled contigs

Place one FASTA file per SAG:

```text
SAGs/
└── SAG_0001.fna
```

The contents of `SAG_0001.fna` look like this (sequences shortened for
display):

```text
>contig_1
ATGCGTACGTTAGCTAGCTAGCTGACTG...
>contig_2
GCTTACGATCGATCGGATCGATGCA...
```

```bash
microsags annotate --input-type contigs SAGs/ -d database -o annotation_output
```

An optional header is accepted:

```text
sag_id	r1	r2
```

or:

```text
sag_id	assembly_fasta
```

Use `-i/--input-type fastq` or `-i/--input-type contigs` to select a route
explicitly. If the option is omitted, `--input-type auto` recognizes
`.fq`, `.fastq`, `.fa`, `.fasta`, and `.fna`, with optional `.gz`. FASTQ files
supplied through a directory or one combined `--file-list` are paired by their
`R1`/`R2` filename suffixes. The C++ reader then validates the decompressed
records, so misleading suffixes, single-end FASTQ, truncated records and
duplicate normalized SAG identifiers fail closed.

```bash
# Explicit paired-read route; a directory may contain many R1/R2 pairs
microsags assemble --input-type fastq SAG_reads/ -d database -o full_output

# Explicit existing-contig route
microsags annotate --input-type contigs SAGs/ -d database -o annotations

# Automatic route selection from standard suffixes
microsags assemble SAGs/ -d database -o full_output
```

## Command 1: species annotation only

Use this command when you only need a species label for each eligible SAG and
do not want Microsags to construct Stage 3A bins:

Paired reads:

```bash
microsags annotate --input-type fastq SAGs/ -d database -o annotation_output
```

Assembled contigs:

```bash
microsags annotate --input-type contigs SAGs/ -d database -o annotation_output
```

The command accepts either paired reads or existing contigs in one invocation.
For reads it first runs fastp and SPAdes; for contigs it skips those two steps.
It stops after the original DNA2bit acceptance rule has produced:

- `annotation_output/02_dna2bit/search_result.csv`: raw search result;
- `annotation_output/02_dna2bit/labels.tsv`: accepted species annotations;
- `annotation_output/03B_unclassified_pending.tsv`: rejected/no-hit SAGs;
- `annotation_output/COMPLETE.json`: annotation-only PASS receipt.

## Command 2: annotation plus Stage 3A and Stage 3B

```bash
microsags assemble SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```

This single command performs input preparation, DNA2bit annotation,
species-guided Stage 3A and DNA2bit-negative Stage 3B. The route is recorded
per SAG in `INPUT_AUDIT.tsv`. CheckM2 and GTDB-Tk executables are found on
`PATH`; their database paths are explicit because those databases are not
stored in the source repository.

## Functional outputs

Microsags v0.2 exposes scientific functions through `annotate`, `assemble`
and `sketch`.

### SAG species annotation

Every SAG with total assembly length at least 1,000 bp is sketched using
DNA2bit parameters `k=17`, `bit_len=55296` and `hash_type=0`. Packed search
uses the strict original `min_ratio=0.01` acceptance rule.

The annotation products are:

- `02_dna2bit/search_result.csv`: compatible raw packed-search output;
- `02_dna2bit/labels.tsv`: accepted reference, taxonomy and species group;
- `03B_unclassified_pending.tsv`: eligible SAGs rejected or without a hit.

The `annotate` subcommand stops here and does not create Stage 3A assemblies.

### Stage 3A: species-guided subassembly

Accepted SAGs are grouped by their DNA2bit species label. Member contigs are
merged with SAG-prefixed record identifiers and passed to the bundled C++
subassembly driver backed by Flye subassembly modules.

Principal outputs:

- `03A_subassemble/groups.tsv`: species group, SAG count, input and bin path;
- `03A_subassemble/<species>/input_subassemblies.fasta`: deterministic merged
  member sequences;
- `03A_subassemble/<species>/run/assembly.fasta`: resulting species bin;
- per-group PASS receipts and no-overlap passthrough evidence where applicable.

### Stage 3B: DNA2bit-negative SAG clustering

The public `assemble` command invokes the Stage 3B engine automatically. It consumes the Stage 3A hand-off
manifest together with a CheckM2-derived quality table and a bac120 nucleotide
marker table. Use Command 3 above.

!!! warning
    The bundled Stage 3B ANI/AF path is experimental and is not a skani clone.
    Stage 3B is documented separately from the stable one-command Stage 1–3A
    workflow. Missing quality or marker evidence causes a hard stop.

## Completion check

A successful stable run ends with `COMPLETE.json` containing `status: PASS`.
Intermediate file counts alone do not establish completion.

```bash
python -m json.tool microsags_output/COMPLETE.json
cat microsags_output/TIMING.tsv
```
