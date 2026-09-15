# Output reference

## Top-level files

| Path | Meaning |
|---|---|
| `INPUT_AUDIT.tsv` | Per-SAG detected type, sources, start stage, skipped stages, assembly length and eligibility |
| `STAGE3B_FASTA_STATS.tsv` | Total bp, maximum contig and GC statistics computed from each assembly |
| `TIMING.tsv` / `TIMING.json` | Monotonic-clock wall time for preflight, sketch, search, Stage 3A and total |
| `COMPLETE.json` | Stable Stage 1–3A terminal receipt; success requires `status: PASS` |

## Stage 1: input preparation

| Path | Contig input | Paired-read input |
|---|---|---|
| `01_assembly/<SAG>/<SAG>.fasta` | Stable link or decompressed view | Link to SPAdes scaffolds |
| `01_assembly/<SAG>/fastp.json` | Not created | fastp machine-readable report |
| `01_assembly/<SAG>/fastp.html` | Not created | fastp human-readable report |
| `01_assembly/<SAG>/spades/` | Not created | SPAdes working directory and scaffolds |
| `01_assembly/excluded_lt1000bp.tsv` | Excluded SAGs | Excluded SAGs |

## Stage 2: DNA2bit annotation

| Path | Meaning |
|---|---|
| `02_dna2bit/bits/<SAG>/` | Per-SAG embedded DNA2bit sketch and receipt |
| `02_dna2bit/search_result.csv` | Raw teacher-compatible packed-search result |
| `02_dna2bit/labels.tsv` | Accepted SAG-to-reference taxonomy and species group |
| `02_dna2bit/SEARCH.PASS` | Search-stage receipt |

## Stage 3A: labelled aggregation

| Path | Meaning |
|---|---|
| `03A_subassemble/groups.tsv` | One row per accepted species group |
| `03A_subassemble/<species>/input_subassemblies.fasta` | Deterministically merged member contigs |
| `03A_subassemble/<species>/run/assembly.fasta` | Final species-guided bin |
| `03A_subassemble/<species>/run/CPP_FULL_PIPELINE_PASS` | Successful C++ subassembly completion marker |
| `03A_subassemble/<species>/run/NO_OVERLAP_PASSTHROUGH.PASS.json` | Audited passthrough receipt for a genuine no-overlap case |

## Stage 3B hand-off

`03B_unclassified_pending.tsv` lists eligible SAGs with no accepted DNA2bit
species label. It is an input to the separate Stage 3B executable, not evidence
that Stage 3B has run.

