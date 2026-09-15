# Usage examples

## Input routes

Microsags reads one tab-separated manifest. Each row represents one SAG. A
single manifest may contain both accepted row forms.

| Input route | Manifest row | Automatic processing |
|---|---|---|
| Paired reads | `SAG_ID<TAB>R1<TAB>R2` | fastp → SPAdes `--sc --careful` → length gate |
| Existing contigs | `SAG_ID<TAB>assembly_fasta` | preserve/decompress FASTA → length gate |

An optional header is accepted:

```text
sag_id	r1	r2
```

or:

```text
sag_id	assembly_fasta
```

Detection uses the decompressed sequence records. File extensions are not used
to decide whether an input is FASTA or FASTQ. Single-end FASTQ, mixed R1/R2
types, truncated records and duplicate normalized SAG identifiers are rejected.

## Command 1: species annotation only

Use this command when you only need a species label for each eligible SAG and
do not want Microsags to construct Stage 3A bins:

```bash
pixi run microsags \
  --manifest SAGs.tsv \
  --out annotation_output \
  --dna-tax /data/GTDB232/genome_taxonomy_1.csv \
  --dna-packed-db /data/GTDB232/dna2bit-packed-index \
  --threads 32 \
  --memory-gb 128 \
  --stop-after annotation
```

The command accepts paired reads, existing contigs, or a mixed manifest. For
reads it first runs fastp and SPAdes; for contigs it skips those two steps.
It stops after the original DNA2bit acceptance rule has produced:

- `annotation_output/02_dna2bit/search_result.csv`: raw search result;
- `annotation_output/02_dna2bit/labels.tsv`: accepted species annotations;
- `annotation_output/03B_unclassified_pending.tsv`: rejected/no-hit SAGs;
- `annotation_output/COMPLETE.json`: annotation-only PASS receipt.

## Command 2: annotation and Stage 3A

```bash
pixi run microsags \
  --manifest SAGs.tsv \
  --out full_output \
  --dna-tax /data/GTDB232/genome_taxonomy_1.csv \
  --dna-packed-db /data/GTDB232/dna2bit-packed-index \
  --threads 32 \
  --memory-gb 128
```

Do not add `--stop-after annotation`. The program performs input preparation,
DNA2bit annotation and species-guided Stage 3A subassembly. The route is
recorded per SAG in `INPUT_AUDIT.tsv`.

## Command 3: continue the same run through Stage 3B

After Command 2 has completed successfully, prepare the required Stage 3B
quality and bac120-marker evidence, then run:

```bash
pixi run sag-stage3b-tractor \
  --manifest full_output/03B_unclassified_pending.tsv \
  --quality-manifest stage3b_quality.tsv \
  --marker-map bac120_marker_nt_map.tsv \
  --ani-engine "$PWD/.pixi/envs/default/bin/gtdb-ani-af" \
  --allow-experimental-ani-engine \
  --leiden-backend "$PWD/python/stage3b_signed_leiden.py" \
  --subass "$PWD/.pixi/envs/default/bin/cpp-subass" \
  --out full_output/03B_unlabelled \
  --threads 32
```

Therefore, the current complete 3A+3B workflow is **Command 2 followed by
Command 3**. Stage 3B is not silently claimed to have run merely because the
pending manifest exists.

## Functional outputs

Microsags v0.1 exposes scientific functions as pipeline stages and files.

### SAG species annotation

Every SAG with total assembly length at least 1,000 bp is sketched using
DNA2bit parameters `k=17`, `bit_len=55296` and `hash_type=0`. Packed search
uses the strict original `min_ratio=0.01` acceptance rule.

The annotation products are:

- `02_dna2bit/search_result.csv`: compatible raw packed-search output;
- `02_dna2bit/labels.tsv`: accepted reference, taxonomy and species group;
- `03B_unclassified_pending.tsv`: eligible SAGs rejected or without a hit.

With `--stop-after annotation`, the command stops here and does not create
Stage 3A group assemblies.

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

Stage 3B is a separate advanced executable. It consumes the Stage 3A hand-off
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
