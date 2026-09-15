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

## Default Stage 1–3A workflow

```bash
pixi run microsags \
  --manifest SAGs.tsv \
  --out microsags_output \
  --dna-tax /data/GTDB232/genome_taxonomy_1.csv \
  --dna-packed-db /data/GTDB232/dna2bit-packed-index \
  --threads 32 \
  --memory-gb 128
```

The same command accepts an all-contigs, all-reads or mixed manifest. The route
is recorded per SAG in `INPUT_AUDIT.tsv`.

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

The v0.1 main command continues from annotation into Stage 3A. It does not yet
provide an annotation-only stop switch; users who only need annotations can
consume the files above and ignore the Stage 3A directory.

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
marker table:

```bash
pixi run sag-stage3b-tractor \
  --manifest microsags_output/03B_unclassified_pending.tsv \
  --quality-manifest stage3b_quality.tsv \
  --marker-map bac120_marker_nt_map.tsv \
  --ani-engine "$PWD/.pixi/envs/default/bin/gtdb-ani-af" \
  --allow-experimental-ani-engine \
  --leiden-backend "$PWD/python/stage3b_signed_leiden.py" \
  --subass "$PWD/.pixi/envs/default/bin/cpp-subass" \
  --out microsags_output/03B_unlabelled \
  --threads 32
```

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

