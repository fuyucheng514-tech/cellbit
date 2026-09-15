# Lake quarter tutorial

This integration example applies the public Microsags v0.1 release to a
deterministic quarter of the 13,742-SAG lake contig collection. It tests
realistic scale and output contracts. The underlying biological data are not
redistributed through GitHub.

## Dataset selection

The authoritative contig manifest contains 13,742 unique SAGs. The example:

1. sorts normalized SAG identifiers lexicographically;
2. selects ranks 1, 5, 9, ...;
3. verifies every selected FASTA exists and is non-empty;
4. records the manifest and receipt checksums.

This selects 3,436 SAGs, totalling 718,690,933 input bytes.

| Evidence | Value |
|---|---|
| Source SAG count | 13,742 |
| Selected SAG count | 3,436 |
| Unique selected IDs | 3,436 |
| Missing/empty FASTA inputs | 0 |
| Source manifest SHA256 | `ef1376f9d4b7ce2c8eee1bea2f2962d6bcf65a74c8c91e73681dc57830ea0f07` |
| Quarter manifest SHA256 | `3c2bb3387e70b22ff4cde5a867fbdff3339960a7dd5146080167f1f6034c5636` |

The repository includes the deterministic selector in
`examples/lake_quarter/select_quarter.py`. Run it against a local copy of the
authoritative manifest:

```bash
python examples/lake_quarter/select_quarter.py \
  --input LAKE_FULL_CONTIGS.tsv \
  --output lake_contigs_quarter.tsv \
  --receipt lake_contigs_quarter.receipt.json
```

## Run the contig route

```bash
pixi run microsags \
  --manifest lake_contigs_quarter.tsv \
  --out lake_quarter_output \
  --dna-tax /data/GTDB232/genome_taxonomy_1.csv \
  --dna-packed-db /data/GTDB232/dna2bit-packed-index \
  --threads 128 \
  --memory-gb 1024
```

Because every selected input is FASTA, the expected route is:

```text
FASTA content detection
  → preserve existing contigs
  → total assembly length >= 1,000 bp
  → embedded DNA2bit sketch and packed search
  → accepted labels to Stage 3A
  → rejected/no-hit SAGs to the Stage 3B hand-off manifest
```

The completed integration run is retained under a write-once server directory:

```text
/home/data/fyc/cellbit_114514/temp/
  Microsags_lake_quarter_example_20260915/attempt_001/
```

## Verified result

The public v0.1 Release completed with exit status 0. `COMPLETE.json`, all 156
Stage 3A group receipts and all 156 non-empty bin FASTA files closed
successfully.

| Result | Value |
|---|---:|
| Input SAGs | 3,436 |
| Detected contig inputs | 3,436 |
| Detected paired-read inputs | 0 |
| Eligible SAGs | 3,436 |
| Excluded below 1,000 bp | 0 |
| DNA2bit-labelled SAGs | 1,873 |
| DNA2bit-negative/no-hit SAGs | 1,563 |
| Stage 3A species groups | 156 |
| Completed Stage 3A bins | 156 |
| Output footprint | 1.3 GiB |

The closure equation is `1,873 labelled + 1,563 pending = 3,436 eligible`.

### Wall-clock timing

| Phase | Seconds |
|---|---:|
| Input preflight and length gate | 37.495 |
| DNA2bit sketch | 2.805 |
| DNA2bit packed search | 30.206 |
| Stage 3A assembly phase | 50.027 |
| Stage 3A finish phase | 203.326 |
| Stage 3A total | 253.353 |
| End-to-end total | 323.859 |

The operating-system measurement was 5 minutes 23.91 seconds wall time,
3,605.07 CPU seconds and 1,474,628 KiB peak resident memory. This is a warm,
uncontrolled shared-filesystem measurement, not a strict cold-cache benchmark.

## Inspect the result

```bash
python -m json.tool lake_quarter_output/COMPLETE.json
column -t -s $'\t' lake_quarter_output/TIMING.tsv
head lake_quarter_output/INPUT_AUDIT.tsv
head lake_quarter_output/02_dna2bit/labels.tsv
head lake_quarter_output/03A_subassemble/groups.tsv
head lake_quarter_output/03B_unclassified_pending.tsv
```

The example validates the contig entrance, annotation contract and Stage 3A.
It does not claim Stage 3B completion merely because a pending manifest exists.
