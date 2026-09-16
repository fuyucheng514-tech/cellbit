# Stage 3B advanced workflow

Stage 3B clusters DNA2bit-negative SAGs. In v0.2 the public `assemble` command
automatically runs its upstream evidence preparation, graph clustering and
subassembly chain. The page documents the advanced scientific boundary.

## Required inputs

- `03B_unclassified_pending.tsv` from the stable Microsags run;
- a quality manifest containing `sag_id`, `max_contig`, `gc_pct` and
  `checkm2_contamination`;
- a bac120 nucleotide marker map generated from the same SAG assemblies;
- explicit paths to the ANI engine, Leiden backend and C++ subassembly driver.

## Main products

- quality-gate audit and graph-node manifest;
- all-pairs ANI/AF statistics and positive edges;
- marker-derived negative-edge evidence;
- signed Leiden membership and unaggregated SAG tables;
- per-cluster merged FASTA, subassembly and PASS receipts;
- a Stage 3B `COMPLETE.json` closure receipt.

See [`README_STAGE3B.md`](https://github.com/fuyucheng514-tech/cellbit/blob/main/README_STAGE3B.md)
for the frozen command contract and current validation boundary.
