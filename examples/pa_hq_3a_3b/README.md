# PA Stage 3A / Stage 3B example

This is a 355-SAG subset of a real 24-sample run: 55 SAGs labelled
`CAG-492_sp000434335` and 300 SAGs without a reads-level DNA2bit label.
It is an **Assembly-mode example**, not a synthetic dataset. The
`annotations.tsv` file is the actual reads-level annotation hand-off from the
complete PA sample. Original FASTQ reads are not in this package; annotating
these contigs again does **not** reproduce that hand-off.

After installing Microsags and configuring the external CheckM2 and GTDB-Tk
databases, run from the repository root:

```bash
microsags assemble examples/pa_hq_3a_3b/SAG_contigs/ \
  -a examples/pa_hq_3a_3b/annotations.tsv \
  -o pa_hq_example_output -t 64
```

The command writes `stage3a.tsv`, `stage3b_clusters.tsv` and the corresponding
FASTA assemblies. The tested subset yielded one Stage 3A HQ bin and two Stage
3B HQ bins. CheckM2 measurements from the independent subset run are in
`EXPECTED_QUALITY.tsv`; HQ means completeness at least 90% and contamination
below 5%. Stage 3B group numbers are assigned anew and are not the original
whole-sample group identifiers.

`SOURCE_INPUTS.tsv` records the SHA-256 and byte size of every packaged SAG
FASTA. `RESULT_PASS.json` records the source run and the independent subset
verification. The full output and logs are retained under
`/home/data/fyc/cellbit_114514/result/hq_3a_3b_tutorial_20260930/attempt_001`
on the source server.
