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
FASTA assemblies. Independent testing with a newly installed user environment
yielded one Stage 3A HQ bin and one Stage 3B HQ bin. Its CheckM2 measurements
are in `NEW_USER_QUALITY.tsv`. The source-server subset run additionally had
a second, borderline Stage 3B HQ bin; those measurements are in
`SOURCE_RUN_QUALITY.tsv`. HQ means completeness at least 90% and contamination
below 5%. Both runs had byte-identical Stage 3B membership tables but different
assembly FASTAs. The borderline bin was 90.39% complete in the source run and
89.84% in the new-user run, so **two Stage 3B HQ bins are not guaranteed**.
Stage 3B group numbers are assigned anew and are not the original whole-sample
group identifiers.

`SOURCE_INPUTS.tsv` records the SHA-256 and byte size of every packaged SAG
FASTA. `RESULT_PASS.json` records the source-server subset run. The full
new-user output and logs, as well as the source output, are retained under
`/home/data/fyc/cellbit_114514/result/hq_3a_3b_tutorial_20260930/attempt_001`
on the source server.
