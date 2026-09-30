# Real-data 3A/3B HQ example

This example is a downloadable subset of one real sample (`PA`) from our
24-sample run. It contains 355 per-SAG contig FASTA files: 55 with a species
label and 300 without one. The reads-derived labels are supplied in
`annotations.tsv`.

## Get the input

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
```

The input is in [`examples/pa_hq_3a_3b/`](https://github.com/fuyucheng514-tech/cellbit/tree/main/examples/pa_hq_3a_3b):

```text
examples/pa_hq_3a_3b/
├── SAG_contigs/           # 355 compressed FASTA files, one per SAG
├── annotations.tsv       # 55 species labels, 300 UNCLASSIFIED
├── SOURCE_INPUTS.tsv      # per-file checksums
├── SOURCE_RUN_QUALITY.tsv # source-server subset result
└── NEW_USER_QUALITY.tsv   # independently reproduced result
```

## Run Assembly mode

First [install Microsags and configure Stage 3B's external programs and
databases](install.md). Then, from the repository root:

```bash
microsags assemble examples/pa_hq_3a_3b/SAG_contigs/ \
  -a examples/pa_hq_3a_3b/annotations.tsv \
  -o pa_hq_example_output -t 64
```

The command writes `stage3a.tsv`, `stage3b_clusters.tsv`, and assembled FASTA
files under `stage3a/` and `stage3b/`. With a separate, newly installed user
environment, Stage 3A produced an HQ `CAG-492_sp000434335` bin (95.24%
completeness, 3.78% contamination). Stage 3B produced four groups, including
one HQ bin (G0002, 92.05%/1.43%). The complete independently reproduced result
is in `NEW_USER_QUALITY.tsv`.

The source-server subset run had a second Stage 3B bin just over the 90%
completeness threshold (90.39%). The new-user run placed the same group just
below it (89.84%). The Stage 3B membership tables were byte-identical, but the
assembled FASTA files differed. Therefore the reproducible target of this
example is **at least one HQ bin on each arm**, not exactly two Stage 3B HQ bins.
The source run is recorded separately in `SOURCE_RUN_QUALITY.tsv`.

This is an **Assembly-mode** example. The bundled annotation table came from
running Annotation mode on the complete sample's FASTQ reads. Those original
reads are not bundled. Running Annotation mode on the packaged contigs changes
which SAGs are labelled, so it is **not** an equivalent way to reproduce the
HQ Assembly result. For an Annotation-mode command using your own reads or
contigs, see [Usage examples](usage.md).

The source and new-user outputs, including assembly FASTAs, logs, and quality
checks, are retained under
`/home/data/fyc/cellbit_114514/result/hq_3a_3b_tutorial_20260930/attempt_001`.
