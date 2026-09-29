# Performance figure sources

`data/` contains the small measurement tables used for the four panels on the
[performance page](../docs/performance.md). The raw sequence datasets and the
full server-side benchmark logs are not part of the Git repository.

The 24-sample Microsags TSVs are the first two complete 50-thread runs. Each
contains 24 reads and 24 contigs jobs. `plot_real_24.py` checks that every
sample/route appears once, input SAG counts match, and both runs produce the
same annotation SHA-256 before drawing. It averages each run's total time and
maximum process RSS. For the time-only plot it subtracts the historical
7.92376-second-per-sample load calibration; the comparator TSV already holds
the archived startup-adjusted estimates. The load calibration is not a fresh
measurement of the new build.

`sim_microsags_first_two.json` identifies the two included 20-thread
simulation runs and retains the excluded third run's measured totals. Other
simulation TSVs are the archived comparator summaries. The plot includes the
same SPAdes stage for each contigs route; direct reads routes do not include
assembly. The simulation values are recorded full wall times, not
startup-adjusted values. Those archived workflows have different measurement
boundaries, so bar length is not a matched benchmark ranking.

To regenerate the two real-data panels with Python and Matplotlib:

```bash
python benchmarks/plot_real_24.py \
  --round1 benchmarks/data/real_round1.tsv \
  --round2 benchmarks/data/real_round2.tsv \
  --comparator-summary benchmarks/data/real_comparator_adjusted_summary.tsv \
  --comparator-samples benchmarks/data/real_comparator_adjusted_by_sample.tsv \
  --microsags-load-seconds 7.92376 \
  --out docs/figures/performance
```

To regenerate the two simulation panels:

```bash
python benchmarks/plot_simulation.py \
  --svg simulation_gh.svg --png simulation_gh.png
python benchmarks/crop_simulation.py simulation_gh.svg docs/figures/performance
```

The simulation script also needs the `data/` directory and Matplotlib. The
crop script edits only SVG viewBoxes; its output remains native vector text
and paths without embedded raster images.
