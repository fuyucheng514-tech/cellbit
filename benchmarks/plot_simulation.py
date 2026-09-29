#!/usr/bin/env python3
"""Plot all simulated g/h input routes from preserved benchmark summaries.

The archived SPAdes stage is exactly the stage that made the 5,000 contig
inputs used by both skani and Microsags. Preserve all original measurements;
do not reinterpret sums of per-SAG runtimes as a single run's wall time.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator


def one_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected one measured row: {path}")
    return rows[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "data")
    parser.add_argument("--optimized-means", type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--png", type=Path)
    args = parser.parse_args()

    data = args.data_dir
    if args.optimized_means is None:
        args.optimized_means = data / "sim_microsags_first_two.json"
    contigs = one_row(data / "sim_archived_microsags_contigs.tsv")
    reads = one_row(data / "sim_archived_microsags_reads.tsv")
    spades = one_row(data / "sim_spades_skani.tsv")
    sourmash = one_row(data / "sim_sourmash_reads.tsv")
    mash = one_row(data / "sim_mash_reads.tsv")
    alternatives = one_row(data / "sim_mash_sourmash_contigs.tsv")

    assert contigs["method"] == reads["method"] == "Microsags v0.5.2"
    assert contigs["input"] == "SPAdes contigs"
    assert reads["input"] == "paired FASTQ reads"
    assert int(contigs["n_sags"]) == int(reads["n_sags"]) == 5000
    assert int(reads["n_exact_truth"]) == int(reads["n_classified"]) == 5000
    assert int(contigs["threads"]) == int(reads["threads_requested"]) == 20
    assert int(contigs["reference_count"]) == int(reads["reference_count"]) == 113104
    assert all(int(row["n_total_samples"]) == 5000 for row in (spades, sourmash, mash))
    assert int(alternatives["n_sags"]) == 5000 and int(alternatives["threads"]) == 20

    spades_seconds = float(spades["spades_total_time_s"])
    spades_rss = float(spades["spades_max_rss_kb"])
    seconds = [
        spades_seconds + float(spades["skani_total_time_s"]),
        spades_seconds + float(alternatives["sourmash_total_recorded_s"]),
        spades_seconds + float(alternatives["mash_total_recorded_s"]),
        spades_seconds + float(contigs["elapsed_s"]),
        float(sourmash["sourmash_total_time_s"]),
        float(mash["mash_total_time_s"]),
        float(reads["elapsed_s"]),
    ]
    rss_gib = [
        max(spades_rss, float(spades["skani_max_rss_kb"])) / 1048576,
        max(spades_rss, float(alternatives["sourmash_stage_max_rss_kb"])) / 1048576,
        max(spades_rss, float(alternatives["mash_stage_max_rss_kb"])) / 1048576,
        max(spades_rss, float(contigs["max_rss_kib"])) / 1048576,
        max(float(sourmash["sourmash_sketch_max_rss_kb"]), float(sourmash["sourmash_search_max_rss_kb"])) / 1048576,
        max(float(mash["mash_sketch_max_rss_kb"]), float(mash["mash_dist_max_rss_kb"])) / 1048576,
        float(reads["max_rss_kib"]) / 1048576,
    ]
    labels = [
        "SPAdes + skani", "SPAdes + sourmash", "SPAdes + Mash",
        "SPAdes + Microsags", "sourmash (reads)", "Mash (reads)",
        "Microsags (reads)",
    ]
    colors = ["#abc5d7", "#789ab0", "#e87359", "#008c88", "#789ab0", "#e87359", "#d59635"]

    if args.optimized_means:
        record = json.loads(args.optimized_means.read_text(encoding="utf-8"))
        runs = record["included_runs"]
        assert len(runs) == 2 and record["threads_requested"] == 20
        assert [run["attempt"] for run in runs] == [
            "attempt_001_candidate3/full5000", "attempt_002_retest_20threads"
        ]
        means = {
            route: {
                metric: sum(float(run[route][metric]) for run in runs) / 2
                for metric in ("wall_seconds", "max_rss_kb")
            }
            for route in ("reads", "contigs")
        }
        seconds[3] = spades_seconds + means["contigs"]["wall_seconds"]
        seconds[6] = means["reads"]["wall_seconds"]
        rss_gib[3] = max(spades_rss, means["contigs"]["max_rss_kb"]) / 1048576
        rss_gib[6] = means["reads"]["max_rss_kb"] / 1048576
        labels[3] = "SPAdes + Microsags"
        labels[6] = "Microsags (reads)"

    plt.rcParams.update({
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "hatch.linewidth": 0.85,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(17.0, 6.4), dpi=160)
    fig.patch.set_facecolor("#ffffff")
    for ax in axes:
        ax.set_facecolor("#fff2e7")
        ax.spines["left"].set_color("#9ca8ae")
        ax.spines["bottom"].set_color("#9ca8ae")
        ax.tick_params(axis="both", length=3, color="#9ca8ae")
        ax.invert_yaxis()

    times_h = [sec / 3600 for sec in seconds]
    time_order = sorted(range(len(labels)), key=lambda i: (times_h[i], labels[i]))
    memory_order = sorted(range(len(labels)), key=lambda i: (rss_gib[i], labels[i]))
    axes[0].barh(
        [labels[i] for i in time_order],
        [times_h[i] for i in time_order],
        color=[colors[i] for i in time_order],
        height=0.62,
    )
    spades_hours = spades_seconds / 3600
    for row, source_idx in enumerate(time_order):
        if source_idx < 4:
            axes[0].barh(
                row,
                spades_hours,
                height=0.62,
                facecolor="none",
                edgecolor=(0.18, 0.25, 0.30, 0.60),
                linewidth=0.25,
                hatch="///",
                zorder=3,
            )
    axes[0].legend(
        handles=[Patch(facecolor="white", edgecolor="#52616a", hatch="///",
                       label=f"SPAdes ({spades_hours:.3f} h)")],
        loc="upper right",
        frameon=False,
        fontsize=9.5,
    )
    axes[0].set_xlim(0, max(times_h) + 6.5)
    axes[0].xaxis.set_major_locator(MultipleLocator(5))
    axes[0].set_xlabel("Elapsed time (hours)")
    axes[0].set_title("g   Total runtime", loc="left", fontweight="bold", pad=12)
    for idx, source_idx in enumerate(time_order):
        seconds_i = seconds[source_idx]
        hours_i = times_h[source_idx]
        if source_idx < 4:
            software_minutes = (seconds_i - spades_seconds) / 60
            label = f"+{software_minutes:.3f} min"
        else:
            label = f"{hours_i:.3f} h" if hours_i >= 1 else f"{seconds_i / 60:.3f} min"
        axes[0].text(hours_i + 0.20, idx, label, va="center", ha="left", fontsize=9.2)

    axes[1].barh(
        [labels[i] for i in memory_order],
        [rss_gib[i] for i in memory_order],
        color=[colors[i] for i in memory_order],
        height=0.62,
    )
    axes[1].set_xlim(0, 66)
    axes[1].xaxis.set_major_locator(MultipleLocator(15))
    axes[1].set_xlabel("Peak process RSS (GiB)")
    axes[1].set_title("h   Peak memory", loc="left", fontweight="bold", pad=12)
    for idx, source_idx in enumerate(memory_order):
        memory = rss_gib[source_idx]
        label = f"{memory:.2f} GiB" if memory < 5 else f"{memory:.1f} GiB"
        if memory > 60:
            axes[1].text(memory - 0.72, idx, label, va="center", ha="right", fontsize=10)
        else:
            axes[1].text(memory + 0.72, idx, label, va="center", ha="left", fontsize=10)

    fig.suptitle("5,000 simulated SAGs · 20 threads/workers", fontsize=15, fontweight="bold", y=0.99)
    fig.text(
        0.055, 0.048,
        "All contig routes include the same archived SPAdes stage; reads routes use FASTQ directly."
        + (" Microsags: first two runs only; third high-load run excluded." if args.optimized_means else ""),
        fontsize=8.4, color="#4f5b61",
    )
    fig.text(
        0.055, 0.026,
        "SPAdes and archived comparators sum per-SAG/stage times; Microsags annotation uses one batch wall-time. Not a matched speed ranking.",
        fontsize=8.4, color="#4f5b61",
    )
    fig.subplots_adjust(left=0.14, right=0.98, top=0.84, bottom=0.19, wspace=0.90)
    args.svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.svg, format="svg", bbox_inches="tight")
    if args.png:
        args.png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.png, format="png", dpi=180, bbox_inches="tight")
    print(f"SVG={args.svg}")
    print(f"reads_seconds={seconds[-1]:.3f} reads_peak_gib={rss_gib[-1]:.6f}")
    print(f"spades_plus_microsags_seconds={seconds[3]:.3f} pipeline_peak_gib={rss_gib[3]:.6f}")


if __name__ == "__main__":
    main()
