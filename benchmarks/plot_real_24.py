#!/usr/bin/env python3
"""Recreate the 24-sample time and memory panels from preserved TSV measurements."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SAMPLES = set("FC GB HA IB JA KA LA LB MA OA OB PA PB SB WC WK WL XA XB XG XH ZA ZB ZF".split())
ROUTES = ("reads", "contigs")
COLORS = {
    "Microsags (reads)": "#d59635",
    "Microsags (contigs)": "#008c88",
    "skani (contigs)": "#abc5d7",
    "Mash (contigs)": "#e87359",
    "sourmash (contigs)": "#789ab0",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def two_run_metrics(paths: list[Path]) -> dict[str, dict[str, float]]:
    if len(paths) != 2:
        raise ValueError("Exactly the first two independent runs are required")
    rounds = []
    for path in paths:
        rows = read_tsv(path)
        by_route: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_route[row["route"]].append(row)
        if len(rows) != 48 or set(by_route) != set(ROUTES):
            raise ValueError(f"Incomplete 24-sample run: {path}")
        for route in ROUTES:
            part = by_route[route]
            if len(part) != 24 or {row["sample"] for row in part} != SAMPLES:
                raise ValueError(f"Missing or duplicate {route} samples: {path}")
            if any(float(row["wall_seconds"]) <= 0 or int(row["max_rss_kb"]) <= 0 for row in part):
                raise ValueError(f"Invalid measurement: {path}")
        rounds.append(by_route)
    for route in ROUTES:
        for sample in SAMPLES:
            first = next(row for row in rounds[0][route] if row["sample"] == sample)
            second = next(row for row in rounds[1][route] if row["sample"] == sample)
            if first["output_sha256"] != second["output_sha256"]:
                raise ValueError(f"Changed annotation output: {route}/{sample}")
            if first["input_sags"] != second["input_sags"]:
                raise ValueError(f"Changed SAG count: {route}/{sample}")
    return {
        route: {
            "mean_total_wall_seconds": mean(sum(float(row["wall_seconds"]) for row in run[route]) for run in rounds),
            "mean_run_peak_rss_gib": mean(max(int(row["max_rss_kb"]) for row in run[route]) for run in rounds) / 1048576,
        }
        for route in ROUTES
    }


def plot_panel(rows: list[tuple[str, float | None, str]], out: Path, title: str,
               xlabel: str, xmax: float, note: str) -> None:
    plt.rcParams.update({"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 12,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9.6, 5.0), dpi=160)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#fff2e7")
    ax.spines["left"].set_color("#9ca8ae")
    ax.spines["bottom"].set_color("#9ca8ae")
    for index, (name, value, label) in enumerate(rows):
        if value is not None:
            ax.barh(index, value, color=COLORS[name], height=0.62)
            ax.text(value + xmax * 0.012, index, label, va="center", ha="left", fontsize=11)
        else:
            ax.text(xmax * 0.012, index, label, va="center", ha="left", fontsize=11)
    ax.set_yticks(range(len(rows)), [row[0] for row in rows])
    ax.invert_yaxis()
    ax.set_ylim(len(rows) - 0.05, -0.75)
    ax.set_xlim(0, xmax)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold", pad=14)
    fig.subplots_adjust(left=0.30, right=0.95, top=0.86, bottom=0.24)
    fig.text(0.04, 0.075, note, fontsize=8.6, color="#4f5b61")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="svg", bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--round1", required=True, type=Path)
    parser.add_argument("--round2", required=True, type=Path)
    parser.add_argument("--comparator-summary", required=True, type=Path)
    parser.add_argument("--comparator-samples", required=True, type=Path)
    parser.add_argument("--microsags-load-seconds", required=True, type=float)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    optimized = two_run_metrics([args.round1, args.round2])
    comparator = {row["route"]: row for row in read_tsv(args.comparator_summary)}
    sample_rows = read_tsv(args.comparator_samples)
    for route in ("skani_contigs", "mash_contigs", "sourmash_contigs"):
        if int(comparator[route]["completed_samples"]) != 24:
            raise ValueError(f"Incomplete archived comparator: {route}")
    load_total = 24 * args.microsags_load_seconds
    times = {
        "Microsags (reads)": (optimized["reads"]["mean_total_wall_seconds"] - load_total) / 3600,
        "Microsags (contigs)": (optimized["contigs"]["mean_total_wall_seconds"] - load_total) / 3600,
        "skani (contigs)": float(comparator["skani_contigs"]["total_hours"]),
        "Mash (contigs)": float(comparator["mash_contigs"]["total_hours"]),
        "sourmash (contigs)": float(comparator["sourmash_contigs"]["total_hours"]),
    }
    if min(times.values()) <= 0:
        raise ValueError("Load correction exceeds measured wall time")
    time_rows = sorted(((name, value, f"{value:.3f} h") for name, value in times.items()), key=lambda row: row[1])
    plot_panel(time_rows, args.out / "real24_time.svg",
               "24 real samples · annotation time", "Estimated time excluding database startup (hours)", 50,
               "Microsags: mean of first two 50-thread runs. Startup correction uses the archived 7.924 s/sample calibration; comparator estimates are archived.")

    memory = {
        "Microsags (reads)": optimized["reads"]["mean_run_peak_rss_gib"],
        "Microsags (contigs)": optimized["contigs"]["mean_run_peak_rss_gib"],
    }
    for route, label in (("skani_contigs", "skani (contigs)"), ("mash_contigs", "Mash (contigs)")):
        part = [row for row in sample_rows if row["route"] == route]
        if len(part) != 24:
            raise ValueError(f"Missing archived memory rows: {route}")
        memory[label] = max(float(row["memory_gib"]) for row in part)
    sour = [row for row in sample_rows if row["route"] == "sourmash_contigs"]
    if len(sour) != 24 or any(row["memory_metric"] != "sampled_aggregate_PSS" for row in sour):
        raise ValueError("Sourmash memory is not the archived PSS metric")
    sour_pss = max(float(row["memory_gib"]) for row in sour)
    memory_rows = sorted(((name, value, f"{value:.2f} GiB") for name, value in memory.items()), key=lambda row: row[1])
    memory_rows.append(("sourmash (contigs)", None, f"PSS only: {sour_pss:.2f} GiB (not RSS)"))
    plot_panel(memory_rows, args.out / "real24_memory.svg",
               "24 real samples · annotation memory", "Peak process RSS (GiB)", 64,
               "Microsags: mean of run-level peaks from the first two 50-thread runs. Other RSS values: archived single runs; sourmash reports PSS only.")
    print(f"reads_time_h={times['Microsags (reads)']:.6f} contigs_time_h={times['Microsags (contigs)']:.6f}")
    print(f"reads_peak_gib={memory['Microsags (reads)']:.6f} contigs_peak_gib={memory['Microsags (contigs)']:.6f}")


if __name__ == "__main__":
    main()
