#!/usr/bin/env python3
"""Read-only duration census for an already completed split Flye finish phase."""

from __future__ import annotations

import argparse
import csv
import re
from datetime import datetime, timezone
from pathlib import Path


STAMP = re.compile(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]")


def first_flye_epoch(log: Path) -> float:
    with log.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            match = STAMP.match(raw)
            if match:
                # Flye's logger uses UTC even though the host displays CST.
                return datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S").replace(
                    tzinfo=timezone.utc
                ).timestamp()
    raise ValueError(f"no Flye timestamp in {log}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--subassemble-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for group_dir in sorted(path for path in args.subassemble_root.iterdir() if path.is_dir()):
        input_fasta = group_dir / "input_subassemblies.fasta"
        run = group_dir / "run"
        log = run / "cpp-subass.log"
        marker = run / "CPP_FULL_PIPELINE_PASS"
        if not all(path.is_file() for path in (input_fasta, log, marker)):
            continue
        start = first_flye_epoch(log)
        finish = marker.stat().st_mtime
        rows.append(
            {
                "group": group_dir.name,
                "input_bytes": input_fasta.stat().st_size,
                "finish_seconds": round(finish - start, 6),
                "start_epoch": round(start, 6),
                "finish_epoch": round(finish, 6),
            }
        )
    if not rows:
        raise SystemExit("no closed group runs found")
    rows.sort(key=lambda row: (-row["finish_seconds"], -row["input_bytes"], row["group"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows[:20]:
        print(f"{row['finish_seconds']:.3f}\t{row['input_bytes']}\t{row['group']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
