#!/usr/bin/env python3
"""Create an exact, write-once FASTA-stat view for the Stage-3B pending SAGs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import uuid


FIELDS = (
    "sag_id", "assembly_fasta", "total_bp", "max_contig", "gc_pct",
    "gc_bases", "acgt_bases",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_rows(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise SystemExit(f"missing {label}: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"missing header in {label}: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pending", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    pending = args.pending.resolve(strict=True)
    stats = args.stats.resolve(strict=True)
    out = args.out_dir.absolute()
    if os.path.lexists(out):
        raise SystemExit(f"write-once destination exists: {out}")

    pending_header, pending_rows = read_rows(pending, "pending manifest")
    if pending_header != ["sag_id", "assembly_fasta", "reason"]:
        raise SystemExit(f"unexpected pending schema: {pending_header}")
    stats_header, stats_rows = read_rows(stats, "FASTA stats")
    if tuple(stats_header) != FIELDS:
        raise SystemExit(f"unexpected stats schema: {stats_header}")

    pending_ids: list[str] = []
    seen_pending: set[str] = set()
    for row in pending_rows:
        sag = row["sag_id"]
        if not sag or sag in seen_pending:
            raise SystemExit(f"empty or duplicate pending SAG: {sag!r}")
        seen_pending.add(sag)
        pending_ids.append(sag)
    if pending_ids != sorted(pending_ids, key=lambda value: value.encode("utf-8")):
        raise SystemExit("pending SAG order is not deterministic bytewise order")

    by_id: dict[str, dict[str, str]] = {}
    for row in stats_rows:
        sag = row["sag_id"]
        if not sag or sag in by_id:
            raise SystemExit(f"empty or duplicate stats SAG: {sag!r}")
        by_id[sag] = row
    missing = seen_pending - set(by_id)
    if missing:
        raise SystemExit(f"stats miss {len(missing)} pending SAGs; first={sorted(missing)[:20]}")

    staging = out.parent / f".{out.name}.incomplete-{os.getpid()}-{uuid.uuid4().hex}"
    staging.mkdir(parents=True, exist_ok=False)
    table = staging / "STAGE3B_FASTA_STATS.tsv"
    with table.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sag in pending_ids:
            writer.writerow(by_id[sag])
        handle.flush()
        os.fsync(handle.fileno())
    receipt = {
        "schema": "stage3b-pending-fasta-stats-view-v1",
        "status": "PASS",
        "pending": {"path": str(pending), "sha256": digest(pending)},
        "source_stats": {"path": str(stats), "sha256": digest(stats), "rows": len(stats_rows)},
        "output_stats": {
            "path": str(out / table.name), "sha256": digest(table), "rows": len(pending_ids)
        },
        "ignored_non_pending_rows": len(stats_rows) - len(pending_ids),
    }
    with (staging / "COMPLETE.json").open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.rename(staging, out)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
