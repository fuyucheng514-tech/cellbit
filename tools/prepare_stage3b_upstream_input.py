#!/usr/bin/env python3
"""Materialize a write-once, auditable input view for CheckM2/GTDB-Tk."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


ALLOWED_REASONS = {
    "dna2bit_rejected_or_no_hit",
    "dna2bit_negative",
    "dna2bit_no_hit",
    "dna2bit_rejected",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pending", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    pending = args.pending.resolve()
    out_dir = args.out_dir.resolve()
    if out_dir.exists():
        raise SystemExit(f"write-once output already exists: {out_dir}")
    if not pending.is_file():
        raise SystemExit(f"pending manifest is missing: {pending}")

    rows: list[tuple[str, Path, str]] = []
    seen: set[str] = set()
    with pending.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sag_id", "assembly_fasta", "reason"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise SystemExit("unexpected pending-manifest schema")
        for row_number, row in enumerate(reader, 2):
            sag = row["sag_id"].strip()
            assembly = Path(row["assembly_fasta"]).resolve()
            reason = row["reason"].strip()
            if not sag or sag in seen or any(c.isspace() for c in sag) or "/" in sag or "\\" in sag:
                raise SystemExit(f"empty, duplicate, or unsafe SAG id at row {row_number}: {sag!r}")
            if reason not in ALLOWED_REASONS:
                raise SystemExit(f"invalid Dna2bit-negative reason at row {row_number}: {reason!r}")
            if not assembly.is_file() or assembly.stat().st_size == 0:
                raise SystemExit(f"missing or empty assembly at row {row_number}: {assembly}")
            seen.add(sag)
            rows.append((sag, assembly, reason))
    if not rows:
        raise SystemExit("pending manifest contains no SAGs")
    rows.sort(key=lambda item: item[0].encode())

    input_dir = out_dir / "input_fna"
    input_dir.mkdir(parents=True)
    canonical_manifest = out_dir / "INPUT_MANIFEST.tsv"
    with canonical_manifest.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sag_id", "source_assembly", "input_link", "reason"])
        for sag, assembly, reason in rows:
            link = input_dir / f"{sag}.fna"
            os.symlink(assembly, link)
            writer.writerow([sag, str(assembly), str(link), reason])

    receipt = {
        "schema": "stage3b-upstream-input-view-v1",
        "status": "PASS",
        "source_pending": {"path": str(pending), "sha256": sha256(pending)},
        "canonical_manifest": {
            "path": str(canonical_manifest),
            "sha256": sha256(canonical_manifest),
        },
        "input_directory": str(input_dir),
        "sag_count": len(rows),
        "unique_sag_ids": len(seen),
        "nonempty_resolved_targets": len(rows),
    }
    temporary = out_dir / ".COMPLETE.json.tmp"
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, out_dir / "COMPLETE.json")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
