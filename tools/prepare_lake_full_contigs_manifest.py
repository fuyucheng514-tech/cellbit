#!/usr/bin/env python3
"""Freeze the complete Lake assembly table as a two-column contig manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_stats(path: Path) -> tuple[int, int, int, float]:
    records = total = maximum = gc_bases = 0
    sequence_length = 0
    seen_header = False
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_length:
                    records += 1
                    total += sequence_length
                    maximum = max(maximum, sequence_length)
                    sequence_length = 0
                if len(line) == 1:
                    raise ValueError(f"empty FASTA header in {path}:{line_number}")
                seen_header = True
                continue
            if not seen_header:
                raise ValueError(f"sequence before FASTA header in {path}:{line_number}")
            if not re.fullmatch(r"[ACGTNacgtn]+", line):
                raise ValueError(f"invalid FASTA sequence in {path}:{line_number}")
            upper = line.upper()
            sequence_length += len(upper)
            gc_bases += upper.count("G") + upper.count("C")
    if sequence_length:
        records += 1
        total += sequence_length
        maximum = max(maximum, sequence_length)
    if records == 0 or total == 0:
        raise ValueError(f"empty FASTA: {path}")
    return records, total, maximum, 100.0 * gc_bases / total


def recover_spades_scaffolds(
    path: Path, expected: tuple[int, int, int, float]
) -> tuple[Path, dict[str, object]]:
    candidates: list[tuple[Path, tuple[int, int, int, float]]] = []
    for candidate in sorted(path.parent.glob("K*/scaffolds.fasta")):
        if not candidate.is_file() or candidate.stat().st_size == 0:
            continue
        observed = fasta_stats(candidate)
        if (
            observed[:3] == expected[:3]
            and round(observed[3], 4) == round(expected[3], 4)
        ):
            candidates.append((candidate.resolve(), observed))
    if len(candidates) != 1:
        raise ValueError(
            f"expected one authenticated SPAdes fallback for {path}, observed {len(candidates)}"
        )
    candidate, observed = candidates[0]
    return candidate, {
        "empty_authority_path": str(path),
        "recovered_path": str(candidate),
        "recovered_sha256": sha256(candidate),
        "verification": {
            "contigs": observed[0],
            "total_bp": observed[1],
            "max_contig": observed[2],
            "gc_pct": observed[3],
            "matches_frozen_table": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-stats", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expected-count", type=int, default=13742)
    args = parser.parse_args()

    source = args.assembly_stats.resolve()
    out_dir = args.out_dir.resolve()
    if out_dir.exists():
        raise SystemExit(f"write-once output already exists: {out_dir}")
    if not source.is_file():
        raise SystemExit(f"assembly stats table is missing: {source}")

    rows: list[tuple[str, Path, int, int, int, float]] = []
    seen: set[str] = set()
    recoveries: list[dict[str, object]] = []
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SAG_id", "assembly_path", "contigs", "total_len", "max_contig", "gc_pct"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise SystemExit("assembly stats table has an unexpected schema")
        for row_number, row in enumerate(reader, 2):
            sag = row["SAG_id"].strip()
            authority_path = Path(row["assembly_path"])
            path = authority_path.resolve()
            if not sag or sag in seen:
                raise SystemExit(f"empty or duplicate SAG_id at row {row_number}: {sag!r}")
            if any(char.isspace() for char in sag) or "/" in sag or "\\" in sag:
                raise SystemExit(f"unsafe SAG_id at row {row_number}: {sag!r}")
            contigs = int(row["contigs"])
            total = int(row["total_len"])
            maximum = int(row["max_contig"])
            gc = float(row["gc_pct"])
            if contigs < 1 or total < 1 or maximum < 1 or not 0.0 <= gc <= 100.0:
                raise SystemExit(f"invalid frozen assembly statistics at row {row_number}")
            if not path.is_file() or path.stat().st_size == 0:
                try:
                    path, recovery = recover_spades_scaffolds(
                        authority_path, (contigs, total, maximum, gc)
                    )
                except ValueError as error:
                    raise SystemExit(
                        f"missing/empty assembly has no unique authenticated recovery at row {row_number}: {error}"
                    ) from error
                recovery["sag_id"] = sag
                recoveries.append(recovery)
            seen.add(sag)
            rows.append((sag, path, contigs, total, maximum, gc))

    if len(rows) != args.expected_count:
        raise SystemExit(f"observed {len(rows)} assemblies; expected {args.expected_count}")
    rows.sort(key=lambda row: row[0].encode())

    out_dir.mkdir(parents=True)
    manifest = out_dir / "LAKE_FULL_CONTIGS.tsv"
    with manifest.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sag_id", "assembly_fasta"])
        for sag, path, *_ in rows:
            writer.writerow([sag, str(path)])

    receipt = {
        "schema": "lake-full-contigs-manifest-v1",
        "status": "PASS",
        "source_table": {"path": str(source), "sha256": sha256(source)},
        "manifest": {"path": str(manifest), "sha256": sha256(manifest)},
        "assembly_count": len(rows),
        "unique_sag_ids": len(seen),
        "existing_nonempty_files": len(rows),
        "frozen_total_contigs": sum(row[2] for row in rows),
        "frozen_total_bp": sum(row[3] for row in rows),
        "frozen_min_total_bp": min(row[3] for row in rows),
        "frozen_max_total_bp": max(row[3] for row in rows),
        "authenticated_spades_fallbacks": recoveries,
        "note": "Sequence-level totals are independently recalculated by dna2bit-sag-pipeline before classification.",
    }
    temporary = out_dir / ".MANIFEST_AUDIT.json.tmp"
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, out_dir / "MANIFEST_AUDIT.json")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
