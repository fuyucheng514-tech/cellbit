#!/usr/bin/env python3
"""Select a deterministic one-in-four SAG subset from a two-column manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    if args.output.exists() or args.receipt.exists():
        raise SystemExit("refusing to overwrite output or receipt")

    rows: list[tuple[str, str]] = []
    with args.input.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            columns = line.rstrip("\r\n").split("\t")
            if line_number == 1 and columns[:2] == ["sag_id", "assembly_fasta"]:
                continue
            if len(columns) != 2 or not columns[0] or not columns[1]:
                raise SystemExit(f"invalid two-column row at line {line_number}")
            rows.append((columns[0], columns[1]))

    rows.sort(key=lambda row: row[0])
    if len({row[0] for row in rows}) != len(rows):
        raise SystemExit("duplicate SAG identifier")
    selected = rows[::4]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write("sag_id\tassembly_fasta\n")
        for sag_id, fasta in selected:
            stream.write(f"{sag_id}\t{fasta}\n")

    receipt = {
        "schema": "microsags-lake-quarter-sample-v1",
        "status": "PASS",
        "method": "sort SAG IDs lexicographically; select ranks 1,5,9,...",
        "source_manifest_sha256": sha256(args.input),
        "source_sags": len(rows),
        "selected_sags": len(selected),
        "unique_sag_ids": len({row[0] for row in selected}),
        "manifest_sha256": sha256(args.output),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    with args.receipt.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()

