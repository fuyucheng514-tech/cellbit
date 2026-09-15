#!/usr/bin/env python3
"""Create a write-once continuation root from an audited failed Lake run.

Only successful Stage 1 and Dna2bit artifacts are copied.  Stage 3A outputs
are deliberately excluded so a corrected subassembly boundary starts cleanly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_regular_files(root: Path) -> tuple[int, int, str]:
    rows: list[str] = []
    total = 0
    count = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        size = path.stat().st_size
        relative = path.relative_to(root).as_posix()
        rows.append(f"{relative}\t{size}\t{sha256(path)}\n")
        total += size
        count += 1
    aggregate = hashlib.sha256("".join(rows).encode()).hexdigest()
    return count, total, aggregate


def copy_required(source: Path, destination: Path) -> None:
    shutil.copytree(source / "01_assembly", destination / "01_assembly", symlinks=True)
    dna_source = source / "02_dna2bit"
    dna_destination = destination / "02_dna2bit"
    dna_destination.mkdir(parents=True)
    shutil.copytree(dna_source / "bits", dna_destination / "bits", symlinks=True)
    for name in ("SKETCH.PASS", "SEARCH.PASS", "search_result.csv"):
        shutil.copy2(dna_source / name, dna_destination / name)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--taxonomy-recovery-receipt", type=Path, required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    destination = args.destination.resolve()
    manifest = args.manifest.resolve()
    recovery = args.taxonomy_recovery_receipt.resolve()
    if destination.exists():
        raise SystemExit(f"destination already exists: {destination}")
    for required in (
        source / "01_assembly",
        source / "02_dna2bit" / "bits",
        source / "02_dna2bit" / "SKETCH.PASS",
        source / "02_dna2bit" / "SEARCH.PASS",
        source / "02_dna2bit" / "search_result.csv",
        manifest,
        recovery,
    ):
        if not required.exists():
            raise SystemExit(f"required source artifact is missing: {required}")

    source_count, source_bytes, source_fingerprint = fingerprint_regular_files(
        source / "02_dna2bit" / "bits"
    )
    if source_count < 1:
        raise SystemExit("source contains no Dna2bit bit files")

    destination.mkdir(parents=True)
    try:
        copy_required(source, destination)
        copied_count, copied_bytes, copied_fingerprint = fingerprint_regular_files(
            destination / "02_dna2bit" / "bits"
        )
        if (copied_count, copied_bytes, copied_fingerprint) != (
            source_count,
            source_bytes,
            source_fingerprint,
        ):
            raise RuntimeError("copied Dna2bit artifacts do not match source fingerprint")

        receipt = {
            "schema": "lake-current-tractor-baseline-continuation-v1",
            "status": "PASS",
            "source_failed_run": str(source),
            "manifest": {"path": str(manifest), "sha256": sha256(manifest)},
            "taxonomy_recovery_receipt": {
                "path": str(recovery),
                "sha256": sha256(recovery),
            },
            "search_result": {
                "source_path": str(source / "02_dna2bit" / "search_result.csv"),
                "sha256": sha256(source / "02_dna2bit" / "search_result.csv"),
            },
            "dna2bit_bits": {
                "regular_file_count": copied_count,
                "bytes": copied_bytes,
                "aggregate_sha256": copied_fingerprint,
            },
            "excluded_from_copy": ["03A_subassemble", "COMPLETE.json", "TIMING.tsv", "TIMING.json"],
            "purpose": "resume frozen Stage1/Dna2bit evidence with corrected no-overlap Stage3A policy",
        }
        temporary = destination / ".CONTINUATION_PROVENANCE.json.tmp"
        temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, destination / "CONTINUATION_PROVENANCE.json")
    except Exception:
        # Keep the incomplete directory as failure evidence; do not delete or
        # silently retry into it.
        raise
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
