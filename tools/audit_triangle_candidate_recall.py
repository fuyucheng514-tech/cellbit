#!/usr/bin/env python3
"""Audit sketch-triangle candidate recall against a frozen pair authority.

This script never recomputes ANI/AF.  It only asks whether every authoritative
undirected pair is present in the deterministic candidate TSV emitted by
``gtdb-ani-af triangle --triangle-mode sketch --triangle-candidates-out``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


class OutsideManifest(ValueError):
    """A pair endpoint is not a member of the explicitly supplied node set."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def choose(header: Iterable[str], aliases: tuple[str, ...], what: str) -> str:
    by_lower = {name.lower(): name for name in header}
    for alias in aliases:
        if alias.lower() in by_lower:
            return by_lower[alias.lower()]
    raise ValueError(f"missing {what} column; accepted aliases: {', '.join(aliases)}")


def optional(header: Iterable[str], aliases: tuple[str, ...]) -> str | None:
    by_lower = {name.lower(): name for name in header}
    for alias in aliases:
        if alias.lower() in by_lower:
            return by_lower[alias.lower()]
    return None


def tsv_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty TSV: {path}")
        return list(reader.fieldnames)


def tsv_rows(path: Path) -> Iterable[dict[str, str]]:
    """Stream data rows so a top-K Lake candidate table need not be materialized twice."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty TSV: {path}")
        yield from reader


def strip_fasta_suffix(name: str) -> str:
    lowered = name.lower()
    for suffix in ("_genomic.fna.gz", "_genomic.fna", ".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa", ".gz"):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]
    return name


class Resolver:
    def __init__(self, manifest: Path | None) -> None:
        self.aliases: dict[str, str] = {}
        self.ids: set[str] = set()
        self.gc: dict[str, float] = {}
        if manifest is None:
            return
        header = tsv_header(manifest)
        id_col = choose(header, ("sag_id", "SAG_id", "id"), "manifest SAG id")
        path_col = choose(header, ("assembly_fasta", "assembly_path", "path", "Ref_file"), "manifest assembly path")
        gc_col = optional(header, ("gc_pct", "GC", "gc"))
        basename_owner: dict[str, str | None] = {}
        for row in tsv_rows(manifest):
            sag = row[id_col].strip()
            raw_path = row[path_col].strip()
            if not sag or not raw_path or sag in self.ids:
                raise ValueError(f"empty/duplicate manifest SAG id: {sag!r}")
            self.ids.add(sag)
            for alias in (raw_path, str(Path(raw_path).resolve(strict=False))):
                old = self.aliases.get(alias)
                if old is not None and old != sag:
                    raise ValueError(f"manifest path maps to multiple SAGs: {alias}")
                self.aliases[alias] = sag
            base = Path(raw_path).name
            basename_owner[base] = sag if base not in basename_owner else None
            if gc_col is not None:
                value = float(row[gc_col])
                if not math.isfinite(value) or not 0.0 <= value <= 100.0:
                    raise ValueError(f"invalid GC for {sag}: {row[gc_col]!r}")
                self.gc[sag] = value
        for base, owner in basename_owner.items():
            if owner is not None:
                self.aliases.setdefault(base, owner)

    def resolve(self, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("empty pair endpoint")
        if value in self.ids:
            return value
        for alias in (value, str(Path(value).resolve(strict=False)), Path(value).name):
            if alias in self.aliases:
                return self.aliases[alias]
        stem = strip_fasta_suffix(Path(value).name)
        if stem in self.ids:
            return stem
        # Without a manifest, exact normalized paths remain a valid authority.
        if not self.ids:
            return str(Path(value).resolve(strict=False))
        raise OutsideManifest(f"cannot map pair endpoint to frozen manifest: {value}")


def pair(a: str, b: str) -> tuple[str, str]:
    if a == b:
        raise ValueError(f"self pair is not valid truth/candidate evidence: {a}")
    return (a, b) if a < b else (b, a)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--missing", type=Path, required=True)
    parser.add_argument("--min-ani", type=float)
    parser.add_argument("--min-max-af", type=float)
    parser.add_argument("--max-gc-diff", type=float)
    parser.add_argument("--require-recall", type=float)
    parser.add_argument(
        "--skip-truth-outside-manifest",
        action="store_true",
        help="explicitly restrict a larger truth TSV to endpoints present in --manifest",
    )
    args = parser.parse_args()

    if args.skip_truth_outside_manifest and args.manifest is None:
        parser.error("--skip-truth-outside-manifest requires --manifest")

    for path in (args.candidates, args.truth, args.manifest):
        if path is not None and not path.is_file():
            raise FileNotFoundError(path)
    for path in (args.summary, args.missing):
        if path.exists():
            raise FileExistsError(f"write-once output exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)

    resolver = Resolver(args.manifest)
    candidate_header = tsv_header(args.candidates)
    ca = choose(candidate_header, ("Ref_file", "Reference", "SAG1"), "candidate endpoint 1")
    cb = choose(candidate_header, ("Query_file", "Query", "SAG2"), "candidate endpoint 2")
    shared_col = optional(candidate_header, ("shared_sketches", "shared"))
    candidates: dict[tuple[str, str], int | None] = {}
    candidate_config_columns = ("sketch_scale", "sketch_size_cap", "max_posting", "min_shared", "max_candidates_per_genome")
    candidate_config: dict[str, str] = {}
    for index, row in enumerate(tsv_rows(args.candidates)):
        key = pair(resolver.resolve(row[ca]), resolver.resolve(row[cb]))
        shared = int(row[shared_col]) if shared_col is not None else None
        if key in candidates:
            raise ValueError(f"duplicate candidate pair: {key}")
        candidates[key] = shared
        for name in candidate_config_columns:
            if name in row:
                if index == 0:
                    candidate_config[name] = row[name]
                elif row[name] != candidate_config[name]:
                    raise ValueError(f"candidate configuration changes between rows: {name}")

    truth_header = tsv_header(args.truth)
    ta = choose(truth_header, ("Ref_file", "Reference", "SAG1"), "truth endpoint 1")
    tb = choose(truth_header, ("Query_file", "Query", "SAG2"), "truth endpoint 2")
    ani_col = optional(truth_header, ("ANI",))
    afr_col = optional(truth_header, ("Align_fraction_ref", "AF_ref"))
    afq_col = optional(truth_header, ("Align_fraction_query", "AF_query"))
    gc_col = optional(truth_header, ("gc_diff", "GC_diff"))
    if args.min_ani is not None and ani_col is None:
        raise ValueError("--min-ani requested but truth TSV has no ANI column")
    if args.min_max_af is not None and (afr_col is None or afq_col is None):
        raise ValueError("--min-max-af requested but truth TSV lacks both AF columns")
    if args.max_gc_diff is not None and gc_col is None and not resolver.gc:
        raise ValueError("--max-gc-diff needs truth gc_diff or manifest GC values")

    truth: dict[tuple[str, str], tuple[str, str]] = {}
    truth_rows_input = 0
    filtered = 0
    outside_manifest = 0
    for row in tsv_rows(args.truth):
        truth_rows_input += 1
        try:
            a, b = resolver.resolve(row[ta]), resolver.resolve(row[tb])
        except OutsideManifest:
            if not args.skip_truth_outside_manifest:
                raise
            outside_manifest += 1
            continue
        if args.min_ani is not None and float(row[ani_col]) < args.min_ani:
            filtered += 1
            continue
        if args.min_max_af is not None and max(float(row[afr_col]), float(row[afq_col])) < args.min_max_af:
            filtered += 1
            continue
        if args.max_gc_diff is not None:
            difference = float(row[gc_col]) if gc_col is not None else abs(resolver.gc[a] - resolver.gc[b])
            if difference > args.max_gc_diff:
                filtered += 1
                continue
        key = pair(a, b)
        if key in truth:
            raise ValueError(f"duplicate retained truth pair: {key}")
        truth[key] = (row[ta], row[tb])
    if not truth:
        raise ValueError("truth filters retained zero pairs")

    missing = sorted(set(truth) - set(candidates))
    recalled = len(truth) - len(missing)
    recall = recalled / len(truth)
    with args.missing.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("node1", "node2", "truth_endpoint1", "truth_endpoint2"))
        for key in missing:
            writer.writerow((*key, *truth[key]))
    summary = {
        "schema": "gtdb-ani-af-triangle-candidate-recall-audit-v1",
        "status": "PASS" if args.require_recall is None or recall + 1e-15 >= args.require_recall else "FAIL",
        "candidate_pairs": len(candidates),
        "truth_rows_input": truth_rows_input,
        "truth_rows_filtered": filtered,
        "truth_rows_outside_manifest": outside_manifest,
        "truth_pairs_retained": len(truth),
        "truth_pairs_recalled": recalled,
        "truth_pairs_missed": len(missing),
        "recall": recall,
        "required_recall": args.require_recall,
        "filters": {"min_ani": args.min_ani, "min_max_af": args.min_max_af, "max_gc_diff": args.max_gc_diff},
        "candidate_configuration": candidate_config,
        "inputs": {
            "candidates": str(args.candidates.resolve()),
            "candidates_sha256": sha256_file(args.candidates),
            "truth": str(args.truth.resolve()),
            "truth_sha256": sha256_file(args.truth),
            "manifest": str(args.manifest.resolve()) if args.manifest else None,
            "manifest_sha256": sha256_file(args.manifest) if args.manifest else None,
        },
        "missing_tsv": str(args.missing.resolve()),
        "missing_tsv_sha256": sha256_file(args.missing),
    }
    with args.summary.open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({key: summary[key] for key in ("status", "candidate_pairs", "truth_pairs_retained", "truth_pairs_recalled", "truth_pairs_missed", "recall")}, sort_keys=True))
    return 0 if summary["status"] == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
