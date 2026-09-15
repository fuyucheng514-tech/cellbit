#!/usr/bin/env python3
"""Fail-closed audit for the large Stage3A Flye finish-thread canary.

The one-thread Flye assemble result is frozen.  The experiment changes only
the finish-phase thread count (2 or 8) and repeats each setting twice.  A PASS
requires exact final FASTA bytes, identical normalized sequence multisets for
every FASTA checkpoint, complete phase markers, and agreement with the
read-only production reference run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


CASES = ((2, 1), (2, 2), (8, 1), (8, 2))
FASTA_CHECKPOINTS = (
    "00-assembly/draft_assembly.fasta",
    "20-repeat/repeat_graph_edges.fasta",
    "30-contigger/contigs.fasta",
    "30-contigger/graph_final.fasta",
    "40-polishing/bubbles_1.fasta",
    "40-polishing/consensus_1.fasta",
    "40-polishing/polished_1.raw.fasta",
    "40-polishing/polished_1.fasta",
    "assembly.fasta",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_fasta(path: Path) -> dict[str, object]:
    sequences: list[str] = []
    current: list[str] = []
    records = 0
    with path.open("rt", encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
                records += 1
            else:
                if records == 0:
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    if records == 0 or records != len(sequences):
        raise ValueError(f"invalid FASTA closure: {path}")
    digest = hashlib.sha256(
        b"\0".join(sequence.encode("ascii") for sequence in sorted(sequences))
    ).hexdigest()
    return {"sha256": digest, "records": records, "bp": sum(map(len, sequences))}


def read_time(path: Path) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        key, value = raw.split("\t", 1)
        if key in {"max_rss_kb", "fs_inputs", "fs_outputs"}:
            values[key] = int(value)
        else:
            values[key] = float(value)
    required = {
        "elapsed_seconds",
        "user_seconds",
        "system_seconds",
        "max_rss_kb",
        "fs_inputs",
        "fs_outputs",
    }
    missing = required - values.keys()
    if missing:
        raise ValueError(f"missing timing fields in {path}: {sorted(missing)}")
    return values


def inspect_run(run: Path, time_path: Path | None) -> dict[str, object]:
    fasta: dict[str, object] = {}
    for relative in FASTA_CHECKPOINTS:
        path = run / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        fasta[relative] = {
            "raw_sha256": sha256(path),
            "normalized": normalized_fasta(path),
        }
    result: dict[str, object] = {"path": str(run), "fasta": fasta}
    if time_path is not None:
        result["time"] = read_time(time_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--reference-run", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    reference_run = args.reference_run.resolve()

    errors: list[str] = []
    cases: dict[str, object] = {}
    try:
        reference = inspect_run(reference_run, None)
    except (OSError, ValueError) as exc:
        errors.append(f"reference: {exc}")
        reference = {}

    for threads, replicate in CASES:
        name = f"t{threads}_r{replicate}"
        run = root / name
        required_markers = (
            run / "CPP_FULL_PIPELINE_PASS",
            run / "40-polishing/CPP_FINISH_PHASE_COMPLETE",
        )
        missing_markers = [str(path) for path in required_markers if not path.is_file()]
        if missing_markers:
            errors.append(f"{name}: missing markers {missing_markers}")
            continue
        try:
            cases[name] = inspect_run(run, root / f"{name}.time.tsv")
        except (OSError, ValueError) as exc:
            errors.append(f"{name}: {exc}")

    all_runs: dict[str, object] = {}
    if reference:
        all_runs["reference"] = reference
    all_runs.update(cases)

    raw_equal: dict[str, bool] = {}
    normalized_equal: dict[str, bool] = {}
    reference_raw_equal: dict[str, bool] = {}
    reference_normalized_equal: dict[str, bool] = {}
    for relative in FASTA_CHECKPOINTS:
        raw_values = {
            entry["fasta"][relative]["raw_sha256"]  # type: ignore[index]
            for entry in all_runs.values()
            if relative in entry.get("fasta", {})  # type: ignore[union-attr]
        }
        normalized_values = {
            entry["fasta"][relative]["normalized"]["sha256"]  # type: ignore[index]
            for entry in all_runs.values()
            if relative in entry.get("fasta", {})  # type: ignore[union-attr]
        }
        raw_equal[relative] = len(raw_values) == 1
        normalized_equal[relative] = len(normalized_values) == 1
        if reference:
            reference_raw = reference["fasta"][relative]["raw_sha256"]  # type: ignore[index]
            reference_normalized = reference["fasta"][relative]["normalized"]["sha256"]  # type: ignore[index]
            reference_raw_equal[relative] = all(
                entry["fasta"][relative]["raw_sha256"] == reference_raw  # type: ignore[index]
                for entry in cases.values()
            )
            reference_normalized_equal[relative] = all(
                entry["fasta"][relative]["normalized"]["sha256"] == reference_normalized  # type: ignore[index]
                for entry in cases.values()
            )

    checks = {
        "four_cases_present": len(cases) == len(CASES),
        "reference_present": bool(reference),
        "final_raw_fasta_byte_identical_including_reference": raw_equal.get("assembly.fasta", False),
        "final_normalized_sequences_identical_including_reference": normalized_equal.get("assembly.fasta", False),
        "all_checkpoint_fastas_normalized_identical_including_reference": bool(normalized_equal)
        and all(normalized_equal.values()),
        "all_checkpoint_fastas_match_reference_normalized": bool(reference_normalized_equal)
        and all(reference_normalized_equal.values()),
    }
    status = "PASS" if not errors and all(checks.values()) else "FAIL"
    result = {
        "schema": "cellbit57.flye-finish-large-group-canary.v1",
        "status": status,
        "group": "Planktophila_sp029977785_",
        "frozen_assemble_threads": 1,
        "changed_variable_only": "finish_threads",
        "thread_counts": [2, 8],
        "replicates_per_count": 2,
        "reference_run": str(reference_run),
        "checks": checks,
        "checkpoint_raw_byte_equal_including_reference": raw_equal,
        "checkpoint_normalized_equal_including_reference": normalized_equal,
        "checkpoint_raw_match_reference": reference_raw_equal,
        "checkpoint_normalized_match_reference": reference_normalized_equal,
        "errors": errors,
        "reference": reference,
        "cases": cases,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    (root / "RESULT.json").write_text(payload, encoding="utf-8")

    with (root / "SUMMARY.tsv").open("wt", encoding="utf-8") as handle:
        handle.write(
            "case\tthreads\treplicate\telapsed_seconds\tuser_seconds\tsystem_seconds\t"
            "max_rss_kb\tfs_inputs\tfs_outputs\traw_sha256\tnormalized_sha256\trecords\tbp\n"
        )
        for threads, replicate in CASES:
            name = f"t{threads}_r{replicate}"
            if name not in cases:
                continue
            entry = cases[name]  # type: ignore[assignment]
            timing = entry["time"]  # type: ignore[index]
            assembly = entry["fasta"]["assembly.fasta"]  # type: ignore[index]
            normalized = assembly["normalized"]
            handle.write(
                f"{name}\t{threads}\t{replicate}\t{timing['elapsed_seconds']}\t"
                f"{timing['user_seconds']}\t{timing['system_seconds']}\t{timing['max_rss_kb']}\t"
                f"{timing['fs_inputs']}\t{timing['fs_outputs']}\t{assembly['raw_sha256']}\t"
                f"{normalized['sha256']}\t{normalized['records']}\t{normalized['bp']}\n"
            )

    if status == "PASS":
        (root / "COMPLETE.json").write_text(payload, encoding="utf-8")
        print("PASS")
        return 0
    print("FAIL")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
