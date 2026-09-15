#!/usr/bin/env python3
"""Audit a write-once Flye finish-phase thread-scaling canary.

The deterministic Flye assemble result is frozen before this experiment.  The
only changed variable is the thread count used by repeat, contigger, alignment,
bubble generation and polisher.  Publication is fail-closed: every run must
finish, the final FASTA must be byte-identical across all thread counts and
replicates, and its sequence multiset must also be identical.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


THREADS = (2, 4, 8)
REPLICATES = (1, 2)
INTERMEDIATE_FASTA = (
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


def fasta_sequences(path: Path) -> tuple[list[str], int]:
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
    return sequences, records


def normalized_sequence_sha(path: Path) -> tuple[str, int, int]:
    sequences, records = fasta_sequences(path)
    normalized = b"\0".join(sequence.encode("ascii") for sequence in sorted(sequences))
    return hashlib.sha256(normalized).hexdigest(), records, sum(map(len, sequences))


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
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()

    cases: dict[str, dict] = {}
    errors: list[str] = []
    for threads in THREADS:
        for replicate in REPLICATES:
            name = f"t{threads}_r{replicate}"
            case = root / name
            required = (
                case / "CPP_FULL_PIPELINE_PASS",
                case / "40-polishing/CPP_FINISH_PHASE_COMPLETE",
                case / "assembly.fasta",
                root / f"{name}.time.tsv",
            )
            missing = [str(path) for path in required if not path.is_file()]
            if missing:
                errors.append(f"{name}: missing {missing}")
                continue
            raw_sha = sha256(case / "assembly.fasta")
            normalized_sha, records, bp = normalized_sequence_sha(case / "assembly.fasta")
            intermediate = {}
            intermediate_normalized = {}
            for relative in INTERMEDIATE_FASTA:
                path = case / relative
                if not path.is_file():
                    errors.append(f"{name}: missing intermediate {relative}")
                    continue
                intermediate[relative] = sha256(path)
                normalized, intermediate_records, intermediate_bp = normalized_sequence_sha(path)
                intermediate_normalized[relative] = {
                    "sha256": normalized,
                    "records": intermediate_records,
                    "bp": intermediate_bp,
                }
            cases[name] = {
                "threads": threads,
                "replicate": replicate,
                "raw_assembly_sha256": raw_sha,
                "normalized_sequence_sha256": normalized_sha,
                "records": records,
                "bp": bp,
                "time": read_time(root / f"{name}.time.tsv"),
                "intermediate_fasta_sha256": intermediate,
                "intermediate_fasta_normalized": intermediate_normalized,
            }

    raw_values = {case["raw_assembly_sha256"] for case in cases.values()}
    normalized_values = {case["normalized_sequence_sha256"] for case in cases.values()}
    intermediate_equal = {
        relative: len(
            {
                case["intermediate_fasta_sha256"].get(relative)
                for case in cases.values()
            }
        ) == 1
        for relative in INTERMEDIATE_FASTA
    }
    intermediate_normalized_equal = {
        relative: len(
            {
                case["intermediate_fasta_normalized"].get(relative, {}).get("sha256")
                for case in cases.values()
            }
        ) == 1
        for relative in INTERMEDIATE_FASTA
    }
    checks = {
        "six_cases_present": len(cases) == len(THREADS) * len(REPLICATES),
        "final_raw_fasta_byte_identical": len(raw_values) == 1,
        "final_normalized_sequences_identical": len(normalized_values) == 1,
        # Flye's polisher may serialize consensus_1 records in worker-completion
        # order.  It is an ephemeral input to deterministic composition, not a
        # released bin.  Require every intermediate sequence multiset to be
        # identical and report (but do not require) raw serialization identity.
        "all_intermediate_fasta_normalized_identical": all(intermediate_normalized_equal.values()),
    }
    status = "PASS" if not errors and all(checks.values()) else "FAIL"
    result = {
        "schema": "cellbit57.flye-finish-thread-scaling-canary.v1",
        "status": status,
        "frozen_assemble_threads": 1,
        "changed_variable_only": "finish_threads",
        "thread_counts": list(THREADS),
        "replicates_per_count": len(REPLICATES),
        "checks": checks,
        "intermediate_fasta_byte_equal": intermediate_equal,
        "intermediate_fasta_normalized_equal": intermediate_normalized_equal,
        "errors": errors,
        "cases": cases,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    (root / "RESULT.json").write_text(payload, encoding="utf-8")

    with (root / "SUMMARY.tsv").open("wt", encoding="utf-8") as handle:
        handle.write("case\tthreads\treplicate\telapsed_seconds\tuser_seconds\tmax_rss_kb\traw_sha256\tnormalized_sha256\trecords\tbp\n")
        for name, case in sorted(cases.items(), key=lambda item: (item[1]["threads"], item[1]["replicate"])):
            timing = case["time"]
            handle.write(
                f"{name}\t{case['threads']}\t{case['replicate']}\t"
                f"{timing['elapsed_seconds']}\t{timing['user_seconds']}\t{timing['max_rss_kb']}\t"
                f"{case['raw_assembly_sha256']}\t{case['normalized_sequence_sha256']}\t"
                f"{case['records']}\t{case['bp']}\n"
            )
    if status == "PASS":
        (root / "COMPLETE.json").write_text(payload, encoding="utf-8")
        print("PASS")
        return 0
    print("FAIL")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
