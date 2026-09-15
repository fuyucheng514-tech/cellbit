#!/usr/bin/env python3
"""Read-only equivalence audit for a 58-group Stage-3B Flye replay.

HISTORICAL ONLY: the 58-group collection came from the superseded
GTDB-reference double-negative entrance.  This audit remains valid only for
that frozen evidence and is not a current Stage3B release gate.

The audit deliberately distinguishes exact reproduction from biological
sequence-level equivalence.  It never edits either assembly collection.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import platform
import shlex
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "cellbit57.stage3b.flye_equivalence_audit.v1"
GROUPS = [f"G{i:04d}" for i in range(1, 59)]
RC_TABLE = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def canonical_sequence_digest(sequences: Iterable[bytes]) -> str:
    per_sequence: list[bytes] = []
    for seq in sequences:
        seq = seq.upper()
        rev = seq.translate(RC_TABLE)[::-1]
        canonical = seq if seq <= rev else rev
        per_sequence.append(hashlib.sha256(canonical).digest())
    per_sequence.sort()
    h = hashlib.sha256()
    for digest in per_sequence:
        h.update(digest)
    return h.hexdigest()


def read_fasta_sequences(path: Path) -> list[bytes]:
    sequences: list[bytes] = []
    chunks: list[bytes] = []
    saw_header = False
    with path.open("rb") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if saw_header:
                    sequences.append(b"".join(chunks))
                saw_header = True
                chunks = []
            else:
                if not saw_header:
                    raise ValueError(f"sequence before FASTA header: {path}")
                chunks.append(line)
    if saw_header:
        sequences.append(b"".join(chunks))
    if not sequences or any(not seq for seq in sequences):
        raise ValueError(f"empty or invalid FASTA: {path}")
    return sequences


def n50(lengths: list[int]) -> int:
    target = (sum(lengths) + 1) // 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    raise AssertionError("unreachable N50")


def fasta_stats(path: Path) -> dict[str, Any]:
    sequences = read_fasta_sequences(path)
    lengths = [len(seq) for seq in sequences]
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "contigs": len(sequences),
        "bp": sum(lengths),
        "n50": n50(lengths),
        "max_contig": max(lengths),
        "canonical_sequence_multiset_sha256": canonical_sequence_digest(sequences),
    }


def write_json(path: Path, obj: Any) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def percentile_summary(values: list[float]) -> dict[str, float]:
    return {
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
    }


def parse_skani(stdout: str, group: str, new_path: Path, historical_path: Path) -> dict[str, Any]:
    rows = list(csv.DictReader(stdout.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"{group}: expected exactly one skani row, observed {len(rows)}")
    row = rows[0]
    if Path(row["Query_file"]).resolve() != new_path.resolve():
        raise RuntimeError(f"{group}: skani query orientation mismatch")
    if Path(row["Ref_file"]).resolve() != historical_path.resolve():
        raise RuntimeError(f"{group}: skani reference orientation mismatch")
    return {
        "group": group,
        "query": str(new_path),
        "reference": str(historical_path),
        "ani": float(row["ANI"]),
        "af_query": float(row["Align_fraction_query"]),
        "af_ref": float(row["Align_fraction_ref"]),
        "query_name": row.get("Query_name", ""),
        "ref_name": row.get("Ref_name", ""),
    }


def run_skani_pair(
    group: str,
    skani: Path,
    new_path: Path,
    historical_path: Path,
    raw_dir: Path,
    threads: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    command = [
        str(skani),
        "dist",
        "--medium",
        "--min-af",
        "0",
        "-t",
        str(threads),
        str(new_path),
        str(historical_path),
    ]
    started = time.monotonic()
    proc = subprocess.run(command, text=True, capture_output=True, check=False)
    elapsed = time.monotonic() - started
    stdout_path = raw_dir / f"{group}.stdout.tsv"
    stderr_path = raw_dir / f"{group}.stderr.log"
    stdout_path.write_text(proc.stdout, encoding="utf-8", newline="\n")
    stderr_path.write_text(proc.stderr, encoding="utf-8", newline="\n")
    command_record = {
        "group": group,
        "command": command,
        "command_shell_escaped": shlex.join(command),
        "returncode": proc.returncode,
        "wall_seconds": elapsed,
        "stdout_path": str(stdout_path),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_path": str(stderr_path),
        "stderr_sha256": sha256_file(stderr_path),
    }
    if proc.returncode != 0:
        raise RuntimeError(f"{group}: skani exited {proc.returncode}; see {stderr_path}")
    result = parse_skani(proc.stdout, group, new_path, historical_path)
    return result, command_record


def count_at_least(rows: list[dict[str, Any]], key: str, threshold: float) -> int:
    return sum(float(row[key]) >= threshold for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-root", required=True, type=Path)
    parser.add_argument("--historical-root", required=True, type=Path)
    parser.add_argument("--historical-input-root", required=True, type=Path)
    parser.add_argument("--skani", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--threads-per-job", type=int, default=2)
    args = parser.parse_args()

    if args.jobs < 1 or args.threads_per_job < 1:
        parser.error("jobs and threads-per-job must be positive")
    if args.jobs * args.threads_per_job > 32:
        parser.error("light audit safety limit: jobs * threads-per-job must be <= 32")
    if args.out.exists():
        raise FileExistsError(f"write-once output already exists: {args.out}")
    for required in (args.new_root, args.historical_root, args.historical_input_root):
        if not required.is_dir():
            raise FileNotFoundError(required)
    if not args.skani.is_file() or not os.access(args.skani, os.X_OK):
        raise FileNotFoundError(f"skani executable unavailable: {args.skani}")

    stage = args.out.parent / f".{args.out.name}.incomplete-{os.getpid()}"
    if stage.exists():
        raise FileExistsError(stage)
    stage.mkdir(parents=False)
    raw_dir = stage / "skani_raw"
    raw_dir.mkdir()

    started_utc = datetime.now(timezone.utc).isoformat()
    preregistration = {
        "schema": SCHEMA,
        "created_before_pairwise_skani": started_utc,
        "scope": "Sequence-level comparison only; CheckM2 quality is intentionally out of scope.",
        "orientation": {
            "query": "new replay assembly",
            "reference": "historical Lake assembly",
        },
        "method": {
            "tool": "skani dist",
            "preset": "--medium",
            "minimum_output_af": 0,
            "learned_ani": "enabled by skani default",
        },
        "claims": {
            "exact_reproduction": "all 58 byte hashes and orientation-invariant sequence multisets equal",
            "near_equivalent": "all 58 ANI >= 99 and both directional AF >= 95",
            "same_species_level": "all 58 ANI >= 95 and both directional AF >= 50",
        },
        "note": "PASS means this preregistered audit completed and the near-equivalent criterion passed; it never means byte-exact when exact_reproduction is false.",
    }
    write_json(stage / "PREREGISTRATION.json", preregistration)

    complete: dict[str, Any]
    try:
        replay_path = args.new_root / "REPLAY_AUDIT.json"
        upstream_complete_path = args.new_root / "COMPLETE.json"
        replay = json.loads(replay_path.read_text(encoding="utf-8"))
        upstream_complete = json.loads(upstream_complete_path.read_text(encoding="utf-8"))
        group_keys = sorted(replay.get("groups", {}))
        failed_upstream_checks = sorted(
            key for key, value in replay.get("checks", {}).items() if value is not True
        )
        permitted_nondeterminism_checks = [
            "historical_assemblies_canonical_exact",
            "historical_output_aggregate_exact",
        ]
        prereq = {
            "group_ids_exact_G0001_G0058": group_keys == GROUPS,
            "replay_all_flye_commands_successful": replay.get("checks", {}).get("all_flye_commands_successful") is True,
            "all_58_returncodes_zero": all(
                replay["groups"][group].get("flye", {}).get("returncode") == 0 for group in GROUPS
            ) if group_keys == GROUPS else False,
            "replay_historical_inputs_byte_exact": replay.get("checks", {}).get("historical_input_aggregate_exact") is True,
            "upstream_failed_checks_only_expected_nondeterminism": failed_upstream_checks == permitted_nondeterminism_checks,
            "upstream_raw_pass": replay.get("checks", {}).get("upstream_raw_pass") is True,
            "upstream_graph_marker_signed_pass": replay.get("checks", {}).get("upstream_graph_marker_signed_pass") is True,
        }
        if not all(prereq.values()):
            raise RuntimeError(f"upstream replay prerequisites failed: {prereq}; failed_checks={failed_upstream_checks}")

        version_proc = subprocess.run([str(args.skani), "--version"], text=True, capture_output=True, check=True)
        versions = {
            "schema": SCHEMA,
            "python": sys.version,
            "platform": platform.platform(),
            "skani_path": str(args.skani.resolve()),
            "skani_sha256": sha256_file(args.skani),
            "skani_version_stdout": version_proc.stdout.strip(),
            "script_path": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
        }
        write_json(stage / "VERSIONS.json", versions)

        paths: dict[str, dict[str, Path]] = {}
        for group in GROUPS:
            paths[group] = {
                "new_input": args.new_root / "inputs" / f"{group}.fasta",
                "historical_input": args.historical_input_root / f"{group}.fasta",
                "new_assembly": args.new_root / "flye" / group / "assembly.fasta",
                "historical_assembly": args.historical_root / group / "assembly.fasta",
            }
            for path in paths[group].values():
                if not path.is_file():
                    raise FileNotFoundError(path)

        # Hash each source directly; do not trust only the upstream receipt.
        hash_tasks: list[tuple[str, str, Path]] = []
        for group in GROUPS:
            for role, path in paths[group].items():
                hash_tasks.append((group, role, path))
        source_rows: list[dict[str, Any]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures = {
                pool.submit(sha256_file, path): (group, role, path)
                for group, role, path in hash_tasks
            }
            for future in concurrent.futures.as_completed(futures):
                group, role, path = futures[future]
                source_rows.append({
                    "group": group,
                    "role": role,
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": future.result(),
                })
        source_rows.sort(key=lambda row: (row["group"], row["role"]))
        source_lookup = {(row["group"], row["role"]): row for row in source_rows}
        direct_input_exact = {
            group: (
                source_lookup[(group, "new_input")]["bytes"] == source_lookup[(group, "historical_input")]["bytes"]
                and source_lookup[(group, "new_input")]["sha256"] == source_lookup[(group, "historical_input")]["sha256"]
            )
            for group in GROUPS
        }
        prereq["direct_all_58_inputs_byte_exact"] = all(direct_input_exact.values())
        if not prereq["direct_all_58_inputs_byte_exact"]:
            raise RuntimeError(
                "direct input hash comparison failed for "
                + ",".join(group for group, ok in direct_input_exact.items() if not ok)
            )
        with (stage / "SOURCE_HASHES.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["group", "role", "path", "bytes", "sha256"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(source_rows)

        # Independently recompute assembly structure and exact-sequence digests.
        structure_rows: list[dict[str, Any]] = []
        for group in GROUPS:
            new_stats = fasta_stats(paths[group]["new_assembly"])
            hist_stats = fasta_stats(paths[group]["historical_assembly"])
            structure_rows.append({
                "group": group,
                "new_contigs": new_stats["contigs"],
                "historical_contigs": hist_stats["contigs"],
                "delta_contigs": new_stats["contigs"] - hist_stats["contigs"],
                "new_bp": new_stats["bp"],
                "historical_bp": hist_stats["bp"],
                "delta_bp": new_stats["bp"] - hist_stats["bp"],
                "delta_bp_percent_of_historical": 100.0 * (new_stats["bp"] - hist_stats["bp"]) / hist_stats["bp"],
                "new_n50": new_stats["n50"],
                "historical_n50": hist_stats["n50"],
                "delta_n50": new_stats["n50"] - hist_stats["n50"],
                "new_max_contig": new_stats["max_contig"],
                "historical_max_contig": hist_stats["max_contig"],
                "delta_max_contig": new_stats["max_contig"] - hist_stats["max_contig"],
                "byte_exact": new_stats["sha256"] == hist_stats["sha256"],
                "canonical_sequence_multiset_exact": (
                    new_stats["canonical_sequence_multiset_sha256"]
                    == hist_stats["canonical_sequence_multiset_sha256"]
                ),
                "new_sha256": new_stats["sha256"],
                "historical_sha256": hist_stats["sha256"],
                "new_canonical_sequence_multiset_sha256": new_stats["canonical_sequence_multiset_sha256"],
                "historical_canonical_sequence_multiset_sha256": hist_stats["canonical_sequence_multiset_sha256"],
            })
        structure_fields = list(structure_rows[0])
        with (stage / "STRUCTURE_DELTAS.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=structure_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(structure_rows)

        pair_rows: list[dict[str, Any]] = []
        command_rows: list[dict[str, Any]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(
                    run_skani_pair,
                    group,
                    args.skani,
                    paths[group]["new_assembly"],
                    paths[group]["historical_assembly"],
                    raw_dir,
                    args.threads_per_job,
                ): group
                for group in GROUPS
            }
            for future in concurrent.futures.as_completed(futures):
                pair, command = future.result()
                pair_rows.append(pair)
                command_rows.append(command)
        pair_rows.sort(key=lambda row: row["group"])
        command_rows.sort(key=lambda row: row["group"])

        with (stage / "PAIRWISE_SKANI.tsv").open("w", encoding="utf-8", newline="") as handle:
            fields = ["group", "query", "reference", "ani", "af_query", "af_ref", "query_name", "ref_name"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(pair_rows)
        write_json(stage / "COMMANDS.json", command_rows)

        exact_reproduction = all(
            row["byte_exact"] and row["canonical_sequence_multiset_exact"]
            for row in structure_rows
        )
        near_equivalent = all(
            row["ani"] >= 99.0 and row["af_query"] >= 95.0 and row["af_ref"] >= 95.0
            for row in pair_rows
        )
        same_species_level = all(
            row["ani"] >= 95.0 and row["af_query"] >= 50.0 and row["af_ref"] >= 50.0
            for row in pair_rows
        )
        new_total_contigs = sum(row["new_contigs"] for row in structure_rows)
        hist_total_contigs = sum(row["historical_contigs"] for row in structure_rows)
        new_total_bp = sum(row["new_bp"] for row in structure_rows)
        hist_total_bp = sum(row["historical_bp"] for row in structure_rows)
        threshold_counts: dict[str, int] = {}
        for threshold in (99.0, 95.0, 50.0):
            label = str(int(threshold))
            threshold_counts[f"ani_ge_{label}"] = count_at_least(pair_rows, "ani", threshold)
            threshold_counts[f"af_query_ge_{label}"] = count_at_least(pair_rows, "af_query", threshold)
            threshold_counts[f"af_ref_ge_{label}"] = count_at_least(pair_rows, "af_ref", threshold)
            threshold_counts[f"both_af_ge_{label}"] = sum(
                row["af_query"] >= threshold and row["af_ref"] >= threshold for row in pair_rows
            )
            threshold_counts[f"ani_and_both_af_ge_{label}"] = sum(
                row["ani"] >= threshold and row["af_query"] >= threshold and row["af_ref"] >= threshold
                for row in pair_rows
            )

        summary = {
            "schema": SCHEMA,
            "status": "PASS" if near_equivalent else "FAIL",
            "started_utc": started_utc,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "prerequisites": prereq,
            "upstream_replay_status": upstream_complete.get("status"),
            "upstream_replay_failed_checks": failed_upstream_checks,
            "groups": len(pair_rows),
            "claims": {
                "exact_reproduction": exact_reproduction,
                "near_equivalent": near_equivalent,
                "same_species_level": same_species_level,
                "wording": (
                    "Biologically near-equivalent but not exact; Flye 2.9.6 was run without a deterministic flag."
                    if near_equivalent and not exact_reproduction
                    else "See thresholds and pair-level results; no exactness is inferred from ANI alone."
                ),
            },
            "skani_percent": {
                "ani": percentile_summary([row["ani"] for row in pair_rows]),
                "af_query_new": percentile_summary([row["af_query"] for row in pair_rows]),
                "af_ref_historical": percentile_summary([row["af_ref"] for row in pair_rows]),
                "threshold_counts_of_58": threshold_counts,
            },
            "structure": {
                "byte_exact_groups": sum(row["byte_exact"] for row in structure_rows),
                "canonical_sequence_multiset_exact_groups": sum(row["canonical_sequence_multiset_exact"] for row in structure_rows),
                "new_total_contigs": new_total_contigs,
                "historical_total_contigs": hist_total_contigs,
                "delta_total_contigs": new_total_contigs - hist_total_contigs,
                "delta_total_contigs_percent": 100.0 * (new_total_contigs - hist_total_contigs) / hist_total_contigs,
                "new_total_bp": new_total_bp,
                "historical_total_bp": hist_total_bp,
                "delta_total_bp": new_total_bp - hist_total_bp,
                "delta_total_bp_percent": 100.0 * (new_total_bp - hist_total_bp) / hist_total_bp,
                "per_group_delta_contigs": percentile_summary([float(row["delta_contigs"]) for row in structure_rows]),
                "per_group_delta_bp": percentile_summary([float(row["delta_bp"]) for row in structure_rows]),
                "per_group_delta_n50": percentile_summary([float(row["delta_n50"]) for row in structure_rows]),
            },
            "outliers": {
                "lowest_ani": sorted(pair_rows, key=lambda row: (row["ani"], row["group"]))[:10],
                "lowest_min_directional_af": sorted(pair_rows, key=lambda row: (min(row["af_query"], row["af_ref"]), row["group"]))[:10],
                "largest_absolute_bp_delta": sorted(structure_rows, key=lambda row: (-abs(row["delta_bp"]), row["group"]))[:10],
                "largest_absolute_contig_delta": sorted(structure_rows, key=lambda row: (-abs(row["delta_contigs"]), row["group"]))[:10],
            },
            "source_bindings": {
                "new_root": str(args.new_root.resolve()),
                "new_replay_audit_sha256": sha256_file(replay_path),
                "new_complete_sha256": sha256_file(upstream_complete_path),
                "historical_root": str(args.historical_root.resolve()),
                "historical_input_root": str(args.historical_input_root.resolve()),
                "source_hashes_tsv_sha256": sha256_file(stage / "SOURCE_HASHES.tsv"),
            },
        }
        write_json(stage / "SUMMARY.json", summary)
        complete = {
            "schema": SCHEMA,
            "status": summary["status"],
            "exact_reproduction": exact_reproduction,
            "near_equivalent": near_equivalent,
            "same_species_level": same_species_level,
            "groups": 58,
            "summary_sha256": sha256_file(stage / "SUMMARY.json"),
            "pairwise_skani_sha256": sha256_file(stage / "PAIRWISE_SKANI.tsv"),
            "structure_deltas_sha256": sha256_file(stage / "STRUCTURE_DELTAS.tsv"),
            "versions_sha256": sha256_file(stage / "VERSIONS.json"),
            "preregistration_sha256": sha256_file(stage / "PREREGISTRATION.json"),
        }
        write_json(stage / "COMPLETE.json", complete)
    except Exception as exc:
        failure = {
            "schema": SCHEMA,
            "status": "FAIL",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json(stage / "FAILURE.json", failure)
        write_json(stage / "COMPLETE.json", failure)
        os.replace(stage, args.out)
        raise

    os.replace(stage, args.out)
    print(json.dumps(complete, indent=2, sort_keys=True))
    return 0 if complete["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
