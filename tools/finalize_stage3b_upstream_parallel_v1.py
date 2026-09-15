#!/usr/bin/env python3
"""Parallel, deterministic Stage-3B upstream finalizer.

This is a conservative acceleration of ``finalize_stage3b_upstream.py``.  A
worker validates one SAG at a time (assembly FASTA, CheckM2 row, marker FASTA,
PFAM hits, and TIGRFAM hits).  The parent process consumes worker results in
the frozen input-manifest order and is the *only* process that writes output.
Consequently both scientific TSVs are byte-for-byte identical to the serial
implementation, independent of worker count.

The legacy module remains the single source of parsing and validation rules;
this file changes scheduling only.  Inputs are read-only and the destination
is still write-once and atomically published.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from dataclasses import dataclass
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
import uuid
from typing import Any, Iterable, Iterator

import finalize_stage3b_upstream as serial


@dataclass(frozen=True)
class SagTask:
    ordinal: int
    sag: str
    pending_assembly: str
    stats_assembly: str
    manifest_assembly: str
    input_link: str
    stats_row: dict[str, str]
    checkm2_row: dict[str, str]
    identify_root: str
    marker_directory: str


@dataclass(frozen=True)
class SagResult:
    ordinal: int
    sag: str
    quality_row: tuple[object, ...]
    fasta_binding_fields: tuple[object, ...]
    input_records: int
    input_bp: int
    protein_records: int
    raw_marker_hits: int
    marker_rows: tuple[tuple[object, ...], ...]
    marker_source_binding_fields: tuple[tuple[object, ...], ...]


def _process_sag(task: SagTask) -> SagResult:
    """Validate and parse one SAG without writing shared output."""
    sag = task.sag
    pending_source = serial.exact_existing_path(
        task.pending_assembly, f"pending assembly for {sag}"
    )
    stats_source = serial.exact_existing_path(
        task.stats_assembly, f"stats assembly for {sag}"
    )
    manifest_source = serial.exact_existing_path(
        task.manifest_assembly, f"manifest source for {sag}"
    )
    if pending_source != stats_source or pending_source != manifest_source:
        raise serial.FinalizeError(f"assembly path binding differs across inputs for {sag}")

    link = Path(task.input_link)
    observed = serial.parse_fasta(link, f"input SAG {sag}", retain_sequences=False)
    if observed.total_bp < 1000:
        raise serial.FinalizeError(f"pending SAG is below the frozen 1000-bp gate: {sag}")

    stats_row = task.stats_row
    declared_total = serial.integer_field(
        stats_row["total_bp"], f"stats total_bp for {sag}", minimum=1
    )
    declared_max = serial.integer_field(
        stats_row["max_contig"], f"stats max_contig for {sag}", minimum=1
    )
    declared_gc_bases = serial.integer_field(
        stats_row["gc_bases"], f"stats gc_bases for {sag}"
    )
    declared_acgt = serial.integer_field(
        stats_row["acgt_bases"], f"stats acgt_bases for {sag}"
    )
    declared_gc_pct = serial.decimal_field(
        stats_row["gc_pct"],
        f"stats gc_pct for {sag}",
        minimum=Decimal("0"),
        maximum=Decimal("100"),
    )
    if (
        declared_gc_bases > declared_acgt
        or declared_acgt > declared_total
        or declared_max > declared_total
    ):
        raise serial.FinalizeError(f"internally inconsistent declared FASTA stats for {sag}")
    if (
        observed.total_bp != declared_total
        or observed.max_contig != declared_max
        or observed.gc_bases != declared_gc_bases
        or observed.acgt_bases != declared_acgt
    ):
        raise serial.FinalizeError(
            f"FASTA/stat count mismatch for {sag}: observed="
            f"({observed.total_bp},{observed.max_contig},{observed.gc_bases},{observed.acgt_bases}) "
            f"declared=({declared_total},{declared_max},{declared_gc_bases},{declared_acgt})"
        )
    calculated_gc_pct = Decimal(observed.gc_bases) * Decimal(100) / Decimal(observed.total_bp)
    if abs(declared_gc_pct - calculated_gc_pct) > Decimal("0.0000051"):
        raise serial.FinalizeError(
            f"FASTA/stat gc_pct mismatch for {sag}: declared={declared_gc_pct}, "
            f"calculated={calculated_gc_pct}"
        )

    completeness = serial.decimal_field(
        task.checkm2_row["Completeness"],
        f"CheckM2 completeness for {sag}",
        minimum=Decimal("0"),
        maximum=Decimal("100"),
    )
    contamination = serial.decimal_field(
        task.checkm2_row["Contamination"],
        f"CheckM2 contamination for {sag}",
        minimum=Decimal("0"),
    )
    quality_row = (
        sag,
        str(pending_source),
        declared_total,
        declared_max,
        serial.decimal_text(declared_gc_pct),
        serial.decimal_text(completeness),
        serial.decimal_text(contamination),
    )
    fasta_binding_fields = (
        sag,
        pending_source,
        observed.file_bytes,
        observed.sha256,
        observed.records,
        observed.total_bp,
    )

    identify_root = Path(task.identify_root)
    directory = Path(task.marker_directory)
    protein_path = directory / f"{sag}_protein.fna"
    protein = serial.parse_fasta(
        protein_path, f"GTDB-Tk nucleotide genes for {sag}", retain_sequences=True
    )
    assert protein.sequences is not None
    source_bindings: list[tuple[object, ...]] = [
        (
            sag,
            protein_path.relative_to(identify_root),
            protein.file_bytes,
            protein.sha256,
        )
    ]
    best: dict[str, tuple[Decimal, int, str]] = {}
    encounter = 0
    raw_marker_hits = 0
    for family in ("pfam", "tigrfam"):
        hit_path = directory / f"{sag}_{family}_tophit.tsv"
        encounter, hit_blob, hit_count = serial.parse_top_hits(
            hit_path, protein.sequences, best, encounter
        )
        raw_marker_hits += hit_count
        source_bindings.append(
            (
                sag,
                hit_path.relative_to(identify_root),
                hit_blob.size,
                hit_blob.sha256,
            )
        )

    marker_rows: list[tuple[object, ...]] = []
    for marker in sorted(best, key=serial.byte_key):
        gene = best[marker][2]
        sequence = protein.sequences[gene]
        if not sequence or set(sequence) - serial.DNA_ALPHABET:
            raise serial.FinalizeError(f"invalid selected marker sequence for {sag}/{marker}")
        marker_rows.append((sag, marker, gene, len(sequence), sequence))

    return SagResult(
        ordinal=task.ordinal,
        sag=sag,
        quality_row=quality_row,
        fasta_binding_fields=fasta_binding_fields,
        input_records=observed.records,
        input_bp=observed.total_bp,
        protein_records=protein.records,
        raw_marker_hits=raw_marker_hits,
        marker_rows=tuple(marker_rows),
        marker_source_binding_fields=tuple(source_bindings),
    )


def _ordered_bounded_results(
    tasks: Iterable[SagTask], workers: int
) -> Iterator[SagResult]:
    """Run a bounded number of processes while yielding strict task order."""
    task_iter = iter(tasks)
    window = max(workers, min(workers * 2, 256))
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        pending: list[tuple[SagTask, concurrent.futures.Future[SagResult]]] = []
        for _ in range(window):
            try:
                task = next(task_iter)
            except StopIteration:
                break
            pending.append((task, executor.submit(_process_sag, task)))
        while pending:
            task, future = pending.pop(0)
            result = future.result()
            if result.ordinal != task.ordinal or result.sag != task.sag:
                raise serial.FinalizeError(
                    f"worker result identity/order mismatch for ordinal {task.ordinal}: {task.sag}"
                )
            yield result
            try:
                next_task = next(task_iter)
            except StopIteration:
                continue
            pending.append((next_task, executor.submit(_process_sag, next_task)))


def run(
    args: argparse.Namespace, staging: Path, final_root: Path
) -> dict[str, Any]:
    pending_path = serial.exact_existing_path(args.pending, "pending manifest")
    stats_path = serial.exact_existing_path(args.stats, "FASTA stats")
    input_view = serial.exact_existing_path(args.input_view, "input view", directory=True)
    checkm2_path = serial.exact_existing_path(args.checkm2_report, "CheckM2 report")
    identify_root = serial.exact_existing_path(
        args.identify_root, "GTDB-Tk identify root", directory=True
    )

    pending_header, pending_rows, pending_blob = serial.strict_tsv(
        pending_path, "pending manifest"
    )
    if tuple(pending_header) != serial.PENDING_FIELDS:
        raise serial.FinalizeError(f"unexpected pending-manifest schema: {pending_header}")
    pending = serial.unique_rows(pending_rows, "sag_id", "pending manifest")
    for sag, row in pending.items():
        if row["reason"] not in serial.ALLOWED_REASONS:
            raise serial.FinalizeError(f"invalid pending reason for {sag}: {row['reason']!r}")
        serial.exact_existing_path(row["assembly_fasta"], f"pending assembly for {sag}")
    expected_ids = set(pending)

    sag_order, input_manifest, input_binding = serial.validate_input_view(
        input_view, pending_path, pending_blob.sha256, pending
    )
    input_dir = Path(input_binding["input_directory"])

    stats_header, stats_rows, stats_blob = serial.strict_tsv(stats_path, "FASTA stats")
    if tuple(stats_header) != serial.STATS_FIELDS:
        raise serial.FinalizeError(f"unexpected FASTA-stats schema: {stats_header}")
    stats = serial.unique_rows(stats_rows, "sag_id", "FASTA stats")
    serial.require_exact_set("FASTA stats", set(stats), expected_ids)

    checkm2_header, checkm2_rows, checkm2_blob = serial.strict_tsv(
        checkm2_path, "CheckM2 report"
    )
    missing_columns = serial.CHECKM2_REQUIRED_FIELDS - set(checkm2_header)
    if missing_columns:
        raise serial.FinalizeError(f"CheckM2 report lacks columns: {sorted(missing_columns)}")
    checkm2 = serial.unique_rows(checkm2_rows, "Name", "CheckM2 report")
    serial.require_exact_set("CheckM2 report", set(checkm2), expected_ids)

    marker_root, gtdb_binding = serial.validate_gtdbtk(
        identify_root, input_dir, expected_ids
    )
    tasks = [
        SagTask(
            ordinal=ordinal,
            sag=sag,
            pending_assembly=pending[sag]["assembly_fasta"],
            stats_assembly=stats[sag]["assembly_fasta"],
            manifest_assembly=input_manifest[sag]["source_assembly"],
            input_link=str(input_dir / f"{sag}.fna"),
            stats_row=stats[sag],
            checkm2_row=checkm2[sag],
            identify_root=str(identify_root),
            marker_directory=str(marker_root / sag),
        )
        for ordinal, sag in enumerate(sag_order)
    ]

    quality_rows: list[tuple[object, ...]] = []
    fasta_binding = hashlib.sha256()
    marker_source_binding = hashlib.sha256()
    total_input_records = 0
    total_input_bp = 0
    protein_records = 0
    raw_marker_hits = 0
    marker_rows = 0
    sags_with_selected_markers = 0
    marker_path = staging / "bac120_marker_nt_map.tsv"
    with marker_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(serial.MARKER_FIELDS)
        for expected_ordinal, result in enumerate(
            _ordered_bounded_results(tasks, args.workers)
        ):
            if result.ordinal != expected_ordinal or result.sag != sag_order[expected_ordinal]:
                raise serial.FinalizeError(
                    f"parent merge order mismatch at ordinal {expected_ordinal}"
                )
            quality_rows.append(result.quality_row)
            serial.update_binding_digest(fasta_binding, result.fasta_binding_fields)
            total_input_records += result.input_records
            total_input_bp += result.input_bp
            protein_records += result.protein_records
            raw_marker_hits += result.raw_marker_hits
            for fields in result.marker_source_binding_fields:
                serial.update_binding_digest(marker_source_binding, fields)
            for row in result.marker_rows:
                writer.writerow(row)
                marker_rows += 1
            if result.marker_rows:
                sags_with_selected_markers += 1
        handle.flush()
        os.fsync(handle.fileno())
    if len(quality_rows) != len(sag_order):
        raise serial.FinalizeError("parallel SAG result cardinality closure failed")
    marker_sha256 = serial.file_sha256(marker_path)

    quality_path = staging / "stage3b_quality.tsv"
    quality_count, quality_sha256 = serial.write_tsv(
        quality_path, serial.QUALITY_FIELDS, quality_rows
    )
    if quality_count != len(expected_ids):
        raise serial.FinalizeError("quality output cardinality closure failed")

    script_path = Path(__file__).resolve(strict=True)
    serial_script_path = Path(serial.__file__).resolve(strict=True)
    output_records = {
        "stage3b_quality": {
            "path": str(final_root / quality_path.name),
            "rows": quality_count,
            "file_bytes": quality_path.stat().st_size,
            "sha256": quality_sha256,
        },
        "bac120_marker_nt_map": {
            "path": str(final_root / marker_path.name),
            "rows": marker_rows,
            "file_bytes": marker_path.stat().st_size,
            "sha256": marker_sha256,
        },
    }
    checks = {
        "pending_nonempty_unique_safe_ids": True,
        "pending_reasons_exact": True,
        "pending_input_view_stats_checkm2_sets_identical": True,
        "input_view_complete_pass_and_hash_bound": True,
        "input_view_exact_symlink_closure": True,
        "all_input_fastas_valid_and_stats_exact": True,
        "all_checkm2_metrics_finite": True,
        "gtdbtk_version_2_7_2": True,
        "gtdbtk_release_r232": True,
        "gtdbtk_identify_completed_on_exact_input_view": True,
        "gtdbtk_failed_genomes_empty": True,
        "gtdbtk_marker_directory_set_exact": True,
        "all_marker_fastas_and_top_hits_valid": True,
        "marker_selection_strict_max_tie_first": True,
        "quality_output_exactly_one_row_per_sag": True,
        "parallel_results_merged_in_manifest_order": True,
    }
    audit: dict[str, Any] = {
        "schema": "stage3b-upstream-finalization-audit-v1",
        "status": "PASS",
        "checks": checks,
        "script": {"path": str(script_path), "sha256": serial.file_sha256(script_path)},
        "validation_rules": {
            "path": str(serial_script_path),
            "sha256": serial.file_sha256(serial_script_path),
        },
        "parallel_execution": {
            "workers": args.workers,
            "unit": "one complete SAG per task",
            "publication": "parent-only strict input-manifest order",
            "bounded_result_window": max(args.workers, min(args.workers * 2, 256)),
        },
        "inputs": {
            "pending": {"path": str(pending_path), "sha256": pending_blob.sha256},
            "stats": {"path": str(stats_path), "sha256": stats_blob.sha256},
            "input_view": input_binding,
            "checkm2_report": {"path": str(checkm2_path), "sha256": checkm2_blob.sha256},
            "gtdbtk_identify": gtdb_binding,
        },
        "counts": {
            "sags": len(expected_ids),
            "input_fasta_records": total_input_records,
            "input_fasta_bp": total_input_bp,
            "protein_fasta_records": protein_records,
            "raw_marker_hits": raw_marker_hits,
            "selected_sag_marker_rows": marker_rows,
            "sags_with_selected_markers": sags_with_selected_markers,
        },
        "aggregate_content_bindings": {
            "input_fastas": {
                "algorithm": "SHA256(length-prefixed sag,path,file_bytes,file_sha256,records,total_bp records in input-manifest order)",
                "sha256": fasta_binding.hexdigest(),
            },
            "marker_sources": {
                "algorithm": "SHA256(length-prefixed sag,relative_path,file_bytes,file_sha256 records; protein then PFAM then TIGRFAM per SAG)",
                "sha256": marker_source_binding.hexdigest(),
            },
        },
        "marker_selection": {
            "family_encounter_order": ["pfam", "tigrfam"],
            "rule": "strictly greatest finite bitscore per (SAG,marker); an exact tie keeps the first encounter",
            "output_order": "input-view bytewise SAG order, then bytewise marker_id order",
        },
        "outputs": output_records,
    }
    audit_path = staging / "AUDIT.json"
    audit_sha256 = serial.write_json(audit_path, audit)
    complete = {
        "schema": "stage3b-upstream-finalization-complete-v1",
        "status": "PASS",
        "sag_count": len(expected_ids),
        "quality_rows": quality_count,
        "marker_rows": marker_rows,
        "audit_sha256": audit_sha256,
        "pending_sha256": pending_blob.sha256,
        "stats_sha256": stats_blob.sha256,
        "input_view_complete_sha256": input_binding["complete"]["sha256"],
        "input_view_manifest_sha256": input_binding["manifest"]["sha256"],
        "checkm2_report_sha256": checkm2_blob.sha256,
        "gtdbtk_json_sha256": gtdb_binding["gtdbtk_json"]["sha256"],
        "failed_genomes_sha256": gtdb_binding["failed_genomes"]["sha256"],
        "stage3b_quality_sha256": quality_sha256,
        "bac120_marker_nt_map_sha256": marker_sha256,
        "input_fasta_binding_sha256": fasta_binding.hexdigest(),
        "marker_source_binding_sha256": marker_source_binding.hexdigest(),
    }
    complete_sha256 = serial.write_json(staging / "COMPLETE.json", complete)
    serial.fsync_directory(staging)
    return {
        "status": "PASS",
        "sags": len(expected_ids),
        "quality_rows": quality_count,
        "marker_rows": marker_rows,
        "workers": args.workers,
        "audit_sha256": audit_sha256,
        "complete_sha256": complete_sha256,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Atomically finalize Stage-3B upstream evidence with ordered SAG-parallel parsing"
    )
    parser.add_argument("--pending", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--input-view", type=Path, required=True)
    parser.add_argument("--checkm2-report", type=Path, required=True)
    parser.add_argument("--identify-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(32, os.cpu_count() or 1),
        help="SAG parser processes (default: min(32, detected CPUs))",
    )
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 256:
        raise serial.FinalizeError("--workers must be in [1,256]")

    final_root = Path(os.path.abspath(os.fspath(args.out_dir)))
    if serial.path_lexists(final_root):
        raise serial.FinalizeError(f"write-once output already exists: {final_root}")
    final_root.parent.mkdir(parents=True, exist_ok=True)
    staging = final_root.parent / (
        f".{final_root.name}.incomplete-{os.getpid()}-{uuid.uuid4().hex}"
    )
    staging.mkdir(exist_ok=False)
    try:
        summary = run(args, staging, final_root)
        if serial.path_lexists(final_root):
            raise serial.FinalizeError(
                f"destination appeared before atomic publication: {final_root}"
            )
        os.rename(staging, final_root)
        serial.fsync_directory(final_root.parent)
    except BaseException as error:
        failure_path = staging / "FAILED.json"
        if staging.is_dir() and not serial.path_lexists(failure_path):
            try:
                serial.write_json(
                    failure_path,
                    {
                        "schema": "stage3b-upstream-finalization-failure-v1",
                        "status": "FAIL",
                        "error": repr(error),
                    },
                )
                serial.fsync_directory(staging)
            except BaseException:
                pass
        raise
    summary["out_dir"] = str(final_root)
    print(json.dumps(summary, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except serial.FinalizeError as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
