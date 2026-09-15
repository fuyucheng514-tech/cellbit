#!/usr/bin/env python3
"""Strictly finalize fresh Stage-3B CheckM2 and GTDB-Tk identify outputs.

The program is intentionally a read-only validator with respect to all inputs.
It publishes one new output directory atomically and refuses to reuse or
overwrite an existing destination.  No expected SAG cardinality is hard-coded:
the exact SAG universe is established by the Dna2bit-negative pending manifest
and must be identical in every other upstream artifact.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path
import shlex
import sys
import uuid
from typing import Any, Iterable


PENDING_FIELDS = ("sag_id", "assembly_fasta", "reason")
INPUT_MANIFEST_FIELDS = ("sag_id", "source_assembly", "input_link", "reason")
STATS_FIELDS = (
    "sag_id",
    "assembly_fasta",
    "total_bp",
    "max_contig",
    "gc_pct",
    "gc_bases",
    "acgt_bases",
)
QUALITY_FIELDS = (
    "sag_id",
    "assembly_fasta",
    "total_bp",
    "max_contig",
    "gc_pct",
    "checkm2_completeness",
    "checkm2_contamination",
)
MARKER_FIELDS = ("SAG_id", "marker_id", "gene_id", "nt_len", "sequence")
CHECKM2_REQUIRED_FIELDS = {"Name", "Completeness", "Contamination"}
ALLOWED_REASONS = {
    "dna2bit_rejected_or_no_hit",
    "dna2bit_negative",
    "dna2bit_no_hit",
    "dna2bit_rejected",
}
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
MAX_METADATA_BYTES = 256 * 1024 * 1024


class FinalizeError(RuntimeError):
    """A fail-closed input, evidence, or publication error."""


@dataclass(frozen=True)
class StableBytes:
    data: bytes
    sha256: str
    size: int


@dataclass(frozen=True)
class FastaStats:
    sha256: str
    file_bytes: int
    records: int
    total_bp: int
    max_contig: int
    gc_bases: int
    acgt_bases: int
    sequences: dict[str, str] | None


def byte_key(value: str) -> bytes:
    return value.encode("utf-8")


def path_lexists(path: Path) -> bool:
    return os.path.lexists(os.fspath(path))


def stable_stat_tuple(path: Path, *, follow_symlinks: bool = True) -> tuple[int, ...]:
    stat = path.stat() if follow_symlinks else path.lstat()
    return (stat.st_dev, stat.st_ino, stat.st_mode, stat.st_size, stat.st_mtime_ns)


def read_stable_bytes(path: Path, label: str, *, max_bytes: int = MAX_METADATA_BYTES) -> StableBytes:
    if not path.is_file():
        raise FinalizeError(f"missing {label}: {path}")
    before = stable_stat_tuple(path)
    if before[3] > max_bytes:
        raise FinalizeError(f"unexpectedly large {label} ({before[3]} bytes): {path}")
    with path.open("rb") as handle:
        data = handle.read(max_bytes + 1)
    after = stable_stat_tuple(path)
    if before != after:
        raise FinalizeError(f"{label} changed while being read: {path}")
    if len(data) > max_bytes or len(data) != after[3]:
        raise FinalizeError(f"unstable or oversized {label}: {path}")
    return StableBytes(data=data, sha256=hashlib.sha256(data).hexdigest(), size=len(data))


def decode_utf8(blob: StableBytes, path: Path, label: str) -> str:
    if b"\0" in blob.data:
        raise FinalizeError(f"NUL byte in {label}: {path}")
    try:
        return blob.data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise FinalizeError(f"non-UTF-8 {label}: {path}") from error


def strict_tsv(path: Path, label: str) -> tuple[list[str], list[dict[str, str]], StableBytes]:
    blob = read_stable_bytes(path, label)
    text = decode_utf8(blob, path, label)
    lines = text.splitlines()
    if not lines:
        raise FinalizeError(f"empty {label}: {path}")
    if any(line == "" for line in lines):
        raise FinalizeError(f"blank row in {label}: {path}")
    header = lines[0].split("\t")
    if not header or any(not field for field in header) or len(header) != len(set(header)):
        raise FinalizeError(f"empty or duplicate column in {label}: {path}")
    rows: list[dict[str, str]] = []
    for row_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise FinalizeError(
                f"{label} row {row_number} has {len(fields)} columns, expected {len(header)}"
            )
        rows.append(dict(zip(header, fields)))
    return header, rows, blob


def reject_json_constant(value: str) -> None:
    raise FinalizeError(f"non-finite JSON number: {value}")


def no_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise FinalizeError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def strict_json(path: Path, label: str) -> tuple[dict[str, Any], StableBytes]:
    blob = read_stable_bytes(path, label)
    text = decode_utf8(blob, path, label)
    try:
        value = json.loads(
            text,
            object_pairs_hook=no_duplicate_json_keys,
            parse_constant=reject_json_constant,
        )
    except (json.JSONDecodeError, TypeError) as error:
        raise FinalizeError(f"invalid {label}: {path}: {error}") from error
    if not isinstance(value, dict):
        raise FinalizeError(f"{label} must be a JSON object: {path}")
    return value, blob


def strict_id(value: str, label: str) -> str:
    if (
        not value
        or value != value.strip()
        or any(character.isspace() for character in value)
        or "/" in value
        or "\\" in value
        or "\0" in value
    ):
        raise FinalizeError(f"empty or unsafe {label}: {value!r}")
    return value


def exact_existing_path(raw: object, label: str, *, directory: bool = False) -> Path:
    if isinstance(raw, os.PathLike):
        raw = os.fspath(raw)
    if not isinstance(raw, str) or not raw or "\0" in raw:
        raise FinalizeError(f"missing or invalid path for {label}")
    path = Path(raw)
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise FinalizeError(f"cannot resolve {label}: {path}") from error
    if directory:
        if not resolved.is_dir():
            raise FinalizeError(f"{label} is not a directory: {resolved}")
    elif not resolved.is_file():
        raise FinalizeError(f"{label} is not a file: {resolved}")
    return resolved


def lexical_absolute(raw: str, label: str) -> Path:
    if not raw or "\0" in raw:
        raise FinalizeError(f"invalid lexical path for {label}: {raw!r}")
    path = Path(raw)
    if not path.is_absolute():
        raise FinalizeError(f"{label} is not absolute: {path}")
    return Path(os.path.abspath(os.fspath(path)))


def integer_field(raw: str, label: str, *, minimum: int = 0) -> int:
    if not raw or raw != raw.strip() or not raw.isascii() or not raw.isdecimal():
        raise FinalizeError(f"invalid integer {label}: {raw!r}")
    value = int(raw)
    if value < minimum:
        raise FinalizeError(f"{label}={value} is below {minimum}")
    return value


def decimal_field(
    raw: str,
    label: str,
    *,
    minimum: Decimal | None = None,
    maximum: Decimal | None = None,
) -> Decimal:
    if not raw or raw != raw.strip():
        raise FinalizeError(f"empty or padded decimal {label}: {raw!r}")
    try:
        value = Decimal(raw)
    except InvalidOperation as error:
        raise FinalizeError(f"invalid decimal {label}: {raw!r}") from error
    if not value.is_finite():
        raise FinalizeError(f"non-finite decimal {label}: {raw!r}")
    if minimum is not None and value < minimum:
        raise FinalizeError(f"{label}={value} is below {minimum}")
    if maximum is not None and value > maximum:
        raise FinalizeError(f"{label}={value} exceeds {maximum}")
    return value


def decimal_text(value: Decimal) -> str:
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text if text not in {"", "-0"} else "0"


def parse_fasta(path: Path, label: str, *, retain_sequences: bool) -> FastaStats:
    if not path.is_file():
        raise FinalizeError(f"missing FASTA for {label}: {path}")
    before = stable_stat_tuple(path)
    link_before = stable_stat_tuple(path, follow_symlinks=False) if path.is_symlink() else None
    digest = hashlib.sha256()
    sequences: dict[str, str] | None = {} if retain_sequences else None
    current_id: str | None = None
    current_length = 0
    current_parts: list[str] = []
    seen_ids: set[str] = set()
    records = 0
    total_bp = 0
    max_contig = 0
    gc_bases = 0
    acgt_bases = 0

    def finish_record() -> None:
        nonlocal current_id, current_length, current_parts, max_contig
        if current_id is None:
            return
        if current_length == 0:
            raise FinalizeError(f"empty FASTA record {current_id!r} in {path}")
        max_contig = max(max_contig, current_length)
        if sequences is not None:
            sequences[current_id] = "".join(current_parts)
        current_id = None
        current_length = 0
        current_parts = []

    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            digest.update(raw_line)
            line_bytes = raw_line.rstrip(b"\r\n")
            if b"\r" in line_bytes or b"\n" in line_bytes or b"\0" in line_bytes:
                raise FinalizeError(f"control byte in FASTA {path} at line {line_number}")
            try:
                line = line_bytes.decode("ascii")
            except UnicodeDecodeError as error:
                raise FinalizeError(f"non-ASCII FASTA {path} at line {line_number}") from error
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                description = line[1:]
                if not description or description != description.strip():
                    raise FinalizeError(f"empty or padded FASTA header in {path} at line {line_number}")
                gene_id = description.split(None, 1)[0]
                if not gene_id or gene_id in seen_ids:
                    raise FinalizeError(f"empty or duplicate FASTA id {gene_id!r} in {path}")
                seen_ids.add(gene_id)
                current_id = gene_id
                records += 1
                continue
            if current_id is None:
                raise FinalizeError(f"sequence before FASTA header in {path} at line {line_number}")
            if any(character.isspace() for character in line):
                raise FinalizeError(f"whitespace inside FASTA sequence in {path} at line {line_number}")
            sequence = line.upper()
            invalid = set(sequence) - DNA_ALPHABET
            if invalid:
                raise FinalizeError(
                    f"invalid nucleotide symbols {sorted(invalid)!r} in {path} at line {line_number}"
                )
            length = len(sequence)
            total_bp += length
            current_length += length
            gc_bases += sequence.count("G") + sequence.count("C")
            acgt_bases += sum(sequence.count(base) for base in "ACGT")
            if sequences is not None:
                current_parts.append(sequence)
    finish_record()
    after = stable_stat_tuple(path)
    link_after = stable_stat_tuple(path, follow_symlinks=False) if path.is_symlink() else None
    if before != after or link_before != link_after:
        raise FinalizeError(f"FASTA changed while being read: {path}")
    if records == 0 or total_bp == 0:
        raise FinalizeError(f"empty FASTA: {path}")
    if sequences is not None and len(sequences) != records:
        raise FinalizeError(f"FASTA record closure failed: {path}")
    return FastaStats(
        sha256=digest.hexdigest(),
        file_bytes=after[3],
        records=records,
        total_bp=total_bp,
        max_contig=max_contig,
        gc_bases=gc_bases,
        acgt_bases=acgt_bases,
        sequences=sequences,
    )


def require_exact_set(label: str, observed: set[str], expected: set[str]) -> None:
    if observed != expected:
        missing = sorted(expected - observed, key=byte_key)[:20]
        extra = sorted(observed - expected, key=byte_key)[:20]
        raise FinalizeError(
            f"{label} SAG set mismatch: observed={len(observed)}, expected={len(expected)}, "
            f"first_missing={missing}, first_extra={extra}"
        )


def unique_rows(
    rows: Iterable[dict[str, str]], key_name: str, label: str
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows, 2):
        key = strict_id(row[key_name], f"{label} {key_name} at row {row_number}")
        if key in result:
            raise FinalizeError(f"duplicate {key_name} in {label}: {key}")
        result[key] = row
    if not result:
        raise FinalizeError(f"{label} contains no data rows")
    return result


def update_binding_digest(
    digest: "hashlib._Hash", fields: Iterable[object]
) -> None:
    encoded = [str(field).encode("utf-8") for field in fields]
    for field in encoded:
        digest.update(len(field).to_bytes(8, "big"))
        digest.update(field)


def write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[Iterable[object]]) -> tuple[int, str]:
    count = 0
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            writer.writerow(row)
            count += 1
        handle.flush()
        os.fsync(handle.fileno())
    return count, file_sha256(path)


def write_json(path: Path, value: object) -> str:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    return file_sha256(path)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(os.fspath(path), os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def validate_input_view(
    input_view: Path,
    pending_path: Path,
    pending_sha256: str,
    pending: dict[str, dict[str, str]],
) -> tuple[list[str], dict[str, dict[str, str]], dict[str, object]]:
    complete_path = input_view / "COMPLETE.json"
    manifest_path = input_view / "INPUT_MANIFEST.tsv"
    input_dir = input_view / "input_fna"
    complete, complete_blob = strict_json(complete_path, "input-view COMPLETE")
    header, manifest_rows, manifest_blob = strict_tsv(manifest_path, "input-view manifest")
    if tuple(header) != INPUT_MANIFEST_FIELDS:
        raise FinalizeError(f"unexpected input-view manifest schema: {header}")
    manifest = unique_rows(manifest_rows, "sag_id", "input-view manifest")
    expected_ids = set(pending)
    require_exact_set("input-view manifest", set(manifest), expected_ids)
    sag_order = [row["sag_id"] for row in manifest_rows]
    if sag_order != sorted(sag_order, key=byte_key):
        raise FinalizeError("input-view manifest is not in deterministic bytewise SAG order")

    if complete.get("schema") != "stage3b-upstream-input-view-v1" or complete.get("status") != "PASS":
        raise FinalizeError("input-view COMPLETE schema/status is not PASS v1")
    source_pending = complete.get("source_pending")
    canonical_manifest = complete.get("canonical_manifest")
    if not isinstance(source_pending, dict) or not isinstance(canonical_manifest, dict):
        raise FinalizeError("input-view COMPLETE lacks source/manifest bindings")
    if exact_existing_path(source_pending.get("path"), "recorded pending") != pending_path:
        raise FinalizeError("input-view COMPLETE pending path binding mismatch")
    if source_pending.get("sha256") != pending_sha256:
        raise FinalizeError("input-view COMPLETE pending SHA256 mismatch")
    if exact_existing_path(canonical_manifest.get("path"), "recorded input manifest") != manifest_path:
        raise FinalizeError("input-view COMPLETE manifest path binding mismatch")
    if canonical_manifest.get("sha256") != manifest_blob.sha256:
        raise FinalizeError("input-view COMPLETE manifest SHA256 mismatch")
    if exact_existing_path(complete.get("input_directory"), "recorded input directory", directory=True) != input_dir:
        raise FinalizeError("input-view COMPLETE input-directory binding mismatch")
    expected_count = len(expected_ids)
    for field in ("sag_count", "unique_sag_ids", "nonempty_resolved_targets"):
        if complete.get(field) != expected_count:
            raise FinalizeError(f"input-view COMPLETE {field} mismatch")

    for sag in sag_order:
        row = manifest[sag]
        pending_row = pending[sag]
        source = exact_existing_path(row["source_assembly"], f"input-view source for {sag}")
        pending_source = exact_existing_path(pending_row["assembly_fasta"], f"pending assembly for {sag}")
        if source != pending_source:
            raise FinalizeError(f"input-view source/pending assembly mismatch for {sag}")
        expected_link = input_dir / f"{sag}.fna"
        recorded_link = lexical_absolute(row["input_link"], f"input link for {sag}")
        if recorded_link != expected_link or not expected_link.is_symlink():
            raise FinalizeError(f"missing or incorrectly named input symlink for {sag}")
        try:
            target = expected_link.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise FinalizeError(f"broken input symlink for {sag}: {expected_link}") from error
        if target != source:
            raise FinalizeError(f"input symlink target mismatch for {sag}")
        if row["reason"] != pending_row["reason"]:
            raise FinalizeError(f"input-view reason mismatch for {sag}")
    expected_names = {f"{sag}.fna" for sag in expected_ids}
    actual_entries = list(input_dir.iterdir())
    actual_names = {entry.name for entry in actual_entries}
    if actual_names != expected_names or len(actual_entries) != expected_count:
        raise FinalizeError("input_fna directory does not contain exactly one link per pending SAG")
    if any(not entry.is_symlink() for entry in actual_entries):
        raise FinalizeError("input_fna contains a non-symlink entry")

    return sag_order, manifest, {
        "complete": {"path": str(complete_path), "sha256": complete_blob.sha256},
        "manifest": {"path": str(manifest_path), "sha256": manifest_blob.sha256},
        "input_directory": str(input_dir),
    }


def validate_gtdbtk(
    identify_root: Path,
    input_dir: Path,
    expected_ids: set[str],
) -> tuple[Path, dict[str, object]]:
    receipt_path = identify_root / "gtdbtk.json"
    failed_path = identify_root / "identify" / "gtdbtk.failed_genomes.tsv"
    marker_root = identify_root / "identify" / "intermediate_results" / "marker_genes"
    receipt, receipt_blob = strict_json(receipt_path, "GTDB-Tk receipt")
    failed_blob = read_stable_bytes(failed_path, "GTDB-Tk failed-genomes table")
    if failed_blob.data.strip():
        raise FinalizeError("GTDB-Tk failed-genomes table is not empty")
    if receipt.get("version") != "2.7.2":
        raise FinalizeError(f"GTDB-Tk version is not 2.7.2: {receipt.get('version')!r}")
    if receipt.get("database_version") != "r232":
        raise FinalizeError(
            f"GTDB-Tk database release is not r232: {receipt.get('database_version')!r}"
        )
    database_path = exact_existing_path(
        receipt.get("database_path"), "GTDB-Tk database", directory=True
    )
    steps = receipt.get("steps")
    if not isinstance(steps, list) or len(steps) != 1 or not isinstance(steps[0], dict):
        raise FinalizeError("GTDB-Tk receipt must contain exactly one identify step")
    step = steps[0]
    checks = {
        "name_identify": step.get("name") == "identify",
        "status_completed": step.get("status") == "completed",
        "extension_fna": step.get("extension") == "fna",
        "write_single_copy_genes": step.get("write_single_copy_genes") is True,
    }
    if not all(checks.values()):
        raise FinalizeError(f"GTDB-Tk identify receipt checks failed: {checks}")
    if exact_existing_path(step.get("genome_dir"), "GTDB-Tk genome_dir", directory=True) != input_dir:
        raise FinalizeError("GTDB-Tk identify genome_dir is not the audited input_fna directory")
    # GTDB-Tk 2.7.2 records the CPU count in ``command_line`` rather than in
    # the per-step JSON object.  Parse that receipt instead of assuming a
    # non-existent ``steps[0].cpus`` field.
    command_line = receipt.get("command_line")
    if not isinstance(command_line, str):
        raise FinalizeError("GTDB-Tk receipt lacks command_line")
    try:
        command = shlex.split(command_line)
    except ValueError as error:
        raise FinalizeError("malformed GTDB-Tk command_line receipt") from error
    if command.count("--cpus") != 1:
        raise FinalizeError("GTDB-Tk command_line must contain exactly one --cpus")
    cpu_index = command.index("--cpus")
    if cpu_index + 1 >= len(command):
        raise FinalizeError("GTDB-Tk command_line has no value after --cpus")
    try:
        cpus = int(command[cpu_index + 1])
    except ValueError as error:
        raise FinalizeError("GTDB-Tk command_line has a non-integer --cpus") from error
    if cpus < 1 or "identify" not in command or "--write_single_copy_genes" not in command:
        raise FinalizeError("GTDB-Tk command_line lacks the required identify contract")
    if not marker_root.is_dir() or marker_root.is_symlink():
        raise FinalizeError(f"missing or symlinked GTDB-Tk marker root: {marker_root}")
    entries = list(marker_root.iterdir())
    if any(not entry.is_dir() or entry.is_symlink() for entry in entries):
        raise FinalizeError("GTDB-Tk marker root contains a non-directory or symlink entry")
    marker_ids = {entry.name for entry in entries}
    if len(entries) != len(marker_ids):
        raise FinalizeError("duplicate marker-directory names are visible")
    require_exact_set("GTDB-Tk marker directories", marker_ids, expected_ids)
    return marker_root, {
        "gtdbtk_json": {"path": str(receipt_path), "sha256": receipt_blob.sha256},
        "failed_genomes": {"path": str(failed_path), "sha256": failed_blob.sha256},
        "version": receipt["version"],
        "database_version": receipt["database_version"],
        "database_path": str(database_path),
        "identify_cpus": cpus,
        "marker_directories": len(marker_ids),
        "receipt_checks": checks,
    }


def parse_top_hits(
    path: Path,
    sequences: dict[str, str],
    best: dict[str, tuple[Decimal, int, str]],
    encounter: int,
) -> tuple[int, StableBytes, int]:
    header, rows, blob = strict_tsv(path, "GTDB-Tk top-hit table")
    expected_header = ["Gene Id", "Top hits (Family id,e-value,bitscore)"]
    if header != expected_header:
        raise FinalizeError(f"unexpected GTDB-Tk top-hit header: {path}: {header}")
    seen_genes: set[str] = set()
    hit_count = 0
    for row_number, row in enumerate(rows, 2):
        gene = row["Gene Id"]
        if not gene or gene != gene.strip() or gene not in sequences:
            raise FinalizeError(f"unknown or padded gene id in {path} at row {row_number}: {gene!r}")
        if gene in seen_genes:
            raise FinalizeError(f"duplicate gene row in GTDB-Tk top-hit table {path}: {gene}")
        seen_genes.add(gene)
        packed_hits = row["Top hits (Family id,e-value,bitscore)"]
        if not packed_hits:
            raise FinalizeError(f"empty top-hit cell in {path} at row {row_number}")
        for raw_hit in packed_hits.split(";"):
            fields = raw_hit.split(",")
            if len(fields) != 3:
                raise FinalizeError(f"malformed marker hit in {path}: {raw_hit!r}")
            marker, raw_evalue, raw_bitscore = fields
            if not marker or marker != marker.strip() or any(char.isspace() for char in marker):
                raise FinalizeError(f"invalid marker id in {path}: {marker!r}")
            decimal_field(raw_evalue, f"e-value in {path}", minimum=Decimal("0"))
            bitscore = decimal_field(
                raw_bitscore, f"bitscore in {path}", minimum=Decimal("0")
            )
            encounter += 1
            hit_count += 1
            previous = best.get(marker)
            # Strict greater-than is deliberate.  Equal scores retain the first
            # encounter, with PFAM traversed before TIGRFAM.
            if previous is None or bitscore > previous[0]:
                best[marker] = (bitscore, encounter, gene)
    return encounter, blob, hit_count


def run(args: argparse.Namespace, staging: Path, final_root: Path) -> dict[str, Any]:
    pending_path = exact_existing_path(args.pending, "pending manifest")
    stats_path = exact_existing_path(args.stats, "FASTA stats")
    input_view = exact_existing_path(args.input_view, "input view", directory=True)
    checkm2_path = exact_existing_path(args.checkm2_report, "CheckM2 report")
    identify_root = exact_existing_path(args.identify_root, "GTDB-Tk identify root", directory=True)

    pending_header, pending_rows, pending_blob = strict_tsv(pending_path, "pending manifest")
    if tuple(pending_header) != PENDING_FIELDS:
        raise FinalizeError(f"unexpected pending-manifest schema: {pending_header}")
    pending = unique_rows(pending_rows, "sag_id", "pending manifest")
    for sag, row in pending.items():
        if row["reason"] not in ALLOWED_REASONS:
            raise FinalizeError(f"invalid pending reason for {sag}: {row['reason']!r}")
        exact_existing_path(row["assembly_fasta"], f"pending assembly for {sag}")
    expected_ids = set(pending)

    sag_order, input_manifest, input_binding = validate_input_view(
        input_view, pending_path, pending_blob.sha256, pending
    )
    input_dir = Path(input_binding["input_directory"])

    stats_header, stats_rows, stats_blob = strict_tsv(stats_path, "FASTA stats")
    if tuple(stats_header) != STATS_FIELDS:
        raise FinalizeError(f"unexpected FASTA-stats schema: {stats_header}")
    stats = unique_rows(stats_rows, "sag_id", "FASTA stats")
    require_exact_set("FASTA stats", set(stats), expected_ids)

    checkm2_header, checkm2_rows, checkm2_blob = strict_tsv(checkm2_path, "CheckM2 report")
    missing_columns = CHECKM2_REQUIRED_FIELDS - set(checkm2_header)
    if missing_columns:
        raise FinalizeError(f"CheckM2 report lacks columns: {sorted(missing_columns)}")
    checkm2 = unique_rows(checkm2_rows, "Name", "CheckM2 report")
    require_exact_set("CheckM2 report", set(checkm2), expected_ids)

    quality_rows: list[tuple[object, ...]] = []
    fasta_binding = hashlib.sha256()
    total_input_records = 0
    total_input_bp = 0
    for sag in sag_order:
        pending_source = exact_existing_path(
            pending[sag]["assembly_fasta"], f"pending assembly for {sag}"
        )
        stats_source = exact_existing_path(
            stats[sag]["assembly_fasta"], f"stats assembly for {sag}"
        )
        manifest_source = exact_existing_path(
            input_manifest[sag]["source_assembly"], f"manifest source for {sag}"
        )
        if pending_source != stats_source or pending_source != manifest_source:
            raise FinalizeError(f"assembly path binding differs across inputs for {sag}")
        link = input_dir / f"{sag}.fna"
        observed = parse_fasta(link, f"input SAG {sag}", retain_sequences=False)
        if observed.total_bp < 1000:
            raise FinalizeError(f"pending SAG is below the frozen 1000-bp gate: {sag}")
        declared_total = integer_field(stats[sag]["total_bp"], f"stats total_bp for {sag}", minimum=1)
        declared_max = integer_field(stats[sag]["max_contig"], f"stats max_contig for {sag}", minimum=1)
        declared_gc_bases = integer_field(stats[sag]["gc_bases"], f"stats gc_bases for {sag}")
        declared_acgt = integer_field(stats[sag]["acgt_bases"], f"stats acgt_bases for {sag}")
        declared_gc_pct = decimal_field(
            stats[sag]["gc_pct"], f"stats gc_pct for {sag}",
            minimum=Decimal("0"), maximum=Decimal("100"),
        )
        if (
            declared_gc_bases > declared_acgt
            or declared_acgt > declared_total
            or declared_max > declared_total
        ):
            raise FinalizeError(f"internally inconsistent declared FASTA stats for {sag}")
        exact_counts = (
            observed.total_bp == declared_total
            and observed.max_contig == declared_max
            and observed.gc_bases == declared_gc_bases
            and observed.acgt_bases == declared_acgt
        )
        if not exact_counts:
            raise FinalizeError(
                f"FASTA/stat count mismatch for {sag}: observed="
                f"({observed.total_bp},{observed.max_contig},{observed.gc_bases},{observed.acgt_bases}) "
                f"declared=({declared_total},{declared_max},{declared_gc_bases},{declared_acgt})"
            )
        calculated_gc_pct = Decimal(observed.gc_bases) * Decimal(100) / Decimal(observed.total_bp)
        if abs(declared_gc_pct - calculated_gc_pct) > Decimal("0.0000051"):
            raise FinalizeError(
                f"FASTA/stat gc_pct mismatch for {sag}: declared={declared_gc_pct}, "
                f"calculated={calculated_gc_pct}"
            )
        completeness = decimal_field(
            checkm2[sag]["Completeness"], f"CheckM2 completeness for {sag}",
            minimum=Decimal("0"), maximum=Decimal("100"),
        )
        contamination = decimal_field(
            checkm2[sag]["Contamination"], f"CheckM2 contamination for {sag}",
            minimum=Decimal("0"),
        )
        quality_rows.append(
            (
                sag,
                str(pending_source),
                declared_total,
                declared_max,
                decimal_text(declared_gc_pct),
                decimal_text(completeness),
                decimal_text(contamination),
            )
        )
        update_binding_digest(
            fasta_binding,
            (sag, pending_source, observed.file_bytes, observed.sha256, observed.records, observed.total_bp),
        )
        total_input_records += observed.records
        total_input_bp += observed.total_bp

    marker_root, gtdb_binding = validate_gtdbtk(identify_root, input_dir, expected_ids)
    marker_path = staging / "bac120_marker_nt_map.tsv"
    marker_rows = 0
    raw_marker_hits = 0
    protein_records = 0
    sags_with_selected_markers = 0
    marker_source_binding = hashlib.sha256()
    with marker_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MARKER_FIELDS)
        for sag in sag_order:
            directory = marker_root / sag
            protein_path = directory / f"{sag}_protein.fna"
            protein = parse_fasta(protein_path, f"GTDB-Tk nucleotide genes for {sag}", retain_sequences=True)
            assert protein.sequences is not None
            protein_records += protein.records
            update_binding_digest(
                marker_source_binding,
                (
                    sag,
                    protein_path.relative_to(identify_root),
                    protein.file_bytes,
                    protein.sha256,
                ),
            )
            best: dict[str, tuple[Decimal, int, str]] = {}
            encounter = 0
            for family in ("pfam", "tigrfam"):
                hit_path = directory / f"{sag}_{family}_tophit.tsv"
                encounter, hit_blob, hit_count = parse_top_hits(
                    hit_path, protein.sequences, best, encounter
                )
                raw_marker_hits += hit_count
                update_binding_digest(
                    marker_source_binding,
                    (sag, hit_path.relative_to(identify_root), hit_blob.size, hit_blob.sha256),
                )
            for marker in sorted(best, key=byte_key):
                gene = best[marker][2]
                sequence = protein.sequences[gene]
                if not sequence or set(sequence) - DNA_ALPHABET:
                    raise FinalizeError(f"invalid selected marker sequence for {sag}/{marker}")
                writer.writerow((sag, marker, gene, len(sequence), sequence))
                marker_rows += 1
            if best:
                sags_with_selected_markers += 1
        handle.flush()
        os.fsync(handle.fileno())
    marker_sha256 = file_sha256(marker_path)

    quality_path = staging / "stage3b_quality.tsv"
    quality_count, quality_sha256 = write_tsv(quality_path, QUALITY_FIELDS, quality_rows)
    if quality_count != len(expected_ids):
        raise FinalizeError("quality output cardinality closure failed")

    script_path = Path(__file__).resolve(strict=True)
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
    }
    audit: dict[str, Any] = {
        "schema": "stage3b-upstream-finalization-audit-v1",
        "status": "PASS",
        "checks": checks,
        "script": {"path": str(script_path), "sha256": file_sha256(script_path)},
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
    audit_sha256 = write_json(audit_path, audit)
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
    complete_sha256 = write_json(staging / "COMPLETE.json", complete)
    fsync_directory(staging)
    return {
        "status": "PASS",
        "sags": len(expected_ids),
        "quality_rows": quality_count,
        "marker_rows": marker_rows,
        "audit_sha256": audit_sha256,
        "complete_sha256": complete_sha256,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Atomically finalize exact Stage-3B CheckM2 and GTDB-Tk upstream evidence"
    )
    parser.add_argument("--pending", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--input-view", type=Path, required=True)
    parser.add_argument("--checkm2-report", type=Path, required=True)
    parser.add_argument("--identify-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    final_root = Path(os.path.abspath(os.fspath(args.out_dir)))
    if path_lexists(final_root):
        raise FinalizeError(f"write-once output already exists: {final_root}")
    final_root.parent.mkdir(parents=True, exist_ok=True)
    staging = final_root.parent / (
        f".{final_root.name}.incomplete-{os.getpid()}-{uuid.uuid4().hex}"
    )
    staging.mkdir(exist_ok=False)
    try:
        summary = run(args, staging, final_root)
        if path_lexists(final_root):
            raise FinalizeError(f"destination appeared before atomic publication: {final_root}")
        os.rename(staging, final_root)
        fsync_directory(final_root.parent)
    except BaseException as error:
        failure_path = staging / "FAILED.json"
        if staging.is_dir() and not path_lexists(failure_path):
            try:
                write_json(
                    failure_path,
                    {
                        "schema": "stage3b-upstream-finalization-failure-v1",
                        "status": "FAIL",
                        "error": repr(error),
                    },
                )
                fsync_directory(staging)
            except BaseException:
                pass
        raise
    summary["out_dir"] = str(final_root)
    print(json.dumps(summary, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FinalizeError as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
