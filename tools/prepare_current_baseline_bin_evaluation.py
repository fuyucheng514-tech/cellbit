#!/usr/bin/env python3
"""Build a write-once CheckM2 input view from current Stage3A/Stage3B bins.

This command does not run CheckM2.  It validates the two producer manifests,
closes their declared SAG membership against the deterministic merged FASTA
headers (and against producer membership sidecars when those are available),
then atomically publishes one directory of non-empty ``*.fna`` symlinks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


STAGE3A_SCHEMA = ("species_group", "sag_count", "input_fasta", "bin_fasta")
STAGE3B_SCHEMA = ("cluster", "n_sags", "n_found", "input_fasta", "bin_fasta")
STAGE3A_LABEL_SCHEMA = ("sag_id", "reference", "taxonomy", "species_group")
STAGE3B_MEMBERSHIP_SCHEMA = ("cluster", "SAG_id", "cluster_size")
STAGE3B_NODE_SCHEMA = ("sag_id", "assembly_fasta", "gc_pct")
STAGE3B_UNAGGREGATED_SCHEMA = ("sag_id", "assembly_fasta", "reason")
STAGE3B_QUALITY_AUDIT_SCHEMA = (
    "sag_id",
    "assembly_fasta",
    "total_bp",
    "max_contig",
    "gc_pct",
    "checkm2_contamination",
    "stage3b_quality_pass",
    "reason",
)
SAFE_ID = re.compile(r"[A-Za-z0-9_.-]+\Z")
POSITIVE_INTEGER = re.compile(r"[1-9][0-9]*\Z")
SEQUENCE_BYTES = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-.*")


class ValidationError(RuntimeError):
    """Raised when the producer outputs cannot be closed safely."""


@dataclass(frozen=True)
class FastaEvidence:
    path: Path
    file_bytes: int
    sha256: str
    records: int
    bases: int
    members: frozenset[str]


@dataclass(frozen=True)
class BinRecord:
    bin_id: str
    module: str
    source_id: str
    declared_sags: int
    members: frozenset[str]
    membership_source: str
    input_fasta: FastaEvidence
    bin_fasta: FastaEvidence
    source_table: Path


@dataclass(frozen=True)
class BinInspectionSpec:
    """Immutable inputs for one independently inspectable output bin."""

    bin_id: str
    module: str
    source_id: str
    declared_sags: int
    input_path: Path
    bin_path: Path
    delimiter: str
    sidecar_members: frozenset[str] | None
    source_table: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_members(members: Iterable[str]) -> str:
    payload = "".join(f"{member}\n" for member in sorted(members)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def require_source_file(path: Path, label: str) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise ValidationError(f"{label} is missing or unresolved: {path}: {error}") from error
    if not resolved.is_file():
        raise ValidationError(f"{label} is not a regular file: {resolved}")
    if resolved.stat().st_size == 0:
        raise ValidationError(f"{label} is empty: {resolved}")
    return resolved


def resolve_row_path(table: Path, raw: str, label: str, row_number: int) -> Path:
    if not raw or raw != raw.strip() or "\x00" in raw:
        raise ValidationError(f"invalid {label} at {table}:{row_number}: {raw!r}")
    candidate = Path(raw)
    if not candidate.is_absolute():
        candidate = table.parent / candidate
    return require_source_file(candidate, f"{label} at row {row_number}")


def read_strict_tsv(path: Path, expected_schema: tuple[str, ...]) -> list[tuple[int, list[str]]]:
    source = require_source_file(path, "TSV")
    rows: list[tuple[int, list[str]]] = []
    try:
        # Read physical LF-delimited records first.  The teacher Dna2bit CSV
        # carries a terminal CR in its taxonomy value; the C++ labels writer
        # preserves that byte immediately before the following TAB.  Feeding
        # ``taxonomy\r\tgroup`` directly to csv.reader makes the bare CR look
        # like a new record.  Normalize only this field-terminal producer
        # artifact (and a conventional final CR); every other CR remains fatal.
        physical_rows: list[list[str]] = []
        with source.open("rb") as handle:
            for line_number, raw in enumerate(handle, 1):
                if raw.endswith(b"\n"):
                    raw = raw[:-1]
                if raw.endswith(b"\r"):
                    raw = raw[:-1]
                raw = raw.replace(b"\r\t", b"\t")
                if b"\r" in raw or b"\n" in raw or b"\x00" in raw:
                    raise ValidationError(f"unexpected control byte at {source}:{line_number}")
                text = raw.decode("utf-8")
                parsed = list(csv.reader([text], delimiter="\t", strict=True))
                if len(parsed) != 1:
                    raise ValidationError(f"ambiguous TSV record at {source}:{line_number}")
                physical_rows.append(parsed[0])
        if not physical_rows:
            raise ValidationError(f"empty TSV: {source}")
        header = physical_rows[0]
        if tuple(header) != expected_schema:
            raise ValidationError(
                f"unexpected schema in {source}; expected {expected_schema!r}, got {tuple(header)!r}"
            )
        for row_number, row in enumerate(physical_rows[1:], 2):
            if len(row) != len(expected_schema):
                raise ValidationError(
                    f"wrong column count at {source}:{row_number}; "
                    f"expected {len(expected_schema)}, got {len(row)}"
                )
            if not any(row):
                raise ValidationError(f"blank row at {source}:{row_number}")
            rows.append((row_number, row))
    except UnicodeDecodeError as error:
        raise ValidationError(f"TSV is not strict UTF-8: {source}: {error}") from error
    return rows


def validate_identifier(value: str, label: str, path: Path, row_number: int) -> str:
    if value != value.strip() or not SAFE_ID.fullmatch(value) or value in {".", ".."}:
        raise ValidationError(f"unsafe {label} at {path}:{row_number}: {value!r}")
    return value


def parse_positive_integer(value: str, label: str, path: Path, row_number: int) -> int:
    if not POSITIVE_INTEGER.fullmatch(value):
        raise ValidationError(f"invalid positive {label} at {path}:{row_number}: {value!r}")
    return int(value)


def membership_from_header(header: str, delimiter: str, expected: frozenset[str] | None) -> str:
    if expected is None:
        if delimiter not in header:
            raise ValidationError(f"merged FASTA header lacks {delimiter!r} SAG delimiter: {header!r}")
        member = header.split(delimiter, 1)[0]
        if not SAFE_ID.fullmatch(member) or member in {".", ".."}:
            raise ValidationError(f"merged FASTA contains an unsafe SAG prefix: {member!r}")
        return member

    # Stage3B permits underscores in a SAG id, including theoretically ``__``.
    # Test every delimiter boundary and take the longest declared prefix so the
    # evidence remains unambiguous even in that uncommon case.
    matches: list[str] = []
    offset = header.find(delimiter)
    while offset >= 0:
        candidate = header[:offset]
        if candidate in expected:
            matches.append(candidate)
        offset = header.find(delimiter, offset + 1)
    if not matches:
        raise ValidationError(f"merged FASTA header has no declared SAG prefix: {header!r}")
    member = max(matches, key=len)
    if sum(1 for candidate in matches if len(candidate) == len(member)) != 1:
        raise ValidationError(f"ambiguous SAG prefix in merged FASTA header: {header!r}")
    return member


def inspect_fasta(
    path: Path,
    label: str,
    *,
    member_delimiter: str | None = None,
    expected_members: frozenset[str] | None = None,
) -> FastaEvidence:
    source = require_source_file(path, label)
    digest = hashlib.sha256()
    records = 0
    bases = 0
    current_bases = 0
    members: set[str] = set()
    saw_header = False
    with source.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            digest.update(raw)
            line = raw.rstrip(b"\r\n")
            if b"\r" in line or b"\n" in line or b"\x00" in line:
                raise ValidationError(f"invalid control byte in {label} at {source}:{line_number}")
            if not line:
                continue
            if line.startswith(b">"):
                if saw_header and current_bases == 0:
                    raise ValidationError(f"zero-length FASTA record in {label} at {source}:{line_number - 1}")
                try:
                    header = line[1:].decode("utf-8")
                except UnicodeDecodeError as error:
                    raise ValidationError(
                        f"non-UTF-8 FASTA header in {label} at {source}:{line_number}"
                    ) from error
                if not header.strip():
                    raise ValidationError(f"empty FASTA header in {label} at {source}:{line_number}")
                records += 1
                saw_header = True
                current_bases = 0
                if member_delimiter is not None:
                    members.add(membership_from_header(header, member_delimiter, expected_members))
                continue
            if not saw_header:
                raise ValidationError(f"sequence precedes first FASTA header in {label}: {source}:{line_number}")
            # Only ASCII FASTA whitespace is ignorable.  Non-ASCII bytes must
            # reach the alphabet check below and fail rather than disappearing.
            sequence = bytes(byte for byte in line if byte not in b" \t\v\f")
            if not sequence:
                continue
            invalid = next((byte for byte in sequence if byte not in SEQUENCE_BYTES), None)
            if invalid is not None:
                raise ValidationError(
                    f"invalid sequence byte 0x{invalid:02x} in {label}: {source}:{line_number}"
                )
            current_bases += len(sequence)
            bases += len(sequence)
    if not saw_header or records == 0 or bases == 0 or current_bases == 0:
        raise ValidationError(f"{label} is not a non-empty FASTA: {source}")
    observed = frozenset(members)
    if member_delimiter is not None and not observed:
        raise ValidationError(f"no SAG membership could be recovered from {label}: {source}")
    if expected_members is not None and observed != expected_members:
        missing = sorted(expected_members - observed)[:5]
        extra = sorted(observed - expected_members)[:5]
        raise ValidationError(
            f"SAG membership mismatch in {label}: missing={missing!r}, extra={extra!r}"
        )
    return FastaEvidence(
        path=source,
        file_bytes=source.stat().st_size,
        sha256=digest.hexdigest(),
        records=records,
        bases=bases,
        members=observed,
    )


def load_stage3a_labels(groups_path: Path) -> tuple[Path | None, dict[str, frozenset[str]]]:
    candidate = groups_path.parent.parent / "02_dna2bit" / "labels.tsv"
    if not candidate.is_file():
        return None, {}
    source = candidate.resolve(strict=True)
    grouped: dict[str, set[str]] = defaultdict(set)
    assigned: set[str] = set()
    for row_number, row in read_strict_tsv(source, STAGE3A_LABEL_SCHEMA):
        sag = validate_identifier(row[0], "Stage3A SAG id", source, row_number)
        group = validate_identifier(row[3], "Stage3A species group", source, row_number)
        if sag in assigned:
            raise ValidationError(f"Stage3A SAG is assigned more than once: {sag}")
        assigned.add(sag)
        grouped[group].add(sag)
    return source, {group: frozenset(members) for group, members in grouped.items()}


def load_stage3b_membership(clusters_path: Path) -> tuple[Path | None, dict[str, frozenset[str]]]:
    current = clusters_path.parent.parent / "05_signed_leiden" / "chosen_membership.tsv"
    legacy = clusters_path.parent.parent / "06_signed_leiden" / "chosen_membership.tsv"
    candidates = [path for path in (current, legacy) if path.is_file()]
    if len(candidates) > 1:
        raise ValidationError(
            "ambiguous Stage3B membership sidecars: both current 05_signed_leiden and "
            "legacy 06_signed_leiden exist"
        )
    if not candidates:
        return None, {}
    source = candidates[0].resolve(strict=True)
    grouped: dict[str, set[str]] = defaultdict(set)
    declared: dict[str, int] = {}
    assigned: set[str] = set()
    for row_number, row in read_strict_tsv(source, STAGE3B_MEMBERSHIP_SCHEMA):
        cluster = validate_identifier(row[0], "Stage3B cluster", source, row_number)
        sag = validate_identifier(row[1], "Stage3B SAG id", source, row_number)
        cluster_size = parse_positive_integer(row[2], "cluster_size", source, row_number)
        if sag in assigned:
            raise ValidationError(f"Stage3B SAG is assigned more than once: {sag}")
        previous = declared.setdefault(cluster, cluster_size)
        if previous != cluster_size:
            raise ValidationError(f"inconsistent Stage3B cluster_size for {cluster}")
        assigned.add(sag)
        grouped[cluster].add(sag)
    for cluster, members in grouped.items():
        if len(members) != declared[cluster]:
            raise ValidationError(
                f"Stage3B membership count mismatch for {cluster}: "
                f"declared={declared[cluster]}, observed={len(members)}"
            )
    return source, {cluster: frozenset(members) for cluster, members in grouped.items()}


def load_current_stage3b_entry_closure(
    clusters_path: Path,
    membership_path: Path | None,
    memberships: dict[str, frozenset[str]],
) -> dict[str, object]:
    root = clusters_path.parent.parent
    is_current_layout = clusters_path.parent.name == "06_subassemble"
    if not is_current_layout:
        return {"mode": "legacy_layout_not_revalidated"}
    if membership_path is None or membership_path.parent.name != "05_signed_leiden":
        raise ValidationError("current Stage3B layout requires 05_signed_leiden/chosen_membership.tsv")

    node_path = require_source_file(
        root / "01_cellbit_negative_quality_pass" / "cellbit_negative_quality_pass.tsv",
        "current Stage3B graph-node table",
    )
    unaggregated_path = require_source_file(
        root / "05_signed_leiden" / "unaggregated_sags.tsv",
        "current Stage3B unaggregated table",
    )
    quality_audit_path = require_source_file(root / "quality_audit.tsv", "Stage3B quality audit")
    triangle_stats_path = require_source_file(
        root / "02_pairwise_ani_af" / "triangle_stats.json",
        "Stage3B exact-triangle stats",
    )
    complete_path = require_source_file(root / "COMPLETE.json", "Stage3B COMPLETE receipt")

    node_paths: dict[str, Path] = {}
    for row_number, row in read_strict_tsv(node_path, STAGE3B_NODE_SCHEMA):
        sag = validate_identifier(row[0], "Stage3B graph-node SAG id", node_path, row_number)
        assembly = resolve_row_path(node_path, row[1], "Stage3B graph-node assembly", row_number)
        if sag in node_paths:
            raise ValidationError(f"duplicate Stage3B graph node: {sag}")
        node_paths[sag] = assembly
    graph_nodes = frozenset(node_paths)

    assigned = frozenset(sag for members in memberships.values() for sag in members)
    unaggregated: set[str] = set()
    for row_number, row in read_strict_tsv(
        unaggregated_path, STAGE3B_UNAGGREGATED_SCHEMA
    ):
        sag = validate_identifier(row[0], "unaggregated Stage3B SAG id", unaggregated_path, row_number)
        assembly = resolve_row_path(
            unaggregated_path, row[1], "unaggregated Stage3B assembly", row_number
        )
        if sag in unaggregated:
            raise ValidationError(f"duplicate unaggregated Stage3B SAG: {sag}")
        if sag not in node_paths or node_paths[sag] != assembly:
            raise ValidationError(f"unaggregated SAG/path is absent from graph nodes: {sag}")
        unaggregated.add(sag)
    unaggregated_set = frozenset(unaggregated)
    if assigned & unaggregated_set:
        raise ValidationError("Stage3B assigned and unaggregated SAG sets overlap")
    if assigned | unaggregated_set != graph_nodes:
        missing = sorted(graph_nodes - assigned - unaggregated_set)[:5]
        extra = sorted((assigned | unaggregated_set) - graph_nodes)[:5]
        raise ValidationError(
            f"Stage3B graph-node closure failed: missing={missing!r}, extra={extra!r}"
        )

    audited: set[str] = set()
    audited_pass: set[str] = set()
    for row_number, row in read_strict_tsv(
        quality_audit_path, STAGE3B_QUALITY_AUDIT_SCHEMA
    ):
        sag = validate_identifier(row[0], "Stage3B quality-audit SAG id", quality_audit_path, row_number)
        assembly = resolve_row_path(
            quality_audit_path, row[1], "Stage3B quality-audit assembly", row_number
        )
        if sag in audited:
            raise ValidationError(f"duplicate Stage3B quality-audit SAG: {sag}")
        audited.add(sag)
        try:
            max_contig = int(row[3])
            contamination = float(row[5])
        except ValueError as error:
            raise ValidationError(
                f"invalid Stage3B quality value at {quality_audit_path}:{row_number}"
            ) from error
        expected_pass = max_contig >= 1000 and contamination < 5.0
        if row[6] not in {"0", "1"} or (row[6] == "1") != expected_pass:
            raise ValidationError(
                f"Stage3B quality decision disagrees with hard gates at {quality_audit_path}:{row_number}"
            )
        if expected_pass:
            if sag not in node_paths or node_paths[sag] != assembly:
                raise ValidationError(f"quality-passing SAG/path is absent from graph nodes: {sag}")
            audited_pass.add(sag)
    if frozenset(audited_pass) != graph_nodes:
        raise ValidationError("Stage3B graph nodes do not equal the quality-pass audit set")

    try:
        complete = json.loads(complete_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValidationError(f"invalid Stage3B COMPLETE receipt: {complete_path}") from error
    expected_fields = {
        "schema": "sag-stage3b-complete-v2",
        "status": "PASS",
        "reference_search": "disabled",
        "reference_positive_exclusion": False,
        "triangle_mode": "exact_all_pairs",
        "candidate_or_sketch_prefilter": False,
        "cellbit_negative_input": len(audited),
        "quality_failed": len(audited) - len(graph_nodes),
        "quality_pass_graph_nodes": len(graph_nodes),
        "pairwise_pairs_requested": len(graph_nodes) * (len(graph_nodes) - 1) // 2,
        "pairwise_pairs_evaluated": len(graph_nodes) * (len(graph_nodes) - 1) // 2,
        "triangle_stats_sha256": sha256_file(triangle_stats_path),
        "clusters_ge10": len(memberships),
        "assigned_sags": len(assigned),
        "unaggregated_sags": len(unaggregated_set),
        "assigned_plus_unaggregated_equals_quality_pass": True,
    }
    mismatches = {
        key: (complete.get(key), expected)
        for key, expected in expected_fields.items()
        if complete.get(key) != expected
    }
    if mismatches:
        raise ValidationError(f"Stage3B COMPLETE contract mismatch: {mismatches!r}")
    return {
        "mode": "current_v2_full_entry_closure",
        "graph_nodes": len(graph_nodes),
        "assigned_sags": len(assigned),
        "unaggregated_sags": len(unaggregated_set),
        "node_table": node_path,
        "unaggregated_table": unaggregated_path,
        "quality_audit": quality_audit_path,
        "triangle_stats": triangle_stats_path,
        "complete": complete_path,
    }


def inspect_bin_spec(spec: BinInspectionSpec) -> tuple[FastaEvidence, FastaEvidence]:
    """Hash and validate one bin's merged input and assembled FASTA."""

    input_evidence = inspect_fasta(
        spec.input_path,
        f"{spec.module} merged input",
        member_delimiter=spec.delimiter,
        expected_members=spec.sidecar_members,
    )
    bin_evidence = inspect_fasta(spec.bin_path, f"{spec.module} bin")
    return input_evidence, bin_evidence


def load_bins(
    stage3a_groups: Path,
    stage3b_clusters: Path,
    *,
    workers: int = 1,
) -> tuple[list[BinRecord], dict[str, object]]:
    stage3a = require_source_file(stage3a_groups, "Stage3A groups manifest")
    stage3b = require_source_file(stage3b_clusters, "Stage3B clusters manifest")
    labels_path, labels = load_stage3a_labels(stage3a)
    membership_path, memberships = load_stage3b_membership(stage3b)
    stage3b_entry_closure = load_current_stage3b_entry_closure(
        stage3b, membership_path, memberships
    )

    specs: list[BinInspectionSpec] = []
    bin_ids: set[str] = set()
    source_ids: dict[str, set[str]] = {"stage3a": set(), "stage3b": set()}
    bin_paths: set[Path] = set()
    input_paths: set[Path] = set()
    all_members: dict[str, str] = {}

    def add_record(
        *,
        module: str,
        source_id: str,
        declared_sags: int,
        input_path: Path,
        bin_path: Path,
        delimiter: str,
        sidecar_members: frozenset[str] | None,
        source_table: Path,
    ) -> None:
        bin_id = f"{module}__{source_id}"
        if bin_id in bin_ids:
            raise ValidationError(f"duplicate evaluation bin id: {bin_id}")
        if source_id in source_ids[module]:
            raise ValidationError(f"duplicate {module} source id: {source_id}")
        source_ids[module].add(source_id)
        if input_path in input_paths:
            raise ValidationError(f"merged input FASTA is reused by multiple bins: {input_path}")
        if bin_path in bin_paths:
            raise ValidationError(f"bin FASTA is reused by multiple bins: {bin_path}")
        input_paths.add(input_path)
        bin_paths.add(bin_path)

        bin_ids.add(bin_id)
        specs.append(
            BinInspectionSpec(
                bin_id=bin_id,
                module=module,
                source_id=source_id,
                declared_sags=declared_sags,
                input_path=input_path,
                bin_path=bin_path,
                delimiter=delimiter,
                sidecar_members=sidecar_members,
                source_table=source_table,
            )
        )

    stage3a_rows = read_strict_tsv(stage3a, STAGE3A_SCHEMA)
    for row_number, row in stage3a_rows:
        group = validate_identifier(row[0], "Stage3A species group", stage3a, row_number)
        count = parse_positive_integer(row[1], "sag_count", stage3a, row_number)
        input_path = resolve_row_path(stage3a, row[2], "Stage3A input_fasta", row_number)
        bin_path = resolve_row_path(stage3a, row[3], "Stage3A bin_fasta", row_number)
        sidecar = labels.get(group) if labels_path is not None else None
        if labels_path is not None and sidecar is None:
            raise ValidationError(f"Stage3A group is absent from labels.tsv: {group}")
        add_record(
            module="stage3a",
            source_id=group,
            declared_sags=count,
            input_path=input_path,
            bin_path=bin_path,
            delimiter="|",
            sidecar_members=sidecar,
            source_table=stage3a,
        )
    if labels_path is not None and set(labels) != source_ids["stage3a"]:
        missing = sorted(set(labels) - source_ids["stage3a"])[:5]
        raise ValidationError(f"labels.tsv contains groups absent from groups.tsv: {missing!r}")

    stage3b_rows = read_strict_tsv(stage3b, STAGE3B_SCHEMA)
    for row_number, row in stage3b_rows:
        cluster = validate_identifier(row[0], "Stage3B cluster", stage3b, row_number)
        count = parse_positive_integer(row[1], "n_sags", stage3b, row_number)
        found = parse_positive_integer(row[2], "n_found", stage3b, row_number)
        if found != count:
            raise ValidationError(
                f"Stage3B n_found/n_sags mismatch at {stage3b}:{row_number}: {found} != {count}"
            )
        input_path = resolve_row_path(stage3b, row[3], "Stage3B input_fasta", row_number)
        bin_path = resolve_row_path(stage3b, row[4], "Stage3B bin_fasta", row_number)
        sidecar = memberships.get(cluster) if membership_path is not None else None
        if membership_path is not None and sidecar is None:
            raise ValidationError(f"Stage3B cluster is absent from chosen_membership.tsv: {cluster}")
        add_record(
            module="stage3b",
            source_id=cluster,
            declared_sags=count,
            input_path=input_path,
            bin_path=bin_path,
            delimiter="__",
            sidecar_members=sidecar,
            source_table=stage3b,
        )
    if membership_path is not None and set(memberships) != source_ids["stage3b"]:
        missing = sorted(set(memberships) - source_ids["stage3b"])[:5]
        raise ValidationError(
            f"chosen_membership.tsv contains clusters absent from clusters.tsv: {missing!r}"
        )
    if not specs:
        raise ValidationError("the two manifests contain no non-empty bins")

    # FASTA validation and hashing dominate this command and every bin is
    # independent.  executor.map preserves manifest order, so the subsequent
    # closure checks and emitted rows are byte-for-byte deterministic.  The
    # default remains one worker for backward-compatible resource behavior.
    effective_workers = min(workers, len(specs))
    if effective_workers == 1:
        inspections = [inspect_bin_spec(spec) for spec in specs]
    else:
        # FASTA alphabet validation is deliberately strict and CPU-heavy.  A
        # process pool gives each independent bin a Python interpreter of its
        # own, avoiding the GIL while retaining the exact same validator.
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            inspections = list(executor.map(inspect_bin_spec, specs))

    records: list[BinRecord] = []
    for spec, (input_evidence, bin_evidence) in zip(specs, inspections, strict=True):
        members = (
            spec.sidecar_members
            if spec.sidecar_members is not None
            else input_evidence.members
        )
        if len(members) != spec.declared_sags:
            raise ValidationError(
                f"{spec.module} membership closure failed for {spec.source_id}: "
                f"declared={spec.declared_sags}, verified={len(members)}"
            )
        for sag in members:
            previous = all_members.setdefault(sag, spec.bin_id)
            if previous != spec.bin_id:
                raise ValidationError(
                    f"SAG {sag} occurs in both {previous} and {spec.bin_id}"
                )
        records.append(
            BinRecord(
                bin_id=spec.bin_id,
                module=spec.module,
                source_id=spec.source_id,
                declared_sags=spec.declared_sags,
                members=members,
                membership_source="producer_sidecar+merged_fasta_headers"
                if spec.sidecar_members is not None
                else "merged_fasta_headers_only",
                input_fasta=input_evidence,
                bin_fasta=bin_evidence,
                source_table=spec.source_table,
            )
        )
    records.sort(key=lambda record: record.bin_id.encode("utf-8"))
    provenance: dict[str, object] = {
        "stage3a_groups": stage3a,
        "stage3b_clusters": stage3b,
        "stage3a_labels": labels_path,
        "stage3b_membership": membership_path,
        "stage3a_membership_mode": "sidecar+headers" if labels_path else "headers_only",
        "stage3b_membership_mode": "sidecar+headers" if membership_path else "headers_only",
        "stage3b_entry_closure": stage3b_entry_closure,
        "unique_sag_members": len(all_members),
    }
    return records, provenance


def write_tsv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def source_receipt(path: Path | None) -> dict[str, object] | None:
    if path is None:
        return None
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage3a-groups", type=Path, required=True)
    parser.add_argument("--stage3b-clusters", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="parallel FASTA inspection workers (default: 1)",
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.workers < 1:
        raise ValidationError(f"--workers must be a positive integer: {args.workers}")
    stage3a_groups = args.stage3a_groups.resolve(strict=True)
    stage3b_clusters = args.stage3b_clusters.resolve(strict=True)
    out_dir = args.out_dir.absolute()
    if out_dir.exists() or out_dir.is_symlink():
        raise ValidationError(f"write-once output already exists: {out_dir}")
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = out_dir.parent / f".{out_dir.name}.incomplete-{os.getpid()}"
    if staging.exists() or staging.is_symlink():
        raise ValidationError(f"staging path already exists: {staging}")

    records, provenance = load_bins(
        stage3a_groups,
        stage3b_clusters,
        workers=args.workers,
    )
    staging.mkdir()
    input_dir = staging / "input_fna"
    input_dir.mkdir()
    manifest_path = staging / "INPUT_MANIFEST.tsv"
    manifest_rows: list[list[object]] = []
    for record in records:
        link_name = f"{record.bin_id}.fna"
        staged_link = input_dir / link_name
        final_link = out_dir / "input_fna" / link_name
        os.symlink(str(record.bin_fasta.path), staged_link)
        if not staged_link.is_symlink():
            raise ValidationError(f"failed to create input symlink: {staged_link}")
        if staged_link.resolve(strict=True) != record.bin_fasta.path:
            raise ValidationError(f"input symlink resolves to the wrong bin: {staged_link}")
        if staged_link.stat().st_size != record.bin_fasta.file_bytes:
            raise ValidationError(f"input symlink size differs from source bin: {staged_link}")
        manifest_rows.append(
            [
                record.bin_id,
                record.module,
                record.source_id,
                record.declared_sags,
                len(record.members),
                record.membership_source,
                sha256_members(record.members),
                str(record.source_table),
                str(record.input_fasta.path),
                record.input_fasta.file_bytes,
                record.input_fasta.sha256,
                record.input_fasta.records,
                record.input_fasta.bases,
                str(record.bin_fasta.path),
                record.bin_fasta.file_bytes,
                record.bin_fasta.sha256,
                record.bin_fasta.records,
                record.bin_fasta.bases,
                str(final_link),
                str(record.bin_fasta.path),
            ]
        )
    write_tsv(
        manifest_path,
        [
            "bin_id",
            "module",
            "source_id",
            "declared_sag_count",
            "verified_sag_count",
            "membership_evidence",
            "member_set_sha256",
            "source_table",
            "source_input_fasta",
            "source_input_bytes",
            "source_input_sha256",
            "source_input_records",
            "source_input_bases",
            "source_bin_fasta",
            "source_bin_bytes",
            "source_bin_sha256",
            "source_bin_records",
            "source_bin_bases",
            "input_link",
            "symlink_target",
        ],
        manifest_rows,
    )

    stage3a_count = sum(record.module == "stage3a" for record in records)
    stage3b_count = sum(record.module == "stage3b" for record in records)
    complete: dict[str, object] = {
        "schema": "current-baseline-bin-evaluation-input-v1",
        "status": "PASS",
        "purpose": "write-once CheckM2 input view; CheckM2 was not run by this command",
        "input_directory": str(out_dir / "input_fna"),
        "input_extension": "fna",
        "bin_count": len(records),
        "stage3a_bin_count": stage3a_count,
        "stage3b_bin_count": stage3b_count,
        "nonempty_bin_count": len(records),
        "unique_bin_id_count": len({record.bin_id for record in records}),
        "symlink_count": len(records),
        "declared_sag_memberships": sum(record.declared_sags for record in records),
        "verified_unique_sag_members": provenance["unique_sag_members"],
        "cross_bin_sag_overlap_count": 0,
        "stage3a_membership_mode": provenance["stage3a_membership_mode"],
        "stage3b_membership_mode": provenance["stage3b_membership_mode"],
        "stage3b_entry_closure_mode": provenance["stage3b_entry_closure"]["mode"],
        "stage3b_graph_nodes": provenance["stage3b_entry_closure"].get("graph_nodes"),
        "stage3b_assigned_sags": provenance["stage3b_entry_closure"].get("assigned_sags"),
        "stage3b_unaggregated_sags": provenance["stage3b_entry_closure"].get(
            "unaggregated_sags"
        ),
        "sources": {
            "stage3a_groups": source_receipt(provenance["stage3a_groups"]),
            "stage3b_clusters": source_receipt(provenance["stage3b_clusters"]),
            "stage3a_labels": source_receipt(provenance["stage3a_labels"]),
            "stage3b_membership": source_receipt(provenance["stage3b_membership"]),
            "stage3b_graph_nodes": source_receipt(
                provenance["stage3b_entry_closure"].get("node_table")
            ),
            "stage3b_unaggregated": source_receipt(
                provenance["stage3b_entry_closure"].get("unaggregated_table")
            ),
            "stage3b_quality_audit": source_receipt(
                provenance["stage3b_entry_closure"].get("quality_audit")
            ),
            "stage3b_triangle_stats": source_receipt(
                provenance["stage3b_entry_closure"].get("triangle_stats")
            ),
            "stage3b_complete": source_receipt(
                provenance["stage3b_entry_closure"].get("complete")
            ),
        },
        "input_manifest": {
            "path": str(out_dir / "INPUT_MANIFEST.tsv"),
            "bytes": manifest_path.stat().st_size,
            "sha256": sha256_file(manifest_path),
        },
        "checks": {
            "exact_manifest_schemas": True,
            "unique_bin_ids": True,
            "unique_resolved_bin_paths": True,
            "all_bins_nonempty_valid_fasta": True,
            "declared_membership_counts_closed": True,
            "merged_fasta_membership_closed": True,
            "available_producer_sidecars_closed": True,
            "current_stage3b_graph_nodes_equal_assigned_union_unaggregated":
                provenance["stage3b_entry_closure"]["mode"]
                == "current_v2_full_entry_closure",
            "sags_disjoint_across_output_bins": True,
            "symlinks_resolve_to_hashed_source_bins": True,
        },
    }
    complete_path = staging / "COMPLETE.json"
    complete_path.write_text(
        json.dumps(complete, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )

    checksum_sources: list[tuple[str, Path]] = [
        ("stage3a_groups", provenance["stage3a_groups"]),
        ("stage3b_clusters", provenance["stage3b_clusters"]),
    ]
    if provenance["stage3a_labels"] is not None:
        checksum_sources.append(("stage3a_labels", provenance["stage3a_labels"]))
    if provenance["stage3b_membership"] is not None:
        checksum_sources.append(("stage3b_membership", provenance["stage3b_membership"]))
    for role, key in (
        ("stage3b_graph_nodes", "node_table"),
        ("stage3b_unaggregated", "unaggregated_table"),
        ("stage3b_quality_audit", "quality_audit"),
        ("stage3b_triangle_stats", "triangle_stats"),
        ("stage3b_complete", "complete"),
    ):
        path = provenance["stage3b_entry_closure"].get(key)
        if path is not None:
            checksum_sources.append((role, path))
    checksum_rows: list[list[object]] = [
        [sha256_file(path), path.stat().st_size, role, str(path)]
        for role, path in checksum_sources
    ]
    checksum_rows.extend(
        [
            [
                sha256_file(manifest_path),
                manifest_path.stat().st_size,
                "generated_input_manifest",
                str(out_dir / "INPUT_MANIFEST.tsv"),
            ],
            [
                sha256_file(complete_path),
                complete_path.stat().st_size,
                "generated_complete_receipt",
                str(out_dir / "COMPLETE.json"),
            ],
        ]
    )
    write_tsv(staging / "SHA256SUMS.tsv", ["sha256", "bytes", "role", "path"], checksum_rows)

    staging.rename(out_dir)
    result = {
        "status": "PASS",
        "out_dir": str(out_dir),
        "bins": len(records),
        "stage3a_bins": stage3a_count,
        "stage3b_bins": stage3b_count,
        "unique_sags": provenance["unique_sag_members"],
        "input_manifest_sha256": complete["input_manifest"]["sha256"],
        "complete_sha256": sha256_file(out_dir / "COMPLETE.json"),
    }
    return result


def main() -> int:
    args = parse_args()
    try:
        result = run(args)
    except (ValidationError, OSError, csv.Error, ValueError) as error:
        raise SystemExit(f"FAIL: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
