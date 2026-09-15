#!/usr/bin/env python3
"""Read-only release audit for the Lake speed-v1 full E2E run.

The auditor never creates, deletes, or modifies a scientific artifact.  It
prints one JSON document to stdout and exits non-zero when any hard gate fails.
Run it only after the E2E runner has finished; ``--deep`` additionally hashes
all paths declared by Stage3B receipts and all merged subassembly inputs.

The old 220-thread Lake run is used as a *scientific authority*, not as a byte
authority for Flye assemblies: Flye's historical multi-thread assemble phase
was proven non-deterministic.  Therefore the audit requires equivalence of the
inputs, labels, quality rows, ANI/AF rows, graph edges, and graph partitions,
while final Flye products are validated against their own write-once manifests
and CheckM2 ID closure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Callable, Iterable


DEFAULT_RUN = Path(
    "/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/"
    "24_lake_full_e2e_220t_speed_v1_20260905/attempt_001"
)
DEFAULT_BASELINE = Path(
    "/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/"
    "11_lake_full_current_tractor_baseline_v1"
)

EXPECTED = {
    "input_sags": 13742,
    "labeled": 7527,
    "pending": 6215,
    "stage3a_groups": 274,
    "quality_nodes": 5834,
    "quality_failed": 381,
    "pairs": 17014861,
    "triangle_rows": 322110,
    "positive_edges": 103555,
    "negative_edges": 48205,
    "stage3b_groups": 62,
    "assigned": 2567,
    "unaggregated": 3267,
    "final_bins": 336,
}

TIME_FILES = [
    "01_STAGE1_3A_TIME.txt",
    "02_INPUT_VIEW_TIME.txt",
    "03_PENDING_STATS_TIME.txt",
    "04_CHECKM2_TIME.txt",
    "05_GTDBTK_TIME.txt",
    "06_FINALIZE_INPUTS_TIME.txt",
    "07_STAGE3B_TIME.txt",
    "08_FINAL_VIEW_TIME.txt",
    "09_FINAL_CHECKM2_TIME.txt",
]


class AuditError(RuntimeError):
    pass


def sha256(path: Path, cache: dict[Path, str] | None = None) -> str:
    path = path.resolve(strict=True)
    if cache is not None and path in cache:
        return cache[path]
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    value = digest.hexdigest()
    if cache is not None:
        cache[path] = value
    return value


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise AuditError(f"missing JSON: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise AuditError(f"invalid JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise AuditError(f"expected JSON object: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def require_fields(value: dict[str, Any], expected: dict[str, Any], label: str) -> None:
    for key, wanted in expected.items():
        observed = value.get(key)
        if observed != wanted:
            raise AuditError(f"{label}: {key}={observed!r}, expected {wanted!r}")


def tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise AuditError(f"missing TSV: {path}")

    # Pipeline TSV records are LF-delimited.  A frozen dna2bit labels file has
    # a bare CR at the end of each taxonomy field (immediately before TAB).
    # ``csv`` otherwise treats that CR as a second record boundary and creates
    # one spurious empty-sag row per real label.  Normalize only CRLF record
    # endings, protect all remaining CR bytes while csv handles quoting/TABs,
    # then restore those CRs as field data.  This keeps the LF record contract
    # without weakening normal CSV quoting semantics.
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            text = handle.read()
    except UnicodeDecodeError as exc:
        raise AuditError(f"TSV is not UTF-8: {path}: {exc}") from exc
    text = text.replace("\r\n", "\n")
    sentinel = "\ue000"
    while sentinel in text:
        codepoint = ord(sentinel) + 1
        if codepoint > 0xF8FF:
            raise AuditError(f"cannot reserve bare-CR sentinel for TSV: {path}")
        sentinel = chr(codepoint)
    protected = text.replace("\r", sentinel)
    try:
        reader = csv.DictReader(io.StringIO(protected, newline=""), delimiter="\t")
        if not reader.fieldnames:
            raise AuditError(f"TSV has no header: {path}")
        rows = list(reader)
    except csv.Error as exc:
        raise AuditError(f"invalid TSV {path}: {exc}") from exc
    if any(None in row for row in rows):
        raise AuditError(f"TSV row has excess columns: {path}")

    def restore(value: str | None) -> str | None:
        return value.replace(sentinel, "\r") if value is not None else None

    header = [restore(name) for name in reader.fieldnames]
    restored_rows = [
        {restore(key): restore(value) for key, value in row.items()}
        for row in rows
    ]
    # Headers are guaranteed non-null by DictReader; keep the runtime check
    # explicit so the return type and downstream dictionary access stay strict.
    if any(name is None for name in header):
        raise AuditError(f"TSV has a null header field: {path}")
    return [str(name) for name in header], [
        {str(key): "" if value is None else value for key, value in row.items()}
        for row in restored_rows
    ]


def decimal_text(raw: str, label: str) -> str:
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise AuditError(f"invalid decimal {label}: {raw!r}") from exc
    if not value.is_finite():
        raise AuditError(f"non-finite decimal {label}: {raw!r}")
    if value == 0:
        return "0"
    return format(value.normalize(), "f")


def digest_records(records: Iterable[tuple[str, ...]]) -> tuple[int, str, int]:
    material = list(records)
    duplicates = len(material) - len(set(material))
    material.sort()
    digest = hashlib.sha256()
    for row in material:
        for item in row:
            data = item.encode("utf-8")
            digest.update(len(data).to_bytes(8, "little"))
            digest.update(data)
    return len(material), digest.hexdigest(), duplicates


def normalized_whole_tsv(path: Path, numeric_columns: set[str] | None = None) -> tuple[int, str, int]:
    header, rows = tsv_rows(path)
    numeric_columns = numeric_columns or set()
    return digest_records(
        tuple(decimal_text(row[name], f"{path}:{name}") if name in numeric_columns else row[name]
              for name in header)
        for row in rows
    )


def normalized_stage1_path_tsv(
    path: Path, stage1_root: Path
) -> tuple[list[str], tuple[int, str, int]]:
    """Normalize the run-root-specific assembly path in a Stage1 table.

    Both the speed run and its authority baseline legitimately bind rows to
    their own write-once ``01_assembly`` trees.  Comparing those absolute path
    strings byte-for-byte is therefore invalid.  We first require every path
    to have the exact canonical location for its SAG in the supplied root,
    then compare all columns after replacing only that root-specific prefix.
    """
    header, rows = tsv_rows(path)
    require("sag_id" in header and "assembly_fasta" in header,
            f"Stage1 authority TSV lacks path-binding columns: {path}")
    records: list[tuple[str, ...]] = []
    for row in rows:
        sag = row.get("sag_id", "")
        require(sag != "", f"Stage1 authority TSV has empty sag_id: {path}")
        observed = Path(row.get("assembly_fasta", ""))
        expected = stage1_root / "01_assembly" / sag / f"{sag}.fasta"
        require(observed == expected,
                f"Stage1 assembly path is not canonically bound: {path}:{sag}: {observed}")
        require(observed.is_file() and observed.stat().st_size > 0,
                f"Stage1 assembly path is missing/empty: {path}:{sag}: {observed}")
        normalized = dict(row)
        normalized["assembly_fasta"] = f"01_assembly/{sag}/{sag}.fasta"
        records.append(tuple(normalized[name] for name in header))
    return header, digest_records(records)


def node_table(path: Path) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    _, rows = tsv_rows(path)
    by_path: dict[str, str] = {}
    nodes: dict[str, tuple[str, str]] = {}
    for row in rows:
        sag = row.get("sag_id", "")
        assembly = row.get("assembly_fasta", "")
        require(sag and assembly and sag not in nodes, f"bad/duplicate node in {path}: {sag!r}")
        require(assembly not in by_path, f"duplicate assembly path in {path}: {assembly}")
        gc = decimal_text(row.get("gc_pct", ""), f"{path}:{sag}:gc_pct")
        nodes[sag] = (assembly, gc)
        by_path[assembly] = sag
    return by_path, nodes


def normalized_triangle(path: Path, path_to_sag: dict[str, str]) -> tuple[int, str, int]:
    _, rows = tsv_rows(path)
    out: list[tuple[str, ...]] = []
    for row in rows:
        try:
            left, right = path_to_sag[row["Ref_file"]], path_to_sag[row["Query_file"]]
        except KeyError as exc:
            raise AuditError(f"triangle path is not bound to a graph node: {exc}") from exc
        af_left = decimal_text(row["Align_fraction_ref"], "triangle AF_ref")
        af_right = decimal_text(row["Align_fraction_query"], "triangle AF_query")
        if right < left:
            left, right, af_left, af_right = right, left, af_right, af_left
        require(left != right, f"triangle self-pair: {left}")
        out.append((left, right, decimal_text(row["ANI"], "triangle ANI"), af_left, af_right))
    return digest_records(out)


def normalized_edge_table(path: Path, kind: str) -> tuple[int, str, int, set[str]]:
    _, rows = tsv_rows(path)
    out: list[tuple[str, ...]] = []
    endpoints: set[str] = set()
    for row in rows:
        left, right = row["SAG1"], row["SAG2"]
        require(left != right, f"{kind} self-edge: {left}")
        endpoints.update((left, right))
        if kind == "positive":
            af_left = decimal_text(row["AF_ref"], "positive AF_ref")
            af_right = decimal_text(row["AF_query"], "positive AF_query")
            if right < left:
                left, right, af_left, af_right = right, left, af_right, af_left
            out.append((
                left, right, decimal_text(row["ANI"], "positive ANI"), af_left, af_right,
                decimal_text(row["gc_diff"], "positive gc_diff"),
                decimal_text(row["weight"], "positive weight"),
            ))
        elif kind == "negative":
            if right < left:
                left, right = right, left
            out.append((
                left, right, str(int(row["n_markers"])),
                decimal_text(row["mean_pident"], "negative mean_pident"),
                decimal_text(row["weight"], "negative weight"),
            ))
        else:
            raise AuditError(f"unknown edge kind: {kind}")
    count, digest, duplicates = digest_records(out)
    return count, digest, duplicates, endpoints


def partition(path: Path) -> tuple[int, int, str, set[str]]:
    _, rows = tsv_rows(path)
    clusters: dict[str, set[str]] = defaultdict(set)
    declared: dict[str, set[int]] = defaultdict(set)
    seen: set[str] = set()
    for row in rows:
        cluster, sag = row["cluster"], row["SAG_id"]
        require(sag not in seen, f"duplicate partition SAG in {path}: {sag}")
        seen.add(sag)
        clusters[cluster].add(sag)
        declared[cluster].add(int(row["cluster_size"]))
    for cluster, members in clusters.items():
        require(declared[cluster] == {len(members)}, f"cluster_size mismatch {path}:{cluster}")
    canonical = [tuple(sorted(members)) for members in clusters.values()]
    count, digest, duplicates = digest_records(canonical)
    require(duplicates == 0, f"duplicate cluster member sets in {path}")
    return count, len(seen), digest, seen


def unaggregated(path: Path) -> set[str]:
    _, rows = tsv_rows(path)
    ids = [row["sag_id"] for row in rows]
    require(len(ids) == len(set(ids)), f"duplicate unaggregated SAG in {path}")
    return set(ids)


def validate_hash_manifest(path: Path, cache: dict[Path, str]) -> int:
    require(path.is_file(), f"missing SHA256 manifest: {path}")
    count = 0
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        match = re.match(r"^([0-9a-f]{64})\s+\*?(.*)$", raw)
        require(match is not None, f"invalid SHA256 manifest line {path}:{line_number}")
        expected, target_text = match.groups()
        target = Path(target_text)
        if not target.is_absolute():
            target = path.parent / target
        require(target.is_file(), f"SHA256 target missing: {target}")
        require(sha256(target, cache) == expected, f"SHA256 mismatch: {target}")
        count += 1
    require(count > 0, f"empty SHA256 manifest: {path}")
    return count


def validate_receipt(path: Path, cache: dict[Path, str], include_inputs: bool) -> int:
    value = load_json(path)
    require(value.get("status") == "PASS", f"non-PASS receipt: {path}")
    checked = 0
    roles = ("inputs", "outputs") if include_inputs else ("outputs",)
    for role in roles:
        entries = value.get(role, [])
        require(isinstance(entries, list), f"receipt {role} is not a list: {path}")
        for entry in entries:
            require(isinstance(entry, dict), f"bad receipt entry: {path}")
            target = Path(entry.get("path", ""))
            require(target.is_file(), f"receipt target missing: {path} -> {target}")
            require(target.stat().st_size == int(entry["bytes"]), f"receipt size mismatch: {target}")
            require(sha256(target, cache) == entry["sha256"], f"receipt SHA mismatch: {target}")
            checked += 1
    return checked


def parse_time_v(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing time file: {path}")
    fields: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if ": " in line:
            key, value = line.rsplit(": ", 1)
            fields[key] = value
    elapsed_text = fields.get("Elapsed (wall clock) time (h:mm:ss or m:ss)")
    require(elapsed_text is not None, f"elapsed wall time absent: {path}")
    pieces = elapsed_text.split(":")
    require(len(pieces) in (2, 3), f"invalid elapsed time {path}: {elapsed_text}")
    if len(pieces) == 2:
        wall = int(pieces[0]) * 60 + float(pieces[1])
    else:
        wall = int(pieces[0]) * 3600 + int(pieces[1]) * 60 + float(pieces[2])
    exit_status = int(fields.get("Exit status", "-1"))
    require(exit_status == 0, f"nonzero /usr/bin/time exit status in {path}: {exit_status}")
    return {
        "wall_seconds": wall,
        "user_seconds": float(fields.get("User time (seconds)", "nan")),
        "system_seconds": float(fields.get("System time (seconds)", "nan")),
        "max_rss_kib": int(fields.get("Maximum resident set size (kbytes)", "-1")),
        "fs_inputs": int(fields.get("File system inputs", "-1")),
        "fs_outputs": int(fields.get("File system outputs", "-1")),
    }


def validate_telemetry(path: Path) -> dict[str, Any]:
    header, rows = tsv_rows(path)
    required = {
        "elapsed_s", "phase", "load1", "mem_available_bytes", "swap_used_bytes",
        "cgroup_memory_current", "cgroup_memory_peak", "cgroup_cpu_usage_usec",
        "cgroup_io_read_bytes", "cgroup_io_write_bytes", "cpu_psi_some_avg10",
        "memory_psi_some_avg10", "memory_psi_full_avg10", "io_psi_some_avg10",
        "io_psi_full_avg10",
    }
    require(required.issubset(header), f"telemetry columns missing: {sorted(required-set(header))}")
    require(len(rows) >= 2, "telemetry has fewer than two samples")
    int_fields = [
        "elapsed_s", "mem_available_bytes", "swap_used_bytes", "cgroup_memory_current",
        "cgroup_memory_peak", "cgroup_cpu_usage_usec", "cgroup_io_read_bytes",
        "cgroup_io_write_bytes",
    ]
    float_fields = [
        "load1", "cpu_psi_some_avg10", "memory_psi_some_avg10",
        "memory_psi_full_avg10", "io_psi_some_avg10", "io_psi_full_avg10",
    ]
    for row in rows:
        for name in int_fields:
            require(int(row[name]) >= 0, f"negative telemetry {name}")
        for name in float_fields:
            value = float(row[name])
            require(math.isfinite(value) and value >= 0, f"invalid telemetry {name}")
    for name in ["elapsed_s", "cgroup_memory_peak", "cgroup_cpu_usage_usec",
                 "cgroup_io_read_bytes", "cgroup_io_write_bytes"]:
        values = [int(row[name]) for row in rows]
        require(values == sorted(values), f"non-monotonic telemetry counter: {name}")
    phases = [row["phase"] for row in rows]
    # Sampling is every 20 seconds, so a short atomic view/finalization phase
    # may legitimately begin and end between samples.  Require the long-lived
    # scientific phases here; the complete ordered phase list is gated below
    # using PHASE_EVENTS.tsv, which records every transition synchronously.
    expected_phases = {
        "stage1_3a", "stage3b_upstream_checkm2",
        "stage3b_upstream_gtdbtk_identify", "stage3b_exact_and_subassemble",
        "final_checkm2_evaluation", "complete",
    }
    require(expected_phases.issubset(set(phases)),
            f"telemetry lacks phases: {sorted(expected_phases-set(phases))}")
    swap = [int(row["swap_used_bytes"]) for row in rows]
    return {
        "samples": len(rows),
        "elapsed_last_s": int(rows[-1]["elapsed_s"]),
        "peak_load1": max(float(row["load1"]) for row in rows),
        "minimum_mem_available_bytes": min(int(row["mem_available_bytes"]) for row in rows),
        "peak_cgroup_memory_bytes": max(int(row["cgroup_memory_peak"]) for row in rows),
        "swap_growth_bytes": max(swap) - min(swap),
        "peak_cpu_psi_some_avg10": max(float(row["cpu_psi_some_avg10"]) for row in rows),
        "peak_memory_psi_some_avg10": max(float(row["memory_psi_some_avg10"]) for row in rows),
        "peak_memory_psi_full_avg10": max(float(row["memory_psi_full_avg10"]) for row in rows),
        "peak_io_psi_some_avg10": max(float(row["io_psi_some_avg10"]) for row in rows),
        "peak_io_psi_full_avg10": max(float(row["io_psi_full_avg10"]) for row in rows),
        "cgroup_cpu_seconds": int(rows[-1]["cgroup_cpu_usage_usec"]) / 1_000_000,
        "cgroup_read_bytes": int(rows[-1]["cgroup_io_read_bytes"]),
        "cgroup_write_bytes": int(rows[-1]["cgroup_io_write_bytes"]),
    }


def group_table(path: Path, group_field: str) -> dict[str, dict[str, str]]:
    _, rows = tsv_rows(path)
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        group = row[group_field]
        require(group and group not in out, f"bad/duplicate group in {path}: {group!r}")
        out[group] = row
    return out


def compare_file_sets_by_sha(
    current: dict[str, Path], baseline: dict[str, Path], cache: dict[Path, str], label: str
) -> None:
    require(current.keys() == baseline.keys(), f"{label} file-key set differs")
    for key in sorted(current):
        require(current[key].is_file() and baseline[key].is_file(), f"{label} missing file: {key}")
        require(sha256(current[key], cache) == sha256(baseline[key], cache),
                f"{label} content differs: {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--baseline-root", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument(
        "--deep", action="store_true",
        help="hash every Stage3B receipt input/output and all merged subassembly inputs",
    )
    args = parser.parse_args()

    result: dict[str, Any] = {
        "schema": "cellbit-lake-speed-release-readonly-audit-v1",
        "status": "FAIL",
        "run_root": str(args.run_root),
        "baseline_root": str(args.baseline_root),
        "mode": "deep" if args.deep else "standard",
        "checks": {},
        "metrics": {},
        "failures": [],
    }
    checks: dict[str, Any] = result["checks"]
    metrics: dict[str, Any] = result["metrics"]
    failures: list[str] = result["failures"]
    cache: dict[Path, str] = {}

    def section(name: str, action: Callable[[], None]) -> None:
        try:
            action()
            checks[name] = "PASS"
        except Exception as exc:
            checks[name] = "FAIL"
            failures.append(f"{name}: {type(exc).__name__}: {exc}")

    run = args.run_root
    base = args.baseline_root
    control = run / "00_control"
    main_out = run / "01_stage1_3a"
    input_view = run / "02_stage3b_input_view"
    pending_stats = run / "03_stage3b_pending_stats"
    upstream = run / "04_stage3b_upstream"
    inputs = run / "05_stage3b_inputs"
    stage3b = run / "06_stage3b_exact"
    final_view = run / "07_final_bin_view"
    final_eval = run / "08_final_checkm2"
    base_main = base / "01_stage1_3a"
    base_inputs = base / "03b_stage3b_inputs_attempt_001"
    base_stage3b = base / "04_stage3b_cellbit_negative_exact_attempt_003"
    base_view = base / "06_attempt003_final_bin_view_v1"

    def audit_run_envelope() -> None:
        run_result = load_json(control / "RUN_RESULT.json")
        require_fields(run_result, {
            "schema": "cellbit-lake-full-e2e-speed-v1", "status": "PASS", "exit_code": 0,
            "threads_requested": 220, "systemd_cpu_quota_cores": 220.0,
            "prior_scientific_artifact_reuse": False,
        }, "RUN_RESULT")
        require(run_result.get("core_pipeline_seconds", 0) > 0, "invalid core runtime")
        require(run_result.get("including_final_evaluation_seconds", 0) >= run_result["core_pipeline_seconds"],
                "total runtime is shorter than core runtime")
        cgroup = load_json(control / "CGROUP_CPU_QUOTA.json")
        require_fields(cgroup, {
            "actual_cpu_quota_cores": 220.0, "required_cpu_quota_cores": 220,
            "process_cpu_affinity": "0-239", "memory_high_bytes": "1288490188800",
            "memory_max_bytes": "1556925644800", "memory_swap_max_bytes": "8589934592",
            "pids_max": "max", "nofile_soft": "1048576", "nofile_hard": "1048576",
        }, "cgroup envelope")
        require(load_json(control / "CACHE_CONDITION.json").get("prior_scientific_artifact_reuse") is False,
                "cache receipt allows prior scientific artifact reuse")
        metrics["run"] = {
            "core_pipeline_seconds": run_result["core_pipeline_seconds"],
            "including_final_evaluation_seconds": run_result["including_final_evaluation_seconds"],
            "threads": run_result["threads_requested"],
        }
    section("run_envelope", audit_run_envelope)

    def audit_stage1() -> None:
        done = load_json(main_out / "COMPLETE.json")
        require_fields(done, {
            "status": "PASS", "dna_search_engine": "packed", "input_sags": 13742,
            "contig_inputs": 13742, "eligible_sags": 13742, "excluded_lt1000bp": 0,
            "labeled": 7527, "groups": 274,
        }, "Stage1-3A COMPLETE")
        _, manifest = tsv_rows(Path(
            "/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/"
            "10_current_tractor_baseline_continuation_v1_code/acceptance/"
            "lake_full_contigs_13742_v1/LAKE_FULL_CONTIGS.tsv"
        ))
        _, labels = tsv_rows(main_out / "02_dna2bit/labels.tsv")
        _, pending = tsv_rows(main_out / "03B_unclassified_pending.tsv")
        manifest_ids = [row["sag_id"] for row in manifest]
        label_ids = [row["sag_id"] for row in labels]
        pending_ids = [row["sag_id"] for row in pending]
        require(len(manifest_ids) == 13742 == len(set(manifest_ids)), "manifest ID closure failed")
        require(len(label_ids) == 7527 == len(set(label_ids)), "label ID closure failed")
        require(len(pending_ids) == 6215 == len(set(pending_ids)), "pending ID closure failed")
        require(not (set(label_ids) & set(pending_ids)), "labels and pending overlap")
        require(set(label_ids) | set(pending_ids) == set(manifest_ids), "labels+pending != manifest")
        groups = group_table(main_out / "03A_subassemble/groups.tsv", "species_group")
        require(len(groups) == 274, "Stage3A group count mismatch")
        observed_counts = Counter(row["species_group"] for row in labels)
        require({key: int(row["sag_count"]) for key, row in groups.items()} == observed_counts,
                "Stage3A group member counts differ from labels")
        for key, row in groups.items():
            for field in ("input_fasta", "bin_fasta"):
                target = Path(row[field])
                require(target.is_file() and target.stat().st_size > 0,
                        f"Stage3A missing/empty {field}: {key}")
            require((Path(row["bin_fasta"]).parent / "CPP_FULL_PIPELINE_PASS").is_file(),
                    f"Stage3A missing CPP_FULL_PIPELINE_PASS: {key}")
        require(sum(1 for _ in (main_out / "01_assembly").glob("*/STAGE1.PASS")) == 13742,
                "STAGE1.PASS marker count mismatch")
        require(sum(1 for _ in (main_out / "02_dna2bit/bits").glob("*/SKETCH.PASS")) == 13742,
                "SKETCH.PASS marker count mismatch")
        for relative in ["02_dna2bit/search_result.csv", "02_dna2bit/labels.tsv"]:
            require((main_out / relative).read_bytes() == (base_main / relative).read_bytes(),
                    f"Stage1 deterministic authority differs: {relative}")
        for relative in ["03B_unclassified_pending.tsv", "STAGE3B_FASTA_STATS.tsv"]:
            current_header, current_normalized = normalized_stage1_path_tsv(
                main_out / relative, main_out
            )
            baseline_header, baseline_normalized = normalized_stage1_path_tsv(
                base_main / relative, base_main
            )
            require(current_header == baseline_header,
                    f"Stage1 authority header differs: {relative}")
            require(current_normalized == baseline_normalized,
                    f"Stage1 normalized authority differs: {relative}")
        current_groups = {key: (int(row["sag_count"]), Path(row["input_fasta"]))
                          for key, row in groups.items()}
        baseline_groups_raw = group_table(base_main / "03A_subassemble/groups.tsv", "species_group")
        baseline_groups = {key: (int(row["sag_count"]), Path(row["input_fasta"]))
                           for key, row in baseline_groups_raw.items()}
        require({k: v[0] for k, v in current_groups.items()} ==
                {k: v[0] for k, v in baseline_groups.items()}, "Stage3A group table differs")
        if args.deep:
            compare_file_sets_by_sha(
                {k: v[1] for k, v in current_groups.items()},
                {k: v[1] for k, v in baseline_groups.items()}, cache,
                "Stage3A merged inputs",
            )
        metrics["stage1_3a"] = {"inputs": 13742, "labeled": 7527, "pending": 6215, "groups": 274}
    section("stage1_3a_closure_and_authority", audit_stage1)

    def audit_upstream() -> None:
        view = load_json(input_view / "COMPLETE.json")
        require_fields(view, {
            "schema": "stage3b-upstream-input-view-v1", "status": "PASS",
            "sag_count": 6215, "unique_sag_ids": 6215, "nonempty_resolved_targets": 6215,
        }, "input-view COMPLETE")
        stats = load_json(pending_stats / "COMPLETE.json")
        require_fields(stats, {
            "schema": "stage3b-pending-fasta-stats-view-v1", "status": "PASS",
            "ignored_non_pending_rows": 7527,
        }, "pending-stats COMPLETE")
        require(stats.get("output_stats", {}).get("rows") == 6215, "pending-stats rows != 6215")
        done = load_json(inputs / "COMPLETE.json")
        require_fields(done, {
            "schema": "stage3b-upstream-finalization-complete-v1", "status": "PASS",
            "sag_count": 6215, "quality_rows": 6215,
        }, "Stage3B inputs COMPLETE")
        audit = load_json(inputs / "AUDIT.json")
        require(audit.get("status") == "PASS" and
                audit.get("schema") == "stage3b-upstream-finalization-audit-v1",
                "Stage3B inputs AUDIT is not PASS v1")
        require(all(value is True for value in audit.get("checks", {}).values()),
                "Stage3B inputs AUDIT contains a false check")
        require(audit.get("counts", {}).get("sags") == 6215 and
                audit.get("counts", {}).get("selected_sag_marker_rows") == 33576,
                "Stage3B inputs AUDIT counts differ")
        bindings = {
            inputs / "AUDIT.json": done["audit_sha256"],
            main_out / "03B_unclassified_pending.tsv": done["pending_sha256"],
            pending_stats / "STAGE3B_FASTA_STATS.tsv": done["stats_sha256"],
            input_view / "COMPLETE.json": done["input_view_complete_sha256"],
            input_view / "INPUT_MANIFEST.tsv": done["input_view_manifest_sha256"],
            upstream / "checkm2_attempt_001/quality_report.tsv": done["checkm2_report_sha256"],
            inputs / "stage3b_quality.tsv": done["stage3b_quality_sha256"],
            inputs / "bac120_marker_nt_map.tsv": done["bac120_marker_nt_map_sha256"],
        }
        for target, wanted in bindings.items():
            require(sha256(target, cache) == wanted, f"upstream receipt SHA mismatch: {target}")
        audit_bindings = [
            audit["inputs"]["pending"], audit["inputs"]["stats"],
            audit["inputs"]["checkm2_report"],
            audit["inputs"]["input_view"]["complete"],
            audit["inputs"]["input_view"]["manifest"],
            audit["inputs"]["gtdbtk_identify"]["gtdbtk_json"],
            audit["inputs"]["gtdbtk_identify"]["failed_genomes"],
            audit["outputs"]["stage3b_quality"],
            audit["outputs"]["bac120_marker_nt_map"],
            audit["script"],
        ]
        for binding in audit_bindings:
            target = Path(binding["path"])
            require(target.is_file(), f"AUDIT binding target missing: {target}")
            if "file_bytes" in binding:
                require(target.stat().st_size == int(binding["file_bytes"]),
                        f"AUDIT binding size mismatch: {target}")
            require(sha256(target, cache) == binding["sha256"],
                    f"AUDIT binding SHA mismatch: {target}")
        require(done["gtdbtk_json_sha256"] ==
                audit["inputs"]["gtdbtk_identify"]["gtdbtk_json"]["sha256"],
                "COMPLETE/AUDIT GTDB-Tk JSON binding differs")
        require(done["failed_genomes_sha256"] ==
                audit["inputs"]["gtdbtk_identify"]["failed_genomes"]["sha256"],
                "COMPLETE/AUDIT failed-genomes binding differs")
        for relative in ["stage3b_quality.tsv", "bac120_marker_nt_map.tsv"]:
            require((inputs / relative).read_bytes() == (base_inputs / relative).read_bytes(),
                    f"upstream scientific authority differs: {relative}")
        current_quality = normalized_whole_tsv(
            upstream / "checkm2_attempt_001/quality_report.tsv",
            {"Completeness", "Contamination", "Coding_Density", "Contig_N50",
             "Average_Gene_Length", "Genome_Size", "GC_Content", "Total_Coding_Sequences"},
        )
        baseline_quality = normalized_whole_tsv(
            base / "03_stage3b_upstream_attempt_002/checkm2_attempt_001/quality_report.tsv",
            {"Completeness", "Contamination", "Coding_Density", "Contig_N50",
             "Average_Gene_Length", "Genome_Size", "GC_Content", "Total_Coding_Sequences"},
        )
        require(current_quality == baseline_quality, "upstream CheckM2 rows differ from baseline")
        metrics["stage3b_upstream"] = {"sags": 6215, "marker_rows": 33576}
    section("stage3b_input_checkm2_gtdb_receipt_closure", audit_upstream)

    def audit_stage3b() -> None:
        done = load_json(stage3b / "COMPLETE.json")
        require_fields(done, {
            "schema": "sag-stage3b-complete-v2", "status": "PASS",
            "cellbit_negative_input": 6215, "quality_failed": 381,
            "quality_pass_graph_nodes": 5834, "pairwise_pairs_requested": 17014861,
            "pairwise_pairs_evaluated": 17014861, "triangle_rows_emitted": 322110,
            "positive_edges": 103555, "negative_edges": 48205, "clusters_ge10": 62,
            "assigned_sags": 2567, "unaggregated_sags": 3267,
            "assigned_plus_unaggregated_equals_quality_pass": True,
        }, "Stage3B COMPLETE")
        triangle_stats = load_json(stage3b / "02_pairwise_ani_af/triangle_stats.json")
        require_fields(triangle_stats, {
            "schema": "gtdb-ani-af-triangle-stats-v1", "status": "PASS",
            "triangle_mode": "exact", "exact_execution": "global-seed-join",
            "nodes": 5834, "pairs_expected": 17014861, "pairs_evaluated": 17014861,
            "rows_emitted": 322110,
        }, "triangle stats")
        require(5834 * 5833 // 2 == 17014861, "nC2 pair arithmetic failed")
        require(sha256(stage3b / "02_pairwise_ani_af/triangle_stats.json", cache) ==
                done["triangle_stats_sha256"], "triangle_stats SHA binding mismatch")

        current_path_map, current_nodes = node_table(
            stage3b / "01_cellbit_negative_quality_pass/cellbit_negative_quality_pass.tsv")
        baseline_path_map, baseline_nodes = node_table(
            base_stage3b / "01_cellbit_negative_quality_pass/cellbit_negative_quality_pass.tsv")
        require(current_nodes == baseline_nodes and len(current_nodes) == 5834,
                "Stage3B graph nodes/GC differ from baseline")
        tri = normalized_triangle(stage3b / "02_pairwise_ani_af/ani_triangle_sparse.tsv", current_path_map)
        base_tri = normalized_triangle(
            base_stage3b / "02_pairwise_ani_af/ani_triangle_sparse.tsv", baseline_path_map)
        require(tri == base_tri and tri[0] == 322110 and tri[2] == 0,
                "ANI/AF sparse triangle differs from baseline or contains duplicate rows")
        pos = normalized_edge_table(stage3b / "02_pairwise_ani_af/positive_edges.tsv", "positive")
        base_pos = normalized_edge_table(base_stage3b / "02_pairwise_ani_af/positive_edges.tsv", "positive")
        require(pos[:3] == base_pos[:3] and pos[0] == 103555 and pos[2] == 0,
                "positive edges differ from baseline or contain duplicates")
        neg = normalized_edge_table(stage3b / "04_marker_blastn/negative_edges.tsv", "negative")
        base_neg = normalized_edge_table(base_stage3b / "04_marker_blastn/negative_edges.tsv", "negative")
        require(neg[:3] == base_neg[:3] and neg[0] == 48205 and neg[2] == 0,
                "negative edges differ from baseline or contain duplicates")
        node_ids = set(current_nodes)
        require(pos[3] <= node_ids and neg[3] <= node_ids, "graph edge endpoint outside node set")

        for relative, numeric in [
            ("04_marker_blastn/marker_pairwise_blastn.tsv", {"n_markers", "mean_pident", "min_pident"}),
            ("quality_audit.tsv", {"total_bp", "max_contig", "gc_pct", "checkm2_contamination", "stage3b_quality_pass"}),
        ]:
            current_digest = normalized_whole_tsv(stage3b / relative, numeric)
            baseline_digest = normalized_whole_tsv(base_stage3b / relative, numeric)
            require(current_digest == baseline_digest, f"Stage3B authority differs: {relative}")

        pre = partition(stage3b / "03_positive_precluster/marker_target_membership.tsv")
        base_pre = partition(base_stage3b / "03_positive_precluster/marker_target_membership.tsv")
        require(pre[:3] == base_pre[:3], "positive precluster partition differs from baseline")
        chosen = partition(stage3b / "05_signed_leiden/chosen_membership.tsv")
        base_chosen = partition(base_stage3b / "05_signed_leiden/chosen_membership.tsv")
        require(chosen[:3] == base_chosen[:3] and chosen[0] == 62 and chosen[1] == 2567,
                "chosen signed-Leiden partition differs from baseline")
        unagg = unaggregated(stage3b / "05_signed_leiden/unaggregated_sags.tsv")
        base_unagg = unaggregated(base_stage3b / "05_signed_leiden/unaggregated_sags.tsv")
        require(unagg == base_unagg and len(unagg) == 3267, "unaggregated set differs from baseline")
        require(not (chosen[3] & unagg) and chosen[3] | unagg == node_ids,
                "assigned/unaggregated do not form exact graph-node partition")

        current_report = load_json(stage3b / "05_signed_leiden/signed_report.json")
        baseline_report = load_json(base_stage3b / "05_signed_leiden/signed_report.json")
        for key in ("chosen", "chosen_score", "clusters_ge10_emitted", "scientific_contract", "sweep"):
            require(current_report.get(key) == baseline_report.get(key),
                    f"signed-Leiden report differs from baseline: {key}")
        require(load_json(stage3b / "03_positive_precluster/precluster_report.json") ==
                load_json(base_stage3b / "03_positive_precluster/precluster_report.json"),
                "positive precluster report differs from baseline")

        clusters = group_table(stage3b / "06_subassemble/clusters.tsv", "cluster")
        base_clusters = group_table(base_stage3b / "06_subassemble/clusters.tsv", "cluster")
        require(len(clusters) == 62 and clusters.keys() == base_clusters.keys(),
                "Stage3B cluster ID set differs")
        member_counts = Counter()
        _, chosen_rows = tsv_rows(stage3b / "05_signed_leiden/chosen_membership.tsv")
        for row in chosen_rows:
            member_counts[row["cluster"]] += 1
        for key, row in clusters.items():
            require(int(row["n_sags"]) == int(row["n_found"]) == member_counts[key],
                    f"Stage3B cluster membership/file count mismatch: {key}")
            require(Path(row["input_fasta"]).is_file() and Path(row["bin_fasta"]).is_file(),
                    f"Stage3B cluster FASTA missing: {key}")
            receipt = load_json(Path(row["input_fasta"]).parent / "MERGED_INPUT.PASS.json")
            require(receipt.get("status") == "PASS" and int(receipt["member_count"]) == member_counts[key],
                    f"Stage3B merged-input receipt mismatch: {key}")
        if args.deep:
            compare_file_sets_by_sha(
                {k: Path(v["input_fasta"]) for k, v in clusters.items()},
                {k: Path(v["input_fasta"]) for k, v in base_clusters.items()}, cache,
                "Stage3B merged inputs",
            )
        metrics["stage3b"] = {
            "nodes": 5834, "pairs_expected_and_evaluated": 17014861,
            "sparse_rows": tri[0], "positive_edges": pos[0], "negative_edges": neg[0],
            "groups": chosen[0], "assigned": chosen[1], "unaggregated": len(unagg),
            "triangle_normalized_sha256": tri[1],
            "positive_edges_normalized_sha256": pos[1],
            "negative_edges_normalized_sha256": neg[1],
            "membership_partition_sha256": chosen[2],
        }
    section("stage3b_exact_graph_and_baseline_equivalence", audit_stage3b)

    def audit_stage3b_receipts() -> None:
        standard = [
            stage3b / "02_pairwise_ani_af/PASS.json",
            stage3b / "03_positive_precluster/PASS.json",
            stage3b / "05_signed_leiden/PASS.json",
        ]
        marker_make = sorted((stage3b / "04_marker_blastn").glob("M*/MAKEBLASTDB.PASS.json"))
        marker_blast = sorted((stage3b / "04_marker_blastn").glob("M*/BLASTN.PASS.json"))
        merged = sorted((stage3b / "06_subassemble").glob("*/MERGED_INPUT.PASS.json"))
        require(len(marker_make) == 123 and len(marker_blast) == 123 and len(merged) == 62,
                "Stage3B receipt file-set count mismatch")
        entries = 0
        for receipt in standard + marker_make + marker_blast:
            entries += validate_receipt(receipt, cache, include_inputs=args.deep)
        for receipt in merged:
            value = load_json(receipt)
            require(value.get("schema") == "sag-stage3b-merged-input-v1" and
                    value.get("status") == "PASS", f"bad merged-input receipt: {receipt}")
        metrics["receipts"] = {
            "standard": len(standard), "makeblastdb": len(marker_make),
            "blastn": len(marker_blast), "merged_inputs": len(merged),
            "individual_path_hashes_checked": entries,
        }
    section("stage3b_receipts", audit_stage3b_receipts)

    def audit_final_quality() -> None:
        view = load_json(final_view / "COMPLETE.json")
        require_fields(view, {
            "schema": "current-baseline-bin-evaluation-input-v1", "status": "PASS",
            "bin_count": 336, "unique_bin_id_count": 336, "nonempty_bin_count": 336,
            "stage3a_bin_count": 274, "stage3b_bin_count": 62,
            "stage3b_graph_nodes": 5834, "stage3b_assigned_sags": 2567,
            "stage3b_unaggregated_sags": 3267,
            "stage3b_entry_closure_mode": "current_v2_full_entry_closure",
            "symlink_count": 336,
        }, "final-view COMPLETE")
        require(all(value is True for value in view.get("checks", {}).values()),
                "final-view COMPLETE contains a false check")
        manifest_binding = view.get("input_manifest", {})
        require(sha256(final_view / "INPUT_MANIFEST.tsv", cache) == manifest_binding.get("sha256") and
                (final_view / "INPUT_MANIFEST.tsv").stat().st_size == int(manifest_binding.get("bytes", -1)),
                "final-view manifest receipt binding mismatch")
        for name, binding in view.get("sources", {}).items():
            target = Path(binding.get("path", ""))
            require(target.is_file(), f"final-view source binding missing ({name}): {target}")
            require(target.stat().st_size == int(binding.get("bytes", -1)),
                    f"final-view source binding size mismatch ({name}): {target}")
            require(sha256(target, cache) == binding.get("sha256"),
                    f"final-view source binding SHA mismatch ({name}): {target}")
        _, manifest_rows = tsv_rows(final_view / "INPUT_MANIFEST.tsv")
        manifest_ids = [row["bin_id"] for row in manifest_rows]
        require(len(manifest_ids) == len(set(manifest_ids)) == 336, "final manifest ID closure failed")
        _, baseline_manifest_rows = tsv_rows(base_view / "INPUT_MANIFEST.tsv")
        require(set(manifest_ids) == {row["bin_id"] for row in baseline_manifest_rows},
                "final 336-bin ID set differs from baseline")
        for row in manifest_rows:
            source = Path(row["source_bin_fasta"])
            link = Path(row["input_link"])
            require(source.is_file() and source.stat().st_size == int(row["source_bin_bytes"]),
                    f"final source bin missing/size mismatch: {row['bin_id']}")
            require(sha256(source, cache) == row["source_bin_sha256"],
                    f"final source bin SHA mismatch: {row['bin_id']}")
            require(link.is_symlink() and link.resolve(strict=True) == source.resolve(strict=True),
                    f"final-view symlink closure failed: {row['bin_id']}")

        report_path = final_eval / "checkm2_raw/quality_report.tsv"
        header, quality_rows = tsv_rows(report_path)
        require({"Name", "Completeness", "Contamination"}.issubset(header),
                "final CheckM2 report lacks required columns")
        names = [row["Name"] for row in quality_rows]
        require(len(names) == len(set(names)) == 336 and set(names) == set(manifest_ids),
                "final CheckM2 IDs do not exactly equal final-view IDs")
        hq: set[str] = set()
        mq_inclusive: set[str] = set()
        for row in quality_rows:
            completeness = float(row["Completeness"])
            contamination = float(row["Contamination"])
            require(math.isfinite(completeness) and math.isfinite(contamination),
                    f"non-finite final CheckM2 metric: {row['Name']}")
            require(0 <= completeness <= 100 and contamination >= 0,
                    f"out-of-range final CheckM2 metric: {row['Name']}")
            if completeness >= 90 and contamination <= 5:
                hq.add(row["Name"])
            if completeness >= 50 and contamination <= 10:
                mq_inclusive.add(row["Name"])
        mq_only = mq_inclusive - hq
        result_json = load_json(final_eval / "RESULT.json")
        require_fields(result_json, {
            "schema": "cellbit-current-tractor-final-quality-v1", "status": "PASS",
            "checkm2_version": "1.0.1", "checkm2_model_mode": "auto",
            "allmodels": False, "minimum_bin_size_filter": None, "total_bins": 336,
            "hq": len(hq), "mq_inclusive": len(mq_inclusive), "mq_only": len(mq_only),
            "lq_or_below": 336 - len(hq) - len(mq_only),
            "exact_bin_id_set_closed": True, "stage3b_full_entry_closure_verified": True,
        }, "final quality RESULT")
        require(sha256(report_path, cache) == result_json["quality_report_sha256"],
                "final CheckM2 report SHA binding mismatch")
        require(sha256(final_view / "INPUT_MANIFEST.tsv", cache) == result_json["input_manifest_sha256"],
                "final manifest SHA binding mismatch")
        metrics["final_quality"] = {
            "total_bins": 336, "hq": len(hq), "mq_inclusive": len(mq_inclusive),
            "mq_only": len(mq_only), "lq_or_below": 336 - len(hq) - len(mq_only),
        }
    section("final_view_336_and_checkm2_id_quality_closure", audit_final_quality)

    def audit_time_and_resources() -> None:
        time_data = {name: parse_time_v(control / name) for name in TIME_FILES}
        telemetry = validate_telemetry(control / "RESOURCE_TELEMETRY.tsv")
        require(telemetry["peak_cgroup_memory_bytes"] <= 1556925644800,
                "telemetry cgroup memory peak exceeds hard memory limit")
        phase_events = (control / "PHASE_EVENTS.tsv").read_text(encoding="utf-8").splitlines()
        expected_event_order = [
            "preflight_authority", "stage1_3a", "stage3b_input_view",
            "stage3b_pending_stats", "stage3b_upstream_checkm2",
            "stage3b_upstream_gtdbtk_identify", "stage3b_finalize_inputs",
            "stage3b_exact_and_subassemble", "final_bin_view",
            "final_checkm2_evaluation", "final_quality_closure", "complete",
        ]
        require(phase_events and phase_events[0] == "timestamp_epoch\tphase",
                "phase-events header mismatch")
        observed_events = [line.split("\t", 1)[1] for line in phase_events[1:]]
        require(observed_events == expected_event_order,
                f"phase-events order differs: {observed_events!r}")
        metrics["times"] = time_data
        metrics["telemetry"] = telemetry
    section("all_times_and_cgroup_telemetry", audit_time_and_resources)

    def audit_final_sha_chains() -> None:
        metrics["sha_manifests"] = {
            "runner_final": validate_hash_manifest(control / "FINAL_SHA256SUMS.txt", cache),
            "control_final": validate_hash_manifest(control / "FINAL_CONTROL_SHA256SUMS.txt", cache),
        }
        run_result = load_json(control / "RUN_RESULT.json")
        final_result = load_json(final_eval / "RESULT.json")
        require(run_result.get("stage3b_complete_sha256") == sha256(stage3b / "COMPLETE.json", cache),
                "RUN_RESULT Stage3B COMPLETE binding mismatch")
        require(run_result.get("final_hq") == final_result.get("hq") and
                run_result.get("final_mq_only") == final_result.get("mq_only"),
                "RUN_RESULT final HQ/MQ values differ from final RESULT")
    section("final_sha_and_result_chains", audit_final_sha_chains)

    result["status"] = "PASS" if not failures else "FAIL"
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
