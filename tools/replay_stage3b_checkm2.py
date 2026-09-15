#!/usr/bin/env python3
"""Write-once CheckM2 validation for a Stage-3B historical Flye replay.

The runner is intentionally server/path specific.  It binds the complete
Flye evidence chain, stages exactly G0001..G0058 as read-only symlinks, runs
the historical CheckM2 1.0.1 command, and reports every numerical difference
from the immutable Lake authority.  It never writes below Lake and has no
resume or overwrite mode.

The historical Flye command did not use ``--deterministic``.  Therefore this
runner can accept a Flye audit whose *only* failed checks are historical
assembly canonical equality and historical output aggregate equality.  All
58 Flye jobs, the raw/signed upstream chain, and the materialised Flye inputs
must still have passed exactly.  This narrowly-scoped exception is recorded
in every downstream receipt.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


LAKE = Path("/home/data/fyc/lake").resolve()
DATABASE = Path(
    "/home/data/fyc/past/root_archive_20260806/checkm2_db/"
    "CheckM2_database/uniref100.KO.1.dmnd"
).resolve()
EXPECTED_DATABASE_BYTES = 3_082_500_605
HISTORICAL_REPORT = (
    LAKE
    / "25_recluster_1530_newmethod/global_signed/checkm2_out/quality_report.tsv"
).resolve()
EXPECTED_GROUPS = tuple(f"G{i:04d}" for i in range(1, 59))
EXPECTED_HISTORICAL_COUNTS = {"HQ": 0, "MQ": 5, "LQ": 53}
EXPECTED_HISTORICAL_MQ = ("G0008", "G0009", "G0015", "G0019", "G0020")
ALLOWED_NONDETERMINISTIC_FLYE_FAILURES = frozenset(
    {
        "historical_assemblies_canonical_exact",
        "historical_output_aggregate_exact",
    }
)
REQUIRED_FLYE_TRUE_CHECKS = frozenset(
    {
        "all_flye_commands_successful",
        "exact_58_groups",
        "expected_effective_sag_count",
        "historical_input_aggregate_exact",
        "historical_inputs_canonical_exact",
        "membership_3759",
        "upstream_graph_marker_signed_pass",
        "upstream_raw_pass",
    }
)


class ValidationError(RuntimeError):
    """A fail-closed evidence or execution error."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def write_json_once(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValidationError(f"refusing to overwrite {path}")
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists():
            raise ValidationError(f"refusing to overwrite {path}")
        os.rename(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_tsv_once(path: Path, header: list[str], rows: list[list[object]]) -> None:
    if path.exists():
        raise ValidationError(f"refusing to overwrite {path}")
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def load_json(path: Path, label: str) -> dict:
    if not path.is_file():
        raise ValidationError(f"missing {label}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"invalid {label}: {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValidationError(f"{label} is not a JSON object: {path}")
    return value


def load_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValidationError(f"missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValidationError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def require_inside(path: Path, root: Path, label: str) -> None:
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValidationError(f"{label} escapes its declared root: {path}") from error


def resolve_executable(specification: str) -> Path:
    if os.sep in specification or (os.altsep and os.altsep in specification):
        candidate = Path(specification).expanduser()
    else:
        located = shutil.which(specification)
        if not located:
            raise ValidationError(f"CheckM2 executable not found on PATH: {specification}")
        candidate = Path(located)
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as error:
        raise ValidationError(f"invalid CheckM2 executable: {candidate}: {error}") from error
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise ValidationError(f"CheckM2 is not an executable file: {resolved}")
    return resolved


def checkm2_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
        check=False,
    )
    version = result.stdout.strip()
    accepted = any(
        re.fullmatch(r"(?:CheckM2(?:\s+version)?\s*)?1\.0\.1", line.strip(), re.I)
        for line in version.splitlines()
    )
    if result.returncode != 0 or not accepted:
        raise ValidationError(
            f"expected CheckM2 1.0.1, rc={result.returncode}, observed={version!r}"
        )
    return version


def bind_flye_replay(root: Path) -> tuple[dict, dict, dict[str, dict], dict]:
    complete_path = root / "COMPLETE.json"
    audit_path = root / "REPLAY_AUDIT.json"
    preregistration_path = root / "PREREGISTRATION.json"
    manifest_path = root / "GROUP_MANIFEST.tsv"
    complete = load_json(complete_path, "Flye COMPLETE")
    audit = load_json(audit_path, "Flye replay audit")
    preregistration = load_json(preregistration_path, "Flye preregistration")

    expected_hashes = {
        "audit_sha256": (audit_path, sha256(audit_path)),
        "preregistration_sha256": (preregistration_path, sha256(preregistration_path)),
        "group_manifest_sha256": (manifest_path, sha256(manifest_path)),
    }
    for field, (path, observed) in expected_hashes.items():
        if complete.get(field) != observed:
            raise ValidationError(f"Flye COMPLETE hash mismatch for {path}")
    if complete.get("config_sha256") != preregistration.get("config_sha256"):
        raise ValidationError("Flye config hash is not closed across COMPLETE/PREREGISTRATION")

    checks = audit.get("checks")
    if not isinstance(checks, dict) or not checks:
        raise ValidationError("Flye audit has no check map")
    missing_required = REQUIRED_FLYE_TRUE_CHECKS - checks.keys()
    if missing_required:
        raise ValidationError(f"Flye audit lacks required checks: {sorted(missing_required)}")
    failed = {name for name, result in checks.items() if result is not True}
    if complete.get("status") == "PASS" and audit.get("status") == "PASS":
        if failed:
            raise ValidationError(f"Flye PASS contains failed checks: {sorted(failed)}")
        exception_applied = False
    elif complete.get("status") == "FAIL" and audit.get("status") == "FAIL":
        if failed != ALLOWED_NONDETERMINISTIC_FLYE_FAILURES:
            raise ValidationError(
                "Flye FAIL is not solely the preregistered non-determinism exception: "
                f"observed={sorted(failed)}"
            )
        exception_applied = True
    else:
        raise ValidationError("Flye COMPLETE and audit status disagree or are invalid")
    for name in REQUIRED_FLYE_TRUE_CHECKS:
        if checks.get(name) is not True:
            raise ValidationError(f"mandatory Flye check failed: {name}")

    observed = audit.get("observed")
    if not isinstance(observed, dict):
        raise ValidationError("Flye audit lacks observed summary")
    required_observed = {
        "groups": 58,
        "membership_sags": 3759,
        "found_sags": 3758,
        "all_input_bytes_equal": True,
        "all_input_sequences_equal": True,
    }
    for field, expected in required_observed.items():
        if observed.get(field) != expected:
            raise ValidationError(
                f"Flye observed closure mismatch: {field}={observed.get(field)!r}, "
                f"expected={expected!r}"
            )

    details = audit.get("groups")
    if not isinstance(details, dict) or tuple(sorted(details)) != EXPECTED_GROUPS:
        raise ValidationError("Flye audit does not contain exactly G0001..G0058")
    _, manifest_rows = load_tsv(manifest_path)
    manifest_by_group: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        group = row.get("cluster", "")
        if group in manifest_by_group:
            raise ValidationError(f"duplicate Flye manifest group: {group}")
        manifest_by_group[group] = row
    if tuple(sorted(manifest_by_group)) != EXPECTED_GROUPS:
        raise ValidationError("Flye GROUP_MANIFEST is not exact G0001..G0058")

    assemblies: dict[str, dict] = {}
    for group in EXPECTED_GROUPS:
        detail = details[group]
        if not isinstance(detail, dict) or detail.get("cluster") != group:
            raise ValidationError(f"invalid Flye group detail: {group}")
        execution = detail.get("flye")
        assembly_record = detail.get("assembly")
        if not isinstance(execution, dict) or execution.get("returncode") != 0:
            raise ValidationError(f"Flye command did not complete successfully: {group}")
        if not isinstance(assembly_record, dict):
            raise ValidationError(f"missing Flye assembly receipt: {group}")
        expected_path = (root / "flye" / group / "assembly.fasta").resolve(strict=True)
        require_inside(expected_path, root, f"Flye assembly {group}")
        recorded_paths = [execution.get("assembly"), assembly_record.get("path")]
        for recorded in recorded_paths:
            if not recorded or Path(recorded).resolve(strict=True) != expected_path:
                raise ValidationError(f"Flye assembly path binding mismatch: {group}")
        if not expected_path.is_file() or expected_path.stat().st_size <= 0:
            raise ValidationError(f"missing/empty Flye assembly: {group}")
        observed_hash = sha256(expected_path)
        if assembly_record.get("sha256") != observed_hash:
            raise ValidationError(f"Flye assembly hash changed after audit: {group}")
        manifest_hash = manifest_by_group[group].get("assembly_sha256")
        if manifest_hash != observed_hash:
            raise ValidationError(f"Flye manifest/audit assembly hash mismatch: {group}")
        assemblies[group] = {
            "path": expected_path,
            "sha256": observed_hash,
            "file_bytes": expected_path.stat().st_size,
            "flye_contigs": assembly_record.get("contigs"),
            "flye_bp": assembly_record.get("bp"),
        }

    binding = {
        "schema": "stage3b-checkm2-flye-binding-v1",
        "flye_root": str(root),
        "complete": {"path": str(complete_path), "sha256": sha256(complete_path)},
        "audit": {"path": str(audit_path), "sha256": sha256(audit_path)},
        "preregistration": {
            "path": str(preregistration_path),
            "sha256": sha256(preregistration_path),
        },
        "group_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        "reported_status": audit.get("status"),
        "failed_checks": sorted(failed),
        "nondeterministic_flye_exception_applied": exception_applied,
        "exception_scope": sorted(ALLOWED_NONDETERMINISTIC_FLYE_FAILURES),
        "mandatory_checks": sorted(REQUIRED_FLYE_TRUE_CHECKS),
    }
    return complete, audit, assemblies, binding


def decimal_field(row: dict[str, str], field: str, group: str) -> Decimal:
    try:
        value = Decimal(row[field])
    except (KeyError, InvalidOperation) as error:
        raise ValidationError(f"invalid {field} for {group}: {row.get(field)!r}") from error
    if not value.is_finite():
        raise ValidationError(f"non-finite {field} for {group}: {value}")
    return value


def classify(completeness: Decimal, contamination: Decimal) -> str:
    if completeness >= Decimal("90") and contamination <= Decimal("5"):
        return "HQ"
    if completeness >= Decimal("50") and contamination <= Decimal("10"):
        return "MQ"
    return "LQ"


def parse_quality_report(path: Path, label: str) -> dict[str, dict[str, object]]:
    header, rows = load_tsv(path)
    required = {"Name", "Completeness", "Contamination"}
    if not required.issubset(header):
        raise ValidationError(f"{label} lacks columns: {sorted(required - set(header))}")
    parsed: dict[str, dict[str, object]] = {}
    for row in rows:
        group = row["Name"].strip()
        if group in parsed:
            raise ValidationError(f"duplicate group in {label}: {group}")
        completeness = decimal_field(row, "Completeness", group)
        contamination = decimal_field(row, "Contamination", group)
        parsed[group] = {
            "completeness": completeness,
            "contamination": contamination,
            "quality": classify(completeness, contamination),
        }
    if tuple(sorted(parsed)) != EXPECTED_GROUPS:
        missing = sorted(set(EXPECTED_GROUPS) - set(parsed))
        extra = sorted(set(parsed) - set(EXPECTED_GROUPS))
        raise ValidationError(
            f"{label} is not exact G0001..G0058: missing={missing}, extra={extra}"
        )
    return parsed


def quality_counts(values: dict[str, dict[str, object]]) -> dict[str, int]:
    return {
        quality: sum(row["quality"] == quality for row in values.values())
        for quality in ("HQ", "MQ", "LQ")
    }


def decimal_text(value: Decimal) -> str:
    return format(value, "f")


def compare_quality(
    observed: dict[str, dict[str, object]], historical: dict[str, dict[str, object]]
) -> tuple[list[list[object]], dict[str, object]]:
    rows: list[list[object]] = []
    changed_metrics: list[str] = []
    changed_classes: list[str] = []
    completeness_deltas: list[Decimal] = []
    contamination_deltas: list[Decimal] = []
    for group in EXPECTED_GROUPS:
        current = observed[group]
        old = historical[group]
        dc = current["completeness"] - old["completeness"]
        dx = current["contamination"] - old["contamination"]
        completeness_deltas.append(dc)
        contamination_deltas.append(dx)
        metric_exact = dc == 0 and dx == 0
        class_exact = current["quality"] == old["quality"]
        if not metric_exact:
            changed_metrics.append(group)
        if not class_exact:
            changed_classes.append(group)
        rows.append(
            [
                group,
                decimal_text(current["completeness"]),
                decimal_text(old["completeness"]),
                decimal_text(dc),
                decimal_text(current["contamination"]),
                decimal_text(old["contamination"]),
                decimal_text(dx),
                current["quality"],
                old["quality"],
                metric_exact,
                class_exact,
            ]
        )
    observed_counts = quality_counts(observed)
    historical_counts = quality_counts(historical)
    count_delta = {
        quality: observed_counts[quality] - historical_counts[quality]
        for quality in ("HQ", "MQ", "LQ")
    }
    n = Decimal(len(EXPECTED_GROUPS))
    summary = {
        "observed_counts": observed_counts,
        "historical_counts": historical_counts,
        "count_delta": count_delta,
        "observed_groups": {
            quality: [g for g in EXPECTED_GROUPS if observed[g]["quality"] == quality]
            for quality in ("HQ", "MQ", "LQ")
        },
        "historical_groups": {
            quality: [g for g in EXPECTED_GROUPS if historical[g]["quality"] == quality]
            for quality in ("HQ", "MQ", "LQ")
        },
        "per_group_metrics_exact": not changed_metrics,
        "per_group_classes_exact": not changed_classes,
        "changed_metric_groups": changed_metrics,
        "changed_class_groups": changed_classes,
        "mean_completeness_delta": decimal_text(sum(completeness_deltas) / n),
        "mean_contamination_delta": decimal_text(sum(contamination_deltas) / n),
        "max_abs_completeness_delta": decimal_text(
            max(abs(value) for value in completeness_deltas)
        ),
        "max_abs_contamination_delta": decimal_text(
            max(abs(value) for value in contamination_deltas)
        ),
    }
    if not changed_metrics:
        summary["equivalence"] = "EXACT"
    elif observed_counts == historical_counts and not changed_classes:
        summary["equivalence"] = "AGGREGATE_AND_CLASSES_EXACT_METRICS_DIFFER"
    elif observed_counts == historical_counts:
        summary["equivalence"] = "AGGREGATE_EXACT_CLASS_ASSIGNMENTS_DIFFER"
    else:
        summary["equivalence"] = "AGGREGATE_DIFFERENT"
    return rows, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flye-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--checkm2", default="checkm2")
    parser.add_argument("--threads", type=int, default=16)
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    flye_root = args.flye_root.resolve(strict=True)
    output = args.out.resolve()
    if output.exists():
        raise ValidationError(f"write-once output already exists: {output}")
    try:
        output.relative_to(LAKE)
    except ValueError:
        pass
    else:
        raise ValidationError("output must not be inside immutable Lake")
    if args.threads < 1:
        raise ValidationError("--threads must be positive")
    if not DATABASE.is_file() or DATABASE.stat().st_size != EXPECTED_DATABASE_BYTES:
        raise ValidationError(
            f"CheckM2 database authority missing/size mismatch: {DATABASE}"
        )
    executable = resolve_executable(args.checkm2)
    version = checkm2_version(executable)
    historical = parse_quality_report(HISTORICAL_REPORT, "historical CheckM2 report")
    historical_counts = quality_counts(historical)
    historical_mq = tuple(
        group for group in EXPECTED_GROUPS if historical[group]["quality"] == "MQ"
    )
    if historical_counts != EXPECTED_HISTORICAL_COUNTS or historical_mq != EXPECTED_HISTORICAL_MQ:
        raise ValidationError(
            "historical authority no longer matches preregistered 0 HQ / 5 MQ closure"
        )
    _, _, assemblies, flye_binding = bind_flye_replay(flye_root)

    output.mkdir(parents=True)
    input_dir = output / "checkm2_in"
    checkm2_out = output / "checkm2_out"
    input_dir.mkdir()
    if checkm2_out.exists():
        raise ValidationError(f"fresh CheckM2 output unexpectedly exists: {checkm2_out}")

    command = [
        str(executable),
        "predict",
        "--database_path",
        str(DATABASE),
        "--input",
        str(input_dir),
        "--output-directory",
        str(checkm2_out),
        "-x",
        "fna",
        "--threads",
        str(args.threads),
        "--force",
    ]
    config = {
        "schema": "stage3b-checkm2-config-v1",
        "groups": list(EXPECTED_GROUPS),
        "thresholds": {
            "HQ": {"completeness_gte": 90, "contamination_lte": 5},
            "MQ": {
                "completeness_gte": 50,
                "contamination_lte": 10,
                "excluding_HQ": True,
            },
        },
        "threads": args.threads,
        "command": command,
        "database": str(DATABASE),
        "historical_report": str(HISTORICAL_REPORT),
        "acceptance": {
            "mandatory": [
                "CheckM2 1.0.1 exits zero",
                "exact G0001..G0058 report closure",
                "all staged Flye assemblies remain byte-identical",
                "observed HQ/MQ/LQ aggregate equals historical 0/5/53",
            ],
            "per_group_metric_exactness": (
                "reported but not mandatory because historical Flye omitted --deterministic"
            ),
        },
    }
    preregistration = {
        "schema": "stage3b-checkm2-preregistration-v1",
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "write_once_root": str(output),
        "config": config,
        "config_sha256": json_sha(config),
        "tool": {
            "path": str(executable),
            "sha256": sha256(executable),
            "version": version,
        },
        "database": {
            "path": str(DATABASE),
            "file_bytes": DATABASE.stat().st_size,
            "mtime_ns": DATABASE.stat().st_mtime_ns,
            "note": "Large immutable authority is bound by exact path and expected size.",
        },
        "historical_authority": {
            "path": str(HISTORICAL_REPORT),
            "sha256": sha256(HISTORICAL_REPORT),
            "counts": historical_counts,
            "mq_groups": list(historical_mq),
        },
        "flye_binding": flye_binding,
    }
    write_json_once(output / "PREREGISTRATION.json", preregistration)

    manifest_rows: list[list[object]] = []
    for group in EXPECTED_GROUPS:
        source = assemblies[group]["path"]
        staged = input_dir / f"{group}.fna"
        if staged.exists() or staged.is_symlink():
            raise ValidationError(f"refusing to overwrite staged input: {staged}")
        staged.symlink_to(source)
        if not staged.is_symlink() or staged.resolve(strict=True) != source:
            raise ValidationError(f"failed to bind staged symlink: {group}")
        if sha256(staged) != assemblies[group]["sha256"]:
            raise ValidationError(f"staged assembly differs from Flye receipt: {group}")
        manifest_rows.append(
            [
                group,
                str(staged),
                os.readlink(staged),
                str(source),
                assemblies[group]["file_bytes"],
                assemblies[group]["sha256"],
                assemblies[group]["flye_contigs"],
                assemblies[group]["flye_bp"],
            ]
        )
    staged_names = tuple(sorted(path.name for path in input_dir.iterdir()))
    if staged_names != tuple(f"{group}.fna" for group in EXPECTED_GROUPS):
        raise ValidationError("staged CheckM2 input is not exact G0001..G0058")
    write_tsv_once(
        output / "INPUT_MANIFEST.tsv",
        [
            "group",
            "staged_path",
            "symlink_target",
            "resolved_source_path",
            "file_bytes",
            "sha256",
            "flye_contigs",
            "flye_bp",
        ],
        manifest_rows,
    )

    log_path = output / "CHECKM2_RUN.log"
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    started = time.monotonic()
    with log_path.open("xb") as log:
        log.write((json.dumps(command) + "\n").encode())
        log.flush()
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        wall_seconds = time.monotonic() - started
        log.write(
            f"\nreturncode={process.returncode}\nwall_seconds={wall_seconds:.6f}\n".encode()
        )
        log.flush()
        os.fsync(log.fileno())
    run_receipt = {
        "schema": "stage3b-checkm2-run-receipt-v1",
        "command": command,
        "returncode": process.returncode,
        "started_utc": started_utc,
        "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "wall_seconds": wall_seconds,
        "log": {"path": str(log_path), "sha256": sha256(log_path)},
    }
    write_json_once(output / "RUN_RECEIPT.json", run_receipt)
    if process.returncode != 0:
        raise ValidationError(f"CheckM2 failed: rc={process.returncode}; see {log_path}")

    for group in EXPECTED_GROUPS:
        staged = input_dir / f"{group}.fna"
        source = assemblies[group]["path"]
        if not staged.is_symlink() or staged.resolve(strict=True) != source:
            raise ValidationError(f"staged symlink changed during CheckM2: {group}")
        if sha256(staged) != assemblies[group]["sha256"]:
            raise ValidationError(f"Flye assembly changed during CheckM2: {group}")
    report = checkm2_out / "quality_report.tsv"
    observed = parse_quality_report(report, "replayed CheckM2 report")
    comparison_rows, comparison = compare_quality(observed, historical)
    write_tsv_once(
        output / "QUALITY_COMPARISON.tsv",
        [
            "group",
            "observed_completeness",
            "historical_completeness",
            "delta_completeness",
            "observed_contamination",
            "historical_contamination",
            "delta_contamination",
            "observed_quality",
            "historical_quality",
            "metrics_exact",
            "class_exact",
        ],
        comparison_rows,
    )

    checks = {
        "flye_evidence_closed_under_declared_exception": True,
        "checkm2_version_exact_1_0_1": True,
        "database_path_and_size_exact": DATABASE.stat().st_size == EXPECTED_DATABASE_BYTES,
        "checkm2_returncode_zero": process.returncode == 0,
        "exact_58_staged_symlinks": staged_names
        == tuple(f"{group}.fna" for group in EXPECTED_GROUPS),
        "staged_sources_unchanged": all(
            sha256(input_dir / f"{group}.fna") == assemblies[group]["sha256"]
            for group in EXPECTED_GROUPS
        ),
        "exact_58_quality_rows": tuple(sorted(observed)) == EXPECTED_GROUPS,
        "historical_authority_is_0_HQ_5_MQ": historical_counts
        == EXPECTED_HISTORICAL_COUNTS,
        "observed_aggregate_quality_matches_historical": comparison["observed_counts"]
        == comparison["historical_counts"],
    }
    audit = {
        "schema": "stage3b-checkm2-audit-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "flye_binding": flye_binding,
        "historical_flye_nondeterminism_policy": {
            "reason": (
                "The historical Flye run omitted --deterministic; only the two declared "
                "assembly/output equivalence checks may be false upstream."
            ),
            "exception_applied": flye_binding[
                "nondeterministic_flye_exception_applied"
            ],
            "allowed_failed_checks": sorted(ALLOWED_NONDETERMINISTIC_FLYE_FAILURES),
            "per_group_CheckM2_metric_exactness_is_mandatory": False,
            "aggregate_HQ_MQ_LQ_exactness_is_mandatory": True,
        },
        "tool": {"path": str(executable), "version": version},
        "database": {"path": str(DATABASE), "file_bytes": DATABASE.stat().st_size},
        "input_manifest": {
            "path": str(output / "INPUT_MANIFEST.tsv"),
            "sha256": sha256(output / "INPUT_MANIFEST.tsv"),
        },
        "observed_report": {"path": str(report), "sha256": sha256(report)},
        "historical_report": {
            "path": str(HISTORICAL_REPORT),
            "sha256": sha256(HISTORICAL_REPORT),
        },
        "comparison": comparison,
        "comparison_table": {
            "path": str(output / "QUALITY_COMPARISON.tsv"),
            "sha256": sha256(output / "QUALITY_COMPARISON.tsv"),
        },
        "interpretation": (
            "EXACT means all 58 completeness/contamination pairs match.  Other equivalence "
            "labels preserve the observed numerical deltas; overall PASS additionally "
            "requires the historical 0 HQ / 5 MQ / 53 LQ aggregate."
        ),
    }
    write_json_once(output / "AUDIT.json", audit)
    complete = {
        "schema": "stage3b-checkm2-complete-v1",
        "status": audit["status"],
        "preregistration_sha256": sha256(output / "PREREGISTRATION.json"),
        "input_manifest_sha256": sha256(output / "INPUT_MANIFEST.tsv"),
        "run_receipt_sha256": sha256(output / "RUN_RECEIPT.json"),
        "quality_report_sha256": sha256(report),
        "quality_comparison_sha256": sha256(output / "QUALITY_COMPARISON.tsv"),
        "audit_sha256": sha256(output / "AUDIT.json"),
        "config_sha256": preregistration["config_sha256"],
        "observed_counts": comparison["observed_counts"],
        "historical_counts": comparison["historical_counts"],
        "equivalence": comparison["equivalence"],
    }
    write_json_once(output / "COMPLETE.json", complete)
    print(
        json.dumps(
            {
                "status": audit["status"],
                "observed_counts": comparison["observed_counts"],
                "historical_counts": comparison["historical_counts"],
                "equivalence": comparison["equivalence"],
                "changed_metric_groups": comparison["changed_metric_groups"],
            },
            indent=2,
        )
    )
    return 0 if audit["status"] == "PASS" else 2


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except (ValidationError, OSError, subprocess.SubprocessError) as error:
        output = args.out.resolve()
        if output.is_dir() and not (output / "FAILURE.json").exists():
            try:
                write_json_once(
                    output / "FAILURE.json",
                    {
                        "schema": "stage3b-checkm2-failure-v1",
                        "status": "FAIL",
                        "error_type": type(error).__name__,
                        "error": str(error),
                    },
                )
            except Exception:
                pass
        print(f"fatal: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
