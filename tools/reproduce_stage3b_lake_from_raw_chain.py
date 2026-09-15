#!/usr/bin/env python3
"""Continue the Stage-3B reproduction from an audited raw-input replay.

HISTORICAL ONLY: this wrapper binds the superseded GTDB-reference
double-negative entrance.  It is not a current production runner, and its
historical checks are intentionally unchanged.

The raw replay regenerates skani search/triangle and the 7,408-SAG entrance.
This wrapper binds those exact files to the historical marker-BLAST and signed
Leiden reproduction engine, so the chain does not silently fall back to the
frozen triangle table.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


LAKE = Path("/home/data/fyc/lake").resolve()


class ChainError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("stage3b_historical_engine", path)
    if spec is None or spec.loader is None:
        raise ChainError(f"cannot import reproduction engine: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--engine-script", required=True, type=Path)
    parser.add_argument("--marker-workers", type=int, default=96)
    parser.add_argument(
        "--marker-map",
        type=Path,
        help="Fresh marker map with sibling PASS audit/COMPLETE receipts",
    )
    args = parser.parse_args()
    if not 1 <= args.marker_workers <= 96:
        raise ChainError("--marker-workers must be in [1,96]")

    raw_root = args.raw_root.resolve()
    output = args.out.resolve()
    engine_path = args.engine_script.resolve()
    if not raw_root.is_dir() or not engine_path.is_file():
        raise ChainError("raw root or engine script is missing")
    if output.exists():
        raise ChainError(f"write-once output exists: {output}")
    if str(output).startswith(str(LAKE) + os.sep):
        raise ChainError("output must not be inside immutable Lake")

    raw_audit_path = raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json"
    if not raw_audit_path.is_file():
        raise ChainError("raw reproduction audit is missing")
    raw_audit = json.loads(raw_audit_path.read_text())
    if raw_audit.get("status") != "PASS" or not all(raw_audit.get("checks", {}).values()):
        raise ChainError("raw reproduction did not pass every check")

    expected_outputs = {
        "neither": raw_root / "quality_filtered_neither_sags.tsv",
        "skani_triangle": raw_root / "skani_triangle_sparse.tsv",
    }
    recorded_paths = raw_audit.get("outputs", {})
    recorded_hashes = raw_audit.get("output_sha256", {})
    for key, path in expected_outputs.items():
        if Path(recorded_paths.get(key, "")).resolve() != path:
            raise ChainError(f"raw audit path binding mismatch for {key}")
        if not path.is_file() or recorded_hashes.get(str(path)) != sha256(path):
            raise ChainError(f"raw output hash mismatch for {key}")

    marker_binding = None
    if args.marker_map is not None:
        marker_map = args.marker_map.resolve()
        marker_audit_path = marker_map.parent / "MARKER_MAP_AUDIT.json"
        marker_complete_path = marker_map.parent / "COMPLETE.json"
        if not marker_map.is_file() or not marker_audit_path.is_file() or not marker_complete_path.is_file():
            raise ChainError("fresh marker-map file or receipt is missing")
        marker_audit = json.loads(marker_audit_path.read_text())
        marker_complete = json.loads(marker_complete_path.read_text())
        marker_output = marker_audit.get("outputs", {}).get("marker_map", {})
        if (
            marker_audit.get("status") != "PASS"
            or not all(marker_audit.get("checks", {}).values())
            or Path(marker_output.get("path", "")).resolve() != marker_map
            or marker_output.get("sha256") != sha256(marker_map)
            or marker_complete.get("status") != "PASS"
            or marker_complete.get("audit_sha256") != sha256(marker_audit_path)
        ):
            raise ChainError("fresh marker-map proof chain did not validate")
        marker_binding = {
            "path": str(marker_map),
            "sha256": sha256(marker_map),
            "audit": {"path": str(marker_audit_path), "sha256": sha256(marker_audit_path)},
            "complete": {"path": str(marker_complete_path), "sha256": sha256(marker_complete_path)},
        }

    engine = load_module(engine_path)
    engine.QUALITY = expected_outputs["neither"]
    engine.SKANI_TRIANGLE = expected_outputs["skani_triangle"]
    if marker_binding is not None:
        engine.MARKER_MAP = Path(marker_binding["path"])
    engine.require_inputs("recompute-markers")
    engine.preregister(output, "recompute-markers-from-raw-replay", args.marker_workers)
    engine.atomic_json(
        output / "UPSTREAM_RAW_BINDING.json",
        {
            "schema": "stage3b-raw-replay-binding-v1",
            "raw_root": str(raw_root),
            "raw_audit": {
                "path": str(raw_audit_path),
                "sha256": sha256(raw_audit_path),
            },
            "engine": {"path": str(engine_path), "sha256": sha256(engine_path)},
            "bound_inputs": {
                key: {"path": str(path), "sha256": sha256(path)}
                for key, path in expected_outputs.items()
            },
            "fresh_marker_map": marker_binding,
        },
    )
    engine.run_positive(output)
    engine.run_marker(output, args.marker_workers)
    engine.run_signed(output)
    result = engine.audit(output, marker_was_recomputed=True)
    chain = {
        "schema": "stage3b-lake-from-raw-chain-audit-v1",
        "status": "PASS" if result.get("status") == "PASS" else "FAIL",
        "raw_audit_status": raw_audit.get("status"),
        "graph_marker_signed_audit_status": result.get("status"),
        "raw_audit_sha256": sha256(raw_audit_path),
        "reproduction_audit_sha256": sha256(output / "REPRODUCTION_AUDIT.json"),
        "fresh_marker_map": marker_binding,
    }
    engine.atomic_json(output / "CHAIN_AUDIT.json", chain)
    print(json.dumps(chain, indent=2, sort_keys=True))
    return 0 if chain["status"] == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ChainError, RuntimeError) as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
