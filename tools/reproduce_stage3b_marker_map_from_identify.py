#!/usr/bin/env python3
"""Write-once Lake Stage-3B GTDB-Tk identify -> nucleotide marker-map replay.

The runner consumes the already-audited ``raw_input_replay_v1`` effective
FASTA manifest.  It does not write under Lake.  It recreates the historical
GTDB-Tk identify call in an isolated directory, then reconstructs
``bac120_marker_nt_map.tsv`` from the identify intermediates.

Selection is the historical rule: for every (SAG, marker), retain the hit with
the strictly greatest bitscore; an exact bitscore tie keeps the first hit seen.
Rows are emitted in effective-QF SAG order and marker-id lexical order.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


LAKE = Path("/home/data/fyc/lake").resolve()
DEFAULT_RAW = Path(
    "/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/"
    "06_stage3b_lake_reproduction_20260904/raw_input_replay_v1"
)
DEFAULT_GTDBTK = Path("/home/data/fyc/biosoft/miniconda3/envs/gtdbtk/bin/gtdbtk")
GTDB_DATA = Path("/home/data/temp/release232")
AUTHORITY_MAP = (
    LAKE
    / "10_cellbit_unknown_sag_cluster_ppt/tables/bac120_marker_nt_map.tsv"
)
AUTHORITY_SHA256 = "c00137b8acbbb89c3c626fceb511a772111447137d768e28f6de2c8d53882ab2"
BROKEN_SAG = "TCACGCGGAGCGAGCAAGGCTCT"
BROKEN_OVERRIDE = (
    LAKE / "07_spades_assemblies" / BROKEN_SAG / "K55/scaffolds.fasta"
).resolve()
EXPECTED_SAGS = 8785


class ReplayError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_exclusive(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ReplayError(f"missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def load_raw_contract(raw_root: Path) -> tuple[list[str], list[dict[str, str]]]:
    audit_path = raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json"
    manifest_path = raw_root / "EFFECTIVE_FASTA_CONTENT_MANIFEST.tsv"
    qf_path = raw_root / "quality_filtered_sags_effective.tsv"
    for path in (audit_path, manifest_path, qf_path):
        if not path.is_file():
            raise ReplayError(f"missing raw-replay authority: {path}")

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "PASS" or not all(audit.get("checks", {}).values()):
        raise ReplayError("upstream raw-input replay is not fully PASS")
    outputs = audit.get("outputs", {})
    recorded = audit.get("output_sha256", {})
    expected = {"content_manifest": manifest_path, "effective_qf": qf_path}
    for name, path in expected.items():
        if Path(outputs.get(name, "")).resolve() != path.resolve():
            raise ReplayError(f"raw audit path binding mismatch: {name}")
        if recorded.get(str(path.resolve())) != sha256(path):
            raise ReplayError(f"raw audit SHA256 binding mismatch: {name}")

    _, manifest = read_tsv(manifest_path)
    qf_fields, qf = read_tsv(qf_path)
    required_manifest = {"SAG_id", "effective_path", "file_bytes", "sha256"}
    if not manifest or not required_manifest <= set(manifest[0]):
        raise ReplayError("bad effective FASTA manifest schema")
    if "SAG_id" not in qf_fields or len(manifest) != EXPECTED_SAGS or len(qf) != EXPECTED_SAGS:
        raise ReplayError("effective input cardinality is not 8,785")
    sag_order = [row["SAG_id"] for row in qf]
    if sag_order != [row["SAG_id"] for row in manifest] or len(set(sag_order)) != EXPECTED_SAGS:
        raise ReplayError("manifest/QF SAG order or uniqueness mismatch")

    overrides = []
    override_owners = []
    for row in manifest:
        path = Path(row["effective_path"]).resolve()
        if not path.is_file() or path.stat().st_size != int(row["file_bytes"]):
            raise ReplayError(f"missing or size-changed effective FASTA: {row['SAG_id']}")
        if row["SAG_id"] == BROKEN_SAG:
            overrides.append(path)
        if path == BROKEN_OVERRIDE:
            override_owners.append(row["SAG_id"])
    if overrides != [BROKEN_OVERRIDE] or override_owners != [BROKEN_SAG]:
        raise ReplayError(
            f"unique K55 override mismatch: paths={overrides}, owners={override_owners}"
        )
    return sag_order, manifest


def parse_fasta(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split(None, 1)[0]
                if not current or current in result:
                    raise ReplayError(f"duplicate/empty FASTA gene id in {path}: {current}")
                result[current] = ""
            elif current is None:
                raise ReplayError(f"sequence before FASTA header: {path}")
            else:
                result[current] += line.upper()
    return result


def extract_marker_map(identify_root: Path, sag_order: list[str], output: Path) -> int:
    marker_root = identify_root / "identify/intermediate_results/marker_genes"
    if not marker_root.is_dir():
        raise ReplayError(f"identify marker directory missing: {marker_root}")
    total = 0
    with output.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SAG_id", "marker_id", "gene_id", "nt_len", "sequence"])
        for sag in sag_order:
            directory = marker_root / sag
            nucleotide = directory / f"{sag}_protein.fna"
            if not nucleotide.is_file():
                raise ReplayError(f"missing identify nucleotide genes for {sag}")
            sequences = parse_fasta(nucleotide)
            best: dict[str, tuple[float, int, str]] = {}
            encounter = 0
            for family in ("pfam", "tigrfam"):
                hits_path = directory / f"{sag}_{family}_tophit.tsv"
                if not hits_path.is_file():
                    raise ReplayError(f"missing identify top-hit table: {hits_path}")
                with hits_path.open(newline="") as hits_handle:
                    rows = csv.reader(hits_handle, delimiter="\t")
                    header = next(rows, None)
                    if header != ["Gene Id", "Top hits (Family id,e-value,bitscore)"]:
                        raise ReplayError(f"unexpected top-hit header: {hits_path}")
                    for row in rows:
                        if len(row) != 2 or row[0] not in sequences:
                            raise ReplayError(f"bad top-hit row in {hits_path}: {row}")
                        for raw_hit in row[1].split(";"):
                            fields = raw_hit.split(",")
                            if len(fields) != 3:
                                raise ReplayError(f"bad marker hit in {hits_path}: {raw_hit}")
                            marker, _evalue, raw_bitscore = fields
                            try:
                                bitscore = float(raw_bitscore)
                            except ValueError as error:
                                raise ReplayError(f"bad bitscore in {hits_path}: {raw_hit}") from error
                            encounter += 1
                            old = best.get(marker)
                            # Strict greater-than is deliberate: ties keep first seen.
                            if old is None or bitscore > old[0]:
                                best[marker] = (bitscore, encounter, row[0])
            for marker in sorted(best):
                gene = best[marker][2]
                sequence = sequences[gene]
                writer.writerow([sag, marker, gene, len(sequence), sequence])
                total += 1
        handle.flush()
        os.fsync(handle.fileno())
    return total


def canonical(path: Path, allowed_sags: set[str]) -> list[tuple[str, str, str, int, str]]:
    _, rows = read_tsv(path)
    required = {"SAG_id", "marker_id", "gene_id", "nt_len", "sequence"}
    if rows and not required <= set(rows[0]):
        raise ReplayError(f"bad marker-map schema: {path}")
    result = []
    keys = set()
    for row in rows:
        if row["SAG_id"] not in allowed_sags:
            continue
        key = (row["SAG_id"], row["marker_id"])
        if key in keys:
            raise ReplayError(f"duplicate (SAG,marker) in {path}: {key}")
        keys.add(key)
        sequence = row["sequence"].upper()
        if int(row["nt_len"]) != len(sequence):
            raise ReplayError(f"nt_len mismatch in {path}: {key}")
        result.append((key[0], key[1], row["gene_id"], len(sequence), sequence))
    result.sort()
    return result


def validate_existing_identify(
    identify_root: Path, raw_root: Path
) -> dict[str, object]:
    """Bind an already completed, exact raw-input identify run to this audit."""
    receipt_path = identify_root / "gtdbtk.json"
    failed_path = identify_root / "identify/gtdbtk.failed_genomes.tsv"
    if not receipt_path.is_file() or not failed_path.is_file():
        raise ReplayError("existing identify receipt or failed-genome table is missing")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    steps = receipt.get("steps", [])
    expected_genomes = (raw_root / "input_fastas").resolve()
    checks = {
        "version": receipt.get("version") == "2.7.2",
        "database_version": receipt.get("database_version") == "r232",
        "database_path": Path(receipt.get("database_path", "")).resolve() == GTDB_DATA,
        "one_completed_identify_step": len(steps) == 1
        and steps[0].get("name") == "identify"
        and steps[0].get("status") == "completed",
        "genome_dir": len(steps) == 1
        and Path(steps[0].get("genome_dir", "")).resolve() == expected_genomes,
        "extension": len(steps) == 1 and steps[0].get("extension") == "fna",
        "write_single_copy_genes": len(steps) == 1
        and steps[0].get("write_single_copy_genes") is True,
        "failed_genomes_empty": failed_path.stat().st_size == 0,
    }
    if not all(checks.values()):
        raise ReplayError(f"existing identify receipt failed validation: {checks}")
    marker_root = identify_root / "identify/intermediate_results/marker_genes"
    marker_dirs = sum(1 for path in marker_root.iterdir() if path.is_dir())
    if marker_dirs != EXPECTED_SAGS:
        raise ReplayError(
            f"existing identify marker directory count {marker_dirs} != {EXPECTED_SAGS}"
        )
    return {
        "path": str(identify_root),
        "gtdbtk_json_sha256": sha256(receipt_path),
        "failed_genomes_sha256": sha256(failed_path),
        "marker_directories": marker_dirs,
        "checks": checks,
    }


def run(args: argparse.Namespace, root: Path) -> int:
    raw_root = args.raw_root.resolve()
    if str(root).startswith(str(LAKE) + os.sep):
        raise ReplayError("output must not be inside immutable Lake")
    if not args.gtdbtk.is_file() or not GTDB_DATA.is_dir():
        raise ReplayError("GTDB-Tk executable or r232 data is missing")
    if not AUTHORITY_MAP.is_file() or sha256(AUTHORITY_MAP) != AUTHORITY_SHA256:
        raise ReplayError("frozen marker-map SHA256 mismatch")
    version = subprocess.check_output([str(args.gtdbtk), "--version"], text=True).strip()
    if "2.7.2" not in version:
        raise ReplayError(f"expected GTDB-Tk 2.7.2, observed: {version}")
    sag_order, manifest = load_raw_contract(raw_root)

    existing_receipt = None
    if args.existing_identify_root is None:
        genomes = root / "genomes"
        genomes.mkdir()
        for row in manifest:
            os.symlink(Path(row["effective_path"]).resolve(), genomes / f"{row['SAG_id']}.fna")
        identify_root = root / "gtdbtk_identify"
        command = [
            str(args.gtdbtk), "identify", "--genome_dir", str(genomes),
            "--out_dir", str(identify_root), "-x", "fna", "--cpus",
            str(args.threads), "--force", "--write_single_copy_genes",
        ]
    else:
        identify_root = args.existing_identify_root.resolve()
        existing_receipt = validate_existing_identify(identify_root, raw_root)
        command = json.loads((identify_root / "gtdbtk.json").read_text(encoding="utf-8"))[
            "command_line"
        ]
    write_json_exclusive(
        root / "PREREGISTRATION.json",
        {
            "schema": "stage3b-identify-marker-map-preregistration-v1",
            "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            "raw_root": str(raw_root),
            "raw_audit_sha256": sha256(raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json"),
            "manifest_sha256": sha256(raw_root / "EFFECTIVE_FASTA_CONTENT_MANIFEST.tsv"),
            "sags": len(sag_order),
            "single_k55_override": {BROKEN_SAG: str(BROKEN_OVERRIDE)},
            "gtdbtk_version": version,
            "gtdb_data": str(GTDB_DATA),
            "gtdb_release": "r232",
            "command": command,
            "existing_identify_receipt": existing_receipt,
            "selection_rule": "strict maximum bitscore per SAG/marker; exact tie keeps first encounter",
            "row_order": "effective-QF SAG order, then marker_id lexical order",
            "authority": {"path": str(AUTHORITY_MAP), "sha256": AUTHORITY_SHA256},
        },
    )
    if existing_receipt is None:
        environment = os.environ.copy()
        environment["GTDBTK_DATA_PATH"] = str(GTDB_DATA)
        environment["PATH"] = str(args.gtdbtk.parent) + os.pathsep + environment.get("PATH", "")
        log_path = root / "gtdbtk_identify.log"
        with log_path.open("x", encoding="utf-8") as log:
            log.write(json.dumps(command) + "\n")
            log.flush()
            result = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
            log.write(f"\nreturncode={result.returncode}\n")
            log.flush()
            os.fsync(log.fileno())
        if result.returncode:
            raise ReplayError(f"GTDB-Tk identify failed with rc={result.returncode}")
    else:
        log_path = identify_root / "gtdbtk.log"
        if not log_path.is_file():
            raise ReplayError("existing identify gtdbtk.log is missing")

    produced = root / "bac120_marker_nt_map.tsv"
    row_count = extract_marker_map(identify_root, sag_order, produced)
    observed = canonical(produced, set(sag_order))
    expected = canonical(AUTHORITY_MAP, set(sag_order))
    observed_set, expected_set = set(observed), set(expected)
    checks = {
        "upstream_raw_replay_pass": True,
        "effective_sag_count_8785": len(sag_order) == EXPECTED_SAGS,
        "single_k55_override_exact": True,
        "gtdbtk_version_2_7_2": True,
        "gtdb_release_r232": True,
        "canonical_marker_map_equal_frozen_qf_subset": observed == expected,
    }
    audit = {
        "schema": "stage3b-identify-marker-map-audit-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed_rows": row_count,
        "authority_rows_for_effective_sags": len(expected),
        "observed_only_rows": len(observed_set - expected_set),
        "authority_only_rows": len(expected_set - observed_set),
        "first_observed_only": [list(row[:4]) for row in sorted(observed_set - expected_set)[:20]],
        "first_authority_only": [list(row[:4]) for row in sorted(expected_set - observed_set)[:20]],
        "outputs": {
            "marker_map": {"path": str(produced), "sha256": sha256(produced)},
            "identify_log": {"path": str(log_path), "sha256": sha256(log_path)},
        },
    }
    write_json_exclusive(root / "MARKER_MAP_AUDIT.json", audit)
    if audit["status"] != "PASS":
        raise ReplayError("canonical marker map differs from frozen authority")
    write_json_exclusive(
        root / "COMPLETE.json",
        {"schema": "stage3b-identify-marker-map-complete-v1", "status": "PASS", "audit_sha256": sha256(root / "MARKER_MAP_AUDIT.json")},
    )
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=96)
    parser.add_argument("--gtdbtk", type=Path, default=DEFAULT_GTDBTK)
    parser.add_argument(
        "--existing-identify-root",
        type=Path,
        help="Reuse only a completed identify run whose gtdbtk.json is strictly validated",
    )
    args = parser.parse_args()
    if not 1 <= args.threads <= 240:
        raise ReplayError("--threads must be in [1,240]")
    root = args.out.resolve()
    if root.exists():
        raise ReplayError(f"write-once output already exists: {root}")
    root.mkdir(parents=True, exist_ok=False)
    try:
        return run(args, root)
    except BaseException as error:
        complete = root / "COMPLETE.json"
        if not complete.exists():
            write_json_exclusive(
                complete,
                {"schema": "stage3b-identify-marker-map-complete-v1", "status": "FAIL", "error": repr(error)},
            )
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReplayError as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
