#!/usr/bin/env python3
"""Rebuild the historical Stage-3B entrance from Lake SAG FASTA files.

HISTORICAL ONLY: this program reproduces the superseded GTDB-reference
double-negative entrance and is not a current production runner.  The old
evidence logic and expected values are intentionally left unchanged.

This program is deliberately separate from the C++ ANI/AF implementation.  It
uses the exact historical skani executable and parameters to establish a
write-once reference run.  Lake is treated as immutable.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


LAKE = Path("/home/data/fyc/lake")
QF = LAKE / "10_cellbit_unknown_sag_cluster_ppt/tables/quality_filtered_sags.tsv"
QF_NEITHER = (
    LAKE
    / "19_cellbit_skani_af50_neither_cluster_flye_checkm2"
    / "tables/quality_filtered_neither_sags.tsv"
)
HIST_LIST = (
    LAKE / "10_cellbit_unknown_sag_cluster_ppt/inputs/quality_filtered_fastas.list"
)
HIST_TRIANGLE = (
    LAKE / "10_cellbit_unknown_sag_cluster_ppt/skani/skani_triangle_sparse.tsv"
)
HIST_BEST = (
    LAKE / "18_skani_gtdb232_sag_taxonomy/tables/skani_best_hit_taxonomy.tsv"
)
SKANI = Path("/home/data/fyc/biosoft/miniconda3/envs/assemble/bin/skani")
SKANI_DB = Path("/home/data/temp/release232/skani/database")
GTDB_META = Path("/home/data/temp/gtdb_metadata_r232/bac120_metadata_r232.tsv.gz")
BROKEN_SAG = "TCACGCGGAGCGAGCAAGGCTCT"
BROKEN_OVERRIDE = (
    LAKE / "07_spades_assemblies" / BROKEN_SAG / "K55/scaffolds.fasta"
)

EXPECTED_SHA256 = {
    str(QF): "2085788b66855fa11b5e4baa1497da1831f8ccd4acf14cc9f9c37224c462c168",
    str(QF_NEITHER): "1c9602c79cb24446c1bd50536ffda898e8fdbde4faeddcf59c1ca29e76e6c021",
    str(HIST_LIST): "64940abbf958fd3ebc6c5baef92a39f1d2a3525fbbab88307c7bdb07550e8505",
    str(HIST_TRIANGLE): "8eef3d81a5219aeb805db56f1cdbde03e46539e62959eee41c0eb42cf292d711",
    str(HIST_BEST): "a2e75b553a7c0b7a1ebfa7694e4113d074170c22034a0a00bc93274b9b753d7d",
    str(SKANI): "b1d20cb7170fe40a964526eadeb6fcdf61eefa748b88ed701ec3ff8dfbc07f5f",
    str(GTDB_META): "1650d3164666e5839c20ee15d82511909a5bbd5269035a7b564b9048dd777893",
}

RANK_ORDER = [
    ("species", "s__"),
    ("genus", "g__"),
    ("family", "f__"),
    ("order", "o__"),
    ("class", "c__"),
    ("phylum", "p__"),
    ("domain", "d__"),
]


class ReproductionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    if path.exists():
        raise ReproductionError(f"refusing to overwrite {path}")
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.rename(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ReproductionError(f"missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def fasta_stats_and_hash(path: Path) -> tuple[int, int, int, int, str]:
    digest = hashlib.sha256()
    contigs = 0
    total = 0
    maximum = 0
    gc = 0
    current = 0
    with path.open("rb") as raw:
        for line in raw:
            digest.update(line)
            if line.startswith(b">"):
                if current:
                    contigs += 1
                    total += current
                    maximum = max(maximum, current)
                current = 0
                continue
            sequence = line.strip().upper()
            current += len(sequence)
            gc += sequence.count(b"G") + sequence.count(b"C")
        if current:
            contigs += 1
            total += current
            maximum = max(maximum, current)
    return contigs, total, maximum, gc, digest.hexdigest()


def norm_accession(value: str) -> str:
    match = re.search(r"(GC[AF]_\d+\.\d+)", value or "")
    if match:
        return match.group(1)
    return (value or "").replace("RS_", "").replace("GB_", "")


def rank_from_taxonomy(taxonomy: str) -> str:
    value = (taxonomy or "").strip()
    if not value or value.startswith("Unclassified"):
        return "unclassified"
    parts = {item[:3]: item[3:] for item in value.split(";") if len(item) >= 3}
    for rank, prefix in RANK_ORDER:
        if parts.get(prefix):
            return rank
    return "classified_empty_rank"


def run_logged(command: list[str], log: Path) -> None:
    if log.exists():
        raise ReproductionError(f"refusing to overwrite log {log}")
    with log.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(command) + "\n")
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        handle.write(f"\nreturncode={result.returncode}\n")
        handle.flush()
        os.fsync(handle.fileno())
    if result.returncode:
        raise ReproductionError(f"command failed with rc={result.returncode}: {command}")


def sag_from_query_path(value: str) -> str:
    path = Path(value)
    if path.suffix in {".fna", ".fa", ".fasta"}:
        return path.stem
    return path.parent.name


def load_metadata() -> dict[str, str]:
    result: dict[str, str] = {}
    with gzip.open(GTDB_META, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"accession", "gtdb_taxonomy"} <= set(reader.fieldnames):
            raise ReproductionError("unexpected GTDB metadata schema")
        for row in reader:
            result[norm_accession(row["accession"])] = row["gtdb_taxonomy"]
    return result


def load_best(search: Path, metadata: dict[str, str]) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    with search.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Ref_file", "Query_file", "ANI", "Align_fraction_ref",
            "Align_fraction_query",
        }
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ReproductionError("unexpected skani search schema")
        for row in reader:
            sag = sag_from_query_path(row["Query_file"])
            accession = norm_accession(row["Ref_file"])
            ani = Decimal(row["ANI"])
            af_ref = Decimal(row["Align_fraction_ref"])
            af_query = Decimal(row["Align_fraction_query"])
            key = (ani, af_query, af_ref)
            if sag not in best or key > best[sag]["key"]:
                taxonomy = metadata.get(accession, "")
                best[sag] = {
                    "key": key,
                    "ref": accession,
                    "ani": ani,
                    "af_ref": af_ref,
                    "af_query": af_query,
                    "taxonomy": taxonomy,
                    "rank": rank_from_taxonomy(taxonomy) if taxonomy else "no_ref_taxonomy",
                }
    return best


def canonical_best_from_authority(qf_ids: set[str]) -> dict[str, tuple[object, ...]]:
    result = {}
    with HIST_BEST.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sag = row["SAG_id"]
            if sag not in qf_ids:
                continue
            result[sag] = (
                row["skani_ref"],
                Decimal(row["ANI"]),
                Decimal(row["AF_ref"]),
                Decimal(row["AF_query"]),
                row["skani_rank"],
                row["skani_taxonomy"],
            )
    return result


def canonical_best_observed(best: dict[str, dict[str, object]]) -> dict[str, tuple[object, ...]]:
    return {
        sag: (
            row["ref"], row["ani"], row["af_ref"], row["af_query"],
            row["rank"], row["taxonomy"],
        )
        for sag, row in best.items()
        if row["taxonomy"]
    }


def canonical_triangle(path: Path) -> list[tuple[object, ...]]:
    rows = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Ref_file", "Query_file", "ANI", "Align_fraction_ref",
            "Align_fraction_query", "Ref_name", "Query_name",
        }
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ReproductionError(f"unexpected skani triangle schema: {path}")
        for row in reader:
            value = (
                sag_from_query_path(row["Ref_file"]),
                sag_from_query_path(row["Query_file"]),
                Decimal(row["ANI"]),
                Decimal(row["Align_fraction_ref"]),
                Decimal(row["Align_fraction_query"]),
                row["Ref_name"],
                row["Query_name"],
            )
            key = value[:2]
            if key in seen:
                raise ReproductionError(f"duplicate triangle pair {key} in {path}")
            seen.add(key)
            rows.append(value)
    rows.sort()
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=96)
    args = parser.parse_args()
    if not 1 <= args.threads <= 240:
        raise ReproductionError("--threads must be in [1,240]")
    root = args.out.resolve()
    if root.exists():
        raise ReproductionError(f"write-once output exists: {root}")
    if str(root).startswith(str(LAKE.resolve()) + os.sep):
        raise ReproductionError("output must not be inside immutable Lake")

    for raw, expected in EXPECTED_SHA256.items():
        path = Path(raw)
        if not path.is_file() or sha256(path) != expected:
            raise ReproductionError(f"authority SHA mismatch: {path}")
    if not (SKANI_DB / "markers.bin").is_file():
        raise ReproductionError(f"missing skani database: {SKANI_DB}")
    version = subprocess.check_output([str(SKANI), "--version"], text=True).strip()
    if version != "skani 0.3.1":
        raise ReproductionError(f"unexpected skani version: {version}")

    root.mkdir(parents=True)
    inputs = root / "input_fastas"
    inputs.mkdir()
    _, qf_rows = read_tsv(QF)
    qf_ids = [row["SAG_id"] for row in qf_rows]
    if len(qf_ids) != 8785 or len(qf_ids) != len(set(qf_ids)):
        raise ReproductionError("QF SAG cardinality/uniqueness mismatch")
    if qf_ids != sorted(qf_ids):
        raise ReproductionError("QF row order is not the historical sorted order")
    historical_order = [Path(line.strip()).stem for line in HIST_LIST.read_text().splitlines() if line.strip()]
    if historical_order != qf_ids:
        raise ReproductionError("historical skani input order does not equal QF order")

    effective_rows = []
    content_manifest = []
    aggregate = {"contigs": 0, "total_bp": 0, "gc_bases": 0, "file_bytes": 0}
    for row in qf_rows:
        sag = row["SAG_id"]
        source = BROKEN_OVERRIDE if sag == BROKEN_SAG else Path(row["assembly_path"])
        if not source.is_file() or source.stat().st_size <= 0:
            raise ReproductionError(f"missing/empty source FASTA: {sag} {source}")
        contigs, total, maximum, gc, digest = fasta_stats_and_hash(source)
        observed_gc = f"{(gc / total * 100.0) if total else 0.0:.4f}"
        expected = (int(row["contigs"]), int(row["total_len"]), int(row["max_contig"]), row["gc_pct"])
        observed = (contigs, total, maximum, observed_gc)
        if observed != expected:
            raise ReproductionError(f"FASTA/QF mismatch for {sag}: {observed} != {expected}")
        link = inputs / f"{sag}.fna"
        os.symlink(source, link)
        copied = dict(row)
        copied["assembly_path"] = str(source)
        effective_rows.append(copied)
        content_manifest.append((sag, str(source), source.stat().st_size, digest))
        aggregate["contigs"] += contigs
        aggregate["total_bp"] += total
        aggregate["gc_bases"] += gc
        aggregate["file_bytes"] += source.stat().st_size

    effective_qf = root / "quality_filtered_sags_effective.tsv"
    with effective_qf.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qf_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(effective_rows)
    manifest_path = root / "EFFECTIVE_FASTA_CONTENT_MANIFEST.tsv"
    with manifest_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["SAG_id", "effective_path", "file_bytes", "sha256"])
        writer.writerows(content_manifest)
    query_list = root / "quality_filtered_fastas.list"
    query_list.write_text("".join(f"{inputs / (sag + '.fna')}\n" for sag in qf_ids))

    prereg = {
        "schema": "stage3b-lake-raw-input-reproduction-preregistration-v1",
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "write_once_root": str(root),
        "threads": args.threads,
        "skani_version": version,
        "skani_sha256": EXPECTED_SHA256[str(SKANI)],
        "skani_search_command": ["search", "-d", str(SKANI_DB), "--ql", str(query_list), "-t", str(args.threads)],
        "skani_triangle_command": ["triangle", "-l", str(query_list), "-t", str(args.threads), "-E", "--medium", "--min-af", "15"],
        "species_positive_rule": "taxonomy nonempty AND ANI>=95 AND AF_query>=50",
        "qf_rows": len(qf_rows),
        "fasta_aggregate": aggregate,
        "single_path_override": {BROKEN_SAG: str(BROKEN_OVERRIDE)},
        "authorities": [{"path": raw, "sha256": digest} for raw, digest in EXPECTED_SHA256.items()],
    }
    atomic_json(root / "PREREGISTRATION.json", prereg)

    search = root / "skani_search.tsv"
    triangle = root / "skani_triangle_sparse.tsv"
    run_logged(
        [str(SKANI), "search", "-d", str(SKANI_DB), "--ql", str(query_list),
         "-t", str(args.threads), "-o", str(search)],
        root / "skani_search.log",
    )
    run_logged(
        [str(SKANI), "triangle", "-l", str(query_list), "-o", str(triangle),
         "-t", str(args.threads), "-E", "--medium", "--min-af", "15"],
        root / "skani_triangle.log",
    )

    metadata = load_metadata()
    best = load_best(search, metadata)
    best_table = root / "skani_best_hit_taxonomy.tsv"
    with best_table.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["SAG_id", "skani_ref", "ANI", "AF_ref", "AF_query", "skani_rank", "skani_taxonomy"])
        for sag in qf_ids:
            if sag not in best or not best[sag]["taxonomy"]:
                continue
            row = best[sag]
            writer.writerow([sag, row["ref"], row["ani"], row["af_ref"], row["af_query"], row["rank"], row["taxonomy"]])

    positive = {
        sag for sag, row in best.items()
        if row["taxonomy"] and row["ani"] >= Decimal("95") and row["af_query"] >= Decimal("50")
    }
    neither_rows = [row for row in effective_rows if row["SAG_id"] not in positive]
    neither_path = root / "quality_filtered_neither_sags.tsv"
    with neither_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qf_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(neither_rows)

    _, authority_neither = read_tsv(QF_NEITHER)
    fields_without_path = [field for field in qf_rows[0] if field != "assembly_path"]
    observed_neither = sorted(tuple(row[field] for field in fields_without_path) for row in neither_rows)
    expected_neither = sorted(tuple(row[field] for field in fields_without_path) for row in authority_neither)
    best_equal = canonical_best_observed(best) == canonical_best_from_authority(set(qf_ids))
    produced_triangle_rows = canonical_triangle(triangle)
    triangle_equal = produced_triangle_rows == canonical_triangle(HIST_TRIANGLE)
    checks = {
        "all_8785_raw_fastas_match_qf": True,
        "one_recovered_path_override_only": sum(row[1] == str(BROKEN_OVERRIDE) for row in content_manifest) == 1,
        "skani_best_hits_equal_historical_qf_subset": best_equal,
        "qf_skani_positive_count_is_1377": len(positive.intersection(qf_ids)) == 1377,
        "neither_count_is_7408": len(neither_rows) == 7408,
        "neither_all_fields_except_recovered_path_equal_authority": observed_neither == expected_neither,
        "triangle_canonical_rows_equal_authority": triangle_equal,
    }
    audit = {
        "schema": "stage3b-lake-raw-input-reproduction-audit-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "observed": {
            "qf_sags": len(qf_rows),
            "skani_hits_with_taxonomy": len(canonical_best_observed(best)),
            "qf_skani_positive": len(positive.intersection(qf_ids)),
            "neither_sags": len(neither_rows),
            "triangle_rows": len(produced_triangle_rows),
            "fasta_aggregate": aggregate,
        },
        "outputs": {
            "effective_qf": str(effective_qf),
            "content_manifest": str(manifest_path),
            "query_list": str(query_list),
            "skani_search": str(search),
            "skani_best": str(best_table),
            "skani_triangle": str(triangle),
            "neither": str(neither_path),
        },
        "output_sha256": {
            str(path): sha256(path)
            for path in (effective_qf, manifest_path, query_list, search,
                         best_table, triangle, neither_path)
        },
    }
    atomic_json(root / "RAW_INPUT_REPRODUCTION_AUDIT.json", audit)
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["status"] == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReproductionError as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
