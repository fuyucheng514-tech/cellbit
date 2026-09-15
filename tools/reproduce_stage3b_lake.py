#!/usr/bin/env python3
"""Re-run and audit the historical Lake Stage-3B graph pipeline.

HISTORICAL ONLY: this preserves the superseded 7,408-node GTDB-reference
double-negative contract.  It is not a current production runner, and its
evidence logic is intentionally unchanged.

This tool deliberately keeps the historical evidence-replay path separate from
the production C++ ANI/AF replacement.  It imports the frozen Lake scripts,
redirects every output to a new write-once directory, and can execute:

1. the positive-only Leiden pre-clustering from the 7,408 quality-filtered SAGs;
2. bac120 marker all-vs-all BLASTN from the original nucleotide marker map;
3. the corrected six-point signed-Leiden sweep;
4. a canonical, order-independent comparison with the historical authority.

It never writes into /home/data/fyc/lake.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path


LAKE = Path("/home/data/fyc/lake")
HIST_ROOT = LAKE / "25_recluster_1530_newmethod"
QUALITY = (
    LAKE
    / "19_cellbit_skani_af50_neither_cluster_flye_checkm2"
    / "tables/quality_filtered_neither_sags.tsv"
)
SKANI_TRIANGLE = (
    LAKE
    / "10_cellbit_unknown_sag_cluster_ppt/skani/skani_triangle_sparse.tsv"
)
MARKER_MAP = (
    LAKE
    / "10_cellbit_unknown_sag_cluster_ppt/tables/bac120_marker_nt_map.tsv"
)
HIST_MARKER_PAIRS = HIST_ROOT / "marker_pairwise_blastn_global.tsv"
HIST_GLOBAL = HIST_ROOT / "global_signed"
HTML_REPORT = HIST_GLOBAL / "report/pipeline_report_unpeeled.html"
PIPELINE_SOURCE = HIST_ROOT / "pipeline_global.py"
MARKER_SOURCE = HIST_ROOT / "marker_blastn_mp.py"
SIGNED_SOURCE = HIST_ROOT / "global_signed.py"
BLAST_BIN = Path("/home/data/fyc/biosoft/miniconda3/envs/assemble/bin")
EXPECTED = {
    "quality_nodes": 7408,
    "positive_edges": 53250,
    "gc_filtered": 10458,
    "marker_targets": 3963,
    "marker_pairs": 624049,
    "marker_pairs_n_ge3": 114648,
    "negative_edges": 60271,
    "chosen": "gs_r2.0_lam30",
    "clusters_ge10": 58,
    "assigned_sags": 3759,
    "largest_cluster": 474,
    "negative_violations": 1208,
}
SWEEP = [
    ("gs_r1.0_lam1", 1.0, 1.0),
    ("gs_r1.0_lam3", 1.0, 3.0),
    ("gs_r1.0_lam10", 1.0, 10.0),
    ("gs_r2.0_lam3", 2.0, 3.0),
    ("gs_r2.0_lam10", 2.0, 10.0),
    ("gs_r2.0_lam30", 2.0, 30.0),
]
EXPECTED_DEPENDENCIES = {
    "networkx": "3.6.1",
    "igraph": "0.11.9",
    "leidenalg": "0.10.2",
}
EXPECTED_SELECTION_SCORES = {
    "gs_r1.0_lam1": (15, 1529),
    "gs_r1.0_lam3": (24, 1893),
    "gs_r1.0_lam10": (28, 2310),
    "gs_r2.0_lam3": (20, 1806),
    "gs_r2.0_lam10": (33, 2219),
    "gs_r2.0_lam30": (50, 3380),
}
EXPECTED_SHA256 = {
    str(QUALITY): "1c9602c79cb24446c1bd50536ffda898e8fdbde4faeddcf59c1ca29e76e6c021",
    str(SKANI_TRIANGLE): "8eef3d81a5219aeb805db56f1cdbde03e46539e62959eee41c0eb42cf292d711",
    str(MARKER_MAP): "c00137b8acbbb89c3c626fceb511a772111447137d768e28f6de2c8d53882ab2",
    str(HIST_MARKER_PAIRS): "6684e513e98c20e1abb1095ef1350f9f763bb19786d8d156d7ed02b85cd2a443",
    str(PIPELINE_SOURCE): "7114ea695996f7953bd4310c86d460603228555fd7e2559b43c51168244e7a45",
    str(MARKER_SOURCE): "08b2cd12d0e787ac9a71f61ca9e0076e25e832d14a640917c69c036802f111f3",
    str(SIGNED_SOURCE): "e60b0afd553e35cf5a2269fda18e7ae8f435a290f3c43db8710bdd9ddc668226",
    str(HTML_REPORT): "093157effeaf76e42f676ecc12a700022009bd592ad4ecff45a49fbf4352da53",
    str(HIST_ROOT / "global_pipeline/graph_report.json"): "cd3a898e25af85ddeec0a03ea64b629a4baf8723a5788c07ef3e1f9eaaee69a8",
    str(HIST_ROOT / "global_pipeline/cluster_membership.tsv"): "cab53c9abef82e42e3da22a7ad28ebdf5dd9ccfa1cf5d2671e14f4d80cb1e645",
    str(HIST_GLOBAL / "signed_report.json"): "f785b5ac093146ec92dd1e941270a7be83bd762b9c1dcbc239b5a6d10344ab9d",
    str(HIST_GLOBAL / "chosen_membership.tsv"): "f00f0fb578ddf786c9ebef8680b36638d765d51773160ed966f7462933a75fd9",
    str(HIST_GLOBAL / "gs_r1.0_lam1_membership.tsv"): "942e9ced4e5c06ead64f5b4194b7a0e2fdb4b50e41e70e0f020f0a2c9a2828ed",
    str(HIST_GLOBAL / "gs_r1.0_lam3_membership.tsv"): "ec00db75fa63bd47718d8e8997cf6057f81f61e1004be1bc902b903b9d781208",
    str(HIST_GLOBAL / "gs_r1.0_lam10_membership.tsv"): "cff1721826f889e664a5538a202255310fae51558e74e49ba2b085d1bb7d11eb",
    str(HIST_GLOBAL / "gs_r2.0_lam3_membership.tsv"): "1e9198addf0754d8e991264846658a268881c96dadcacf342203c8d89b2d4ca2",
    str(HIST_GLOBAL / "gs_r2.0_lam10_membership.tsv"): "59db95a165b38b3714b7aae1a3ac8a16fe977149e0b22624e218c4e87c77ad64",
    str(HIST_GLOBAL / "gs_r2.0_lam30_membership.tsv"): "f00f0fb578ddf786c9ebef8680b36638d765d51773160ed966f7462933a75fd9",
    str(BLAST_BIN / "makeblastdb"): "407bf59ec36ee6444104a90f28eab282e374a0dfe4c7ffee89a65d82dc51c4bf",
    str(BLAST_BIN / "blastn"): "0b88d4a00cb7fa579c203653151175b24c83d27fe240d420abe7b5261a3083d1",
}


class ReproductionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(block.count(b"\n") for block in iter(
            lambda: handle.read(8 * 1024 * 1024), b""))


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ReproductionError(f"refusing to overwrite: {path}")
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


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ReproductionError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    # ProcessPoolExecutor serialises historical marker workers by qualified
    # module name.  Register the redirected module before executing it so the
    # worker is importable/pickleable; the historical scientific code itself
    # is unchanged.
    previous = sys.modules.get(name)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        if previous is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = previous
        raise
    return module


def require_inputs(mode: str) -> None:
    required = [
        QUALITY,
        SKANI_TRIANGLE,
        MARKER_MAP,
        HIST_MARKER_PAIRS,
        PIPELINE_SOURCE,
        MARKER_SOURCE,
        SIGNED_SOURCE,
        HTML_REPORT,
        HIST_ROOT / "global_pipeline/graph_report.json",
        HIST_ROOT / "global_pipeline/cluster_membership.tsv",
        HIST_GLOBAL / "chosen_membership.tsv",
        HIST_GLOBAL / "signed_report.json",
        BLAST_BIN / "makeblastdb",
        BLAST_BIN / "blastn",
    ]
    required.extend(HIST_GLOBAL / f"{tag}_membership.tsv" for tag, _, _ in SWEEP)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ReproductionError("missing historical input(s): " + ", ".join(missing))
    mismatched = []
    for raw_path, expected in EXPECTED_SHA256.items():
        path = Path(raw_path)
        if mode == "frozen-evidence" and path in {
            MARKER_MAP, MARKER_SOURCE, BLAST_BIN / "makeblastdb", BLAST_BIN / "blastn"
        }:
            continue
        if path.is_file() and sha256(path) != expected:
            mismatched.append(str(path))
    if mismatched:
        raise ReproductionError(
            "historical authority SHA256 mismatch: " + ", ".join(mismatched)
        )
    observed_dependencies = {
        package: importlib.metadata.version(package)
        for package in EXPECTED_DEPENDENCIES
    }
    if observed_dependencies != EXPECTED_DEPENDENCIES:
        raise ReproductionError(
            "historical Python dependency mismatch: "
            f"expected {EXPECTED_DEPENDENCIES}, observed {observed_dependencies}"
        )
    if sys.version_info[:3] != (3, 12, 9):
        raise ReproductionError(
            "historical Python mismatch: expected 3.12.9, "
            f"observed {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        )


def preregister(root: Path, mode: str, marker_workers: int) -> None:
    if root.exists():
        raise ReproductionError(f"write-once output already exists: {root}")
    root.mkdir(parents=True)
    files = [
        Path(__file__).resolve(),
        QUALITY,
        SKANI_TRIANGLE,
        MARKER_MAP,
        HIST_MARKER_PAIRS,
        PIPELINE_SOURCE,
        MARKER_SOURCE,
        SIGNED_SOURCE,
        HTML_REPORT,
        HIST_GLOBAL / "chosen_membership.tsv",
        HIST_GLOBAL / "signed_report.json",
        HIST_ROOT / "global_pipeline/graph_report.json",
        HIST_ROOT / "global_pipeline/cluster_membership.tsv",
        BLAST_BIN / "makeblastdb",
        BLAST_BIN / "blastn",
    ]
    files.extend(HIST_GLOBAL / f"{tag}_membership.tsv" for tag, _, _ in SWEEP)
    inventory = []
    for path in files:
        record = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        if path.suffix in {".tsv", ".py", ".json"}:
            record["lines"] = line_count(path)
        inventory.append(record)
    atomic_json(
        root / "PREREGISTRATION.json",
        {
            "schema": "stage3b-lake-reproduction-preregistration-v1",
            "mode": mode,
            "write_once_root": str(root),
            "python": sys.version,
            "python_executable": sys.executable,
            "pythonhashseed": os.environ.get("PYTHONHASHSEED"),
            "dependency_versions": EXPECTED_DEPENDENCIES,
            "marker_workers": marker_workers,
            "threads_per_marker": 2,
            "seed": 20260811,
            "parameter_sweep": [[res, lam] for _, res, lam in SWEEP],
            "selection": (
                "lexicographic maximum of (number of size>=10 clusters with "
                "rounded marker purity>=0.9, SAG count in those clusters)"
            ),
            "expected_html_authority": EXPECTED,
            "inputs": inventory,
        },
    )


def run_positive(root: Path) -> None:
    target = root / "global_pipeline"
    if target.exists():
        raise ReproductionError(f"refusing to overwrite positive stage: {target}")
    module = load_module("lake_pipeline_global_reproduction", PIPELINE_SOURCE)
    module.QF = str(QUALITY)
    module.SKANI = str(SKANI_TRIANGLE)
    module.OUT = str(target)
    module.main()


def run_marker(root: Path, workers: int) -> None:
    output = root / "marker_pairwise_blastn_global.tsv"
    if output.exists():
        raise ReproductionError(f"refusing to overwrite marker result: {output}")
    membership = root / "global_pipeline/cluster_membership.tsv"
    if not membership.is_file():
        raise ReproductionError(f"positive-stage membership is missing: {membership}")
    module = load_module("lake_marker_blastn_reproduction", MARKER_SOURCE)
    module.BASE = str(root)
    module.GP = str(root / "global_pipeline")
    module.MAP = str(MARKER_MAP)
    module.BIN = str(BLAST_BIN)
    module.WORK = str(root / "global_pipeline/marker_blast")
    module.WORKERS = workers
    module.THREADS_PER_JOB = 2
    keep = module.target_sags()
    per = module.load_markers(keep)
    expected_jobs = sum(len(sequences) >= 2 for sequences in per.values())
    if len(keep) != EXPECTED["marker_targets"]:
        raise ReproductionError(
            f"expected {EXPECTED['marker_targets']} marker target SAGs, observed {len(keep)}"
        )
    if expected_jobs != 123:
        raise ReproductionError(
            f"expected 123 marker jobs from frozen inputs, observed {expected_jobs}"
        )
    module.main()
    residual = sorted(
        str(path.relative_to(Path(module.WORK)))
        for path in Path(module.WORK).rglob("*")
        if path.is_file()
    )
    if residual:
        raise ReproductionError(
            "marker worker left residual files, indicating a swallowed BLAST failure: "
            + ", ".join(residual[:20])
        )
    atomic_json(
        root / "MARKER_STAGE_AUDIT.json",
        {
            "schema": "stage3b-marker-stage-audit-v1",
            "status": "PASS",
            "target_sags": len(keep),
            "marker_jobs": expected_jobs,
            "workers": workers,
            "threads_per_job": 2,
            "residual_work_files": 0,
            "output_rows": line_count(output) - 1,
        },
    )


def copy_frozen_marker_evidence(root: Path) -> None:
    output = root / "marker_pairwise_blastn_global.tsv"
    if output.exists():
        raise ReproductionError(f"refusing to overwrite marker result: {output}")
    shutil.copyfile(HIST_MARKER_PAIRS, output)


def run_signed(root: Path) -> None:
    target = root / "global_signed"
    if target.exists():
        raise ReproductionError(f"refusing to overwrite signed stage: {target}")
    marker_pairs = root / "marker_pairwise_blastn_global.tsv"
    if not marker_pairs.is_file():
        raise ReproductionError(f"marker evidence is missing: {marker_pairs}")
    module = load_module("lake_global_signed_reproduction", SIGNED_SOURCE)
    module.QF = str(QUALITY)
    module.SKANI = str(SKANI_TRIANGLE)
    module.MK = str(marker_pairs)
    module.GP = str(root / "global_pipeline")
    module.OUT = str(target)
    module.main()


def table_rows(path: Path, fields: tuple[str, ...]) -> list[tuple[str, ...]]:
    values = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(fields) <= set(reader.fieldnames):
            raise ReproductionError(f"bad table schema: {path}")
        for row in reader:
            value = tuple(row[name] for name in fields)
            if value in seen:
                raise ReproductionError(f"duplicate canonical row in {path}: {value[:2]}")
            seen.add(value)
            values.append(value)
    values.sort()
    return values


def canonical_marker_rows(
    path: Path,
) -> list[tuple[str, str, int, Decimal, Decimal]]:
    values = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SAG1", "SAG2", "n_markers", "mean_pident", "min_pident"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ReproductionError(f"bad marker-pair schema: {path}")
        for row in reader:
            a, b = sorted((row["SAG1"], row["SAG2"]))
            if not a or not b or a == b:
                raise ReproductionError(f"invalid marker pair in {path}: {(a, b)}")
            key = (a, b)
            if key in seen:
                raise ReproductionError(f"duplicate undirected marker pair in {path}: {key}")
            seen.add(key)
            values.append(
                (
                    a,
                    b,
                    int(row["n_markers"]),
                    Decimal(row["mean_pident"]),
                    Decimal(row["min_pident"]),
                )
            )
    values.sort()
    return values


def canonical_clusters(path: Path) -> list[tuple[str, ...]]:
    clusters: dict[str, list[str]] = {}
    sizes: dict[str, int] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cluster = row["cluster"]
            sag = row["SAG_id"]
            declared = int(row["cluster_size"])
            clusters.setdefault(cluster, []).append(sag)
            previous = sizes.setdefault(cluster, declared)
            if previous != declared:
                raise ReproductionError(f"inconsistent cluster_size for {cluster}")
    result = []
    all_sags = set()
    for cluster, members in clusters.items():
        if len(members) != sizes[cluster] or len(members) != len(set(members)):
            raise ReproductionError(f"membership closure failed for {cluster}")
        overlap = all_sags.intersection(members)
        if overlap:
            raise ReproductionError(f"SAG assigned more than once: {min(overlap)}")
        all_sags.update(members)
        result.append(tuple(sorted(members)))
    result.sort(key=lambda members: (-len(members), members))
    return result


def membership_map(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    """Return SAG->cluster and cluster->size after enforcing table closure."""
    labels: dict[str, str] = {}
    sizes: dict[str, int] = {}
    counts: dict[str, int] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cluster", "SAG_id", "cluster_size"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ReproductionError(f"bad membership schema: {path}")
        for row in reader:
            cluster, sag = row["cluster"], row["SAG_id"]
            declared = int(row["cluster_size"])
            if sag in labels:
                raise ReproductionError(f"SAG assigned more than once in {path}: {sag}")
            labels[sag] = cluster
            old = sizes.setdefault(cluster, declared)
            if old != declared:
                raise ReproductionError(f"inconsistent cluster_size in {path}: {cluster}")
            counts[cluster] = counts.get(cluster, 0) + 1
    for cluster, declared in sizes.items():
        if counts.get(cluster) != declared or declared < 10:
            raise ReproductionError(f"membership closure failed in {path}: {cluster}")
    return labels, sizes


def independent_selection_score(
    membership: Path, marker_pairs: Path
) -> tuple[int, int]:
    """Recompute the historical rounded-purity objective independently."""
    labels, sizes = membership_map(membership)
    ok: dict[str, int] = {}
    bad: dict[str, int] = {}
    with marker_pairs.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SAG1", "SAG2", "n_markers", "mean_pident"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ReproductionError(f"bad marker-pair schema: {marker_pairs}")
        for row in reader:
            if int(row["n_markers"]) < 3:
                continue
            cluster = labels.get(row["SAG1"])
            if cluster is None or cluster != labels.get(row["SAG2"]):
                continue
            target = ok if float(row["mean_pident"]) >= 97.0 else bad
            target[cluster] = target.get(cluster, 0) + 1
    pure = []
    for cluster, size in sizes.items():
        total = ok.get(cluster, 0) + bad.get(cluster, 0)
        purity = round(ok.get(cluster, 0) / total, 3) if total else None
        if purity is not None and purity >= 0.9:
            pure.append(size)
    return len(pure), sum(pure)


def audit(root: Path, marker_was_recomputed: bool) -> dict[str, object]:
    produced_report_path = root / "global_signed/signed_report.json"
    produced_membership = root / "global_signed/chosen_membership.tsv"
    if not produced_report_path.is_file() or not produced_membership.is_file():
        raise ReproductionError("signed outputs are incomplete")
    with produced_report_path.open() as handle:
        produced_report = json.load(handle)
    with (HIST_GLOBAL / "signed_report.json").open() as handle:
        authority_report = json.load(handle)
    with (root / "global_pipeline/graph_report.json").open() as handle:
        graph_report = json.load(handle)
    with (HIST_ROOT / "global_pipeline/graph_report.json").open() as handle:
        authority_graph_report = json.load(handle)

    produced_clusters = canonical_clusters(produced_membership)
    authority_clusters = canonical_clusters(HIST_GLOBAL / "chosen_membership.tsv")
    produced_positive = canonical_clusters(root / "global_pipeline/cluster_membership.tsv")
    authority_positive = canonical_clusters(
        HIST_ROOT / "global_pipeline/cluster_membership.tsv"
    )
    marker_fields = ("SAG1", "SAG2", "n_markers", "mean_pident", "min_pident")
    marker_equal = None
    if marker_was_recomputed:
        marker_equal = canonical_marker_rows(
            root / "marker_pairwise_blastn_global.tsv"
        ) == canonical_marker_rows(
            HIST_MARKER_PAIRS
        )

    chosen_sweep = next(
        row for row in produced_report["sweep"] if row["tag"] == produced_report["chosen"]
    )
    marker_target_rows = line_count(root / "global_pipeline/cluster_membership.tsv") - 1
    observed = {
        "quality_nodes": graph_report["n_sags"],
        "positive_edges": graph_report["edges"],
        "gc_filtered": graph_report["gc_filtered"],
        "marker_targets": marker_target_rows,
        "marker_pairs": line_count(root / "marker_pairwise_blastn_global.tsv") - 1,
        "marker_pairs_n_ge3": sum(
            row[2] >= 3
            for row in canonical_marker_rows(
                root / "marker_pairwise_blastn_global.tsv"
            )
        ),
        "negative_edges": chosen_sweep["neg_edges"],
        "chosen": produced_report["chosen"],
        "clusters_ge10": len(produced_clusters),
        "assigned_sags": sum(map(len, produced_clusters)),
        "largest_cluster": max(map(len, produced_clusters), default=0),
        "negative_violations": chosen_sweep["neg_violated"],
    }
    checks = {
        key: observed[key] == expected
        for key, expected in EXPECTED.items()
        if key in observed
    }
    expected_tags = [tag for tag, _, _ in SWEEP]
    produced_tags = [row.get("tag") for row in produced_report.get("sweep", [])]
    checks["sweep_tags_and_order_exact"] = produced_tags == expected_tags
    checks["complete_signed_report_equals_authority"] = (
        produced_report == authority_report
    )
    checks["complete_positive_graph_report_equals_authority"] = (
        graph_report == authority_graph_report
    )
    checks["canonical_positive_membership_equals_authority"] = (
        produced_positive == authority_positive
    )
    sweep_membership_checks = {}
    independent_scores = []
    for tag, res, lam in SWEEP:
        produced_path = root / f"global_signed/{tag}_membership.tsv"
        authority_path = HIST_GLOBAL / f"{tag}_membership.tsv"
        equal = (
            produced_path.is_file()
            and authority_path.is_file()
            and canonical_clusters(produced_path) == canonical_clusters(authority_path)
        )
        sweep_membership_checks[tag] = equal
        checks[f"canonical_membership_equals_authority:{tag}"] = equal
        if produced_path.is_file():
            independent_scores.append(
                {
                    "tag": tag,
                    "resolution": res,
                    "lambda": lam,
                    "score": list(independent_selection_score(
                        produced_path,
                        root / "marker_pairwise_blastn_global.tsv",
                    )),
                }
            )
    if len(independent_scores) == len(SWEEP):
        best_score = max(tuple(row["score"]) for row in independent_scores)
        independently_chosen = next(
            row["tag"]
            for row in independent_scores
            if tuple(row["score"]) == best_score
        )
    else:
        independently_chosen = None
    checks["chosen_equals_independent_lexicographic_selection"] = (
        produced_report.get("chosen") == independently_chosen
    )
    checks["independent_selection_scores_equal_authority"] = (
        len(independent_scores) == len(SWEEP)
        and all(
            tuple(row["score"]) == EXPECTED_SELECTION_SCORES[row["tag"]]
            for row in independent_scores
        )
    )
    checks["canonical_membership_equals_authority"] = produced_clusters == authority_clusters
    checks["chosen_report_equals_authority"] = produced_report["chosen"] == authority_report["chosen"]
    if marker_equal is not None:
        checks["canonical_marker_table_equals_authority"] = marker_equal
    result = {
        "schema": "stage3b-lake-reproduction-audit-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "observed": observed,
        "expected": EXPECTED,
        "checks": checks,
        "independent_selection": {
            "scores_in_frozen_order": independent_scores,
            "chosen": independently_chosen,
        },
        "sweep_membership_checks": sweep_membership_checks,
        "produced": {
            "signed_report": str(produced_report_path),
            "membership": str(produced_membership),
            "marker_pairs": str(root / "marker_pairwise_blastn_global.tsv"),
        },
        "authority": {
            "signed_report": str(HIST_GLOBAL / "signed_report.json"),
            "membership": str(HIST_GLOBAL / "chosen_membership.tsv"),
            "marker_pairs": str(HIST_MARKER_PAIRS),
        },
    }
    atomic_json(root / "REPRODUCTION_AUDIT.json", result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--mode", choices=("frozen-evidence", "recompute-markers"),
        default="recompute-markers",
        help=(
            "frozen-evidence validates the parameter/Leiden chain using the original "
            "marker-pair table; recompute-markers also regenerates marker BLASTN evidence"
        ),
    )
    parser.add_argument("--marker-workers", type=int, default=96)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 1 <= args.marker_workers <= 96:
        raise ReproductionError("--marker-workers must be in [1, 96]")
    require_inputs(args.mode)
    root = args.out.resolve()
    if str(root).startswith(str(LAKE.resolve()) + os.sep):
        raise ReproductionError("output must not be inside the immutable Lake authority")
    preregister(root, args.mode, args.marker_workers)
    run_positive(root)
    marker_was_recomputed = args.mode == "recompute-markers"
    if marker_was_recomputed:
        run_marker(root, args.marker_workers)
    else:
        copy_frozen_marker_evidence(root)
    run_signed(root)
    result = audit(root, marker_was_recomputed)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReproductionError as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2)
