#!/usr/bin/env python3
"""Parameter-frozen signed Leiden backend for SAG Stage 3B.

Scientific implementation is intentionally delegated to the official Python
packages ``igraph`` and ``leidenalg``.  The C++ executable validates and builds
the two evidence layers; this backend validates them again, performs the exact
six-run sweep, and writes the selected size>=10 membership once.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import multiprocessing
import os
import tempfile
from collections import defaultdict
from pathlib import Path

SEED = 20260811
PARAMETERS = ((1.0, 1.0), (1.0, 3.0), (1.0, 10.0),
              (2.0, 3.0), (2.0, 10.0), (2.0, 30.0))

# Populated once in each spawned parameter worker.  Process isolation keeps the
# native igraph/leidenalg optimiser (and its seeded RNG) independent for every
# concurrent run.  The parent remains solely responsible for score/tie order.
_SIGNED_WORKER_INPUTS = None


class ContractError(RuntimeError):
    pass


def rows(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ContractError(f"missing TSV header: {path}")
        yield from reader


def finite(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"invalid {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise ContractError(f"non-finite {label}: {value!r}")
    return result


def pair(a: str, b: str) -> tuple[str, str]:
    if not a or not b or a == b:
        raise ContractError(f"invalid SAG pair: {a!r}, {b!r}")
    return (a, b) if a < b else (b, a)


def load_inputs(node_path: Path, positive_path: Path, marker_path: Path,
                negative_path: Path):
    nodes, positive = load_nodes_positive(node_path, positive_path)
    node_set = set(nodes)

    marker_pairs = {}
    for row in rows(marker_path):
        key = pair(row.get("SAG1", ""), row.get("SAG2", ""))
        if not set(key) <= node_set or key in marker_pairs:
            raise ContractError(f"unknown or duplicate marker pair: {key}")
        count = int(row.get("n_markers", ""))
        pid = finite(row.get("mean_pident", ""), "mean_pident")
        if count < 1 or not 0 <= pid <= 100:
            raise ContractError(f"invalid marker evidence: {key}")
        marker_pairs[key] = (count, pid)

    negative = {}
    for row in rows(negative_path):
        key = pair(row.get("SAG1", ""), row.get("SAG2", ""))
        if not set(key) <= node_set or key in negative:
            raise ContractError(f"unknown or duplicate negative pair: {key}")
        count = int(row.get("n_markers", ""))
        pid = finite(row.get("mean_pident", ""), "mean_pident")
        weight = finite(row.get("weight", ""), "negative weight")
        expected = (97.0 - pid) / 97.0
        if count < 3 or pid >= 97.0 or weight <= 0 or not math.isclose(
                weight, expected, rel_tol=1e-8, abs_tol=1e-10):
            raise ContractError(f"invalid negative-edge contract: {key}")
        if marker_pairs.get(key) != (count, pid):
            raise ContractError(f"negative edge disagrees with marker table: {key}")
        negative[key] = weight
    expected_negative = {key for key, (n, pid) in marker_pairs.items()
                         if n >= 3 and pid < 97.0}
    if set(negative) != expected_negative:
        raise ContractError("negative edge table is not the exact n>=3, mean<97 subset")
    return nodes, positive, marker_pairs, negative


def load_nodes_positive(node_path: Path, positive_path: Path):
    node_rows = list(rows(node_path))
    if not node_rows or "sag_id" not in node_rows[0]:
        raise ContractError("nodes TSV requires sag_id and at least one row")
    # Preserve the audited node-manifest order for the historical positive-only
    # find_partition call.  The C++ writer currently emits a deterministic order.
    nodes = [r["sag_id"] for r in node_rows]
    if any(not n for n in nodes) or len(nodes) != len(set(nodes)):
        raise ContractError("nodes must be nonempty and unique")
    node_set = set(nodes)

    positive = {}
    for row in rows(positive_path):
        key = pair(row.get("SAG1", ""), row.get("SAG2", ""))
        if not set(key) <= node_set or key in positive:
            raise ContractError(f"unknown or duplicate positive pair: {key}")
        ani = finite(row.get("ANI", ""), "ANI")
        afr = finite(row.get("AF_ref", ""), "AF_ref")
        afq = finite(row.get("AF_query", ""), "AF_query")
        gcd = finite(row.get("gc_diff", ""), "gc_diff")
        weight = finite(row.get("weight", ""), "positive weight")
        expected = (ani - 95.0) / 5.0 * (max(afr, afq) / 100.0) + 0.01
        if ani < 95.0 or gcd > 2.0 or weight <= 0 or not math.isclose(
                weight, expected, rel_tol=1e-8, abs_tol=1e-10):
            raise ContractError(f"invalid positive-edge contract: {key}")
        positive[key] = weight

    return nodes, positive


def positive_precluster(nodes, positive):
    """Exact marker-target precluster from pipeline_global.py:70-80."""
    if not positive:
        return [{node} for node in nodes]
    try:
        import igraph as ig
        import leidenalg as la
    except ImportError as exc:
        raise ContractError(
            "official python-igraph and leidenalg packages are required") from exc
    index = {node: i for i, node in enumerate(nodes)}
    graph = ig.Graph(n=len(nodes),
                     edges=[(index[a], index[b]) for a, b in positive])
    graph.es["weight"] = [positive[key] for key in positive]
    partition = la.find_partition(
        graph, la.RBConfigurationVertexPartition, weights="weight",
        resolution_parameter=1.0, seed=SEED, n_iterations=-1)
    return [{nodes[i] for i in community} for community in partition]


def signed_leiden(nodes, positive, negative, resolution, lam):
    try:
        import igraph as ig
        import leidenalg as la
    except ImportError as exc:
        raise ContractError(
            "official python-igraph and leidenalg packages are required") from exc

    index = {node: i for i, node in enumerate(nodes)}

    def layer(edges, scale=1.0):
        graph = ig.Graph(n=len(nodes),
                         edges=[(index[a], index[b]) for a, b in edges])
        graph.es["weight"] = [edges[key] * scale for key in edges]
        return graph

    gp = layer(positive)
    gn = layer(negative, lam)
    pp = la.RBConfigurationVertexPartition(
        gp, weights="weight", resolution_parameter=resolution)
    # Frozen correction from global_signed.py: CPM resolution 0 is a pure
    # internal-repulsion penalty when the multiplex layer weight is -1.
    pn = la.CPMVertexPartition(
        gn, weights="weight", resolution_parameter=0.0)
    optimiser = la.Optimiser()
    optimiser.set_rng_seed(SEED)
    optimiser.consider_comms = la.ALL_COMMS
    optimiser.optimise_partition_multiplex(
        [pp, pn], layer_weights=[1, -1], n_iterations=-1)
    return [{nodes[i] for i in community} for community in pp if community]


def _init_signed_worker(nodes, positive, negative):
    global _SIGNED_WORKER_INPUTS
    _SIGNED_WORKER_INPUTS = (nodes, positive, negative)


def _run_signed_parameter(index_and_parameter):
    if _SIGNED_WORKER_INPUTS is None:
        raise ContractError("signed-Leiden worker was not initialised")
    index, (resolution, lam) = index_and_parameter
    nodes, positive, negative = _SIGNED_WORKER_INPUTS
    parts = signed_leiden(nodes, positive, negative, resolution, lam)
    return index, parts


def evaluate_parameters(nodes, positive, negative, workers=1, parameters=PARAMETERS):
    """Evaluate the requested sweep, returning partitions in parameter order."""
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ContractError("workers must be a positive integer")
    if not positive:
        # With no attraction edge, the only valid maximum is singleton output;
        # avoiding RB's zero-total-weight division is deterministic.
        return [[{node} for node in nodes] for _ in parameters]
    if workers == 1:
        # Preserve the historical execution path exactly by default.
        return [signed_leiden(nodes, positive, negative, resolution, lam)
                for resolution, lam in parameters]

    worker_count = min(workers, len(parameters))
    indexed_parameters = tuple(enumerate(parameters))
    # Always spawn: unlike threads, each optimiser has an isolated native RNG;
    # unlike fork, this has the same semantics on Linux, macOS and Windows.
    context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=context,
            initializer=_init_signed_worker,
            initargs=(nodes, positive, negative)) as executor:
        completed = list(executor.map(
            _run_signed_parameter, indexed_parameters, chunksize=1))

    ordered = [None] * len(PARAMETERS)
    for index, parts in completed:
        if not 0 <= index < len(PARAMETERS) or ordered[index] is not None:
            raise ContractError("invalid or duplicate signed-Leiden worker result")
        ordered[index] = parts
    if any(parts is None for parts in ordered):
        raise ContractError("incomplete signed-Leiden parameter sweep")
    return ordered


def summarise(parts, marker_pairs, negative, minimum_size=10):
    label = {node: i for i, community in enumerate(parts)
             for node in community}
    ok = defaultdict(int)
    bad = defaultdict(int)
    worst = {}
    for (a, b), (count, pid) in marker_pairs.items():
        if count < 3 or label.get(a) != label.get(b):
            continue
        cluster = label.get(a)
        if cluster is None:
            continue
        if pid >= 97.0:
            ok[cluster] += 1
        else:
            bad[cluster] += 1
        worst[cluster] = min(worst.get(cluster, 101.0), pid)
    violated = sum(label.get(a) == label.get(b) for a, b in negative)
    detail = []
    for i, community in enumerate(parts):
        if len(community) < minimum_size:
            continue
        total = ok[i] + bad[i]
        detail.append({
            "partition_index": i,
            "size": len(community),
            # Authority rounds to three decimals before applying >=0.9.
            "purity": round(ok[i] / total, 3) if total else None,
            "marker_pairs_n_ge3": total,
            "min_marker_pident": worst.get(i),
        })
    detail.sort(key=lambda r: (-r["size"], r["partition_index"]))
    pure = [r for r in detail
            if r["purity"] is not None and r["purity"] >= 0.9]
    return detail, (len(pure), sum(r["size"] for r in pure)), violated


def write_once(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ContractError(f"refusing to overwrite: {path}")
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp",
                                     dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        # O_EXCL above and the same-filesystem replace make partial output
        # impossible.  The C++ orchestrator rejects unreceipted output.
        if path.exists():
            raise ContractError(f"refusing to overwrite: {path}")
        os.rename(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def membership_text(parts, minimum_size=10, prefix="G"):
    output = ["cluster\tSAG_id\tcluster_size"]
    emitted = 0
    ordered = sorted(parts, key=lambda p: (-len(p), tuple(sorted(p))))
    for community in ordered:
        if len(community) < minimum_size:
            continue
        emitted += 1
        for node in sorted(community):
            output.append(f"{prefix}{emitted:04d}\t{node}\t{len(community)}")
    return "\n".join(output) + "\n", emitted


def run(args):
    nodes, positive, marker_pairs, negative = load_inputs(
        args.nodes, args.positive, args.marker_pairs, args.negative)
    reports = []
    best = None
    workers = getattr(args, "workers", 1)
    parameters = PARAMETERS
    if getattr(args, "leiden_resolution", None) is not None:
        resolution = args.leiden_resolution
        parameters = tuple((resolution, lam) for lam in (1.0, 3.0, 10.0, 30.0))
    partitions = evaluate_parameters(
        nodes, positive, negative, workers=workers, parameters=parameters)
    for (resolution, lam), parts in zip(parameters, partitions):
        detail, score, violated = summarise(parts, marker_pairs, negative)
        tag = f"gs_r{resolution:.1f}_lam{int(lam)}"
        record = {
            "tag": tag,
            "res_pos": resolution,
            "lambda": lam,
            "clusters_ge10": len(detail),
            "pure_clusters_ge10": score[0],
            "sags_in_pure_clusters": score[1],
            "negative_edges": len(negative),
            "negative_edges_violated": violated,
            "clusters": detail,
        }
        reports.append(record)
        if best is None or score > best[0]:
            best = (score, tag, parts)
    assert best is not None
    membership, emitted = membership_text(best[2])
    report = {
        "schema": "sag-stage3b-signed-leiden-v1",
        "chosen": best[1],
        "chosen_score": {
            "pure_clusters_ge10": best[0][0],
            "sags_in_pure_clusters": best[0][1],
        },
        "clusters_ge10_emitted": emitted,
        "nodes": len(nodes),
        "positive_edges": len(positive),
        "negative_edges": len(negative),
        "scientific_contract": {
            "seed": SEED,
            "n_iterations": -1,
            "positive_partition": "RBConfigurationVertexPartition",
            "negative_partition": "CPMVertexPartition",
            "negative_resolution": 0.0,
            "multiplex_layer_weights": [1, -1],
            "selection": "lexicographic (pure clusters size>=10 and purity>=0.9, SAGs therein)",
            "parameter_order": [list(x) for x in parameters],
            "user_resolution": getattr(args, "leiden_resolution", None),
        },
        "sweep": reports,
    }
    write_once(args.membership, membership)
    write_once(args.report, json.dumps(report, indent=2, sort_keys=True) + "\n")


def run_precluster(args):
    nodes, positive = load_nodes_positive(args.nodes, args.positive)
    parts = positive_precluster(nodes, positive)
    membership, emitted = membership_text(parts, prefix="C")
    targeted = sum(len(p) for p in parts if len(p) >= 10)
    report = {
        "schema": "sag-stage3b-positive-precluster-v1",
        "method": "leidenalg.find_partition",
        "partition": "RBConfigurationVertexPartition",
        "resolution_parameter": 1.0,
        "seed": SEED,
        "n_iterations": -1,
        "nodes": len(nodes),
        "positive_edges": len(positive),
        "communities_ge10": emitted,
        "marker_target_sags": targeted,
    }
    write_once(args.membership, membership)
    write_once(args.report, json.dumps(report, indent=2, sort_keys=True) + "\n")


def self_test():
    parts = [set(f"S{i:02d}" for i in range(10)),
             set(f"T{i:02d}" for i in range(9))]
    marker = {(f"S{i:02d}", f"S{i+1:02d}"): (3, 99.0)
              for i in range(8)}
    marker[("S08", "S09")] = (3, 96.0)
    detail, score, violated = summarise(parts, marker, {("S08", "S09"): 1/97})
    assert detail[0]["purity"] == round(8/9, 3)
    assert score == (0, 0) and violated == 1
    text, count = membership_text(parts)
    assert count == 1 and "G0001\tS00\t10" in text and "T00" not in text
    pre_text, pre_count = membership_text(
        positive_precluster(["A", "B"], {}), prefix="C")
    assert pre_count == 0 and pre_text == "cluster\tSAG_id\tcluster_size\n"
    # Exercise the four TSV contracts and the frozen edge formula without
    # importing or running Leiden.
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "nodes.tsv").write_text(
            "sag_id\tassembly_fasta\tgc_pct\nA\ta.fna\t40\nB\tb.fna\t40\n")
        (root / "positive.tsv").write_text(
            "SAG1\tSAG2\tANI\tAF_ref\tAF_query\tgc_diff\tweight\n"
            "A\tB\t96\t20\t30\t0\t0.07\n")
        (root / "markers.tsv").write_text(
            "SAG1\tSAG2\tn_markers\tmean_pident\tmin_pident\n"
            "A\tB\t3\t96.000\t95.000\n")
        (root / "negative.tsv").write_text(
            "SAG1\tSAG2\tn_markers\tmean_pident\tweight\n"
            f"A\tB\t3\t96.000\t{1/97:.12f}\n")
        loaded = load_inputs(root / "nodes.tsv", root / "positive.tsv",
                             root / "markers.tsv", root / "negative.tsv")
        assert loaded[0] == ["A", "B"] and len(loaded[1]) == len(loaded[3]) == 1
        zero = root / "zero"
        zero.mkdir()
        (zero / "nodes.tsv").write_text(
            "sag_id\tassembly_fasta\tgc_pct\nA\ta.fna\t40\nB\tb.fna\t40\n")
        (zero / "positive.tsv").write_text(
            "SAG1\tSAG2\tANI\tAF_ref\tAF_query\tgc_diff\tweight\n")
        (zero / "markers.tsv").write_text(
            "SAG1\tSAG2\tn_markers\tmean_pident\tmin_pident\n")
        (zero / "negative.tsv").write_text(
            "SAG1\tSAG2\tn_markers\tmean_pident\tweight\n")
        zero_args = argparse.Namespace(
            nodes=zero / "nodes.tsv", positive=zero / "positive.tsv",
            marker_pairs=zero / "markers.tsv", negative=zero / "negative.tsv",
            membership=zero / "membership.tsv", report=zero / "report.json",
            workers=2)
        run(zero_args)
        assert zero_args.membership.read_text() == "cluster\tSAG_id\tcluster_size\n"
        zero_report = json.loads(zero_args.report.read_text())
        assert zero_report["nodes"] == 2 and zero_report["negative_edges"] == 0
    assert round(0.8996, 3) == 0.9  # authority scores the rounded purity
    assert PARAMETERS[-1] == (2.0, 30.0)
    print("PASS stage3b Leiden-backend synthetic self-test")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--positive", type=Path)
    parser.add_argument("--marker-pairs", type=Path)
    parser.add_argument("--negative", type=Path)
    parser.add_argument("--membership", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument(
        "--workers", type=int, default=1,
        help=("independent signed-Leiden parameter processes (default: 1; "
              "values above 6 are capped at the six frozen parameter sets)"))
    parser.add_argument(
        "--leiden-resolution", type=float,
        help="override positive-layer Leiden resolution; all other rules remain unchanged")
    parser.add_argument("--precluster", action="store_true",
                        help="positive-only Leiden marker-target stage")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be >=1")
    if args.leiden_resolution is not None and (not math.isfinite(args.leiden_resolution) or args.leiden_resolution <= 0):
        parser.error("--leiden-resolution must be a finite number greater than zero")
    if args.self_test:
        self_test()
        return
    if args.precluster:
        required = (args.nodes, args.positive, args.membership, args.report)
        if any(x is None for x in required):
            parser.error("precluster requires --nodes --positive --membership --report")
        run_precluster(args)
        return
    required = (args.nodes, args.positive, args.marker_pairs, args.negative,
                args.membership, args.report)
    if any(x is None for x in required):
        parser.error("all input/output arguments are required unless --self-test")
    run(args)


if __name__ == "__main__":
    try:
        cli()
    except (ContractError, OSError, ValueError) as exc:
        raise SystemExit(f"fatal: {exc}")
