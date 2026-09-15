#!/usr/bin/env python3
"""Write-once Stage-3B Flye replay and historical-equivalence audit.

This runner deliberately starts only after the raw-input and signed-Leiden
proof chains have passed.  It materialises exactly G0001..G0058, invokes the
historical Flye command, and compares both inputs and assemblies with the Lake
authority.  It never writes below Lake and has no resume/overwrite mode.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import os
from pathlib import Path
import subprocess
import sys
import time


LAKE = Path("/home/data/fyc/lake").resolve()
HISTORICAL_QF = (
    LAKE / "10_cellbit_unknown_sag_cluster_ppt/tables/quality_filtered_sags.tsv"
)
BROKEN_SAG = "TCACGCGGAGCGAGCAAGGCTCT"
EXPECTED_GROUPS = tuple(f"G{i:04d}" for i in range(1, 59))
EXPECTED_MEMBERS = 3759
EXPECTED_HISTORICAL_FOUND = 3758
EXPECTED_HISTORICAL_INPUT = {"contigs": 1_306_173, "bp": 728_088_150}
EXPECTED_HISTORICAL_OUTPUT = {"contigs": 3262, "bp": 29_260_880}
CHUNK_RECORDS = 100_000


class ReplayError(RuntimeError):
    pass


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
        raise ReplayError(f"refusing to overwrite {path}")
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists():
            raise ReplayError(f"refusing to overwrite {path}")
        os.rename(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_tsv_once(path: Path, header: list[str], rows: list[list[object]]) -> None:
    if path.exists():
        raise ReplayError(f"refusing to overwrite {path}")
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def load_json(path: Path, label: str) -> dict:
    if not path.is_file():
        raise ReplayError(f"missing {label}: {path}")
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ReplayError(f"invalid {label}: {path}: {error}") from error
    if not isinstance(value, dict):
        raise ReplayError(f"{label} is not a JSON object: {path}")
    return value


def all_checks_true(value: dict) -> bool:
    checks = value.get("checks")
    return isinstance(checks, dict) and bool(checks) and all(checks.values())


def load_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ReplayError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def resolve_record(row: dict[str, str], names: tuple[str, ...], label: str) -> str:
    for name in names:
        if name in row:
            return row[name]
    raise ReplayError(f"missing {label} column; accepted={names}")


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def fasta_records(path: Path):
    header = None
    pieces: list[bytes] = []
    with open_fasta(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    if not pieces:
                        raise ReplayError(f"empty FASTA record in {path}: {header!r}")
                    yield header, b"".join(pieces).upper()
                header = line[1:].strip()
                if not header:
                    raise ReplayError(f"empty FASTA header in {path}")
                pieces = []
            else:
                if header is None:
                    raise ReplayError(f"sequence before FASTA header in {path}")
                sequence = b"".join(line.split())
                if any(base not in b"ACGTURYSWKMBDHVNacgturyswkmbdhvn.-" for base in sequence):
                    raise ReplayError(f"invalid FASTA sequence byte in {path}")
                pieces.append(sequence)
    if header is not None:
        if not pieces:
            raise ReplayError(f"empty FASTA record in {path}: {header!r}")
        yield header, b"".join(pieces).upper()


def flush_digest_chunk(values: list[bytes], scratch: Path, number: int) -> Path:
    values.sort()
    target = scratch / f"chunk_{number:06d}.bin"
    with target.open("xb") as handle:
        handle.write(b"".join(values))
    return target


def read_digest(handle):
    value = handle.read(32)
    if not value:
        return None
    if len(value) != 32:
        raise ReplayError("truncated canonical digest chunk")
    return value


def canonical_sequence_digest(path: Path, scratch: Path) -> tuple[str, list[int]]:
    scratch.mkdir(parents=True)
    chunks: list[Path] = []
    values: list[bytes] = []
    lengths: list[int] = []
    for _, sequence in fasta_records(path):
        lengths.append(len(sequence))
        values.append(hashlib.sha256(sequence).digest())
        if len(values) >= CHUNK_RECORDS:
            chunks.append(flush_digest_chunk(values, scratch, len(chunks)))
            values = []
    if values:
        chunks.append(flush_digest_chunk(values, scratch, len(chunks)))
    if not lengths:
        raise ReplayError(f"empty FASTA: {path}")
    handles = [chunk.open("rb") for chunk in chunks]
    heap = []
    for index, handle in enumerate(handles):
        value = read_digest(handle)
        if value is not None:
            heapq.heappush(heap, (value, index))
    digest = hashlib.sha256()
    seen = 0
    try:
        while heap:
            value, index = heapq.heappop(heap)
            digest.update(value)
            seen += 1
            nxt = read_digest(handles[index])
            if nxt is not None:
                heapq.heappush(heap, (nxt, index))
    finally:
        for handle in handles:
            handle.close()
        for chunk in chunks:
            chunk.unlink()
        scratch.rmdir()
    if seen != len(lengths):
        raise ReplayError(f"canonical digest count mismatch for {path}")
    return digest.hexdigest(), lengths


def fasta_stats(path: Path, scratch: Path) -> dict[str, object]:
    canonical, lengths = canonical_sequence_digest(path, scratch)
    total = sum(lengths)
    threshold = (total + 1) // 2
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            n50 = length
            break
    return {
        "path": str(path.resolve()),
        "file_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "canonical_sequence_multiset_sha256": canonical,
        "contigs": len(lengths),
        "bp": total,
        "n50": n50,
        "max_contig": max(lengths),
    }


def load_membership(path: Path) -> dict[str, list[str]]:
    _, rows = load_tsv(path)
    groups: dict[str, list[str]] = {}
    seen: set[str] = set()
    declared: dict[str, int] = {}
    for row in rows:
        cluster = resolve_record(row, ("cluster",), "cluster")
        sag = resolve_record(row, ("SAG_id", "sag_id"), "SAG")
        size = int(resolve_record(row, ("cluster_size",), "cluster_size"))
        if sag in seen:
            raise ReplayError(f"SAG occurs more than once in membership: {sag}")
        seen.add(sag)
        groups.setdefault(cluster, []).append(sag)
        if cluster in declared and declared[cluster] != size:
            raise ReplayError(f"inconsistent declared size for {cluster}")
        declared[cluster] = size
    if tuple(sorted(groups)) != EXPECTED_GROUPS or len(seen) != EXPECTED_MEMBERS:
        raise ReplayError(
            f"membership is not exact 58/3759 closure: groups={len(groups)} members={len(seen)}"
        )
    for cluster in EXPECTED_GROUPS:
        groups[cluster].sort()
        if len(groups[cluster]) != declared[cluster] or declared[cluster] < 10:
            raise ReplayError(f"membership closure failed: {cluster}")
    if BROKEN_SAG not in groups["G0001"]:
        raise ReplayError(f"historical empty SAG is not in G0001: {BROKEN_SAG}")
    return groups


def bind_upstream(raw_root: Path, chain_root: Path) -> tuple[Path, Path, Path, dict, dict]:
    raw_audit_path = raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json"
    raw_audit = load_json(raw_audit_path, "raw-input audit")
    if raw_audit.get("status") != "PASS" or not all_checks_true(raw_audit):
        raise ReplayError("raw-input proof chain is not PASS")
    effective_qf = raw_root / "quality_filtered_sags_effective.tsv"
    content_manifest = raw_root / "EFFECTIVE_FASTA_CONTENT_MANIFEST.tsv"
    for key, path in (("effective_qf", effective_qf), ("content_manifest", content_manifest)):
        recorded = raw_audit.get("outputs", {}).get(key)
        recorded_hash = raw_audit.get("output_sha256", {}).get(str(path))
        if Path(recorded or "").resolve() != path.resolve() or recorded_hash != sha256(path):
            raise ReplayError(f"raw-input proof binding mismatch: {key}")

    chain_audit_path = chain_root / "CHAIN_AUDIT.json"
    chain_audit = load_json(chain_audit_path, "graph/marker/signed chain audit")
    reproduction_path = chain_root / "REPRODUCTION_AUDIT.json"
    reproduction = load_json(reproduction_path, "signed reproduction audit")
    if chain_audit.get("status") != "PASS" or reproduction.get("status") != "PASS":
        raise ReplayError("graph/marker/signed proof chain is not PASS")
    if not all_checks_true(reproduction):
        raise ReplayError("signed reproduction has a failed check")
    if chain_audit.get("raw_audit_sha256") != sha256(raw_audit_path):
        raise ReplayError("signed chain is not bound to this raw replay")
    if chain_audit.get("reproduction_audit_sha256") != sha256(reproduction_path):
        raise ReplayError("signed-chain reproduction audit hash mismatch")
    membership = chain_root / "global_signed/chosen_membership.tsv"
    if not membership.is_file():
        raise ReplayError(f"missing chosen membership: {membership}")
    return effective_qf, content_manifest, membership, raw_audit, chain_audit


def require_historical_qf_authority(raw_root: Path, historical_qf: Path) -> None:
    prereg = load_json(raw_root / "PREREGISTRATION.json", "raw-input preregistration")
    authorities = prereg.get("authorities")
    if not isinstance(authorities, list):
        raise ReplayError("raw-input preregistration has no authority list")
    observed = sha256(historical_qf)
    matches = [
        row for row in authorities
        if isinstance(row, dict)
        and Path(row.get("path", "")).resolve() == historical_qf
        and row.get("sha256") == observed
    ]
    if len(matches) != 1:
        raise ReplayError("historical QF is not uniquely bound by raw-input preregistration")


def load_sources(
    effective_qf: Path, content_manifest: Path, historical_qf: Path, mode: str
) -> tuple[dict[str, Path], dict[str, dict[str, object]]]:
    _, qf_rows = load_tsv(effective_qf)
    qf = {
        resolve_record(row, ("SAG_id", "sag_id"), "SAG"): Path(
            resolve_record(row, ("assembly_path", "assembly_fasta"), "assembly path")
        ).resolve()
        for row in qf_rows
    }
    _, manifest_rows = load_tsv(content_manifest)
    manifest = {
        resolve_record(row, ("SAG_id", "sag_id"), "SAG"): row
        for row in manifest_rows
    }
    if set(qf) != set(manifest):
        raise ReplayError("effective QF and content-manifest SAG sets differ")
    evidence: dict[str, dict[str, object]] = {}
    for sag, path in qf.items():
        row = manifest[sag]
        recorded = Path(resolve_record(row, ("effective_path",), "effective_path")).resolve()
        if recorded != path or not path.is_file():
            raise ReplayError(f"effective source path mismatch: {sag}")
        if int(resolve_record(row, ("file_bytes",), "file_bytes")) != path.stat().st_size:
            raise ReplayError(f"effective source size mismatch: {sag}")
        if resolve_record(row, ("sha256",), "sha256") != sha256(path):
            raise ReplayError(f"effective source hash mismatch: {sag}")
        evidence[sag] = {"path": str(path), "sha256": row["sha256"]}

    if mode == "historical":
        _, old_rows = load_tsv(historical_qf)
        old = {
            resolve_record(row, ("SAG_id", "sag_id"), "SAG"): Path(
                resolve_record(row, ("assembly_path", "assembly_fasta"), "assembly path")
            ).resolve()
            for row in old_rows
        }
        broken = old.get(BROKEN_SAG)
        if broken is None or not broken.is_file() or broken.stat().st_size != 0:
            raise ReplayError("historical branch requires the exact top-level zero-byte SAG")
        qf[BROKEN_SAG] = broken
        evidence[BROKEN_SAG] = {
            "path": str(broken), "sha256": sha256(broken), "historical_zero_contribution": True
        }
    return qf, evidence


def materialise_group(target: Path, members: list[str], sources: dict[str, Path]) -> tuple[int, list[str]]:
    if target.exists():
        raise ReplayError(f"refusing to overwrite {target}")
    found = 0
    skipped: list[str] = []
    with target.open("xb") as output:
        for sag in members:
            source = sources.get(sag)
            if source is None or not source.is_file():
                raise ReplayError(f"missing source binding for {sag}")
            if source.stat().st_size == 0:
                if sag != BROKEN_SAG:
                    raise ReplayError(f"unexpected empty SAG: {sag}")
                skipped.append(sag)
                continue
            records = 0
            with open_fasta(source) as handle:
                saw_header = False
                for raw in handle:
                    if raw.startswith(b">"):
                        original = raw[1:].strip()
                        if not original:
                            raise ReplayError(f"empty header in {source}")
                        output.write(b">" + sag.encode("ascii") + b"__" + original + b"\n")
                        saw_header = True
                        records += 1
                    else:
                        if raw.strip() and not saw_header:
                            raise ReplayError(f"sequence before header in {source}")
                        output.write(raw)
            if records == 0:
                raise ReplayError(f"non-empty source has no FASTA records: {source}")
            found += 1
        output.flush()
        os.fsync(output.fileno())
    return found, skipped


def resolve_authority(root: Path, cluster: str, kind: str) -> Path:
    if kind == "input":
        candidates = [
            root / f"{cluster}.fasta", root / f"{cluster}.fa", root / f"{cluster}.fna",
            root / cluster / "input_subassemblies.fasta", root / cluster / "input.fasta",
        ]
    else:
        candidates = [
            root / cluster / "assembly.fasta", root / cluster / "run/assembly.fasta",
            root / f"{cluster}.fasta", root / f"{cluster}.fa", root / f"{cluster}.fna",
        ]
    found = [path.resolve() for path in candidates if path.is_file()]
    unique = list(dict.fromkeys(found))
    if len(unique) != 1:
        raise ReplayError(
            f"expected one historical {kind} for {cluster}, found {len(unique)} under {root}"
        )
    return unique[0]


def flye_version(flye: Path) -> str:
    result = subprocess.run(
        [str(flye), "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    version = result.stdout.strip()
    if result.returncode != 0 or "2.9.6-b1802" not in version:
        raise ReplayError(f"expected Flye 2.9.6-b1802, observed: {version!r}")
    return version


def launch_flye(
    flye: Path, jobs: list[tuple[str, Path, Path, Path]], concurrency: int
) -> dict[str, dict[str, object]]:
    pending = list(jobs)
    running: dict[subprocess.Popen, tuple[str, Path, object, float, list[str]]] = {}
    result: dict[str, dict[str, object]] = {}
    try:
        while pending or running:
            while pending and len(running) < concurrency:
                cluster, input_path, run_dir, log_path = pending.pop(0)
                run_dir.parent.mkdir(parents=True, exist_ok=True)
                command = [
                    str(flye), "--subassemblies", str(input_path),
                    "--out-dir", str(run_dir), "--threads", "16",
                ]
                log = log_path.open("xb")
                log.write((json.dumps(command) + "\n").encode())
                log.flush()
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                running[process] = (cluster, run_dir, log, time.monotonic(), command)
            finished = None
            while finished is None:
                for process in running:
                    if process.poll() is not None:
                        finished = process
                        break
                if finished is None:
                    time.sleep(1)
            cluster, run_dir, log, started, command = running.pop(finished)
            elapsed = time.monotonic() - started
            log.write(f"\nreturncode={finished.returncode}\nwall_seconds={elapsed:.6f}\n".encode())
            log.flush()
            os.fsync(log.fileno())
            log.close()
            assembly = run_dir / "assembly.fasta"
            if finished.returncode != 0 or not assembly.is_file() or assembly.stat().st_size == 0:
                raise ReplayError(f"Flye failed for {cluster}: rc={finished.returncode}")
            result[cluster] = {
                "command": command, "returncode": finished.returncode,
                "wall_seconds": elapsed, "assembly": str(assembly.resolve()),
            }
    except BaseException:
        for process, (_, _, _, _, _) in running.items():
            if process.poll() is None:
                process.terminate()
        for process, (_, _, log, _, _) in running.items():
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            if not log.closed:
                log.write(f"\nterminated_after_peer_failure_rc={process.returncode}\n".encode())
                log.close()
        raise
    return result


def compare_stats(observed: dict, authority: dict) -> dict[str, object]:
    return {
        "authority": authority,
        "byte_exact": observed["sha256"] == authority["sha256"],
        "canonical_sequence_exact": (
            observed["canonical_sequence_multiset_sha256"]
            == authority["canonical_sequence_multiset_sha256"]
        ),
        "structure_exact": all(
            observed[key] == authority[key] for key in ("contigs", "bp", "n50", "max_contig")
        ),
        "delta": {
            key: observed[key] - authority[key] for key in ("contigs", "bp", "n50", "max_contig")
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--chain-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--flye", required=True, type=Path)
    parser.add_argument("--historical-input-root", required=True, type=Path)
    parser.add_argument("--historical-assembly-root", required=True, type=Path)
    parser.add_argument("--historical-qf", type=Path, default=HISTORICAL_QF)
    parser.add_argument("--mode", choices=("historical", "corrected"), default="historical")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--max-total-threads", type=int, default=192)
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    raw_root = args.raw_root.resolve()
    chain_root = args.chain_root.resolve()
    output = args.out.resolve()
    flye = args.flye.resolve()
    history_inputs = args.historical_input_root.resolve()
    history_assemblies = args.historical_assembly_root.resolve()
    historical_qf = args.historical_qf.resolve()
    if output.exists():
        raise ReplayError(f"write-once output exists: {output}")
    if str(output).startswith(str(LAKE) + os.sep):
        raise ReplayError("output must not be inside immutable Lake")
    if args.jobs < 1 or args.jobs * 16 > args.max_total_threads:
        raise ReplayError("jobs must be positive and jobs*16 must not exceed max-total-threads")
    for path, label in ((flye, "Flye"), (history_inputs, "historical inputs"),
                        (history_assemblies, "historical assemblies"),
                        (historical_qf, "historical QF")):
        if not path.exists():
            raise ReplayError(f"missing {label}: {path}")

    effective_qf, content_manifest, membership, raw_audit, chain_audit = bind_upstream(
        raw_root, chain_root
    )
    require_historical_qf_authority(raw_root, historical_qf)
    groups = load_membership(membership)
    sources, source_evidence = load_sources(
        effective_qf, content_manifest, historical_qf, args.mode
    )
    version = flye_version(flye)
    historical_input_paths = {
        cluster: resolve_authority(history_inputs, cluster, "input") for cluster in EXPECTED_GROUPS
    }
    historical_assembly_paths = {
        cluster: resolve_authority(history_assemblies, cluster, "assembly")
        for cluster in EXPECTED_GROUPS
    }

    output.mkdir(parents=True)
    config = {
        "schema": "stage3b-flye-replay-config-v1",
        "mode": args.mode,
        "groups": list(EXPECTED_GROUPS),
        "jobs": args.jobs,
        "threads_per_flye": 16,
        "max_total_threads": args.max_total_threads,
        "command_template": [
            str(flye), "--subassemblies", "{INPUT}", "--out-dir", "{OUT}", "--threads", "16"
        ],
        "historical_zero_contribution_sag": BROKEN_SAG if args.mode == "historical" else None,
    }
    prereg = {
        "schema": "stage3b-flye-replay-preregistration-v1",
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "config": config,
        "config_sha256": json_sha(config),
        "write_once_root": str(output),
        "tool": {"path": str(flye), "sha256": sha256(flye), "version": version},
        "upstream": {
            "raw_audit": {"path": str(raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json"),
                          "sha256": sha256(raw_root / "RAW_INPUT_REPRODUCTION_AUDIT.json")},
            "chain_audit": {"path": str(chain_root / "CHAIN_AUDIT.json"),
                            "sha256": sha256(chain_root / "CHAIN_AUDIT.json")},
            "membership": {"path": str(membership), "sha256": sha256(membership)},
            "effective_qf": {"path": str(effective_qf), "sha256": sha256(effective_qf)},
            "content_manifest": {"path": str(content_manifest), "sha256": sha256(content_manifest)},
        },
        "authorities": {
            "historical_qf": {"path": str(historical_qf), "sha256": sha256(historical_qf)},
            "historical_inputs_root": str(history_inputs),
            "historical_assemblies_root": str(history_assemblies),
        },
    }
    write_json_once(output / "PREREGISTRATION.json", prereg)
    write_json_once(output / "SOURCE_BINDING.json", {
        "schema": "stage3b-flye-source-binding-v1", "mode": args.mode,
        "source_count": len(source_evidence), "sources": source_evidence,
        "raw_status": raw_audit.get("status"), "chain_status": chain_audit.get("status"),
    })

    input_dir = output / "inputs"
    run_root = output / "flye"
    scratch_root = output / ".canonical_scratch"
    input_dir.mkdir()
    run_root.mkdir()
    scratch_root.mkdir()
    group_rows = []
    details: dict[str, dict[str, object]] = {}
    total_found = 0
    generated_jobs = []
    for cluster in EXPECTED_GROUPS:
        target = input_dir / f"{cluster}.fasta"
        found, skipped = materialise_group(target, groups[cluster], sources)
        total_found += found
        generated = fasta_stats(target, scratch_root / f"input_{cluster}_new")
        authority = fasta_stats(
            historical_input_paths[cluster], scratch_root / f"input_{cluster}_old"
        )
        input_comparison = compare_stats(generated, authority)
        details[cluster] = {
            "cluster": cluster, "n_sags": len(groups[cluster]), "n_found": found,
            "skipped_sags": skipped, "input": generated,
            "historical_input_comparison": input_comparison,
        }
        generated_jobs.append((cluster, target, run_root / cluster, run_root / f"{cluster}.log"))

    historical_expected_found = (
        EXPECTED_HISTORICAL_FOUND if args.mode == "historical" else EXPECTED_MEMBERS
    )
    if total_found != historical_expected_found:
        raise ReplayError(f"unexpected n_found: {total_found} != {historical_expected_found}")
    aggregate_input = {
        "contigs": sum(int(details[c]["input"]["contigs"]) for c in EXPECTED_GROUPS),
        "bp": sum(int(details[c]["input"]["bp"]) for c in EXPECTED_GROUPS),
    }
    if args.mode == "historical" and aggregate_input != EXPECTED_HISTORICAL_INPUT:
        raise ReplayError(f"historical input aggregate mismatch: {aggregate_input}")

    executions = launch_flye(flye, generated_jobs, args.jobs)
    for cluster in EXPECTED_GROUPS:
        observed = fasta_stats(
            Path(executions[cluster]["assembly"]), scratch_root / f"assembly_{cluster}_new"
        )
        authority = fasta_stats(
            historical_assembly_paths[cluster], scratch_root / f"assembly_{cluster}_old"
        )
        comparison = compare_stats(observed, authority)
        details[cluster]["flye"] = executions[cluster]
        details[cluster]["assembly"] = observed
        details[cluster]["historical_assembly_comparison"] = comparison
        group_rows.append([
            cluster, len(groups[cluster]), details[cluster]["n_found"],
            details[cluster]["input"]["contigs"], details[cluster]["input"]["bp"],
            details[cluster]["input"]["n50"], details[cluster]["input"]["sha256"],
            observed["contigs"], observed["bp"], observed["n50"], observed["sha256"],
            comparison["byte_exact"], comparison["canonical_sequence_exact"],
            comparison["structure_exact"], f"{executions[cluster]['wall_seconds']:.6f}",
        ])
    scratch_root.rmdir()

    aggregate_output = {
        "contigs": sum(int(details[c]["assembly"]["contigs"]) for c in EXPECTED_GROUPS),
        "bp": sum(int(details[c]["assembly"]["bp"]) for c in EXPECTED_GROUPS),
    }
    input_canonical_exact = all(
        bool(details[c]["historical_input_comparison"]["canonical_sequence_exact"])
        for c in EXPECTED_GROUPS
    )
    input_byte_exact = all(
        bool(details[c]["historical_input_comparison"]["byte_exact"]) for c in EXPECTED_GROUPS
    )
    assembly_canonical_exact = all(
        bool(details[c]["historical_assembly_comparison"]["canonical_sequence_exact"])
        for c in EXPECTED_GROUPS
    )
    assembly_structure_exact = all(
        bool(details[c]["historical_assembly_comparison"]["structure_exact"])
        for c in EXPECTED_GROUPS
    )
    checks = {
        "upstream_raw_pass": raw_audit.get("status") == "PASS",
        "upstream_graph_marker_signed_pass": chain_audit.get("status") == "PASS",
        "exact_58_groups": tuple(details) == EXPECTED_GROUPS,
        "membership_3759": sum(len(v) for v in groups.values()) == EXPECTED_MEMBERS,
        "expected_effective_sag_count": total_found == historical_expected_found,
        "historical_input_aggregate_exact": (
            args.mode != "historical" or aggregate_input == EXPECTED_HISTORICAL_INPUT
        ),
        "all_flye_commands_successful": len(executions) == 58,
        "historical_inputs_canonical_exact": (
            input_canonical_exact if args.mode == "historical" else True
        ),
        "historical_output_aggregate_exact": (
            args.mode != "historical" or aggregate_output == EXPECTED_HISTORICAL_OUTPUT
        ),
        "historical_assemblies_canonical_exact": (
            assembly_canonical_exact if args.mode == "historical" else True
        ),
    }
    audit = {
        "schema": "stage3b-flye-replay-audit-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "mode": args.mode,
        "note": (
            "Flye was historically run without --deterministic; byte equality is reported "
            "but is not an acceptance condition. In historical mode canonical sequence "
            "equality is required; corrected mode reports the expected scientific delta."
        ),
        "observed": {
            "groups": len(details), "membership_sags": EXPECTED_MEMBERS,
            "found_sags": total_found, "input": aggregate_input, "output": aggregate_output,
            "all_input_bytes_equal": input_byte_exact,
            "all_input_sequences_equal": input_canonical_exact,
            "all_assembly_sequences_equal": assembly_canonical_exact,
            "all_assembly_structures_equal": assembly_structure_exact,
        },
        "historical_expected": {
            "found_sags": EXPECTED_HISTORICAL_FOUND,
            "input": EXPECTED_HISTORICAL_INPUT,
            "output": EXPECTED_HISTORICAL_OUTPUT,
        },
        "groups": details,
    }
    write_tsv_once(
        output / "GROUP_MANIFEST.tsv",
        ["cluster", "n_sags", "n_found", "input_contigs", "input_bp", "input_n50",
         "input_sha256", "assembly_contigs", "assembly_bp", "assembly_n50",
         "assembly_sha256", "historical_byte_exact", "historical_sequence_exact",
         "historical_structure_exact", "wall_seconds"],
        group_rows,
    )
    write_json_once(output / "REPLAY_AUDIT.json", audit)
    complete = {
        "schema": "stage3b-flye-replay-complete-v1",
        "status": audit["status"],
        "preregistration_sha256": sha256(output / "PREREGISTRATION.json"),
        "group_manifest_sha256": sha256(output / "GROUP_MANIFEST.tsv"),
        "audit_sha256": sha256(output / "REPLAY_AUDIT.json"),
        "config_sha256": prereg["config_sha256"],
    }
    write_json_once(output / "COMPLETE.json", complete)
    print(json.dumps({"status": audit["status"], "observed": audit["observed"]}, indent=2))
    return 0 if audit["status"] == "PASS" else 2


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except ReplayError as error:
        output = args.out.resolve()
        if output.is_dir() and not (output / "FAILURE.json").exists():
            try:
                write_json_once(output / "FAILURE.json", {
                    "schema": "stage3b-flye-replay-failure-v1",
                    "status": "FAIL", "error": str(error),
                })
            except Exception:
                pass
        print(f"fatal: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
