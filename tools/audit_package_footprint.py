#!/usr/bin/env python3
"""Fail-closed byte audit for the small Cellbit SAG code package.

Reference databases, input data, and third-party runtime prefixes are measured
as external components but are never charged to (or followed from) the package
root.  The final archive can be checked independently by passing --archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

GIB = 1024**3
REFERENCE_BASENAMES = {
    "REFS.tsv",
    "SKETCHES.bin",
    "POSTINGS.bin",
    "POSTINGS_LOOKUP.bin",
    "uniref100.KO.1.dmnd",
}
REFERENCE_SUFFIXES = {".bit", ".dmnd", ".mmi", ".msh"}


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected NAME=PATH")
    name, raw = value.split("=", 1)
    if not name or not raw:
        raise argparse.ArgumentTypeError("expected non-empty NAME=PATH")
    return name, Path(os.path.expandvars(os.path.expanduser(raw)))


def parse_named_bytes(value: str) -> tuple[str, int]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected NAME=INTEGER_BYTES")
    name, raw = value.split("=", 1)
    try:
        size = int(raw)
    except ValueError as error:
        raise argparse.ArgumentTypeError("byte count must be an integer") from error
    if not name or size < 0:
        raise argparse.ArgumentTypeError("expected non-empty NAME and non-negative bytes")
    return name, size


def allocated_bytes(stat_result: os.stat_result) -> int:
    blocks = getattr(stat_result, "st_blocks", None)
    return int(blocks) * 512 if blocks is not None else int(stat_result.st_size)


def stable_sha256(path: Path) -> str:
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise RuntimeError(f"file changed while hashing: {path}")
    return digest.hexdigest()


def inventory_tree(root: Path, hash_files: bool) -> dict:
    totals = {
        "logical_file_bytes": 0,
        "allocated_file_bytes": 0,
        "unique_inode_logical_bytes": 0,
        "unique_inode_allocated_bytes": 0,
        "regular_files": 0,
        "directories": 0,
        "symlinks": 0,
        "other_entries": 0,
    }
    files: list[dict] = []
    links: list[dict] = []
    seen_inodes: set[tuple[int, int]] = set()

    if root.is_file():
        entries = [(root.parent, root.name)]
        display_root = root.parent
    else:
        entries = [(root, "")]
        display_root = root

    stack = entries[:]
    while stack:
        directory, only_name = stack.pop()
        scan = [next(e for e in os.scandir(directory) if e.name == only_name)] if only_name else list(os.scandir(directory))
        for entry in sorted(scan, key=lambda item: item.name, reverse=True):
            path = Path(entry.path)
            relative = path.relative_to(display_root).as_posix()
            if entry.is_symlink():
                totals["symlinks"] += 1
                target = os.readlink(path)
                links.append({"path": relative, "target": target})
            elif entry.is_dir(follow_symlinks=False):
                totals["directories"] += 1
                stack.append((path, ""))
            elif entry.is_file(follow_symlinks=False):
                stat_result = entry.stat(follow_symlinks=False)
                logical = int(stat_result.st_size)
                allocated = allocated_bytes(stat_result)
                totals["regular_files"] += 1
                totals["logical_file_bytes"] += logical
                totals["allocated_file_bytes"] += allocated
                inode = (int(stat_result.st_dev), int(stat_result.st_ino))
                if inode not in seen_inodes:
                    seen_inodes.add(inode)
                    totals["unique_inode_logical_bytes"] += logical
                    totals["unique_inode_allocated_bytes"] += allocated
                record = {"path": relative, "bytes": logical}
                if hash_files:
                    record["sha256"] = stable_sha256(path)
                files.append(record)
            else:
                totals["other_entries"] += 1
    files.sort(key=lambda item: item["path"])
    links.sort(key=lambda item: item["path"])
    return {**totals, "files": files if hash_files else None, "symlinks_detail": links}


def is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def inspect_component(name: str, component_class: str, configured: Path, known: int | None) -> dict:
    result = {
        "name": name,
        "class": component_class,
        "configured_path": str(configured),
        "expected_logical_bytes": known,
        "exists": configured.exists(),
    }
    if not configured.exists():
        return result
    resolved = configured.resolve(strict=True)
    result["resolved_path"] = str(resolved)
    result["inventory"] = inventory_tree(resolved, hash_files=False)
    return result


LDD_RESOLVED = re.compile(r"^\s*(\S+)\s+=>\s+(\S+)\s+\(")
LDD_DIRECT = re.compile(r"^\s*(/\S+)\s+\(")


def inspect_ldd(binary: Path, package_root: Path) -> dict:
    configured = binary
    if not binary.is_absolute():
        binary = package_root / binary
    result = {"configured_path": str(configured), "path": str(binary), "exists": binary.is_file()}
    if not binary.is_file():
        return result
    ldd = shutil.which("ldd")
    if not ldd:
        result["error"] = "ldd not found"
        return result
    process = subprocess.run([ldd, str(binary)], text=True, capture_output=True, check=False)
    # Loader addresses in raw ldd output vary between runs; retain only the
    # stable structured resolution below (and stderr when there is a failure).
    result.update({"returncode": process.returncode, "stderr": process.stderr})
    dependencies = []
    for line in process.stdout.splitlines():
        if "=> not found" in line:
            dependencies.append({"needed": line.split("=>", 1)[0].strip(), "resolved": None})
            continue
        match = LDD_RESOLVED.match(line)
        if match:
            needed, resolved_text = match.groups()
        else:
            direct = LDD_DIRECT.match(line)
            if not direct:
                continue
            needed = Path(direct.group(1)).name
            resolved_text = direct.group(1)
        resolved = Path(resolved_text).resolve(strict=False)
        record = {"needed": needed, "resolved": str(resolved), "inside_package": is_within(resolved, package_root)}
        if resolved.is_file():
            record["bytes"] = resolved.stat().st_size
        dependencies.append(record)
    result["dependencies"] = dependencies
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--archive", type=Path, help="final distributable archive; exact st_size is gated")
    parser.add_argument("--require-archive", action="store_true")
    parser.add_argument("--limit-bytes", type=int, default=GIB)
    parser.add_argument("--external-runtime", action="append", default=[], type=parse_named_path)
    parser.add_argument("--external-reference", action="append", default=[], type=parse_named_path)
    parser.add_argument("--external-input", action="append", default=[], type=parse_named_path)
    parser.add_argument("--known-bytes", action="append", default=[], type=parse_named_bytes)
    parser.add_argument("--ldd-binary", action="append", default=[], type=Path)
    parser.add_argument("--allow-missing-external", action="store_true")
    parser.add_argument("--json-out", type=Path, help="write report with exclusive create; default is stdout")
    args = parser.parse_args()

    violations: list[str] = []
    if args.limit_bytes <= 0:
        parser.error("--limit-bytes must be positive")
    if not args.package_root.is_dir():
        parser.error("--package-root must be an existing directory")
    package_root = args.package_root.resolve(strict=True)
    package = inventory_tree(package_root, hash_files=True)
    if package["logical_file_bytes"] > args.limit_bytes:
        violations.append("package root logical bytes exceed limit")

    forbidden_payload = []
    for item in package["files"]:
        packaged_path = Path(item["path"])
        if packaged_path.name in REFERENCE_BASENAMES or packaged_path.suffix.lower() in REFERENCE_SUFFIXES:
            forbidden_payload.append(item["path"])
            violations.append(f"reference-database payload found in package: {item['path']}")

    # External symlinks are forbidden in the release tree: a later tar command
    # using --dereference could otherwise ingest a multi-GiB database silently.
    for link in package["symlinks_detail"]:
        link_path = package_root / link["path"]
        resolved = link_path.resolve(strict=False)
        if not is_within(resolved, package_root):
            violations.append(f"package symlink escapes root: {link['path']} -> {link['target']}")

    known_bytes = dict(args.known_bytes)
    components = []
    names: set[str] = set()
    for component_class, values in (
        ("external_runtime", args.external_runtime),
        ("external_reference", args.external_reference),
        ("external_input", args.external_input),
    ):
        for name, path in values:
            if name in names:
                parser.error(f"duplicate external component name: {name}")
            names.add(name)
            component = inspect_component(name, component_class, path, known_bytes.get(name))
            components.append(component)
            if not component["exists"]:
                if not args.allow_missing_external:
                    violations.append(f"missing {component_class}: {name}={path}")
                continue
            resolved = Path(component["resolved_path"])
            if is_within(resolved, package_root):
                violations.append(f"external component is inside package root: {name}")
            expected = component.get("expected_logical_bytes")
            observed = component["inventory"]["logical_file_bytes"]
            if expected is not None and observed != expected:
                violations.append(f"known byte mismatch for {name}: expected {expected}, observed {observed}")

    archive = None
    if args.archive:
        archive_path = args.archive.resolve(strict=False)
        archive = {"path": str(archive_path), "exists": archive_path.is_file()}
        if archive_path.is_file():
            archive["bytes"] = archive_path.stat().st_size
            archive["sha256"] = stable_sha256(archive_path)
            if archive["bytes"] > args.limit_bytes:
                violations.append("final archive bytes exceed limit")
        else:
            violations.append("--archive is not an existing regular file")
    elif args.require_archive:
        violations.append("final archive was required but --archive was omitted")

    dynamic = [inspect_ldd(path, package_root) for path in args.ldd_binary]
    for item in dynamic:
        if not item.get("exists"):
            violations.append(f"ldd target missing: {item['path']}")
        elif item.get("error"):
            violations.append(f"ldd audit failed for {item['path']}: {item['error']}")
        elif item.get("returncode") != 0:
            violations.append(f"ldd returned nonzero for {item['path']}")
        elif any(dep.get("resolved") is None for dep in item.get("dependencies", [])):
            violations.append(f"unresolved dynamic dependency for {item['path']}")

    report = {
        "schema": "cellbit-sag-package-footprint-v1",
        "status": "PASS" if not violations else "FAIL",
        "limit_bytes": args.limit_bytes,
        "package_root": str(package_root),
        "package": package,
        "forbidden_reference_payload_matches": forbidden_payload,
        "archive": archive,
        "external_components": sorted(components, key=lambda item: (item["class"], item["name"])),
        "dynamic_link_audit": dynamic,
        "violations": violations,
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        with args.json_out.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(payload)
    else:
        sys.stdout.write(payload)
    return 0 if not violations else 2


if __name__ == "__main__":
    raise SystemExit(main())
