#!/usr/bin/env python3
"""Synthetic CLI test for cpp-subass' explicit no-overlap policy."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile


def fasta_records(data: bytes) -> dict[str, bytes]:
    records: dict[str, bytearray] = {}
    current: str | None = None
    for raw in data.splitlines():
        if not raw:
            continue
        if raw.startswith(b">"):
            current = raw[1:].split(None, 1)[0].decode("ascii")
            if not current or current in records:
                raise AssertionError("empty/duplicate synthetic FASTA ID")
            records[current] = bytearray()
        else:
            if current is None:
                raise AssertionError("sequence before header")
            records[current].extend(raw.strip().upper())
    return {name: bytes(sequence) for name, sequence in records.items()}


def invoke(binary: Path, root: Path, input_bytes: bytes, mode: str, policy: bool,
           phase: str = "all") -> subprocess.CompletedProcess[str]:
    root.mkdir(parents=True, exist_ok=True)
    reads = root / "input.fasta"
    reads.write_bytes(input_bytes)
    config = root / "asm_subasm.cfg"
    config.write_text("fixture\n", encoding="utf-8")
    command = [
        str(binary),
        "--reads", str(reads),
        "--out-dir", str(root / "out"),
        "--flye-modules", str(root.parent / "fake-flye-modules"),
        "--minimap2", str(root.parent / "unused-minimap2"),
        "--samtools", str(root.parent / "unused-samtools"),
        "--package-root", str(root),
        "--config", str(config),
        "--threads", "2",
        "--assemble-threads", "1",
        "--phase", phase,
    ]
    if policy:
        command.extend(["--no-overlap-policy", "passthrough"])
    environment = dict(os.environ)
    environment["FAKE_FLYE_MODE"] = mode
    return subprocess.run(command, text=True, capture_output=True, env=environment, check=False)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_no_overlap_cli.py /path/to/cpp-subass")
    binary = Path(sys.argv[1]).resolve(strict=True)
    single = b">SAG01__contig_1 description\nACGTNN\n>SAG01__contig_2\nGGCC\n"
    multiple = b">SAG01__a\nAAAA\n>SAG02__b\nCCCC\n>SAG03__c\nNNNN\n"

    with tempfile.TemporaryDirectory(prefix="cpp-subass-no-overlap-cli-") as temporary_text:
        temporary = Path(temporary_text)
        fake = temporary / "fake-flye-modules"
        fake.write_text(
            """#!/usr/bin/env python3
import os
from pathlib import Path
import sys

args = sys.argv[1:]
stage = args[0]
def value(flag):
    return Path(args[args.index(flag) + 1])

mode = os.environ.get("FAKE_FLYE_MODE", "no_overlap")
if stage == "assemble":
    draft, log = value("--out-asm"), value("--log")
    draft.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    if mode == "normal":
        draft.write_text(">disjointig_1\\nACGTACGT\\n", encoding="ascii")
        log.write_text("Assembled 1 disjointig\\n", encoding="ascii")
        raise SystemExit(0)
    draft.write_bytes(b"")
    if mode == "no_evidence":
        log.write_text("assembly ended quietly\\n", encoding="ascii")
    else:
        log.write_text("No overlaps found!\\n", encoding="ascii")
    raise SystemExit(17 if mode == "bad_rc" else 0)
if stage == "repeat":
    output = value("--out-dir")
    output.mkdir(parents=True, exist_ok=True)
    (output / "FAKE_REPEAT_REACHED").write_text("yes\\n", encoding="ascii")
    raise SystemExit(23)
raise SystemExit(29)
""",
            encoding="utf-8",
        )
        fake.chmod(fake.stat().st_mode | stat.S_IXUSR)

        strict = invoke(binary, temporary / "strict", single, "no_overlap", policy=False)
        assert strict.returncode != 0
        assert not (temporary / "strict/out/assembly.fasta").exists()

        for name, fasta, expected_records in (("single", single, 2), ("multiple", multiple, 3)):
            result = invoke(binary, temporary / name, fasta, "no_overlap", policy=True)
            assert result.returncode == 0, result.stderr
            assembly = temporary / name / "out/assembly.fasta"
            assert assembly.read_bytes() == fasta
            assert fasta_records(assembly.read_bytes()) == fasta_records(fasta)
            receipt = json.loads((temporary / name / "out/NO_OVERLAP_PASSTHROUGH.PASS.json").read_text())
            assert receipt["status"] == "PASS"
            assert receipt["mode"] == "byte_exact_input_passthrough"
            assert receipt["fasta_stats"]["records"] == expected_records
            assert receipt["output"]["sha256"] == hashlib.sha256(fasta).hexdigest()
            assert all(receipt["closure"].values())
            marker = (temporary / name / "out/CPP_FULL_PIPELINE_PASS").read_text()
            assert "mode=byte_exact_input_passthrough" in marker

        split_root = temporary / "split-no-overlap"
        split_assemble = invoke(binary, split_root, multiple, "no_overlap", policy=True, phase="assemble")
        assert split_assemble.returncode == 0, split_assemble.stderr
        assert (split_root / "out/00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE").is_file()
        split_finish = invoke(binary, split_root, multiple, "no_overlap", policy=True, phase="finish")
        assert split_finish.returncode == 0, split_finish.stderr
        assert (split_root / "out/40-polishing/CPP_FINISH_PHASE_COMPLETE").is_file()
        assert (split_root / "out/assembly.fasta").read_bytes() == multiple

        normal = invoke(binary, temporary / "normal", single, "normal", policy=True)
        assert normal.returncode != 0  # fake repeat exits 23 after proving normal continuation
        assert (temporary / "normal/out/20-repeat/FAKE_REPEAT_REACHED").is_file()
        assert not (temporary / "normal/out/NO_OVERLAP_PASSTHROUGH.PASS.json").exists()
        assert not (temporary / "normal/out/assembly.fasta").exists()

        for name, mode in (("no-evidence", "no_evidence"), ("bad-rc", "bad_rc")):
            result = invoke(binary, temporary / name, single, mode, policy=True)
            assert result.returncode != 0
            assert not (temporary / name / "out/assembly.fasta").exists()

    print("PASS cpp-subass CLI: strict, single/multi/split passthrough, normal continuation, fail-closed errors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
