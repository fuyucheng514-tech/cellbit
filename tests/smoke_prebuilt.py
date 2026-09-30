#!/usr/bin/env python3
"""Exercise a relocated Linux archive without an activated Conda environment."""

import csv
import gzip
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    program = Path(sys.argv[1]).resolve()
    assert program.is_file(), program
    clean_env = os.environ.copy()
    for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "LD_LIBRARY_PATH", "PYTHONPATH"):
        clean_env.pop(key, None)
    clean_env["PATH"] = "/usr/bin:/bin"

    def run(*args):
        command = [str(program), *map(str, args)]
        result = subprocess.run(command, env=clean_env, capture_output=True,
                                text=True)
        if result.returncode:
            raise RuntimeError(f"{command!r} exited {result.returncode}\n"
                               f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    run("--help")
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        refs = root / "references"
        refs.mkdir()
        contigs = root / "contigs"
        contigs.mkdir()
        reads = root / "reads"
        reads.mkdir()
        rng = random.Random(20260930)
        genomes = {}
        for number in (1, 2):
            accession = f"GCF_{number:09d}.1"
            sequence = "".join(rng.choices("ACGT", k=200000))
            genomes[number] = sequence
            (refs / f"{accession}_genomic.fna").write_text(
                f">{accession}\n{sequence}\n")
        (root / "taxonomy.csv").write_text(
            "GCF_000000001.1,d__Bacteria;p__Example;c__Example;o__Example;"
            "f__Example;g__Example;s__Example_A\n"
            "GCF_000000002.1,d__Bacteria;p__Example;c__Example;o__Example;"
            "f__Example;g__Example;s__Example_B\n")
        for number in (1, 2):
            sid = f"SAG_{number:04d}"
            seq = genomes[number]
            (contigs / f"{sid}.fna").write_text(f">{sid}\n{seq}\n")
            for mate in (1, 2):
                with gzip.open(reads / f"{sid}_R{mate}.fastq.gz", "wt") as out:
                    for index in range(400):
                        start = (index * 479 + mate * 37) % (len(seq) - 151)
                        read = seq[start:start + 150]
                        out.write(f"@{sid}_{index}/{mate}\n{read}\n+\n{'I'*150}\n")
        db = root / "database"
        run("sketch", refs, "-x", root / "taxonomy.csv", "-o", db, "-t", "2")
        assert (db / "COMPLETE.json").is_file()
        for mode, inputs in (("contigs", contigs), ("fastq", reads)):
            output = root / f"annotation_{mode}"
            run("annotate", "--input-type", mode, inputs,
                "-d", db, "-o", output, "-t", "2")
            assert {path.name for path in output.iterdir()} == {"annotations.tsv"}
            with (output / "annotations.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            assert rows == [
                {"sag_id": "SAG_0001", "species": "Example_A"},
                {"sag_id": "SAG_0002", "species": "Example_B"},
            ], (mode, rows)
    print("PASS relocated prebuilt sketch and annotation (contigs and FASTQ)")


if __name__ == "__main__":
    main()
