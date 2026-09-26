"""Byte-for-byte scheduler regression tests; no change to DNA2bit semantics."""
import argparse
import gzip
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

parser = argparse.ArgumentParser()
parser.add_argument("--binary", required=True)
parser.add_argument("--legacy-binary")
options, unittest_args = parser.parse_known_args()
BINARY = str(Path(options.binary).resolve())
LEGACY = str(Path(options.legacy_binary or options.binary).resolve())
sys.argv = [sys.argv[0], *unittest_args]


class BatchSketch(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_tool(self, args, ok=True, binary=BINARY):
        result = subprocess.run([binary, *map(str, args)], capture_output=True, text=True, timeout=45)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        return result

    def references(self):
        rng = random.Random(20260927)
        refs = []
        for i in range(24):
            seq = "".join(rng.choice("ACGTNacgt") for _ in range(2000 + i * 73))
            text = ">first\n" + seq + "\n>short\nACGT\n>second\n" + seq[::-1] + "\n"
            p = self.root / f"GCF_{i:09d}.1_genomic.fna{'.gz' if i % 2 else ''}"
            if i % 2:
                with gzip.open(p, "wt") as f:
                    f.write(text)
            else:
                p.write_text(text)
            refs.append(p)
        return refs

    def test_bytes_match_legacy_for_one_and_many_threads(self):
        refs = self.references()
        listing = self.root / "references.list"
        # Relative paths, spaces, and CRLF are accepted without changing names.
        spaced = self.root / "GCF_999999999.1_genomic with spaces.fna"
        spaced.write_bytes(refs[0].read_bytes()); refs.append(spaced)
        listing.write_bytes(("\r\n".join(p.name for p in refs) + "\r\n").encode())
        legacy = self.root / "legacy"; legacy.mkdir()
        for ref in refs:
            self.run_tool([ref, legacy / (ref.name + ".k.17.l.55296.bit")], binary=LEGACY)
        for workers in (1, 4, 50):
            out = self.root / f"batch{workers}"
            self.run_tool(["--batch-list", listing, "--output-dir", out, "--threads", workers])
            self.assertEqual(len(list(out.iterdir())), len(refs))
            for expected in legacy.iterdir():
                self.assertEqual(expected.stat().st_size, 6912)
                self.assertEqual((out / expected.name).read_bytes(), expected.read_bytes())

    def test_paired_input_legacy_interface_is_unchanged(self):
        refs = self.references()[:2]
        a, b = self.root / "old.bit", self.root / "new.bit"
        self.run_tool([*refs, a], binary=LEGACY)
        self.run_tool([*refs, b])
        self.assertEqual(a.read_bytes(), b.read_bytes())

    def test_duplicate_basename_is_rejected_before_writing(self):
        refs=[]
        for name in ("one", "two"):
            d=self.root/name; d.mkdir()
            p=d/"same.fna"; p.write_text(">a\nACGT\n"); refs.append(p)
        listing=self.root/"refs"; listing.write_text("\n".join(map(str,refs))+"\n")
        out=self.root/"bits"
        result=self.run_tool(["--batch-list",listing,"--output-dir",out],ok=False)
        self.assertIn("duplicate",result.stderr)
        self.assertFalse(out.exists())

    def test_empty_missing_and_invalid_options_fail(self):
        listing=self.root/"refs"; listing.write_text("")
        out=self.root/"bits"
        self.run_tool(["--batch-list",listing,"--output-dir",out],ok=False)
        listing.write_text("missing.fna\n")
        self.run_tool(["--batch-list",listing,"--output-dir",out],ok=False)
        for value in ("0","-1","two","1.2","99999999999999999999999"):
            self.run_tool(["--batch-list",listing,"--output-dir",out,"--threads",value],ok=False)
        self.run_tool(["--batch-list",listing,"--output-dir"],ok=False)
        self.run_tool(["--batch-list",listing,"--unknown",out],ok=False)
        self.run_tool(["--batch-list",listing,"--threads","1","--threads","2"],ok=False)
        self.assertFalse(out.exists())

    def test_existing_output_is_never_overwritten(self):
        ref=self.root/"one.fna"; ref.write_text(">a\nACGT\n")
        listing=self.root/"refs"; listing.write_text(str(ref)+"\n")
        out=self.root/"bits"; out.mkdir()
        sentinel=out/"one.fna.k.17.l.55296.bit"; sentinel.write_bytes(b"keep me")
        self.run_tool(["--batch-list",listing,"--output-dir",out],ok=False)
        self.assertEqual(sentinel.read_bytes(),b"keep me")


if __name__ == "__main__":
    unittest.main()
