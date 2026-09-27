"""Byte-for-byte scheduler regression tests; no change to DNA2bit semantics."""
import argparse
import gzip
from pathlib import Path
import random
import struct
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


def reference_bits17(sequences):
    """Small independent integer oracle for the archived wyhash-17 contract."""
    mask = (1 << 64) - 1
    secret0, secret1 = 0xa0761d6478bd642f, 0xe7037ed1a0b428db

    def mix(a, b):
        product = a * b
        return (product & mask) ^ (product >> 64)

    seed0 = mix(secret0, secret1)

    def hash17(window):
        seed = mix(int.from_bytes(window[:8], 'little') ^ secret1,
                   int.from_bytes(window[8:16], 'little') ^ seed0)
        product = ((int.from_bytes(window[1:9], 'little') ^ secret1) *
                   (int.from_bytes(window[9:17], 'little') ^ seed))
        return mix((product & mask) ^ secret0 ^ 17,
                   (product >> 64) ^ secret1)

    counts = [0] * 55296
    complement = {ord('A'): ord('T'), ord('T'): ord('A'),
                  ord('G'): ord('C'), ord('C'): ord('G')}
    for sequence in sequences:
        if len(sequence) <= 17:
            continue
        reverse = bytes(complement.get(c, ord('N')) for c in sequence[-2::-1]) + b'\n'
        for strand in (sequence, reverse):
            for i in range(len(strand) - 17):
                value = hash17(strand[i:i+17])
                counts[(value & ((1 << 63) - 1)) % len(counts)] += 2 * (value >> 63) - 1
    packed = bytearray()
    for i in range(0, len(counts), 64):
        word = 0
        for count in counts[i:i+64]:
            word = (word << 1) | int(count >= 0)
        packed.extend(struct.pack('<Q', word))
    return bytes(packed)


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

    def test_gzip_container_variants_preserve_all_records(self):
        # The optional bounded decoder must fall back for concatenated members
        # and must not discard the second half of a FASTA/FASTQ stream.
        rng = random.Random(20260928)
        seqs = ["".join(rng.choice("ACGTNacgt") for _ in range(n))
                for n in (16, 17, 18, 65, 4097, 131073)]
        first = "".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs[:3])).encode()
        second = "".join(f">s{i+3}\n{s}\n" for i, s in enumerate(seqs[3:])).encode()
        plain = self.root / "reference.fa"; plain.write_bytes(first + second)
        expected = self.root / "expected.bit"
        self.run_tool([plain, expected], binary=LEGACY)
        variants = {
            "single.fa.gz": gzip.compress(first + second),
            "multi.fa.gz": gzip.compress(first) + gzip.compress(second),
            "padded.fa.gz": gzip.compress(first + second) + b"\0" * 64,
            "uncompressed.fa.gz": first + second,
            "records.fastq.gz": gzip.compress("".join(
                f"@s{i}\r\n{s}\r\n+\r\n{'I' * len(s)}\r\n"
                for i, s in enumerate(seqs)).encode()),
        }
        for name, content in variants.items():
            with self.subTest(container=name):
                source = self.root / name; source.write_bytes(content)
                output = self.root / (name + ".bit")
                self.run_tool([source, output])
                self.assertEqual(output.read_bytes(), expected.read_bytes())

    def test_integer_oracle_preserves_hash_and_legacy_window_rules(self):
        rng = random.Random(20260928)
        sequences = [bytes(rng.choice(b'ACGTNacgtRYSWKMBDHV.-') for _ in range(n))
                     for n in (0, 1, 16, 17, 18, 19, 31, 32, 65, 1001, 4097)]
        sequences.extend([b'A' * 10000, b'C' * 10000, b'ACGT' * 2048])
        expected = reference_bits17(sequences)
        source = self.root / 'oracle.fna.gz'
        source.write_bytes(gzip.compress(b''.join(b'>s\n' + s + b'\n' for s in sequences)))
        output = self.root / 'oracle.bit'
        self.run_tool([source, output])
        self.assertEqual(output.read_bytes(), expected)

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
