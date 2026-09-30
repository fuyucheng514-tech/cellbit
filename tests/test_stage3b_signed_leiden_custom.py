"""Regression checks for an explicit Stage 3B Leiden resolution."""

import importlib.util
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch


BACKEND_PATH = Path(__file__).resolve().parents[1] / "python" / "stage3b_signed_leiden.py"
SPEC = importlib.util.spec_from_file_location("stage3b_signed_leiden", BACKEND_PATH)
backend = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(backend)


class FakeExecutor:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def map(self, function, indexed_parameters, chunksize):
        assert chunksize == 1
        return [(index, [{"A", "B"}]) for index, _ in indexed_parameters]


class CustomResolution(unittest.TestCase):
    def test_parallel_sweep_accepts_four_custom_parameters(self):
        parameters = tuple((0.18, lam) for lam in (1.0, 3.0, 10.0, 30.0))
        with patch.object(backend.concurrent.futures, "ProcessPoolExecutor", FakeExecutor):
            parts = backend.evaluate_parameters(
                ["A", "B"], {("A", "B"): 1.0}, {},
                workers=6, parameters=parameters)
        self.assertEqual(len(parts), 4)
        self.assertTrue(all(part == [{"A", "B"}] for part in parts))

    def test_explicit_resolution_keeps_exact_value_in_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "nodes.tsv").write_text(
                "sag_id\tassembly_fasta\tgc_pct\nA\ta.fna\t40\nB\tb.fna\t40\n")
            (root / "positive.tsv").write_text(
                "SAG1\tSAG2\tANI\tAF_ref\tAF_query\tgc_diff\tweight\n")
            (root / "markers.tsv").write_text(
                "SAG1\tSAG2\tn_markers\tmean_pident\tmin_pident\n")
            (root / "negative.tsv").write_text(
                "SAG1\tSAG2\tn_markers\tmean_pident\tweight\n")
            args = Namespace(
                nodes=root / "nodes.tsv", positive=root / "positive.tsv",
                marker_pairs=root / "markers.tsv", negative=root / "negative.tsv",
                membership=root / "membership.tsv", report=root / "report.json",
                workers=6, leiden_resolution=0.18)
            backend.run(args)
            report = json.loads(args.report.read_text())
            self.assertEqual(len(report["sweep"]), 4)
            self.assertEqual(report["chosen"], "gs_r0.18_lam1")
            self.assertEqual(report["scientific_contract"]["user_resolution"], 0.18)
            self.assertEqual(report["sweep"][-1]["tag"], "gs_r0.18_lam30")


if __name__ == "__main__":
    unittest.main()
