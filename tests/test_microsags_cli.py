import importlib.machinery, importlib.util, tempfile, unittest
from argparse import Namespace
from pathlib import Path

path = Path(__file__).parents[1] / "tools" / "microsags"
loader = importlib.machinery.SourceFileLoader("microsags_cli", str(path))
spec = importlib.util.spec_from_loader(loader.name, loader)
cli = importlib.util.module_from_spec(spec); loader.exec_module(cli)

class Inputs(unittest.TestCase):
    def test_fasta_basename_ids(self):
        with tempfile.TemporaryDirectory() as d:
            ps=[Path(d)/"SAG_A.fna",Path(d)/"SAG-B.fa.gz"]
            for p in ps: p.write_text(">x\nACGT\n")
            a=Namespace(file_list=None,inputs=list(map(str,ps)))
            self.assertEqual([cli.sag_id(p) for p in cli.fasta_paths(a)],["SAG_A","SAG-B"])
    def test_list_relative_to_itself(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"one.fa"; p.write_text(">x\nA\n")
            l=Path(d)/"list.txt"; l.write_text("one.fa\n")
            self.assertEqual(cli.list_paths(l),[p.resolve()])
    def test_read_pair_id(self):
        with tempfile.TemporaryDirectory() as d:
            a,b=Path(d)/"sample_R1.fastq",Path(d)/"sample_R2.fastq"
            a.write_text("@x\nA\n+\nI\n"); b.write_text("@x\nT\n+\nI\n")
            n=Namespace(inputs=[],file_list=None,reads1=str(a),reads2=str(b),r1_list=None,r2_list=None)
            self.assertEqual(cli.rows(n)[1][0][0],"sample")
if __name__ == "__main__": unittest.main()
