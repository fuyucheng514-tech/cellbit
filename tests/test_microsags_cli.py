import importlib.machinery, importlib.util, tempfile, unittest
from unittest.mock import patch
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
            n=Namespace(inputs=[],file_list=None,reads1=str(a),reads2=str(b),r1_list=None,r2_list=None,input_type="auto")
            self.assertEqual(cli.rows(n)[1][0][0],"sample")
    def test_generic_scaffolds_uses_parent_sag_id(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"SAG_007"/"scaffolds.fasta"; p.parent.mkdir(); p.write_text(">x\nA\n")
            self.assertEqual(cli.sag_id(p),"SAG_007")
    def test_directory_input_discovers_sags(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/"SAGs"
            for name in ("A","B"):
                p=root/name/"scaffolds.fasta"; p.parent.mkdir(parents=True); p.write_text(">x\nA\n")
            a=Namespace(file_list=None,inputs=[str(root)])
            self.assertEqual([cli.sag_id(p) for p in cli.fasta_paths(a)],["A","B"])
    def test_auto_detects_and_pairs_fastq_directory(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/"reads"; root.mkdir()
            for name in ("SAG_A_R1.fastq.gz","SAG_A_R2.fastq.gz"):
                (root/name).write_text("@x\nA\n+\nI\n")
            a=Namespace(inputs=[str(root)],file_list=None,reads1=None,reads2=None,
                        r1_list=None,r2_list=None,input_type="auto")
            self.assertEqual(cli.rows(a)[1][0][0],"SAG_A")
    def test_explicit_type_rejects_suffix_conflict(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"SAG_A.fna"; p.write_text(">x\nA\n")
            a=Namespace(inputs=[str(p)],file_list=None,reads1=None,reads2=None,
                        r1_list=None,r2_list=None,input_type="reads")
            with self.assertRaises(cli.Error): cli.rows(a)
    def test_assemble_passes_flye_root(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); db=root/"db"; db.mkdir()
            for name in ("references.pack","references.tsv","genome_taxonomy.csv"):
                (db/name).write_text("x\n")
            args=Namespace(database=str(db),output=str(root/"out"),threads=2)
            prefix=root/"microsags"
            binaries={"dna2bit-sag-pipeline":str(prefix/"bin/dna2bit-sag-pipeline"),
                      "cpp-subass":str(prefix/"bin/cpp-subass"),
                      "flye":str(prefix/"bin/flye")}
            with patch.object(cli,"exe",side_effect=lambda name:binaries[name]), \
                 patch.object(cli,"run") as invoked:
                cli.pipeline(args,root/"manifest.tsv",False)
            command=invoked.call_args.args[0]
            self.assertEqual(command[command.index("--flye-root")+1],prefix.resolve())
if __name__ == "__main__": unittest.main()
