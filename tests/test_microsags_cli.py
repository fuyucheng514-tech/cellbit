import importlib.machinery, importlib.util, os, tempfile, unittest
from unittest.mock import patch
from argparse import Namespace
from pathlib import Path

path = Path(__file__).parents[1] / "tools" / "microsags"
loader = importlib.machinery.SourceFileLoader("microsags_cli", str(path))
spec = importlib.util.spec_from_loader(loader.name, loader)
cli = importlib.util.module_from_spec(spec); loader.exec_module(cli)

class Inputs(unittest.TestCase):
    def test_stage3b_database_paths_load_from_user_config(self):
        with tempfile.TemporaryDirectory() as d:
            config=Path(d)/"microsags/paths.env"; config.parent.mkdir()
            config.write_text("CHECKM2DB=/db/checkm2.dmnd\nGTDBTK_DATA_PATH=/db/gtdbtk\n")
            with patch.dict(os.environ,{"XDG_CONFIG_HOME":d},clear=False):
                values=cli.configured_paths()
            self.assertEqual(values["CHECKM2DB"],"/db/checkm2.dmnd")
            self.assertEqual(values["GTDBTK_DATA_PATH"],"/db/gtdbtk")
    def test_sketch_mode_has_standard_database_builder_options(self):
        args=cli.make_parser().parse_args(["sketch","references","-x","taxonomy.csv","-o","database","-t","8"])
        self.assertEqual(args.command,"sketch")
        self.assertEqual(args.inputs,["references"])
        self.assertEqual(args.taxonomy,"taxonomy.csv")
        self.assertEqual(args.output,"database")
        self.assertEqual(args.threads,8)
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
    def test_fastq_input_type_alias(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/"reads"; root.mkdir()
            for name in ("SAG_A_R1.fastq","SAG_A_R2.fastq"):
                (root/name).write_text("@x\nA\n+\nI\n")
            a=Namespace(inputs=[str(root)],file_list=None,reads1=None,reads2=None,
                        r1_list=None,r2_list=None,input_type="fastq")
            self.assertEqual(cli.rows(a)[1][0][0],"SAG_A")
    def test_singleton_fastq_input(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"SAG_A.fastq.gz"; p.write_text("@x\nA\n+\nI\n")
            a=Namespace(inputs=[str(p)],file_list=None,reads1=None,reads2=None,
                        r1_list=None,r2_list=None,input_type="singleton")
            header, data = cli.rows(a)
            self.assertEqual(header,["sag_id","singleton_fastq"])
            self.assertEqual(data[0][0],"SAG_A")
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
            annotations=root/"annotations"; (annotations/"02_dna2bit").mkdir(parents=True)
            (annotations/"02_dna2bit/labels.tsv").write_text("sag_id\treference\ttaxonomy\tspecies_group\n")
            (annotations/"03B_unclassified_pending.tsv").write_text("sag_id\tassembly_fasta\treason\n")
            args=Namespace(annotations=str(annotations),output=str(root/"out"),threads=2)
            prefix=root/"microsags"
            binaries={"dna2bit-sag-pipeline":str(prefix/"bin/dna2bit-sag-pipeline"),
                      "cpp-subass":str(prefix/"bin/cpp-subass"),
                      "flye":str(prefix/"bin/flye")}
            with patch.object(cli,"exe",side_effect=lambda name:binaries[name]), \
                 patch.object(cli,"run") as invoked:
                cli.pipeline(args,root/"manifest.tsv",False)
            command=invoked.call_args.args[0]
            self.assertEqual(command[command.index("--flye-root")+1],prefix.resolve())
    def test_annotation_publication_is_two_columns_and_complete(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); work=root/"work"; (work/"02_dna2bit").mkdir(parents=True)
            (work/"02_dna2bit/labels.tsv").write_text(
                "sag_id\treference\ttaxonomy\tspecies_group\nA\tGCF_1\td__Bacteria;s__Species_alpha\tSpecies_alpha\n")
            cli._publish_annotations(work,root/"result",[["A","a.fna"],["B","b.fna"]])
            self.assertEqual((root/"result/annotations.tsv").read_text(),
                             "sag_id\tspecies\nA\tSpecies_alpha\nB\tUNCLASSIFIED\n")
            self.assertEqual([p.name for p in (root/"result").iterdir()],["annotations.tsv"])
    def test_two_column_annotations_rebuild_internal_handoff(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); table=root/"annotations.tsv"
            table.write_text("sag_id\tspecies\nA\tSpecies_alpha\nB\tUNCLASSIFIED\n")
            cli._compat_annotations(table,[["A","a.fna"],["B","b.fna"]],root/"compat")
            self.assertIn("A\tNA\ts__Species_alpha\tSpecies_alpha",
                          (root/"compat/02_dna2bit/labels.tsv").read_text())
            self.assertIn("B\tb.fna\tdna2bit_rejected_or_no_hit",
                          (root/"compat/03B_unclassified_pending.tsv").read_text())
    def test_final_assembly_publication_is_minimal(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); work=root/"work"
            (work/"03A_subassemble/species/run").mkdir(parents=True)
            a=work/"03A_subassemble/species/run/assembly.fasta"; a.write_text(">a\nACGT\n")
            (work/"03A_subassemble/groups.tsv").write_text(
                f"species_group\tsag_count\tinput_fasta\tbin_fasta\nSpecies_alpha\t2\tin.fa\t{a}\n")
            (work/"stage3b/result/05_signed_leiden").mkdir(parents=True)
            (work/"stage3b/result/05_signed_leiden/chosen_membership.tsv").write_text(
                "cluster\tSAG_id\tcluster_size\nG0001\tB\t10\n")
            (work/"stage3b/result/06_subassemble/G0001/run").mkdir(parents=True)
            b=work/"stage3b/result/06_subassemble/G0001/run/assembly.fasta"; b.write_text(">b\nTGCA\n")
            (work/"stage3b/result/06_subassemble/clusters.tsv").write_text(
                f"cluster\tn_sags\tn_found\tinput_fasta\tbin_fasta\nG0001\t10\t10\tin.fa\t{b}\n")
            cli._publish_assemblies(work,root/"result")
            self.assertEqual(sorted(p.name for p in (root/"result").iterdir()),
                             ["stage3a","stage3a.tsv","stage3b","stage3b_clusters.tsv"])
            self.assertTrue((root/"result/stage3a/Species_alpha.fasta").is_file())
            self.assertTrue((root/"result/stage3b/G0001.fasta").is_file())
if __name__ == "__main__": unittest.main()
