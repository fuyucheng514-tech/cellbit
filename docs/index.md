# Microsags

**Species labels first. SAG assembly second.**

Microsags connects two tasks for single-amplified genomes (SAGs): identify the
species of each SAG, then assemble groups of SAGs from their contigs. You can
stop after annotation or continue to assembly.

## The workflow

1. **Annotate** FASTQ reads or per-SAG contigs. The result is
   `annotations.tsv`: one species call or `UNCLASSIFIED` for each SAG.
2. **Assemble** per-SAG contigs using that annotation table. Labelled SAGs enter
   species-based assembly; SAGs without a species label enter label-free
   clustering. The results are group membership tables and assembled FASTA
   files.

```bash
microsags annotate SAGs/ -d database -o annotation_output
microsags assemble SAG_contigs/ -a annotation_output/annotations.tsv -o assembly_output
```

Assembly requires contigs. If you annotated FASTQ reads, assemble each SAG
individually first, keeping the same SAG identifiers. Annotation itself does not
assemble reads.

**Start here:** [Install Microsags](install.md) ·
[Run your own SAGs](usage.md) ·
[Try the tested PA example](pa-hq-tutorial.md)
