# Microsags

Microsags provides DNA2bit species annotation and species-guided subassembly
for single-amplified genomes (SAGs).

## Annotation mode

Annotation mode accepts paired FASTQ reads, singleton FASTQ reads, or assembled
contigs. FASTQ reads are quality-controlled with fastp and passed directly to
DNA2bit; they are not assembled first.

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

The principal result maps each SAG identifier to an accepted species label.

## Assembly mode

Assembly mode is the next step. It accepts one contig FASTA per SAG and a
completed Annotation-mode output directory.

```bash
microsags assemble SAG_contigs/ \
  -a annotation_output/annotations.tsv \
  -o assembly_output
```

Labelled SAGs enter Stage 3A; unclassified SAGs enter Stage 3B. Both routes use
the existing `cpp-subass`/Flye subassembly workflow. If annotation was run
from FASTQ reads, users must first assemble every SAG independently and preserve
the same SAG identifiers.

Continue with [Install](install.md) or [Usage examples](usage.md).
