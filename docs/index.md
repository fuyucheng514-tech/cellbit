# Microsags

Microsags is a Linux command-line workflow for single-amplified genome (SAG)
species annotation and species-guided subassembly. It accepts assembled SAG
contigs or paired short reads through the same manifest and determines the input
route from sequence content rather than filename extensions.

The default v0.1 workflow performs four auditable operations:

1. validate and route every SAG input;
2. assemble paired reads with fastp and SPAdes, while preserving supplied
   contigs without reassembly;
3. classify eligible SAGs with the embedded, teacher-compatible DNA2bit
   `k=17` packed search;
4. subassemble DNA2bit-labelled SAGs by species in Stage 3A.

DNA2bit-negative SAGs are written to a stable hand-off manifest. Stage 3B is
available as a separate advanced executable and is not represented as part of
the one-command stable workflow.

## Start here

- [Install](install.md) describes Pixi, Conda and source builds.
- [Usage examples](usage.md) explains the input routes and functional outputs.
- [Lake quarter tutorial](lake-quarter-tutorial.md) records a reproducible
  3,436-SAG integration example.
- [Output reference](outputs.md) defines the files emitted by each stage.

## Scope and scientific identity

Microsags embeds the compatible C++ implementation of the original DNA2bit
sketch and packed-search semantics. It does not use Cellbit57 ALC, Top16,
skani-based GTDB classification, or a trained Cellbit57 model. Mature external
algorithms such as fastp, SPAdes and Flye remain declared, versioned
dependencies; orchestration and scientific receipts do not rename those tools
as native Microsags algorithms.

