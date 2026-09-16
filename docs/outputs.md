# Output reference

Microsags publishes a deliberately small result surface. Intermediate files
are written to a sibling work directory. Failed work is retained for diagnosis;
successful work is removed after the final files have been published.

## Annotation mode

```text
annotation_output/
└── annotations.tsv
```

`annotations.tsv` has exactly two columns, `sag_id` and `species`, and exactly
one data row per input SAG. `UNCLASSIFIED` means DNA2bit did not accept a label.

## Assembly mode

```text
assembly_output/
├── stage3a.tsv
├── stage3a/
│   └── <species>.fasta
├── stage3b_clusters.tsv
└── stage3b/
    └── <group_id>.fasta
```

`stage3a.tsv` maps each accepted species label to its final assembled FASTA.
`stage3b_clusters.tsv` lists each Stage 3B group and every SAG assigned to it.
The two FASTA directories contain only final cluster assemblies.
