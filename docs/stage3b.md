# Adjust Leiden clustering

Leiden controls how the DNA2bit-unclassified SAGs are divided into Stage 3B
clusters before subassembly. Most users can keep the default setting.

Use `--leiden-resolution` only when you want to change the cluster granularity:

- a lower value usually produces fewer, larger clusters;
- a higher value usually produces more, smaller clusters;
- changing this value can change both cluster membership and the resulting
  Stage 3B assemblies.

For example, set the resolution to `0.18` with:

```bash
microsags assemble SAG_contigs/ \
  -a annotation_output/annotations.tsv \
  -o assembly_output \
  --leiden-resolution 0.18
```

Replace `0.18` with the value you want to test. The value must be greater than
zero. If the option is omitted, Microsags uses its default Leiden parameter
selection. This option changes Stage 3B clustering only; it does not change
DNA2bit annotation or Stage 3A grouping.
