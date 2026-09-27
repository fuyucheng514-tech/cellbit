# Sketch optimization compatibility (v0.5.2)

## Scope

This release retains the v0.5.1 persistent-worker scheduler. It accelerates
bounded, single-member gzip decoding with optional libdeflate and uses narrower
working counters with an exact, overflow-guarded flush into the legacy wide
accumulator. The original kseq parser and DNA2bit traversal, reverse-complement,
ambiguous-base and bit-serialization rules are preserved.

Plain input, large gzip streams, concatenated members and unsupported containers
use the existing zlib streaming path. Builds without libdeflate are supported.
The scientific parameters, database format and public workflow are unchanged.

## Full-data compatibility evidence

The validated engine source SHA256 is
`2801405c6139aa5436651ff82a7b4b0c75f1de530c8b990779b1e981ab92332e`.
It was checked against the v0.5.1 baseline on 199,923 compressed GTDB232 reference
genomes. Every byte of `references.pack`, `references.tsv` and
`genome_taxonomy.csv` matched. Fifty reference-derived annotation canaries also
produced byte-identical `annotations.tsv` files. These are compatibility checks,
not an independent assessment of species-assignment accuracy.

Separate differential API checks passed all 56 cases in each configuration:

- normal libdeflate-accelerated build;
- forced small counter budget, exercising flush and wide-counter fallback;
- build without libdeflate.

Cases include multiple hash algorithms, sketch dimensions, paired inputs,
short records, ambiguous and lowercase bases, and compressed input. The
scientific parameters of the complete database check were not changed.

## Reproducible synthetic regression tests

```bash
ctest --test-dir build --output-on-failure
python3 tests/test_batch_sketch.py \
  --binary build/dna2bit-embedded-sketch \
  --legacy-binary /path/to/v0.5.1/dna2bit-embedded-sketch
```

The suite checks exact sketch bytes, thread-count independence, paired input,
path validation, existing-output protection, gzip container variants, and an
independent integer oracle for the archived wyhash-17 contract.

To check the optional-dependency fallback in a separate build directory:

```bash
cmake -S . -B build-no-libdeflate -DCMAKE_BUILD_TYPE=Release \
  -DDNA2BIT_ENABLE_LIBDEFLATE=OFF
cmake --build build-no-libdeflate --target dna2bit-embedded-sketch -j2
python3 tests/test_batch_sketch.py \
  --binary build-no-libdeflate/dna2bit-embedded-sketch \
  --legacy-binary build/dna2bit-embedded-sketch
```

No timing assertion is part of these tests. Full benchmark receipts and data
remain outside this source-only repository.
