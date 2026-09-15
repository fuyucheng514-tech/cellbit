# cpp-subass-full

Independent C++17 orchestration and glue implementation of the Flye 2.9.6
`--subassemblies` path.  The scientific assembly kernels are the official
`flye-modules` C++ executable; this project replaces the Python driver and the
subassembly-specific intermediate processing with C++.

The release is not considered validated until the XA cold comparison reaches
the final `assembly.fasta` on both implementations and every stage receipt is
closed.

Pipeline: configure -> assemble -> repeat -> contigger -> align -> bubbles ->
polisher -> coverage filter -> polished graph -> finalize.

The default `--no-overlap-policy strict` preserves Flye's normal failure when
assembly returns an empty `draft_assembly.fasta`.  The production Stage 3A
launcher opts in explicitly with `--no-overlap-policy passthrough`.  That mode
is deliberately narrow: passthrough is allowed only when `flye-modules
assemble` returned zero, wrote a real zero-byte draft, and the dedicated log
from that exact invocation contains Flye's explicit no-overlap/zero-disjointig
message.  The input FASTA is validated and copied byte-for-byte to
`assembly.fasta`; `NO_OVERLAP_PASSTHROUGH.PASS.json` binds input, output, empty
draft, evidence log, command, SHA-256 values, record/base/GC/length statistics,
and ID/record-set closure.  `CPP_FULL_PIPELINE_PASS` records the mode.  A
non-zero exit, missing/malformed input, missing draft, stale/missing evidence,
or any other error remains fail-closed.  A non-empty draft always follows the
unchanged full Flye path.

## Build

`htslib` may be discovered through `pkg-config`, or from a self-contained
installation prefix with `-DHTS_ROOT=...`; `zlib` is resolved with CMake's
standard `FindZLIB` module.  A lean production build is:

```bash
export PKG_CONFIG_PATH=/path/to/environment/lib/pkgconfig
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --parallel
```

For Conda-style environments that provide `include/htslib` and `lib/libhts.so`
but no `htslib.pc`, use:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
  -DHTS_ROOT=/path/to/environment
cmake --build build --parallel
```

Release builds use function/data sections and linker garbage collection.  IPO/LTO
is enabled only when CMake verifies that the compiler and linker support it; use
`-DSUBASS_ENABLE_IPO=OFF` for toolchains where LTO is undesirable.

With `BUILD_TESTING=ON`, CMake also builds focused `cpp-subass-bubble-test` and
`cpp-subass-polish-glue-test` executables.  They are intended for fixed-input,
byte-for-byte regression checks and are not part of the production executable.

The optimized implementation preserves the original ordering, random seed,
coverage rules, bubble boundaries, branch selection, and output formatting.  It
reduces memory by allocating insertion maps only where insertions exist, reuses a
single indexed BAM reader across contigs, moves large alignment/sequence objects
instead of copying them, and releases per-contig state before writing the next
contig.
