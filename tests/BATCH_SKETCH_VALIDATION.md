# Native batch sketch validation (v0.5.1)

## Change under test

`microsags sketch -t 50` now starts one native sketch process containing up to
50 persistent C++ workers, rather than starting a separate executable for each
reference. The per-reference DNA2bit implementation, sketch parameters, packed
database builder, annotation rules, and assembly workflow are unchanged.

## Full database check

Measured on 2026-09-27 (server time, UTC+08:00), using the same frozen list of
199,923 GTDB232 references and the same taxonomy as the previous measurement.
Both runs requested 50 workers. Wall time covers the complete `microsags sketch`
command, including reference sketching and database packing.

| Measurement | Previous scheduler | Native batch scheduler |
| --- | ---: | ---: |
| Wall time (seconds) | 589.58 | 330.26 |
| User CPU time (seconds) | 12,634.15 | 10,963.19 |
| System CPU time (seconds) | 1,328.16 | 166.59 |
| Sampled peak process-group PSS (KiB) | 1,257,992 | 2,022,960 |
| Sampled peak process-group RSS (KiB) | 1,338,512 | 2,034,756 |

The observed wall-time reduction was 44.0% (1.79x faster). These are single,
non-exclusive server measurements; concurrent workload and filesystem caching
were not controlled. They are not a claim of a controlled or universal speedup.
Memory samples include the complete benchmark process group and were collected
every approximately two seconds; sampled peaks can miss shorter spikes.

The new and previous `references.pack`, `references.tsv`, and
`genome_taxonomy.csv` were each byte-identical, including reference row order.
The packed payload is 1,381,867,776 bytes (199,923 rows of 6,912 bytes).

| File | SHA256 (both builds) |
| --- | --- |
| `references.pack` | `ebaaca5b152a3475cddbda18b98c572e00b6cfd1d5d621926d1ddfa9db31c6df` |
| `references.tsv` | `ab6a2557161921b589c6bb7d3b8831168880e4bbdc41b99ea30b0ea73b727abe` |
| `genome_taxonomy.csv` | `b10a9a80d42455f7fdc419afcbf231fb06119ee5de4346d5cb5bf8ba9749b78b` |

Three reference FASTA files were also submitted through the old and new public
Annotation commands with their respective databases. The resulting
`annotations.tsv` files were byte-identical. This is a compatibility smoke test,
not an independent assessment of annotation accuracy.

## Regression checks

All seven CTest targets passed. The native batch suite additionally passed
against the previously compiled v0.5.0 sketch executable:

- Exact sketch bytes for plain and compressed FASTA, multiple records, ambiguous
  bases, lowercase bases, short sequences, paths with spaces, and CRLF lists.
- Worker counts of 1, 4, and 50 (capped at the number of references).
- The legacy paired-input interface.
- Failures for empty or missing input lists, invalid arguments, and colliding
  output basenames.
- Preservation of an existing output file rather than overwriting it.

The public CLI suite (19 tests) also verifies that sketch construction invokes
one batch command and then the existing pack builder.

To repeat the synthetic checks after building:

```bash
ctest --test-dir build --output-on-failure
python3 tests/test_batch_sketch.py --binary build/dna2bit-embedded-sketch
```

For cross-version byte checking, append
`--legacy-binary /path/to/previous/dna2bit-embedded-sketch` to the second command.

The full benchmark retains its exact command, source snapshot, binary hashes,
input hashes, build logs, resource samples, old/new annotation outputs, and
database equivalence receipt separately from the source-only repository.
