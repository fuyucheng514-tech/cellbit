# Cellbit SAG Pipeline speed-v1

This branch changes execution strategy only. Scientific inputs, seeds, ANI/AF
formulae, thresholds, CheckM2 models, Flye configuration, and output acceptance
rules remain unchanged.

## Frozen lake baselines

- 220-thread core pipeline: 6,825.91 s.
- 50-thread core pipeline: 13,782 s.
- The old 220-thread Stage3B run requested and evaluated all 17,014,861 pairs.
- Old 220-thread long stages: Dna2bit search 315.23 s, Stage3A subassembly
  509.86 s, upstream CheckM2 894.35 s, Stage3B 4,772 s.

## Accepted optimizations

| Area | Frozen comparison | speed-v1 result | Verification |
|---|---:|---:|---|
| Dna2bit reference search | 1,450.65 s at 50 threads | 27.60 s | byte-identical teacher output |
| Stage3B exact ANI/AF, 1,024 nodes | 306.75 s | 57.76–57.94 s | all rows byte-identical |
| Signed Leiden, full 5,834 nodes | 616.83 s | 191.22 s | membership and report byte-identical |
| cpp-subass G0001 | 226.52 s | 136.10 s | core assembly artifacts byte-identical |
| CheckM2, 6,215 bins | 1,082.31 s | 807.37 s | all 13 original DIAMOND chunks byte-identical; normalized per-bin report exact |
| Final 336-bin evaluation view | 265.94 s | 27.23 s | IDs, order, source FASTA hashes, links, and manifest semantics exact |

## Complete 220-thread lake acceptance

The fresh-output, write-once acceptance run completed successfully at:

`/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/24_lake_full_e2e_220t_speed_v1_20260905/attempt_001`

| Scope | Frozen 220-thread baseline | speed-v1 | Speedup |
|---|---:|---:|---:|
| Core pipeline | 6,825.91 s | 2,786 s | 2.450x |
| Core plus final evaluation | about 7,238.50 s | 2,926 s | 2.474x |
| Stage3B | 4,772 s | 937.68 s | 5.089x |
| Upstream CheckM2 | 894.35 s | 425.19 s | 2.103x |
| Final-bin view | 265.94 s | 26.83 s | 9.912x |
| Final CheckM2 | 146.65 s | 111.43 s | 1.316x |

The acceptance run closed all 13,742 inputs, all 5,834 Stage3B graph nodes,
and all 17,014,861 exact ANI/AF pairs.  It produced the same 274 Stage3A plus
62 Stage3B bins.  Peak cgroup memory was 126,201,655,296 bytes, with no OOM.

Final quality was HQ=0 and MQ-only=22.  This matches the deterministic
50-thread Flye branch.  It is a known scientific-quality limitation of the
current grouping rules, not a speed-pipeline counting error.

An independent deep release audit passed all eight sections.  It verified the
13,742-input authority closure, Stage3B exact graph and receipt closure, final
336-bin identity closure, CheckM2 ID/quality closure, result/SHA chains, timing,
and cgroup telemetry.  The formal audit record is
`SPEED_RELEASE_DEEP_AUDIT_V2.json` beside the acceptance root and has SHA256
`9da70fefc8a913361546344751a1a3fb989e98b778f06a3295ac262080810539`.

## Accepted follow-on canaries

Flye's deterministic finishing phase was tested on frozen G0001 with 2, 4,
and 8 threads, twice at every setting.  The final FASTA was byte-identical in
all six runs, and all intermediate FASTAs were sequence-normalized identical.
Only an internal consensus record order differed and did not propagate to the
final assembly.  Mean wall times were 89.985 s, 59.02 s, and 40.43 s,
respectively.  Eight finish threads therefore gave a 2.226x speedup over two
threads while preserving the final raw FASTA exactly.  A second gate on the
largest Stage3A group is required before this adaptive policy is released.

The Stage3B upstream finalizer was also changed from serial SAG parsing to a
bounded 64-process ordered parse.  On all 6,215 SAGs it fell from 216.63 s to
9.58 s (22.61x).  Both scientific outputs were byte-identical to the serial
result, and every COMPLETE binding field except the independently generated
audit-file digest was identical.  This canary is preserved at
`/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/29_finalize_parallel_full_canary_20260905/attempt_001`.

## What changed

1. Dna2bit stores the frozen reference sketches in a compact packed database
   and scans that database reference-major in cache-friendly query blocks.
2. Exact ANI/AF still covers every unordered pair, but one global seed index
   dispatches candidate seed hits instead of rebuilding the same reference
   lookup for each pair. The DP trace uses a compact byte layout and stable row
   pointers.
3. The six independent signed-Leiden parameter runs execute concurrently and
   retain their original fixed seeds and deterministic selection order.
4. Flye's non-deterministic `assemble` phase is fixed at one thread per group;
   many groups run concurrently. Later deterministic finishing phases use two
   threads per group. Longest-processing-time scheduling reduces tail latency.
5. CheckM2 preserves its original 500-bin DIAMOND partitions and models, but
   independent partitions execute concurrently under one strict total thread
   budget. Model reference data is loaded once rather than once per batch.
6. Final-bin validation reads independent bins in a 16-process pool while
   preserving the original deterministic manifest order and strict hashes.

## Rejected experiments

- Merging or repartitioning CheckM2 DIAMOND queries changed hits and was rejected.
- Running all 6,215 CheckM2 model rows as one prediction batch changed floating
  output and was rejected; the proven-equivalent 1,000-row model batching is
  retained.
- ANI/AF bitmap-rank and alternative map/layout experiments that did not improve
  wall time were reverted.
- No comparison pair, quality threshold, model, or biological rule was removed.
