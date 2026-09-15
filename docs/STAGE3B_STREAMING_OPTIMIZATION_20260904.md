# Stage3B C++ streaming/orchestration optimization (2026-09-04)

> **Historical execution-equivalence record; superseded for current production.**
> This document audits the former GTDB-reference double-negative contract.  Its
> evidence and expected values remain immutable; the current contract is defined
> in `README_STAGE3B.md` and sends every quality-passing Dna2bit-negative SAG to
> exact SAG–SAG clustering without reference filtering.

This change is an execution-only refactor of `src/stage3b.cpp`.  It does not
change the Stage3B scientific contract.

## Frozen scientific contract

- quality: `contamination < 5`, `max_contig >= 1000`;
- label-positive: non-empty GTDB taxonomy, `ANI >= 95`, `AF_query >= 50`;
- graph edge: `ANI >= 95`, `abs(GC1-GC2) <= 2` and the original weight;
- marker sequence length `>=150`, longest sequence per SAG/marker;
- BLASTN `evalue=1e-5`, `max_target_seqs=10000`, two threads, original four
  output fields, longest-alignment selection, `n_markers>=3`, mean PID `<97`;
- positive precluster: RB, resolution 1, seed 20260811, iterations -1;
- signed sweep, in declared order: `(1,1), (1,3), (1,10), (2,3), (2,10),
  (2,30)`; lexicographic score selection; expected Lake winner `r=2, lambda=30`;
- emitted clusters have size `>=10`;
- the subassembly command and every Flye executable/config/threads argument are
  unchanged.

## Execution changes

1. TSVs are consumed one row at a time.  The former
   `vector<unordered_map<string,string>>` representation no longer retains the
   entire taxonomy, search, triangle, marker-map or membership table.
2. FASTA validation, total/max/GC statistics and SHA-256 are computed in one
   byte pass.  A size/mtime and (on Linux) device/inode/mtime/ctime-bound digest
   cache reuses that digest.  A changed file is rejected rather than silently
   re-hashed into the same run.
3. Stage output SHA-256 values are computed once and reused in receipts instead
   of reading every output again to render its receipt.
4. Large derived TSVs are atomically streamed to a temporary file instead of
   being retained in `ostringstream` memory.  Resume still compares size and
   SHA-256 before discarding the recomputation.
5. Marker pair evidence stores `(count, sum, minimum)` rather than one heap
   vector of PID values per pair.  Marker sequence maps are released as soon as
   their deterministic marker FASTA has been materialized.
6. Independent marker BLAST jobs run with a hard bound of
   `max(1, floor(--threads/2))`.  Each BLAST still uses exactly two threads, so
   active BLAST threads never exceed the requested total thread budget when
   `--threads >= 2` (`--threads=1` retains the historical two-thread BLAST
   command rather than changing a scientific execution parameter).
   Results are parsed later in the original lexicographic marker order; worker
   completion order cannot change marker means, edge rows or Leiden input.
7. Linux tools are invoked with `fork/execvp` argument vectors rather than a
   shell.  The human-readable command string and its receipt hash remain in the
   original form.  Windows retains the quoted shell fallback for build/self-test
   compatibility; production remains Linux.
8. Membership cluster-size validation is linear rather than
   `groups x membership_rows`.
9. Merged group FASTAs gain `MERGED_INPUT.PASS.json`, binding ordered source
   fingerprints, member count and merged-output fingerprint.  New-run content
   is unchanged.  A future resume avoids reconstructing a second full merged
   FASTA; an old receipt-less run is verified once by deterministic replay and
   then receives the new receipt.
10. Group membership vectors and the positive/marker edge maps are released as
    soon as their last consumer completes.  Runtime I/O/cache/concurrency
    counters are printed to stderr for benchmarking without adding
    non-deterministic fields to scientific artifacts.

For the frozen raw Lake source set (1,621,776,239 input-file bytes), the new
single-pass stat+SHA path removes at least one complete duplicate source read
and also avoids the later double-negative subset re-read.  Exact wall/RSS gains
must be reported only after a separate v0/v1 cold benchmark; they are not
invented here.  Marker throughput can now scale across independent markers,
subject to filesystem contention, while the thread ceiling remains fixed.

The stripped `-O2` orchestration executable compiled on `bio-codex` is 289,080
bytes (0.276 MiB), versus 227,640 bytes for v0.  The 61,440-byte increase buys
streaming, immutable digest caching, safe argv execution and bounded workers;
it is negligible beside the scientific databases.  A later package-size pass
may use `-Os -flto -s`, but changing build flags was intentionally outside this
file-scoped task.

## Frozen Lake equivalence gate

Do not run v1 against or resume into the frozen v0 directory.  Produce two
separate write-once roots with identical frozen inputs/tools, then run:

```bash
python3 tests/audit_stage3b_v0_v1_frozen_lake_20260904.py \
  --v0 /path/to/frozen_v0_output \
  --v1 /path/to/new_v1_output \
  --json-out /new/write_once/path/STAGE3B_V0_V1_AUDIT.json
```

The audit fails closed on:

- any deterministic file-set, byte or line-order difference (run-root prefixes
  alone are normalized);
- any marker FASTA or BLAST row difference;
- any merged input difference;
- any canonical Flye FASTA-record difference;
- broken PASS receipt output size/SHA bindings;
- Lake count drift (`7408` nodes, `53250` positive edges, `3963` marker targets,
  `624049` marker pairs, `60271` negative edges);
- any change in the six sweep records or their full cluster details;
- a winner other than `gs_r2.0_lam30`, score other than `(50,3380)`, chosen
  membership other than 58 clusters/3759 SAGs/largest 474, or the frozen
  membership/report SHA-256 values.

The auditor is read-only unless a new `--json-out` is supplied, and that output
is opened with exclusive-create semantics.
