# Lake speed-v1 full-run release audit

This checklist is for the completed write-once run only.  Every command is
read-only with respect to the scientific run and the old baseline.

## 1. Confirm the runner is finished

```bash
systemctl --user show cellbit-lake-full-e2e-220t-speed-v2-20260905.service \
  -p ActiveState -p SubState -p Result -p ExecMainStatus -p NRestarts \
  -p MemoryPeak -p CPUUsageNSec
journalctl --user -u cellbit-lake-full-e2e-220t-speed-v2-20260905.service \
  --no-pager -n 120
```

Run the scientific audit only when the unit is inactive, `Result=success`,
`ExecMainStatus=0`, and no child process remains.

## 2. Run the strict read-only audit

Standard audit (all closures, normalized scientific comparisons, final bins,
quality, times, telemetry, and receipt outputs):

```bash
/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/tools/audit_lake_speed_release_v1.py
```

Deep release audit (recommended before release; also hashes every path declared
as a Stage3B receipt input and compares all 336 merged subassembly inputs):

```bash
/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/tools/audit_lake_speed_release_v1.py --deep
```

The script writes nothing.  It emits one JSON object on stdout and exits 0 only
when all hard gates pass.  If a permanent audit record is wanted, redirect it
to a separately claimed, new write-once audit directory; never redirect into
the scientific run root.

## 3. Required hard gates

- `RUN_RESULT.json`: PASS, exit 0, 220-thread quota, no prior scientific reuse.
- Cgroup receipt: 220 cores, CPU affinity 0-239, 1.2 TiB high water, 1.45 TiB
  hard memory ceiling, 8 GiB swap ceiling, unlimited tasks, expected nofile.
- Stage1-3A: 13,742 inputs = 7,527 labeled + 6,215 pending; 274 groups;
  13,742 Stage1 and sketch markers; packed search/labels/pending/statistics equal
  to the old teacher baseline.  In deep mode, all 274 merged inputs are hashed
  against the baseline.
- Stage3B upstream: input-view, pending-statistics, CheckM2, GTDB-Tk finalizer,
  AUDIT and COMPLETE receipts close on exactly 6,215 SAGs; 33,576 selected
  marker rows; CheckM2 scientific rows, quality table, and marker map equal the
  old baseline.
- Stage3B exact graph: 5,834 nodes; `5834 choose 2 = 17,014,861` requested and
  evaluated pairs; 322,110 sparse ANI/AF rows; 103,555 positive and 48,205
  negative edges; 62 final partitions; 2,567 assigned + 3,267 unaggregated =
  5,834.  ANI/AF is normalized by unordered SAG pair with AF direction retained;
  graph partitions are normalized as sets, so harmless cluster-label ordering
  cannot mask or create a scientific difference.
- Stage receipts: three main PASS receipts, 123 makeblastdb receipts, 123 BLASTN
  receipts, and 62 merged-input receipts.  Standard mode hashes receipt outputs;
  deep mode hashes inputs too.
- Final view: exactly 336 unique IDs (274 Stage3A + 62 Stage3B), same ID set as
  baseline, every symlink bound to the manifest source, every source size/SHA
  correct.
- Final CheckM2: exactly the same 336 IDs as final view, finite/bounded metrics,
  CheckM2 1.0.1 auto mode, and independently recalculated HQ/MQ counts equal
  `RESULT.json` and top-level `RUN_RESULT.json`.
- Every one of the nine `/usr/bin/time -v` records exits 0; all synchronous
  phase events occur once in order; telemetry counters are valid/monotonic and
  peak memory stays under the hard limit.  Peak load, PSI, memory, I/O, CPU time,
  and global swap growth are reported in the audit JSON.
- Both final SHA256 manifests and all cross-file receipt hashes verify.

## 4. Important interpretation

Do not require final Flye assembly bytes or HQ/MQ counts to equal the old
multi-thread baseline.  The old Flye assemble step was proven non-deterministic.
The hard equivalence boundary is its input plus all upstream scientific graph
artifacts; the new final FASTAs are instead closed against their own immutable
manifest and exact CheckM2 ID set.  Any upstream ANI/AF, edge, membership, or
input-subassembly difference is a release blocker.

If any gate fails, preserve the entire write-once run unchanged and diagnose
the first failed artifact.  Do not edit a receipt, clear an incomplete path, or
rerun into the same output root.
