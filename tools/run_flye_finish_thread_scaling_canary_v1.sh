#!/usr/bin/env bash
set -euo pipefail

BASE=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass
RELEASE="$BASE/17_cellbit_sag_pipeline_cpp_speed_v1_20260905"
ROOT="$BASE/28_flye_finish_thread_scaling_canary_20260905/attempt_001"
CONTROL="$BASE/28_flye_finish_thread_scaling_canary_20260905/control_attempt_003"
INPUT="$BASE/11_lake_full_current_tractor_baseline_v1/04_stage3b_cellbit_negative_exact_attempt_003/06_subassemble/G0001/input_subassemblies.fasta"
# This phase-split canary was produced by the current deterministic release and
# contains both the frozen one-thread draft and its phase-completion marker.
ASSEMBLY_SOURCE="$BASE/18_speed_canary_20260905/subass_g0001_phase_split_v2/run/00-assembly"
SUBASS="$RELEASE/build_speed_v1/subass/cpp-subass"
AUDITOR="$RELEASE/tools/audit_flye_finish_thread_scaling_v1.py"
FLYE_ROOT=/home/data/fyc/biosoft/miniconda3/envs/assemble
MAIN_UNIT=cellbit-lake-full-e2e-220t-speed-v2-20260905.service
LOCK="$BASE/.flye_finish_thread_scaling_canary_v1.lock"

mkdir -p "$CONTROL"
exec 9>"$LOCK"
if ! flock -n 9; then
  printf '%s\n' 'another Flye finish thread-scaling canary owns the lock' >&2
  exit 73
fi
if [[ -e "$ROOT" ]]; then
  printf '%s\n' "refusing to reuse write-once root: $ROOT" >&2
  exit 74
fi

for path in "$INPUT" "$ASSEMBLY_SOURCE/draft_assembly.fasta" \
  "$ASSEMBLY_SOURCE/CPP_ASSEMBLY_PHASE_COMPLETE" "$SUBASS" "$AUDITOR" \
  "$FLYE_ROOT/bin/flye-modules" "$FLYE_ROOT/bin/flye-minimap2" \
  "$FLYE_ROOT/bin/flye-samtools" \
  "$FLYE_ROOT/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg"; do
  [[ -s "$path" ]] || { printf 'missing required authority: %s\n' "$path" >&2; exit 66; }
done

# Never compete with the formal end-to-end acceptance run.  Afterwards require
# three consecutive low-pressure samples before creating the scientific root.
stable=0
while (( stable < 3 )); do
  active=$(systemctl --user is-active "$MAIN_UNIT" 2>/dev/null || true)
  load1=$(awk '{print $1}' /proc/loadavg)
  mem_kib=$(awk '/^MemAvailable:/{print $2}' /proc/meminfo)
  io_some=$(awk '$1=="some"{for(i=1;i<=NF;i++)if($i~/^avg10=/){split($i,a,"=");print a[2]}}' /proc/pressure/io)
  mem_some=$(awk '$1=="some"{for(i=1;i<=NF;i++)if($i~/^avg10=/){split($i,a,"=");print a[2]}}' /proc/pressure/memory)
  pass=$(awk -v a="$active" -v l="$load1" -v m="$mem_kib" -v i="$io_some" -v p="$mem_some" \
    'BEGIN{print (a!="active" && l<80 && m>838860800 && i<2 && p<1)?1:0}')
  if [[ "$pass" == 1 ]]; then stable=$((stable+1)); else stable=0; fi
  printf '%s\tmain=%s\tload1=%s\tmem_available_kib=%s\tio_some_avg10=%s\tmem_some_avg10=%s\tstable=%s/3\n' \
    "$(date --iso-8601=seconds)" "$active" "$load1" "$mem_kib" "$io_some" "$mem_some" "$stable" \
    | tee -a "$CONTROL/GATE.tsv"
  (( stable == 3 )) || sleep 30
done

mkdir -p "$ROOT"
{
  sha256sum "$INPUT"
  sha256sum "$ASSEMBLY_SOURCE/draft_assembly.fasta"
  sha256sum "$ASSEMBLY_SOURCE/CPP_ASSEMBLY_PHASE_COMPLETE"
  sha256sum "$SUBASS"
  sha256sum "$AUDITOR"
  sha256sum "$FLYE_ROOT/bin/flye-modules"
  sha256sum "$FLYE_ROOT/bin/flye-minimap2"
  sha256sum "$FLYE_ROOT/bin/flye-samtools"
  sha256sum "$FLYE_ROOT/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg"
} > "$ROOT/AUTHORITY_SHA256.txt"
printf 'assemble_threads\t1\ninput\t%s\nassembly_source\t%s\n' "$INPUT" "$ASSEMBLY_SOURCE" > "$ROOT/DESIGN.tsv"

for threads in 2 4 8; do
  for replicate in 1 2; do
    name="t${threads}_r${replicate}"
    out="$ROOT/$name"
    mkdir -p "$out/00-assembly"
    cp --reflink=auto "$ASSEMBLY_SOURCE/draft_assembly.fasta" "$out/00-assembly/draft_assembly.fasta"
    cp --reflink=auto "$ASSEMBLY_SOURCE/CPP_ASSEMBLY_PHASE_COMPLETE" "$out/00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE"
    /usr/bin/time -f $'elapsed_seconds\t%e\nuser_seconds\t%U\nsystem_seconds\t%S\nmax_rss_kb\t%M\nfs_inputs\t%I\nfs_outputs\t%O' \
      -o "$ROOT/$name.time.tsv" \
      "$SUBASS" --reads "$INPUT" --out-dir "$out" \
      --flye-modules "$FLYE_ROOT/bin/flye-modules" \
      --minimap2 "$FLYE_ROOT/bin/flye-minimap2" \
      --samtools "$FLYE_ROOT/bin/flye-samtools" \
      --package-root "$FLYE_ROOT/lib/python3.9/site-packages/flye" \
      --config "$FLYE_ROOT/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg" \
      --threads "$threads" --assemble-threads 1 --phase finish \
      --no-overlap-policy passthrough \
      > "$ROOT/$name.stdout.log" 2> "$ROOT/$name.stderr.log"
  done
done

python3 "$AUDITOR" --root "$ROOT"
sha256sum "$ROOT/AUTHORITY_SHA256.txt" "$ROOT/DESIGN.tsv" "$ROOT/SUMMARY.tsv" "$ROOT/RESULT.json" \
  > "$ROOT/CONTROL_SHA256.txt"
