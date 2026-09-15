#!/usr/bin/env bash
set -Eeuo pipefail

# Write-once full Lake Stage3B run for the corrected entry rule:
# Cellbit/Dna2bit negative + historical quality gates.  GTDB reference status
# is intentionally absent from the decision path.

readonly CODE_ROOT="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/12_cellbit_negative_stage3b_code_v1"
readonly BASE_ROOT="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1"
readonly OUT_ROOT="${BASE_ROOT}/04_stage3b_cellbit_negative_exact_attempt_002"
readonly CONTROL_ROOT="${BASE_ROOT}/05_stage3b_cellbit_negative_exact_control_attempt_002"
readonly THREADS=220

readonly MANIFEST="${BASE_ROOT}/01_stage1_3a/03B_unclassified_pending.tsv"
readonly QUALITY="${BASE_ROOT}/03b_stage3b_inputs_attempt_001/stage3b_quality.tsv"
readonly MARKERS="${BASE_ROOT}/03b_stage3b_inputs_attempt_001/bac120_marker_nt_map.tsv"
readonly TRACTOR="${CODE_ROOT}/build/sag-stage3b-tractor"
readonly ANI_ENGINE="${CODE_ROOT}/build/gtdb-ani-af"
readonly SUBASS="${CODE_ROOT}/build/subass/cpp-subass"
readonly LEIDEN="${CODE_ROOT}/python/stage3b_signed_leiden.py"
readonly PYTHON="/home/data/fyc/biosoft/miniconda3/bin/python"
readonly MAKEBLASTDB="/home/data/fyc/biosoft/miniconda3/envs/decon_bench/bin/makeblastdb"
readonly BLASTN="/home/data/fyc/biosoft/miniconda3/envs/decon_bench/bin/blastn"
readonly FLYE_ROOT="/home/data/fyc/biosoft/miniconda3/envs/assemble"

if [[ -e "${OUT_ROOT}" || -e "${CONTROL_ROOT}" ]]; then
  printf 'write-once target already exists; refusing to start\n' >&2
  exit 73
fi

for required in "${MANIFEST}" "${QUALITY}" "${MARKERS}" "${TRACTOR}" \
                "${ANI_ENGINE}" "${SUBASS}" "${LEIDEN}" "${PYTHON}" \
                "${MAKEBLASTDB}" "${BLASTN}"; do
  [[ -e "${required}" ]] || { printf 'missing required path: %s\n' "${required}" >&2; exit 66; }
done

mkdir -p "${CONTROL_ROOT}"
umask 022

readonly TELEMETRY="${CONTROL_ROOT}/RESOURCE_TELEMETRY.tsv"
readonly TIME_FILE="${CONTROL_ROOT}/STAGE3B_TIME.txt"
readonly STDOUT_LOG="${CONTROL_ROOT}/STAGE3B_STDOUT.log"
readonly STDERR_LOG="${CONTROL_ROOT}/STAGE3B_STDERR.log"
readonly FOOTPRINT="${CONTROL_ROOT}/SOFTWARE_FOOTPRINT.tsv"

exec > >(tee -a "${STDOUT_LOG}") 2> >(tee -a "${STDERR_LOG}" >&2)

# Record the exact code and runtime payload actually used by this run.  Shared
# operating-system libraries and unrelated files in the Conda environments are
# deliberately not charged to the program; every row is a concrete bound path.
record_path_size() {
  local component="$1" path="$2"
  [[ -e "${path}" ]] || { printf 'missing footprint path: %s\n' "${path}" >&2; exit 66; }
  printf '%s\t%s\t%s\n' "${component}" "$(du -sb "${path}" | awk '{print $1}')" "${path}" >> "${FOOTPRINT}"
}

printf 'component\tbytes\tpath\n' > "${FOOTPRINT}"
record_path_size tractor_source_and_build "${CODE_ROOT}"
record_path_size tractor_binary "${TRACTOR}"
record_path_size ani_engine_binary "${ANI_ENGINE}"
record_path_size subass_binary "${SUBASS}"
record_path_size blast_makeblastdb_binary "${MAKEBLASTDB}"
record_path_size blastn_binary "${BLASTN}"
record_path_size flye_modules "${FLYE_ROOT}/bin/flye-modules"
record_path_size flye_minimap2 "${FLYE_ROOT}/bin/flye-minimap2"
record_path_size flye_samtools "${FLYE_ROOT}/bin/flye-samtools"
record_path_size flye_python_package "${FLYE_ROOT}/lib/python3.9/site-packages/flye"
while IFS=$'\t' read -r component path; do
  record_path_size "${component}" "${path}"
done < <("${PYTHON}" - <<'PY'
import pathlib
import igraph
import leidenalg
print("python_igraph_package\t" + str(pathlib.Path(igraph.__file__).resolve().parent))
print("python_leidenalg_package\t" + str(pathlib.Path(leidenalg.__file__).resolve().parent))
PY
)

start_epoch="$(date +%s)"
printf 'timestamp_epoch\telapsed_s\tload1\tmem_available_bytes\tswap_used_bytes\tcgroup_memory_current\tcgroup_memory_peak\tcgroup_cpu_usage_usec\tcgroup_io_read_bytes\tcgroup_io_write_bytes\tcpu_psi_some_avg10\tmemory_psi_some_avg10\tmemory_psi_full_avg10\tio_psi_some_avg10\tio_psi_full_avg10\n' > "${TELEMETRY}"

psi_value() {
  local file="$1" row="$2"
  awk -v wanted="${row}" '$1 == wanted { for (i=2; i<=NF; ++i) if ($i ~ /^avg10=/) { split($i, a, "="); print a[2]; exit } }' "${file}"
}

sample_resources() {
  local now elapsed cg_rel cg mem_available swap_total swap_free swap_used
  local mem_current mem_peak cpu_usage io_read io_write load1
  now="$(date +%s)"
  elapsed="$((now - start_epoch))"
  cg_rel="$(awk -F: '$1 == "0" {print $3; exit}' /proc/self/cgroup)"
  cg="/sys/fs/cgroup${cg_rel}"
  load1="$(awk '{print $1}' /proc/loadavg)"
  mem_available="$(awk '$1 == "MemAvailable:" {print $2 * 1024}' /proc/meminfo)"
  swap_total="$(awk '$1 == "SwapTotal:" {print $2 * 1024}' /proc/meminfo)"
  swap_free="$(awk '$1 == "SwapFree:" {print $2 * 1024}' /proc/meminfo)"
  swap_used="$((swap_total - swap_free))"
  mem_current="$(cat "${cg}/memory.current" 2>/dev/null || printf '0')"
  mem_peak="$(cat "${cg}/memory.peak" 2>/dev/null || printf '0')"
  cpu_usage="$(awk '$1 == "usage_usec" {print $2}' "${cg}/cpu.stat" 2>/dev/null || printf '0')"
  io_read="$(awk '{for(i=1;i<=NF;i++) if($i ~ /^rbytes=/){split($i,a,"="); s+=a[2]}} END{printf "%.0f",s+0}' "${cg}/io.stat" 2>/dev/null || printf '0')"
  io_write="$(awk '{for(i=1;i<=NF;i++) if($i ~ /^wbytes=/){split($i,a,"="); s+=a[2]}} END{printf "%.0f",s+0}' "${cg}/io.stat" 2>/dev/null || printf '0')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${now}" "${elapsed}" "${load1}" "${mem_available}" "${swap_used}" \
    "${mem_current}" "${mem_peak}" "${cpu_usage}" "${io_read}" "${io_write}" \
    "$(psi_value /proc/pressure/cpu some)" \
    "$(psi_value /proc/pressure/memory some)" \
    "$(psi_value /proc/pressure/memory full)" \
    "$(psi_value /proc/pressure/io some)" \
    "$(psi_value /proc/pressure/io full)" >> "${TELEMETRY}"
}

monitor_resources() {
  while true; do
    sample_resources
    sleep 20
  done
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

printf 'Starting corrected Cellbit-negative Stage3B at %s\n' "$(date --iso-8601=seconds)"
printf 'threads=%s output=%s\n' "${THREADS}" "${OUT_ROOT}"

monitor_resources &
monitor_pid="$!"
cleanup_monitor() {
  kill "${monitor_pid}" 2>/dev/null || true
  wait "${monitor_pid}" 2>/dev/null || true
}
trap cleanup_monitor EXIT

set +e
/usr/bin/time -v -o "${TIME_FILE}" \
  "${TRACTOR}" \
    --manifest "${MANIFEST}" \
    --quality-manifest "${QUALITY}" \
    --marker-map "${MARKERS}" \
    --ani-engine "${ANI_ENGINE}" \
    --allow-experimental-ani-engine \
    --makeblastdb "${MAKEBLASTDB}" \
    --blastn "${BLASTN}" \
    --python "${PYTHON}" \
    --leiden-backend "${LEIDEN}" \
    --subass "${SUBASS}" \
    --flye-root "${FLYE_ROOT}" \
    --threads "${THREADS}" \
    --out "${OUT_ROOT}"
rc="$?"
set -e

sample_resources
cleanup_monitor
trap - EXIT

end_epoch="$(date +%s)"
elapsed="$((end_epoch - start_epoch))"

{
  printf 'component\tbytes\n'
  du -sb "${CODE_ROOT}" | awk '{print "code_source_build\t" $1}'
  for binary in "${TRACTOR}" "${ANI_ENGINE}" "${SUBASS}"; do
    printf 'binary:%s\t%s\n' "${binary}" "$(stat -c '%s' "${binary}")"
  done
  if [[ -e "${OUT_ROOT}" ]]; then
    du -sb "${OUT_ROOT}" | awk '{print "stage3b_output_total\t" $1}'
    for stage in "${OUT_ROOT}"/*; do
      [[ -e "${stage}" ]] || continue
      du -sb "${stage}" | awk -F '\t' '{n=$2; sub(/^.*\//,"",n); print "output:" n "\t" $1}'
    done
  fi
} > "${CONTROL_ROOT}/DISK_USAGE.tsv"

complete_sha=""
if [[ -f "${OUT_ROOT}/COMPLETE.json" ]]; then
  complete_sha="$(sha256sum "${OUT_ROOT}/COMPLETE.json" | awk '{print $1}')"
fi

cat > "${CONTROL_ROOT}/RUN_RESULT.json" <<EOF
{
  "schema": "cellbit-stage3b-resource-run-v1",
  "status": "$([[ ${rc} -eq 0 ]] && printf 'PASS' || printf 'FAIL')",
  "entry_rule": "Cellbit/Dna2bit negative && max_contig>=1000 && CheckM2 contamination<5",
  "reference_positive_exclusion": false,
  "threads": ${THREADS},
  "start_epoch": ${start_epoch},
  "end_epoch": ${end_epoch},
  "elapsed_seconds": ${elapsed},
  "exit_code": ${rc},
  "output_root": "${OUT_ROOT}",
  "complete_sha256": "${complete_sha}"
}
EOF

printf 'Finished corrected Stage3B rc=%s elapsed=%ss at %s\n' "${rc}" "${elapsed}" "$(date --iso-8601=seconds)"
exit "${rc}"
