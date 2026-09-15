#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C TZ=UTC

# Fresh, write-once DNA2bit-SAG Original Embedded Lake benchmark from the frozen raw-contig manifest
# through Stage3B bins.  The final CheckM2 pass is an external quality
# evaluation and is timed separately from the core pipeline.  The systemd unit
# must enforce exactly the same CPU quota as CELLBIT_THREADS; every explicit
# tool budget below is bounded by that value.

readonly THREADS="${CELLBIT_THREADS:-50}"
readonly CODE_ROOT="/home/data/fyc/dna2bit_sag_original_integrated_20260909"
readonly MANIFEST="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/10_current_tractor_baseline_continuation_v1_code/acceptance/lake_full_contigs_13742_v1/LAKE_FULL_CONTIGS.tsv"
readonly RUN_ROOT="${CELLBIT_RUN_ROOT:-/home/data/fyc/dna2bit_sag_original_integrated_runs_20260909/attempt_001}"
readonly CONTROL="${RUN_ROOT}/00_control"
readonly MAIN_OUT="${RUN_ROOT}/01_stage1_3a"
readonly INPUT_VIEW="${RUN_ROOT}/02_stage3b_input_view"
readonly PENDING_STATS="${RUN_ROOT}/03_stage3b_pending_stats"
readonly UPSTREAM="${RUN_ROOT}/04_stage3b_upstream"
readonly STAGE3B_INPUTS="${RUN_ROOT}/05_stage3b_inputs"
readonly STAGE3B_OUT="${RUN_ROOT}/06_stage3b_exact"
readonly FINAL_VIEW="${RUN_ROOT}/07_final_bin_view"
readonly FINAL_EVAL="${RUN_ROOT}/08_final_checkm2"

readonly PYTHON="/home/data/fyc/biosoft/miniconda3/bin/python"
readonly MAIN="${CODE_ROOT}/build_embedded_original/dna2bit-sag-pipeline"
readonly TRACTOR="${CODE_ROOT}/build_embedded_original/sag-stage3b-tractor"
readonly ANI_ENGINE="${CODE_ROOT}/build_embedded_original/gtdb-ani-af"
readonly SUBASS="${CODE_ROOT}/build_embedded_original/subass/cpp-subass"
readonly LEIDEN="${CODE_ROOT}/python/stage3b_signed_leiden.py"
readonly DNA_TAX="/home/data/shared/software/Dna2bit/GTDB232/genome_taxonomy_1.csv"
readonly DNA_PACKED_DB="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/dna2bit_gtdb232_packed_v1"
readonly CHECKM2_ENV="/home/data/fyc/biosoft/miniconda3/envs/checkm2_env"
readonly CHECKM2_SPEED_ROOT="/home/data/fyc/cellbit_114514/result/20_checkm2_speed_v2_balanced_20260905"
readonly CHECKM2_SPEED_V1_ROOT="/home/data/fyc/cellbit_114514/result/19_checkm2_speed_v1_20260905"
readonly CHECKM2="${CHECKM2_SPEED_ROOT}/checkm2-balanced-v2"
readonly CHECKM2_DB="/home/data/shared/software/CheckM2_database/uniref100.KO.1.dmnd"
readonly GTDBTK="/home/data/fyc/biosoft/miniconda3/envs/gtdbtk/bin/gtdbtk"
readonly GTDBTK_DATA="/home/data/temp/release232"
readonly MAKEBLASTDB="/home/data/fyc/biosoft/miniconda3/envs/decon_bench/bin/makeblastdb"
readonly BLASTN="/home/data/fyc/biosoft/miniconda3/envs/decon_bench/bin/blastn"
readonly FLYE_ROOT="/home/data/fyc/biosoft/miniconda3/envs/assemble"
readonly ORIGINAL_PATH="${PATH}"
readonly GLOBAL_RUN_LOCK="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/.dna2bit-sag-original-embedded-v1.lock"

readonly EXPECTED_MANIFEST_SHA256="ef1376f9d4b7ce2c8eee1bea2f2962d6bcf65a74c8c91e73681dc57830ea0f07"
readonly EXPECTED_MAIN_SHA256="735f66cf351ef80a461012fbbfed5ff306d6075303645d467919b14157fad5fa"
readonly EXPECTED_TRACTOR_SHA256="b28bb6ab405d31de8322191bc8450095fa5e8c4217169b042c0f1ad0570e910e"
readonly EXPECTED_ANI_ENGINE_SHA256="c5455124238369a45ac7b478f27f0856de166f8c02d8b9f3bef2670e9e2ba9f0"
readonly EXPECTED_SUBASS_SHA256="fe860c40972aa8c04fbda78aa32252fcdd2a90d367b2796d3ff74359519941d6"
readonly EXPECTED_LEIDEN_SHA256="3592e4e6ceabd3bbd4f0ff9a8bd05e7b714d85731aa54d2adc185b31254709b1"
readonly EXPECTED_DNA_PACKED_COMPLETE_SHA256="46bdd12d83d3552e5c82ea15d0c82416088dc7452f2edb2ba053341a13696ffd"
readonly EXPECTED_CHECKM2_WRAPPER_SHA256="e8d470d1fa195f40899c71a2269a28ac39ea116a359b168fbc600edbdece7239"
readonly EXPECTED_CHECKM2_INIT_OVERLAY_SHA256="eb3841498d9784791f14981cd9260ef592ac1607957a539279fae51a80b708c0"
readonly EXPECTED_CHECKM2_DIAMOND_OVERLAY_SHA256="13f714ea8d5dbf1a5bbcd09a53f9971932ed58a9b95fd733f5e54a04b6851cbb"
readonly EXPECTED_CHECKM2_V1_DIAMOND_OVERLAY_SHA256="e8b86e9a87b0fc5c004ad97f2dac091a438ec5624a3c598c1d6bae0c2535e0d7"
readonly EXPECTED_CHECKM2_PREDICT_OVERLAY_SHA256="892540d09db8663472345623c0e6ef3830a676a874e42a7d6f61b7404ea70617"
readonly EXPECTED_PREPARE_UPSTREAM_SHA256="82104547e71b3a533b01020f9561572e9956b512b25a5b4f7a432e159bbc6f4a"
readonly EXPECTED_PENDING_STATS_SHA256="770637b20f5b50dc269efe04ea6578138ee167c328a3ef9861d66b90c6f5bccb"
readonly EXPECTED_FINALIZE_UPSTREAM_SHA256="fa8ccf6d92e018206aac8263690fe9cbfc2eda56a810b135ded9fb17b16fe5dd"
readonly EXPECTED_FINAL_VIEW_SHA256="a974bd0f79a641e162142bca0a95bff45d7e4e01217e7ce35adade806c274445"
readonly EXPECTED_CHECKM2_DB_SHA256="aff3badcd2bd389539e833b5a0d446c13fb257d71abcb6c73fac38e8ab070acb"
readonly EXPECTED_CHECKM2_ENTRYPOINT_SHA256="e9409a71dea68fb43c7bea3337caec1d224610197f374980d8b358376d15318d"
readonly EXPECTED_GTDBTK_ENTRYPOINT_SHA256="0baa6b1732c2689932ea9dfb2a3342c7c318aaf9bafc2e0f35f5c6bb7c9e3369"
readonly EXPECTED_GTDBTK_METADATA_SHA256="028481e18d6e7dddd11c126047010cfd1157f5c148277bf898c51801fc5373b8"

[[ "${THREADS}" =~ ^[1-9][0-9]*$ && "${THREADS}" -le 220 ]] || {
  printf 'CELLBIT_THREADS must be an integer in [1,220], got %s\n' "${THREADS}" >&2
  exit 64
}
[[ "${RUN_ROOT}" == /* ]] || { printf 'CELLBIT_RUN_ROOT must be absolute\n' >&2; exit 64; }

for required in \
  "${MANIFEST}" "${PYTHON}" "${MAIN}" "${TRACTOR}" "${ANI_ENGINE}" \
  "${SUBASS}" "${LEIDEN}" "${DNA_TAX}" \
  "${DNA_PACKED_DB}/COMPLETE.json" \
  "${CHECKM2_SPEED_ROOT}/checkm2/__init__.py" \
  "${CHECKM2_SPEED_ROOT}/checkm2/diamond.py" "${CHECKM2_SPEED_ROOT}/checkm2/predictQuality.py" \
  "${CHECKM2_SPEED_V1_ROOT}/checkm2/diamond.py" \
  "${CHECKM2}" "${CHECKM2_ENV}/bin/checkm2" "${CHECKM2_DB}" \
  "${GTDBTK}" "${GTDBTK_DATA}" "${GTDBTK_DATA}/metadata/metadata.txt" \
  "${MAKEBLASTDB}" "${BLASTN}" \
  "${CODE_ROOT}/tools/prepare_stage3b_upstream_input.py" \
  "${CODE_ROOT}/tools/prepare_stage3b_pending_stats.py" \
  "${CODE_ROOT}/tools/finalize_stage3b_upstream.py" \
  "${CODE_ROOT}/tools/prepare_current_baseline_bin_evaluation.py"; do
  [[ -e "${required}" ]] || { printf 'missing required path: %s\n' "${required}" >&2; exit 66; }
done

require_sha256() {
  local path="$1" expected="$2" role="$3" actual
  actual="$(sha256sum "${path}" | awk '{print $1}')"
  [[ "${actual}" == "${expected}" ]] || {
    printf 'frozen authority mismatch for %s: expected=%s actual=%s path=%s\n' \
      "${role}" "${expected}" "${actual}" "${path}" >&2
    exit 65
  }
}

assert_frozen_authority() {
  local mode="${1:-fast}"
  require_sha256 "${MANIFEST}" "${EXPECTED_MANIFEST_SHA256}" manifest
  require_sha256 "${MAIN}" "${EXPECTED_MAIN_SHA256}" main_binary
  require_sha256 "${TRACTOR}" "${EXPECTED_TRACTOR_SHA256}" stage3b_tractor
  require_sha256 "${ANI_ENGINE}" "${EXPECTED_ANI_ENGINE_SHA256}" ani_engine
  require_sha256 "${SUBASS}" "${EXPECTED_SUBASS_SHA256}" subassemble_binary
  require_sha256 "${LEIDEN}" "${EXPECTED_LEIDEN_SHA256}" leiden_backend
  require_sha256 "${DNA_PACKED_DB}/COMPLETE.json" "${EXPECTED_DNA_PACKED_COMPLETE_SHA256}" dna2bit_packed_database_receipt
  require_sha256 "${CHECKM2}" "${EXPECTED_CHECKM2_WRAPPER_SHA256}" checkm2_speed_wrapper
  require_sha256 "${CHECKM2_SPEED_ROOT}/checkm2/__init__.py" "${EXPECTED_CHECKM2_INIT_OVERLAY_SHA256}" checkm2_namespace_overlay
  require_sha256 "${CHECKM2_SPEED_ROOT}/checkm2/diamond.py" "${EXPECTED_CHECKM2_DIAMOND_OVERLAY_SHA256}" checkm2_diamond_overlay
  require_sha256 "${CHECKM2_SPEED_V1_ROOT}/checkm2/diamond.py" "${EXPECTED_CHECKM2_V1_DIAMOND_OVERLAY_SHA256}" checkm2_base_diamond_overlay
  require_sha256 "${CHECKM2_SPEED_ROOT}/checkm2/predictQuality.py" "${EXPECTED_CHECKM2_PREDICT_OVERLAY_SHA256}" checkm2_prediction_overlay
  require_sha256 "${CODE_ROOT}/tools/prepare_stage3b_upstream_input.py" "${EXPECTED_PREPARE_UPSTREAM_SHA256}" stage3b_input_view_helper
  require_sha256 "${CODE_ROOT}/tools/prepare_stage3b_pending_stats.py" "${EXPECTED_PENDING_STATS_SHA256}" stage3b_pending_stats_helper
  require_sha256 "${CODE_ROOT}/tools/finalize_stage3b_upstream.py" "${EXPECTED_FINALIZE_UPSTREAM_SHA256}" stage3b_finalize_helper
  require_sha256 "${CODE_ROOT}/tools/prepare_current_baseline_bin_evaluation.py" "${EXPECTED_FINAL_VIEW_SHA256}" final_view_helper
  require_sha256 "${CHECKM2_ENV}/bin/checkm2" "${EXPECTED_CHECKM2_ENTRYPOINT_SHA256}" checkm2_entrypoint
  require_sha256 "${GTDBTK}" "${EXPECTED_GTDBTK_ENTRYPOINT_SHA256}" gtdbtk_entrypoint
  require_sha256 "${GTDBTK_DATA}/metadata/metadata.txt" "${EXPECTED_GTDBTK_METADATA_SHA256}" gtdbtk_release_metadata
  if [[ "${mode}" == "full" ]]; then
    require_sha256 "${CHECKM2_DB}" "${EXPECTED_CHECKM2_DB_SHA256}" checkm2_database
  elif [[ "${mode}" != "fast" ]]; then
    printf 'unknown authority-check mode: %s\n' "${mode}" >&2
    exit 64
  fi
}

# Validate the frozen scientific authorities, evaluator version, and the real
# cgroup limit before claiming the write-once output path.
assert_frozen_authority full
readonly CHECKM2_VERSION="$({ PATH="${CHECKM2_ENV}/bin:${ORIGINAL_PATH}" "${CHECKM2}" --version; } | tr -d '\r\n')"
[[ "${CHECKM2_VERSION}" == "1.0.1" ]] || {
  printf 'unexpected CheckM2 version before science: %s\n' "${CHECKM2_VERSION}" >&2
  exit 65
}
readonly CGROUP_REL="$(awk -F: '$1 == "0" {print $3; exit}' /proc/self/cgroup)"
readonly CGROUP_ROOT="/sys/fs/cgroup${CGROUP_REL}"
[[ -r "${CGROUP_ROOT}/cpu.max" ]] || { printf 'cannot read cgroup cpu.max\n' >&2; exit 65; }
read -r CPU_QUOTA_US CPU_PERIOD_US < "${CGROUP_ROOT}/cpu.max"
readonly CPU_QUOTA_US CPU_PERIOD_US
[[ "${CPU_QUOTA_US}" =~ ^[0-9]+$ && "${CPU_PERIOD_US}" =~ ^[0-9]+$ && "${CPU_PERIOD_US}" -gt 0 ]] || {
  printf 'finite CPU quota is not enforced: cpu.max=%s %s\n' "${CPU_QUOTA_US}" "${CPU_PERIOD_US}" >&2
  exit 65
}
(( CPU_QUOTA_US == THREADS * CPU_PERIOD_US )) || {
  printf 'cgroup quota is not exactly %s cores: cpu.max=%s %s\n' \
    "${THREADS}" "${CPU_QUOTA_US}" "${CPU_PERIOD_US}" >&2
  exit 65
}
readonly ACTUAL_CPU_QUOTA_CORES="$(awk -v q="${CPU_QUOTA_US}" -v p="${CPU_PERIOD_US}" 'BEGIN{printf "%.6f",q/p}')"
CPUSET_EFFECTIVE="unavailable_not_delegated"
if [[ -r "${CGROUP_ROOT}/cpuset.cpus.effective" ]]; then
  observed_cpuset="$(tr -d '\r\n' < "${CGROUP_ROOT}/cpuset.cpus.effective")"
  [[ -z "${observed_cpuset}" ]] || CPUSET_EFFECTIVE="${observed_cpuset}"
fi
readonly CPUSET_EFFECTIVE
readonly PROCESS_CPU_AFFINITY="$(awk '/^Cpus_allowed_list:/ {sub(/^[^:]*:[[:space:]]*/, ""); print; exit}' /proc/self/status)"
readonly MEMORY_HIGH_BYTES="$(tr -d '\r\n' < "${CGROUP_ROOT}/memory.high")"
readonly MEMORY_MAX_BYTES="$(tr -d '\r\n' < "${CGROUP_ROOT}/memory.max")"
readonly MEMORY_SWAP_MAX_BYTES="$(tr -d '\r\n' < "${CGROUP_ROOT}/memory.swap.max")"
readonly PIDS_MAX="$(tr -d '\r\n' < "${CGROUP_ROOT}/pids.max")"
readonly NOFILE_SOFT="$(ulimit -Sn)"
readonly NOFILE_HARD="$(ulimit -Hn)"
[[ "${PROCESS_CPU_AFFINITY}" == "0-239" ]] || {
  printf 'process CPU affinity is not exactly 0-239: %s\n' "${PROCESS_CPU_AFFINITY}" >&2
  exit 65
}
[[ "${CPUSET_EFFECTIVE}" == "unavailable_not_delegated" || "${CPUSET_EFFECTIVE}" == "0-239" ]] || {
  printf 'unexpected delegated cpuset: %s\n' "${CPUSET_EFFECTIVE}" >&2
  exit 65
}
[[ "${MEMORY_HIGH_BYTES}" == "1288490188800" ]] || { printf 'unexpected memory.high: %s\n' "${MEMORY_HIGH_BYTES}" >&2; exit 65; }
[[ "${MEMORY_MAX_BYTES}" == "1556925644800" ]] || { printf 'unexpected memory.max: %s\n' "${MEMORY_MAX_BYTES}" >&2; exit 65; }
[[ "${MEMORY_SWAP_MAX_BYTES}" == "8589934592" ]] || { printf 'unexpected memory.swap.max: %s\n' "${MEMORY_SWAP_MAX_BYTES}" >&2; exit 65; }
[[ "${PIDS_MAX}" == "max" ]] || { printf 'unexpected pids.max: %s\n' "${PIDS_MAX}" >&2; exit 65; }
[[ "${NOFILE_SOFT}" == "1048576" && "${NOFILE_HARD}" == "1048576" ]] || {
  printf 'unexpected RLIMIT_NOFILE soft=%s hard=%s\n' "${NOFILE_SOFT}" "${NOFILE_HARD}" >&2
  exit 65
}

exec 9> "${GLOBAL_RUN_LOCK}"
flock -n 9 || { printf 'another full E2E runner holds %s\n' "${GLOBAL_RUN_LOCK}" >&2; exit 75; }

[[ ! -e "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  printf 'write-once run root exists: %s\n' "${RUN_ROOT}" >&2; exit 73;
}
mkdir -p "$(dirname "${RUN_ROOT}")"
mkdir "${RUN_ROOT}" || { printf 'atomic write-once claim failed: %s\n' "${RUN_ROOT}" >&2; exit 73; }
mkdir "${CONTROL}"
umask 022
readonly STDOUT_LOG="${CONTROL}/RUN_STDOUT.log"
readonly STDERR_LOG="${CONTROL}/RUN_STDERR.log"
readonly TELEMETRY="${CONTROL}/RESOURCE_TELEMETRY.tsv"
readonly PHASE_FILE="${CONTROL}/CURRENT_PHASE.txt"
readonly PHASE_EVENTS="${CONTROL}/PHASE_EVENTS.tsv"
readonly RESULT="${CONTROL}/RUN_RESULT.json"
exec > >(tee -a "${STDOUT_LOG}") 2> >(tee -a "${STDERR_LOG}" >&2)

readonly START_EPOCH="$(date +%s)"
core_start_epoch=0
pipeline_end_epoch=0
monitor_pid=""
current_phase="preflight"

set_phase() {
  current_phase="$1"
  printf '%s\n' "${current_phase}" > "${PHASE_FILE}.tmp"
  mv "${PHASE_FILE}.tmp" "${PHASE_FILE}"
  printf '%s\t%s\n' "$(date +%s)" "${current_phase}" >> "${PHASE_EVENTS}"
}

psi_value() {
  local file="$1" row="$2"
  awk -v wanted="${row}" '$1 == wanted {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]; exit}}' "${file}" 2>/dev/null || printf 'NA'
}

sample_resources() {
  local now elapsed cg_rel cg mem_available swap_total swap_free swap_used
  local mem_current mem_peak cpu_usage io_read io_write phase
  now="$(date +%s)"; elapsed="$((now - START_EPOCH))"
  phase="$(cat "${PHASE_FILE}" 2>/dev/null || printf '%s' "${current_phase}")"
  cg_rel="$(awk -F: '$1 == "0" {print $3; exit}' /proc/self/cgroup)"; cg="/sys/fs/cgroup${cg_rel}"
  mem_available="$(awk '$1 == "MemAvailable:" {print $2 * 1024}' /proc/meminfo)"
  swap_total="$(awk '$1 == "SwapTotal:" {print $2 * 1024}' /proc/meminfo)"
  swap_free="$(awk '$1 == "SwapFree:" {print $2 * 1024}' /proc/meminfo)"; swap_used="$((swap_total - swap_free))"
  mem_current="$(cat "${cg}/memory.current" 2>/dev/null || printf '0')"
  mem_peak="$(cat "${cg}/memory.peak" 2>/dev/null || printf '0')"
  cpu_usage="$(awk '$1 == "usage_usec" {print $2}' "${cg}/cpu.stat" 2>/dev/null || printf '0')"
  io_read="$(awk '{for(i=1;i<=NF;i++) if($i~/^rbytes=/){split($i,a,"=");s+=a[2]}} END{printf "%.0f",s+0}' "${cg}/io.stat" 2>/dev/null || printf '0')"
  io_write="$(awk '{for(i=1;i<=NF;i++) if($i~/^wbytes=/){split($i,a,"=");s+=a[2]}} END{printf "%.0f",s+0}' "${cg}/io.stat" 2>/dev/null || printf '0')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${now}" "${elapsed}" "${phase}" "$(awk '{print $1}' /proc/loadavg)" \
    "${mem_available}" "${swap_used}" "${mem_current}" "${mem_peak}" "${cpu_usage}" \
    "${io_read}" "${io_write}" "$(psi_value /proc/pressure/cpu some)" \
    "$(psi_value /proc/pressure/memory some)" "$(psi_value /proc/pressure/memory full)" \
    "$(psi_value /proc/pressure/io some)" "$(psi_value /proc/pressure/io full)" >> "${TELEMETRY}"
}

monitor_resources() {
  while true; do sample_resources; sleep 20; done
}

finish() {
  local rc="$?" end_epoch elapsed core_elapsed complete_sha final_hq="null" final_mq="null"
  if [[ -n "${monitor_pid}" ]]; then kill "${monitor_pid}" 2>/dev/null || true; wait "${monitor_pid}" 2>/dev/null || true; fi
  sample_resources || true
  end_epoch="$(date +%s)"; elapsed="$((end_epoch - START_EPOCH))"
  core_elapsed="null"
  [[ "${pipeline_end_epoch}" -gt 0 && "${core_start_epoch}" -gt 0 ]] && \
    core_elapsed="$((pipeline_end_epoch - core_start_epoch))"
  complete_sha=""; [[ -f "${STAGE3B_OUT}/COMPLETE.json" ]] && complete_sha="$(sha256sum "${STAGE3B_OUT}/COMPLETE.json" | awk '{print $1}')"
  if [[ -f "${FINAL_EVAL}/RESULT.json" ]]; then
    final_hq="$(${PYTHON} -c 'import json,sys; print(json.load(open(sys.argv[1]))["hq"])' "${FINAL_EVAL}/RESULT.json" 2>/dev/null || printf 'null')"
    final_mq="$(${PYTHON} -c 'import json,sys; print(json.load(open(sys.argv[1]))["mq_only"])' "${FINAL_EVAL}/RESULT.json" 2>/dev/null || printf 'null')"
  fi
  du -sb "${RUN_ROOT}"/* 2>/dev/null | sort -k2,2 > "${CONTROL}/DISK_USAGE.tsv" || true
  cat > "${RESULT}.tmp" <<EOF
{
  "schema": "cellbit-lake-full-e2e-speed-v1",
  "status": "$([[ ${rc} -eq 0 ]] && printf 'PASS' || printf 'FAIL')",
  "exit_code": ${rc},
  "threads_requested": ${THREADS},
  "systemd_cpu_quota_cores": ${ACTUAL_CPU_QUOTA_CORES},
  "cgroup_cpu_max_quota_us": ${CPU_QUOTA_US},
  "cgroup_cpu_max_period_us": ${CPU_PERIOD_US},
  "application_cache": "cold_new_write_once_outputs_no_resume",
  "filesystem_page_cache": "natural_warm_uncontrolled_no_drop_caches",
  "prior_scientific_artifact_reuse": false,
  "start_epoch": ${START_EPOCH},
  "core_start_epoch": ${core_start_epoch},
  "pipeline_end_epoch": ${pipeline_end_epoch},
  "end_epoch": ${end_epoch},
  "core_pipeline_seconds": ${core_elapsed},
  "including_final_evaluation_seconds": ${elapsed},
  "stage3b_complete_sha256": "${complete_sha}",
  "final_hq": ${final_hq},
  "final_mq_only": ${final_mq},
  "run_root": "${RUN_ROOT}"
}
EOF
  mv "${RESULT}.tmp" "${RESULT}"
  {
    for evidence in \
      "${RESULT}" "${PHASE_EVENTS}" "${TELEMETRY}" \
      "${CONTROL}/CODE_AND_INPUT_AUTHORITY.tsv" \
      "${CONTROL}/CACHE_CONDITION.json" "${CONTROL}/CGROUP_CPU_QUOTA.json" \
      "${CONTROL}/FINAL_SHA256SUMS.txt" \
      "${MAIN_OUT}/COMPLETE.json" "${STAGE3B_INPUTS}/COMPLETE.json" \
      "${STAGE3B_OUT}/COMPLETE.json" "${FINAL_VIEW}/COMPLETE.json" \
      "${FINAL_EVAL}/RESULT.json"; do
      [[ -f "${evidence}" ]] && sha256sum "${evidence}"
    done
    find "${CONTROL}" -maxdepth 1 -type f -name '*_TIME.txt' -print0 2>/dev/null | \
      sort -z | xargs -0 -r sha256sum
  } | sort -k2,2 > "${CONTROL}/FINAL_CONTROL_SHA256SUMS.txt.tmp"
  mv "${CONTROL}/FINAL_CONTROL_SHA256SUMS.txt.tmp" "${CONTROL}/FINAL_CONTROL_SHA256SUMS.txt"
  printf 'Finished %s-thread speed-v1 E2E rc=%s core=%s total=%ss at %s\n' "${THREADS}" "${rc}" "${core_elapsed}" "${elapsed}" "$(date --iso-8601=seconds)"
}
trap finish EXIT

printf 'timestamp_epoch\telapsed_s\tphase\tload1\tmem_available_bytes\tswap_used_bytes\tcgroup_memory_current\tcgroup_memory_peak\tcgroup_cpu_usage_usec\tcgroup_io_read_bytes\tcgroup_io_write_bytes\tcpu_psi_some_avg10\tmemory_psi_some_avg10\tmemory_psi_full_avg10\tio_psi_some_avg10\tio_psi_full_avg10\n' > "${TELEMETRY}"
printf 'timestamp_epoch\tphase\n' > "${PHASE_EVENTS}"
set_phase preflight_authority

cat > "${CONTROL}/CACHE_CONDITION.json" <<EOF
{
  "schema": "cellbit-benchmark-cache-condition-v1",
  "application_cache": "cold_new_write_once_outputs_no_resume",
  "filesystem_page_cache": "natural_warm_uncontrolled_no_drop_caches",
  "reason_no_drop_caches": "avoid perturbing other server workloads",
  "prior_scientific_artifact_reuse": false,
  "within_run_stage_to_stage_cache": "natural_and_inherent",
  "comparison_note": "attempt_003 also used fresh output after an earlier run and was likely page-cache warm"
}
EOF
printf '%s\n' "${CHECKM2_VERSION}" > "${CONTROL}/CHECKM2_VERSION_PREFLIGHT.txt"
cat > "${CONTROL}/CGROUP_CPU_QUOTA.json" <<EOF
{
  "cpu_max_quota_us": ${CPU_QUOTA_US},
  "cpu_max_period_us": ${CPU_PERIOD_US},
  "actual_cpu_quota_cores": ${ACTUAL_CPU_QUOTA_CORES},
  "required_cpu_quota_cores": ${THREADS},
  "cgroup_relative_path": "${CGROUP_REL}",
  "cpuset_cpus_effective": "${CPUSET_EFFECTIVE}",
  "process_cpu_affinity": "${PROCESS_CPU_AFFINITY}",
  "cpu_affinity_enforcement": "taskset_0-239_with_process_status_verification",
  "memory_high_bytes": "${MEMORY_HIGH_BYTES}",
  "memory_max_bytes": "${MEMORY_MAX_BYTES}",
  "memory_swap_max_bytes": "${MEMORY_SWAP_MAX_BYTES}",
  "pids_max": "${PIDS_MAX}",
  "nofile_soft": "${NOFILE_SOFT}",
  "nofile_hard": "${NOFILE_HARD}",
  "global_run_lock": "${GLOBAL_RUN_LOCK}"
}
EOF

{
  printf 'artifact\tbytes\tsha256\tpath\n'
  for path in \
    "${MANIFEST}" "${CODE_ROOT}/src/main.cpp" "${CODE_ROOT}/src/dna2bit_embedded.cpp" "${CODE_ROOT}/src/dna2bit_packed_api.cpp" "${CODE_ROOT}/src/stage3b.cpp" \
    "${CODE_ROOT}/python/stage3b_signed_leiden.py" \
    "${CODE_ROOT}/tools/run_lake_full_e2e_50t_embedded_original.sh" \
    "${CODE_ROOT}/tools/prepare_stage3b_upstream_input.py" \
    "${CODE_ROOT}/tools/prepare_stage3b_pending_stats.py" \
    "${CODE_ROOT}/tools/finalize_stage3b_upstream.py" \
    "${CODE_ROOT}/tools/prepare_current_baseline_bin_evaluation.py" \
    "${MAIN}" "${TRACTOR}" "${ANI_ENGINE}" "${SUBASS}" "${LEIDEN}" \
    "${DNA_TAX}" "${DNA_PACKED_DB}/COMPLETE.json" \
    "${CHECKM2}" "${CHECKM2_SPEED_ROOT}/checkm2/__init__.py" \
    "${CHECKM2_SPEED_ROOT}/checkm2/diamond.py" \
    "${CHECKM2_SPEED_V1_ROOT}/checkm2/diamond.py" \
    "${CHECKM2_SPEED_ROOT}/checkm2/predictQuality.py" \
    "${GTDBTK}" "${MAKEBLASTDB}" "${BLASTN}"; do
    printf '%s\t%s\t%s\t%s\n' "$(basename "${path}")" "$(stat -c '%s' "${path}")" "$(sha256sum "${path}" | awk '{print $1}')" "${path}"
  done
} > "${CONTROL}/CODE_AND_INPUT_AUTHORITY.tsv"
{
  printf 'component\tbytes\tmeasurement\tpath\n'
  printf 'code_root\t%s\tlive_du_sb\t%s\n' "$(du -sb "${CODE_ROOT}" | awk '{print $1}')" "${CODE_ROOT}"
  printf 'dna2bit_packed_database\t%s\tlive_du_sb\t%s\n' "$(du -sb "${DNA_PACKED_DB}" | awk '{print $1}')" "${DNA_PACKED_DB}"
  printf 'checkm2_database\t3082500605\tprior_audited_du_sb\t%s\n' "${CHECKM2_DB}"
  printf 'gtdbtk_r232_database\t100482230310\tprior_audited_du_sb\t%s\n' "${GTDBTK_DATA}"
} > "${CONTROL}/SOFTWARE_AND_DATABASE_FOOTPRINT.tsv"
cat /proc/meminfo > "${CONTROL}/MEMINFO_AT_START.txt"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CHECKM2_SPEED_MODEL_CHUNK=1000
monitor_resources & monitor_pid="$!"

assert_frozen_authority
core_start_epoch="$(date +%s)"
printf '%s\n' "${core_start_epoch}" > "${CONTROL}/CORE_PIPELINE_START_EPOCH.txt"
set_phase stage1_3a
/usr/bin/time -v -o "${CONTROL}/01_STAGE1_3A_TIME.txt" \
  "${MAIN}" --manifest "${MANIFEST}" --out "${MAIN_OUT}" --threads "${THREADS}" --memory-gb 1200 \
  --dna-search-engine packed --dna-packed-db "${DNA_PACKED_DB}" \
  --subass "${SUBASS}" --flye-root "${FLYE_ROOT}"
"${PYTHON}" - "${MANIFEST}" "${MAIN_OUT}" <<'PY'
import csv, json, pathlib, sys
manifest, root = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
complete = json.loads((root / "COMPLETE.json").read_text(encoding="utf-8"))
expected = {
    "status": "PASS", "dna_search_engine": "packed", "input_sags": 13742,
    "contig_inputs": 13742, "eligible_sags": 13742,
    "excluded_lt1000bp": 0, "labeled": 7527, "groups": 274,
}
for key, value in expected.items():
    if complete.get(key) != value:
        raise SystemExit(f"Stage1-3A COMPLETE mismatch: {key}={complete.get(key)!r} != {value!r}")
def byte_rows(path, columns):
    lines = [line[:-1] if line.endswith(b"\r") else line
             for line in path.read_bytes().split(b"\n")]
    if not lines or lines[0].split(b"\t") != columns:
        raise SystemExit(f"unexpected TSV header: {path}")
    rows = [line.split(b"\t") for line in lines[1:] if line]
    if any(len(row) != len(columns) for row in rows):
        raise SystemExit(f"wrong TSV width: {path}")
    return rows
manifest_rows = byte_rows(manifest, [b"sag_id", b"assembly_fasta"])
label_rows = byte_rows(root / "02_dna2bit/labels.tsv", [b"sag_id", b"reference", b"taxonomy", b"species_group"])
pending_rows = byte_rows(root / "03B_unclassified_pending.tsv", [b"sag_id", b"assembly_fasta", b"reason"])
manifest_ids = {row[0] for row in manifest_rows}
label_ids = {row[0] for row in label_rows}
pending_ids = {row[0] for row in pending_rows}
if len(manifest_rows) != 13742 or len(manifest_ids) != 13742:
    raise SystemExit("manifest count/uniqueness mismatch")
if len(label_rows) != 7527 or len(label_ids) != 7527:
    raise SystemExit("labels count/uniqueness mismatch")
if len(pending_rows) != 6215 or len(pending_ids) != 6215:
    raise SystemExit("pending count/uniqueness mismatch")
if label_ids & pending_ids or label_ids | pending_ids != manifest_ids:
    raise SystemExit("labels/pending do not form an exact manifest partition")
with (root / "03A_subassemble/groups.tsv").open(newline="", encoding="utf-8") as handle:
    groups = list(csv.DictReader(handle, delimiter="\t"))
if len(groups) != 274 or sum(int(row["sag_count"]) for row in groups) != 7527:
    raise SystemExit("Stage3A group count/member total mismatch")
if {row["species_group"].encode() for row in groups} != {row[3].replace(b"\r", b"") for row in label_rows}:
    raise SystemExit("Stage3A group IDs differ from label groups")
for row in groups:
    output = pathlib.Path(row["bin_fasta"])
    if not output.is_file() or output.stat().st_size == 0:
        raise SystemExit(f"missing/empty Stage3A bin: {output}")
print("PASS Stage1-3A fixed lake closure 13742=7527+6215 groups=274")
PY

set_phase stage3b_input_view
/usr/bin/time -v -o "${CONTROL}/02_INPUT_VIEW_TIME.txt" \
  "${PYTHON}" "${CODE_ROOT}/tools/prepare_stage3b_upstream_input.py" \
  --pending "${MAIN_OUT}/03B_unclassified_pending.tsv" --out-dir "${INPUT_VIEW}"
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("schema")=="stage3b-upstream-input-view-v1" and d.get("status")=="PASS" and d.get("sag_count")==6215 and d.get("unique_sag_ids")==6215 and d.get("nonempty_resolved_targets")==6215' "${INPUT_VIEW}/COMPLETE.json"

set_phase stage3b_pending_stats
/usr/bin/time -v -o "${CONTROL}/03_PENDING_STATS_TIME.txt" \
  "${PYTHON}" "${CODE_ROOT}/tools/prepare_stage3b_pending_stats.py" \
  --pending "${MAIN_OUT}/03B_unclassified_pending.tsv" \
  --stats "${MAIN_OUT}/STAGE3B_FASTA_STATS.tsv" --out-dir "${PENDING_STATS}"
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("schema")=="stage3b-pending-fasta-stats-view-v1" and d.get("status")=="PASS" and d.get("output_stats",{}).get("rows")==6215 and d.get("ignored_non_pending_rows")==7527' "${PENDING_STATS}/COMPLETE.json"

mkdir "${UPSTREAM}"
export PATH="${CHECKM2_ENV}/bin:${ORIGINAL_PATH}"
set_phase stage3b_upstream_checkm2
/usr/bin/time -v -o "${CONTROL}/04_CHECKM2_TIME.txt" \
  "${CHECKM2}" predict --input "${INPUT_VIEW}/input_fna" \
  --output-directory "${UPSTREAM}/checkm2_attempt_001" --database_path "${CHECKM2_DB}" \
  --extension fna --threads "${THREADS}" --force \
  >"${UPSTREAM}/CHECKM2_STDOUT.log" 2>"${UPSTREAM}/CHECKM2_STDERR.log"
[[ -s "${UPSTREAM}/checkm2_attempt_001/quality_report.tsv" ]]
"${PYTHON}" - "${INPUT_VIEW}/INPUT_MANIFEST.tsv" "${UPSTREAM}/checkm2_attempt_001/quality_report.tsv" <<'PY'
import csv, pathlib, sys
manifest, report = map(pathlib.Path, sys.argv[1:])
with manifest.open(newline="", encoding="utf-8") as handle:
    input_ids = [row["sag_id"] for row in csv.DictReader(handle, delimiter="\t")]
with report.open(newline="", encoding="utf-8") as handle:
    quality_ids = [row["Name"] for row in csv.DictReader(handle, delimiter="\t")]
if len(input_ids) != 6215 or len(set(input_ids)) != 6215:
    raise SystemExit("upstream input manifest count/ID mismatch")
if len(quality_ids) != 6215 or len(set(quality_ids)) != 6215 or set(quality_ids) != set(input_ids):
    raise SystemExit("CheckM2 report does not exactly cover 6215 input SAG IDs")
PY

set_phase stage3b_upstream_gtdbtk_identify
export GTDBTK_DATA_PATH="${GTDBTK_DATA}"
export PATH="/home/data/fyc/biosoft/miniconda3/envs/gtdbtk/bin:/home/data/fyc/biosoft/miniconda3/envs/eggnog/bin:${ORIGINAL_PATH}"
/usr/bin/time -v -o "${CONTROL}/05_GTDBTK_TIME.txt" \
  "${GTDBTK}" identify --genome_dir "${INPUT_VIEW}/input_fna" \
  --out_dir "${UPSTREAM}/gtdbtk_identify_attempt_001" -x fna --cpus "${THREADS}" \
  --force --write_single_copy_genes \
  >"${UPSTREAM}/GTDBTK_STDOUT.log" 2>"${UPSTREAM}/GTDBTK_STDERR.log"
[[ -s "${UPSTREAM}/gtdbtk_identify_attempt_001/identify/gtdbtk.bac120.markers_summary.tsv" ]]
printf 'PASS\ncheckm2_threads=%s\ngtdbtk_cpus=%s\n' "${THREADS}" "${THREADS}" > "${UPSTREAM}/UPSTREAM.PASS"

set_phase stage3b_finalize_inputs
/usr/bin/time -v -o "${CONTROL}/06_FINALIZE_INPUTS_TIME.txt" \
  "${PYTHON}" "${CODE_ROOT}/tools/finalize_stage3b_upstream.py" \
  --pending "${MAIN_OUT}/03B_unclassified_pending.tsv" \
  --stats "${PENDING_STATS}/STAGE3B_FASTA_STATS.tsv" --input-view "${INPUT_VIEW}" \
  --checkm2-report "${UPSTREAM}/checkm2_attempt_001/quality_report.tsv" \
  --identify-root "${UPSTREAM}/gtdbtk_identify_attempt_001" --out-dir "${STAGE3B_INPUTS}"
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("schema")=="stage3b-upstream-finalization-complete-v1" and d.get("status")=="PASS" and d.get("sag_count")==6215 and d.get("quality_rows")==6215' "${STAGE3B_INPUTS}/COMPLETE.json"

set_phase stage3b_exact_and_subassemble
assert_frozen_authority
/usr/bin/time -v -o "${CONTROL}/07_STAGE3B_TIME.txt" \
  "${TRACTOR}" --manifest "${MAIN_OUT}/03B_unclassified_pending.tsv" \
  --quality-manifest "${STAGE3B_INPUTS}/stage3b_quality.tsv" \
  --marker-map "${STAGE3B_INPUTS}/bac120_marker_nt_map.tsv" \
  --ani-engine "${ANI_ENGINE}" --allow-experimental-ani-engine \
  --makeblastdb "${MAKEBLASTDB}" --blastn "${BLASTN}" --python "${PYTHON}" \
  --leiden-backend "${LEIDEN}" --subass "${SUBASS}" --flye-root "${FLYE_ROOT}" \
  --threads "${THREADS}" --out "${STAGE3B_OUT}"
"${PYTHON}" - "${STAGE3B_OUT}" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
done = json.loads((root / "COMPLETE.json").read_text(encoding="utf-8"))
expected = {
    "schema": "sag-stage3b-complete-v2", "status": "PASS",
    "quality_pass_graph_nodes": 5834,
    "pairwise_pairs_requested": 17014861,
    "pairwise_pairs_evaluated": 17014861,
    "triangle_rows_emitted": 322110,
    "positive_edges": 103555, "negative_edges": 48205,
    "clusters_ge10": 62, "assigned_sags": 2567,
    "unaggregated_sags": 3267,
    "assigned_plus_unaggregated_equals_quality_pass": True,
}
for key, value in expected.items():
    if done.get(key) != value:
        raise SystemExit(f"Stage3B closure mismatch: {key}={done.get(key)!r} != {value!r}")
stats = json.loads((root / "02_pairwise_ani_af/triangle_stats.json").read_text(encoding="utf-8"))
if not (stats.get("status") == "PASS" and stats.get("triangle_mode") == "exact" and
        stats.get("exact_execution") == "global-seed-join" and
        stats.get("pairs_expected") == 17014861 and stats.get("pairs_evaluated") == 17014861 and
        stats.get("rows_emitted") == 322110):
    raise SystemExit("Stage3B exact-triangle receipt mismatch")
print("PASS Stage3B fixed lake closure 5834 nodes 17014861 exact pairs 62 groups")
PY
pipeline_end_epoch="$(date +%s)"
printf '%s\n' "${pipeline_end_epoch}" > "${CONTROL}/CORE_PIPELINE_END_EPOCH.txt"

# External benchmark-only quality evaluation; this is not part of the software's
# core contigs-to-bins runtime above.
set_phase final_bin_view
/usr/bin/time -v -o "${CONTROL}/08_FINAL_VIEW_TIME.txt" \
  "${PYTHON}" "${CODE_ROOT}/tools/prepare_current_baseline_bin_evaluation.py" \
  --stage3a-groups "${MAIN_OUT}/03A_subassemble/groups.tsv" \
  --stage3b-clusters "${STAGE3B_OUT}/06_subassemble/clusters.tsv" \
  --workers "$(( THREADS < 16 ? THREADS : 16 ))" --out-dir "${FINAL_VIEW}"

mkdir "${FINAL_EVAL}"
set_phase final_checkm2_evaluation
export PATH="${CHECKM2_ENV}/bin:${ORIGINAL_PATH}"
"${CHECKM2}" --version > "${FINAL_EVAL}/CHECKM2_VERSION.txt"
[[ "$(tr -d '\r\n' < "${FINAL_EVAL}/CHECKM2_VERSION.txt")" == "${CHECKM2_VERSION}" ]]
/usr/bin/time -v -o "${CONTROL}/09_FINAL_CHECKM2_TIME.txt" \
  "${CHECKM2}" predict --database_path "${CHECKM2_DB}" --input "${FINAL_VIEW}/input_fna" \
  --output-directory "${FINAL_EVAL}/checkm2_raw" -x fna --threads "${THREADS}" \
  >"${FINAL_EVAL}/CHECKM2_STDOUT.log" 2>"${FINAL_EVAL}/CHECKM2_STDERR.log"

set_phase final_quality_closure
"${PYTHON}" - "${FINAL_VIEW}/COMPLETE.json" "${FINAL_VIEW}/INPUT_MANIFEST.tsv" \
  "${FINAL_EVAL}/checkm2_raw/quality_report.tsv" "${FINAL_EVAL}/RESULT.json" <<'PY'
import csv, hashlib, json, math, pathlib, sys
view_receipt, input_manifest, report, result_path = map(pathlib.Path, sys.argv[1:])
view = json.loads(view_receipt.read_text(encoding="utf-8"))
if view.get("status") != "PASS" or view.get("schema") != "current-baseline-bin-evaluation-input-v1":
    raise SystemExit("final-bin input view receipt is not an expected PASS")
if view.get("stage3b_entry_closure_mode") != "current_v2_full_entry_closure":
    raise SystemExit("final-bin view lacks current Stage3B full-entry closure")
if view.get("checks", {}).get("current_stage3b_graph_nodes_equal_assigned_union_unaggregated") is not True:
    raise SystemExit("final-bin view Stage3B node closure check is not true")
expected = int(view["bin_count"])
if expected <= 0 or int(view.get("unique_bin_id_count", -1)) != expected:
    raise SystemExit("invalid final-bin count/unique-ID closure")
with input_manifest.open(newline="", encoding="utf-8") as handle:
    manifest_reader = csv.DictReader(handle, delimiter="\t")
    if not manifest_reader.fieldnames or "bin_id" not in manifest_reader.fieldnames:
        raise SystemExit("final-bin INPUT_MANIFEST lacks bin_id")
    manifest_rows = list(manifest_reader)
manifest_ids = [row.get("bin_id", "") for row in manifest_rows]
if len(manifest_ids) != expected or any(not name for name in manifest_ids) or len(set(manifest_ids)) != expected:
    raise SystemExit("final-bin manifest count/ID uniqueness mismatch")
with report.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required_columns = {"Name", "Completeness", "Contamination"}
    if not reader.fieldnames or not required_columns.issubset(reader.fieldnames):
        raise SystemExit("CheckM2 report lacks required columns")
    rows = list(reader)
if len(rows) != expected: raise SystemExit(f"row mismatch: {len(rows)} != {expected}")
names = [row.get("Name", "") for row in rows]
if len(names) != len(set(names)) or any(not name for name in names): raise SystemExit("empty/duplicate bin names")
if set(names) != set(manifest_ids):
    raise SystemExit("CheckM2 Name set differs from final-bin manifest bin_id set")
def number(row, key):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(f"invalid {key} for {row.get('Name','?')}: {exc}")
    if not math.isfinite(value):
        raise SystemExit(f"non-finite {key} for {row.get('Name','?')}")
    if key == "Completeness" and not (0.0 <= value <= 100.0):
        raise SystemExit(f"out-of-range completeness for {row.get('Name','?')}: {value}")
    if key == "Contamination" and value < 0.0:
        raise SystemExit(f"negative contamination for {row.get('Name','?')}: {value}")
    return value
for row in rows:
    number(row, "Completeness"); number(row, "Contamination")
hq = [r for r in rows if number(r,"Completeness") >= 90 and number(r,"Contamination") <= 5]
mq_inc = [r for r in rows if number(r,"Completeness") >= 50 and number(r,"Contamination") <= 10]
hq_names = {r["Name"] for r in hq}
mq_only = [r for r in mq_inc if r["Name"] not in hq_names]
payload = {
    "schema":"cellbit-current-tractor-final-quality-v1", "status":"PASS",
    "checkm2_version":"1.0.1", "checkm2_model_mode":"auto", "allmodels":False,
    "minimum_bin_size_filter":None, "total_bins":len(rows), "hq":len(hq),
    "mq_inclusive":len(mq_inc), "mq_only":len(mq_only),
    "lq_or_below":len(rows)-len(hq)-len(mq_only),
    "hq_rule":"Completeness >= 90 and Contamination <= 5",
    "mq_rule":"Completeness >= 50 and Contamination <= 10",
    "quality_report":str(report.resolve()), "quality_report_sha256":hashlib.sha256(report.read_bytes()).hexdigest(),
    "input_manifest":str(input_manifest.resolve()), "input_manifest_sha256":hashlib.sha256(input_manifest.read_bytes()).hexdigest(),
    "exact_bin_id_set_closed":True,
    "stage3b_full_entry_closure_verified":True,
    "input_view_complete_sha256":hashlib.sha256(view_receipt.read_bytes()).hexdigest(),
}
result_path.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n", encoding="utf-8")
print(json.dumps(payload, sort_keys=True))
PY

sha256sum "${MAIN_OUT}/COMPLETE.json" "${STAGE3B_INPUTS}/COMPLETE.json" \
  "${STAGE3B_OUT}/COMPLETE.json" "${FINAL_VIEW}/COMPLETE.json" \
  "${FINAL_VIEW}/INPUT_MANIFEST.tsv" "${CONTROL}/CODE_AND_INPUT_AUTHORITY.tsv" \
  "${CONTROL}/CACHE_CONDITION.json" "${CONTROL}/CGROUP_CPU_QUOTA.json" \
  "${FINAL_EVAL}/checkm2_raw/quality_report.tsv" "${FINAL_EVAL}/RESULT.json" \
  > "${CONTROL}/FINAL_SHA256SUMS.txt"
set_phase complete
