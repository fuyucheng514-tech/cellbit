#!/usr/bin/env bash
set -Eeuo pipefail

readonly UNIT="cellbit-lake-full-e2e-220t-speed-v2-20260905.service"
readonly RUNNER="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/tools/run_lake_full_e2e_50t_benchmark_v1.sh"
readonly EXPECTED_RUNNER_SHA256="e9f089344010c3f0dda738e84a6bf3c8e245937e9f234dce76d936aa5de70e08"
readonly RUN_ROOT="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/24_lake_full_e2e_220t_speed_v1_20260905/attempt_001"

[[ -x "${RUNNER}" ]] || { printf 'runner is missing or not executable: %s\n' "${RUNNER}" >&2; exit 66; }
[[ "$(sha256sum "${RUNNER}" | awk '{print $1}')" == "${EXPECTED_RUNNER_SHA256}" ]] || {
  printf 'runner SHA256 no longer matches the frozen speed-v1 release\n' >&2
  exit 65
}
[[ ! -e "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  printf 'write-once run root already exists: %s\n' "${RUN_ROOT}" >&2
  exit 73
}
if systemctl --user is-active --quiet "${UNIT}"; then
  printf 'unit is already active: %s\n' "${UNIT}" >&2
  exit 75
fi
running_conflicts="$(ps -eo pid=,comm=,args= | awk '
  $2 == "diamond" || $2 == "gtdb-ani-af" || $2 == "cpp-subass" ||
  $2 == "flye-modules" || $2 == "flye-minimap2" || $2 == "flye-samtools" ||
  $2 == "dna2bit-sag-pip" || $2 == "sag-stage3b-tra" ||
  ($2 ~ /^python/ && ($0 ~ /\/checkm2 predict/ || $0 ~ /gtdbtk identify/ ||
                      $0 ~ /stage3b_signed_leiden/)) {print}
')"
[[ -z "${running_conflicts}" ]] || {
  printf 'refusing overlap with active scientific processes:\n%s\n' "${running_conflicts}" >&2
  exit 75
}

read -r load1 _ < /proc/loadavg
mem_available_kib="$(awk '$1 == "MemAvailable:" {print $2}' /proc/meminfo)"
io_some="$(awk '$1 == "some" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' /proc/pressure/io)"
io_full="$(awk '$1 == "full" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' /proc/pressure/io)"
memory_some="$(awk '$1 == "some" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' /proc/pressure/memory)"
memory_full="$(awk '$1 == "full" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' /proc/pressure/memory)"
readonly USER_MEMORY_PRESSURE="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/memory.pressure"
[[ -r "${USER_MEMORY_PRESSURE}" ]] || {
  printf 'cannot read ancestor user-cgroup memory PSI: %s\n' "${USER_MEMORY_PRESSURE}" >&2
  exit 75
}
user_memory_some="$(awk '$1 == "some" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' "${USER_MEMORY_PRESSURE}")"
user_memory_full="$(awk '$1 == "full" {for(i=2;i<=NF;i++) if($i ~ /^avg10=/){split($i,a,"="); print a[2]}}' "${USER_MEMORY_PRESSURE}")"
awk -v load="${load1}" -v mem="${mem_available_kib}" -v ios="${io_some}" -v iof="${io_full}" \
  -v ms="${memory_some}" -v mf="${memory_full}" \
  -v ums="${user_memory_some}" -v umf="${user_memory_full}" \
  'BEGIN {exit !(load < 100 && mem >= 838860800 && ios < 2 && iof < 2 && ms < 1 && mf < 1 && ums < 1 && umf < 1)}' || {
  printf 'resource gate failed: load1=%s MemAvailableKiB=%s io.some/full=%s/%s memory.some/full=%s/%s user-memory.some/full=%s/%s\n' \
    "${load1}" "${mem_available_kib}" "${io_some}" "${io_full}" \
    "${memory_some}" "${memory_full}" "${user_memory_some}" "${user_memory_full}" >&2
  exit 75
}

systemd-run --user --unit="${UNIT%.service}" \
  --property=AllowedCPUs=0-239 \
  --property=CPUQuota=22000% \
  --property=CPUWeight=20 \
  --property=IOWeight=20 \
  --property=MemoryHigh=1200G \
  --property=MemoryMax=1450G \
  --property=MemorySwapMax=8G \
  --property=ManagedOOMPreference=omit \
  --property=Restart=no \
  --property=KillMode=control-group \
  --property=TasksMax=infinity \
  --property=LimitNOFILE=1048576 \
  --setenv=CELLBIT_THREADS=220 \
  --setenv=CELLBIT_RUN_ROOT="${RUN_ROOT}" \
  /usr/bin/taskset --cpu-list 0-239 "${RUNNER}"

systemctl --user show "${UNIT}" \
  -p ActiveState -p SubState -p CPUQuotaPerSecUSec -p MemoryHigh -p MemoryMax -p MemorySwapMax
