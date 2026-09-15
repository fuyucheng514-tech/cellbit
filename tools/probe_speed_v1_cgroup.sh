#!/usr/bin/env bash
set -euo pipefail

rel="$(awk -F: '$1 == "0" {print $3; exit}' /proc/self/cgroup)"
root="/sys/fs/cgroup${rel}"
printf 'self_cgroup=%s\n' "${rel}"
printf 'root=%s\n' "${root}"
for name in cpuset.cpus cpuset.cpus.effective cpu.max memory.high memory.max memory.swap.max pids.max; do
  printf '%s=' "${name}"
  if [[ -r "${root}/${name}" ]]; then
    tr -d '\r\n' < "${root}/${name}"
  else
    printf '<unreadable>'
  fi
  printf '\n'
done
printf 'status_allowed_list='
awk -F: '$1 == "Cpus_allowed_list" {gsub(/^[[:space:]]+/, "", $2); print $2}' /proc/self/status
printf 'affinity_taskset='
taskset -pc $$ | sed 's/^.*: //'
printf 'root_listing:\n'
ls -la "${root}" || true
printf 'status_excerpt:\n'
grep -E 'Cpus_allowed|Mems_allowed' /proc/self/status || true
sleep "${CGROUP_PROBE_SLEEP_SECONDS:-0}"
