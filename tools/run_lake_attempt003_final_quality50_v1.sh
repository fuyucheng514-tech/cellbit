#!/usr/bin/env bash
set -Eeuo pipefail

readonly THREADS=50
readonly CODE_ROOT="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/13_cellbit_sag_pipeline_cpp_v1_positive_edge_precision_only_v1"
readonly BASE="/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1"
readonly STAGE3A="${BASE}/01_stage1_3a"
readonly STAGE3B="${BASE}/04_stage3b_cellbit_negative_exact_attempt_003"
readonly VIEW="${BASE}/06_attempt003_final_bin_view_v1"
readonly EVAL="${BASE}/07_attempt003_final_checkm2_50t_v1"
readonly CONTROL="${BASE}/08_attempt003_final_quality_control_v1"
readonly PYTHON="/home/data/fyc/biosoft/miniconda3/bin/python"
readonly CHECKM2_ENV="/home/data/fyc/biosoft/miniconda3/envs/checkm2_env"
readonly CHECKM2="${CHECKM2_ENV}/bin/checkm2"
readonly CHECKM2_DB="/home/data/shared/software/CheckM2_database/uniref100.KO.1.dmnd"
readonly PREPARE="${CODE_ROOT}/tools/prepare_current_baseline_bin_evaluation.py"

# Everything before the CONTROL mkdir is read-only preflight.  CONTROL is the
# single atomic lease for this write-once evaluation; VIEW and EVAL are then
# created by commands which also fail closed if either path appears.
for path in "${VIEW}" "${EVAL}" "${CONTROL}"; do
  [[ ! -e "${path}" && ! -L "${path}" ]] || {
    printf 'write-once target exists (including dangling symlink): %s\n' "${path}" >&2
    exit 73
  }
done
for path in "${STAGE3A}/03A_subassemble/groups.tsv" "${STAGE3B}/06_subassemble/clusters.tsv" \
            "${STAGE3B}/COMPLETE.json" "${PYTHON}" "${CHECKM2}" "${CHECKM2_DB}" "${PREPARE}"; do
  [[ -e "${path}" ]] || { printf 'missing required path: %s\n' "${path}" >&2; exit 66; }
done
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("status")=="PASS" and d.get("schema")=="sag-stage3b-complete-v2"' "${STAGE3B}/COMPLETE.json"

cgroup_rel="$(awk -F: '$1 == "0" { print $3; exit }' /proc/self/cgroup)"
[[ -n "${cgroup_rel}" ]] || { printf 'cannot resolve unified cgroup path\n' >&2; exit 69; }
cpu_max_path="/sys/fs/cgroup${cgroup_rel%/}/cpu.max"
[[ -r "${cpu_max_path}" ]] || { printf 'cannot read cgroup CPU quota: %s\n' "${cpu_max_path}" >&2; exit 69; }
read -r cpu_quota cpu_period cpu_extra < "${cpu_max_path}"
[[ -z "${cpu_extra:-}" && "${cpu_quota}" != "max" \
   && "${cpu_quota}" =~ ^[0-9]+$ && "${cpu_period}" =~ ^[0-9]+$ \
   && "${cpu_period}" -gt 0 ]] || {
  printf 'invalid or unlimited cgroup cpu.max: %s\n' "$(<"${cpu_max_path}")" >&2
  exit 69
}
expected_cpu_quota=$((THREADS * cpu_period))
[[ "${cpu_quota}" -eq "${expected_cpu_quota}" ]] || {
  printf 'cgroup CPU quota is not exactly %s cores: quota=%s period=%s\n' \
    "${THREADS}" "${cpu_quota}" "${cpu_period}" >&2
  exit 69
}
actual_cpu_quota_cores="$(awk -v q="${cpu_quota}" -v p="${cpu_period}" 'BEGIN { printf "%.6f", q / p }')"

checkm2_version="$(PATH="${CHECKM2_ENV}/bin:${PATH}" "${CHECKM2}" --version 2>&1 | tr -d '\r\n')"
[[ "${checkm2_version}" == "1.0.1" ]] || {
  printf 'unexpected CheckM2 version before science: %s\n' "${checkm2_version}" >&2
  exit 69
}

mkdir -- "${CONTROL}" || { printf 'failed to atomically claim write-once control path: %s\n' "${CONTROL}" >&2; exit 73; }
exec > >(tee -a "${CONTROL}/STDOUT.log") 2> >(tee -a "${CONTROL}/STDERR.log" >&2)
start_epoch="$(date +%s)"
{
  printf 'key\tvalue\n'
  printf 'write_once_control_claim\ttrue\n'
  printf 'cgroup_relative_path\t%s\n' "${cgroup_rel}"
  printf 'cpu_max_path\t%s\n' "${cpu_max_path}"
  printf 'cpu_quota_us\t%s\n' "${cpu_quota}"
  printf 'cpu_period_us\t%s\n' "${cpu_period}"
  printf 'actual_cpu_quota_cores\t%s\n' "${actual_cpu_quota_cores}"
  printf 'required_thread_budget\t%s\n' "${THREADS}"
  printf 'checkm2_version_preflight\t%s\n' "${checkm2_version}"
} > "${CONTROL}/PREFLIGHT.tsv"

/usr/bin/time -v -o "${CONTROL}/FINAL_VIEW_TIME.txt" \
  "${PYTHON}" "${PREPARE}" --stage3a-groups "${STAGE3A}/03A_subassemble/groups.tsv" \
  --stage3b-clusters "${STAGE3B}/06_subassemble/clusters.tsv" --out-dir "${VIEW}"

mkdir "${EVAL}"
export PATH="${CHECKM2_ENV}/bin:${PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '%s\n' "${checkm2_version}" > "${EVAL}/CHECKM2_VERSION.txt"
/usr/bin/time -v -o "${CONTROL}/FINAL_CHECKM2_TIME.txt" \
  "${CHECKM2}" predict --database_path "${CHECKM2_DB}" --input "${VIEW}/input_fna" \
  --output-directory "${EVAL}/checkm2_raw" -x fna --threads "${THREADS}" \
  >"${EVAL}/CHECKM2_STDOUT.log" 2>"${EVAL}/CHECKM2_STDERR.log"

"${PYTHON}" - "${VIEW}/COMPLETE.json" "${VIEW}/INPUT_MANIFEST.tsv" \
  "${EVAL}/checkm2_raw/quality_report.tsv" "${EVAL}/RESULT.json" <<'PY'
import csv, hashlib, json, math, pathlib, sys

view_receipt, manifest, report, result_path = map(pathlib.Path, sys.argv[1:])
view = json.loads(view_receipt.read_text(encoding="utf-8"))
if view.get("schema") != "current-baseline-bin-evaluation-input-v1":
    raise SystemExit("unexpected final-bin view schema")
if view.get("status") != "PASS":
    raise SystemExit("final-bin view not PASS")
if view.get("stage3b_entry_closure_mode") != "current_v2_full_entry_closure":
    raise SystemExit("final-bin view lacks current-v2 full entry closure")
checks = view.get("checks")
if not isinstance(checks, dict):
    raise SystemExit("final-bin view checks object missing")
if checks.get("current_stage3b_graph_nodes_equal_assigned_union_unaggregated") is not True:
    raise SystemExit("Stage3B graph/assigned/unaggregated closure is not true")

expected = view.get("bin_count")
if type(expected) is not int or expected <= 0:
    raise SystemExit("invalid final-bin bin_count")
for key in ("nonempty_bin_count", "unique_bin_id_count", "symlink_count"):
    if view.get(key) != expected:
        raise SystemExit(f"final-bin {key} does not equal bin_count")
stage3a_count = view.get("stage3a_bin_count")
stage3b_count = view.get("stage3b_bin_count")
if type(stage3a_count) is not int or type(stage3b_count) is not int or stage3a_count + stage3b_count != expected:
    raise SystemExit("stage3a/stage3b bin counts do not close to bin_count")
graph_nodes = view.get("stage3b_graph_nodes")
assigned = view.get("stage3b_assigned_sags")
unaggregated = view.get("stage3b_unaggregated_sags")
if not all(type(value) is int and value >= 0 for value in (graph_nodes, assigned, unaggregated)):
    raise SystemExit("invalid Stage3B graph closure counts")
if graph_nodes != assigned + unaggregated:
    raise SystemExit("Stage3B graph nodes do not equal assigned plus unaggregated")

manifest_receipt = view.get("input_manifest")
if not isinstance(manifest_receipt, dict):
    raise SystemExit("input-manifest receipt missing")
if pathlib.Path(manifest_receipt.get("path", "")).resolve() != manifest.resolve():
    raise SystemExit("input-manifest receipt path mismatch")
manifest_bytes = manifest.read_bytes()
manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
if manifest_receipt.get("bytes") != len(manifest_bytes) or manifest_receipt.get("sha256") != manifest_sha256:
    raise SystemExit("input-manifest receipt size/hash mismatch")
with manifest.open(newline="", encoding="utf-8") as handle:
    manifest_reader = csv.DictReader(handle, delimiter="\t")
    if not manifest_reader.fieldnames or "bin_id" not in manifest_reader.fieldnames:
        raise SystemExit("input manifest lacks bin_id column")
    manifest_rows = list(manifest_reader)
manifest_names = [row.get("bin_id", "") for row in manifest_rows]
if len(manifest_rows) != expected:
    raise SystemExit("input-manifest row-count mismatch")
if any(not name for name in manifest_names) or len(manifest_names) != len(set(manifest_names)):
    raise SystemExit("empty/duplicate input-manifest bin IDs")
manifest_name_set = set(manifest_names)

with report.open(newline="", encoding="utf-8") as handle:
    report_reader = csv.DictReader(handle, delimiter="\t")
    required_columns = {"Name", "Completeness", "Contamination"}
    if not report_reader.fieldnames or not required_columns.issubset(report_reader.fieldnames):
        raise SystemExit("CheckM2 report lacks required columns")
    rows = list(report_reader)
if len(rows) != expected:
    raise SystemExit("CheckM2 row-count mismatch")
names = [row.get("Name", "") for row in rows]
if len(names) != len(set(names)) or any(not name for name in names):
    raise SystemExit("empty/duplicate CheckM2 bin names")
if set(names) != manifest_name_set:
    missing = sorted(manifest_name_set - set(names))
    extra = sorted(set(names) - manifest_name_set)
    raise SystemExit(f"CheckM2/input-manifest bin-ID set mismatch: missing={missing[:5]} extra={extra[:5]}")

def n(row, key):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(f"invalid CheckM2 {key} for {row.get('Name', '<unknown>')}: {exc}")
    if not math.isfinite(value):
        raise SystemExit(f"non-finite CheckM2 {key} for {row.get('Name', '<unknown>')}")
    return value

numeric_rows = []
for row in rows:
    completeness = n(row, "Completeness")
    contamination = n(row, "Contamination")
    if not 0.0 <= completeness <= 100.0:
        raise SystemExit(f"out-of-range completeness for {row['Name']}: {completeness}")
    if contamination < 0.0:
        raise SystemExit(f"negative contamination for {row['Name']}: {contamination}")
    numeric_rows.append((row, completeness, contamination))

hq = [r for r, completeness, contamination in numeric_rows if completeness >= 90 and contamination <= 5]
mq_inc = [r for r, completeness, contamination in numeric_rows if completeness >= 50 and contamination <= 10]
hq_names = {r["Name"] for r in hq}
mq_only = [r for r in mq_inc if r["Name"] not in hq_names]
id_set_sha256 = hashlib.sha256(("\n".join(sorted(manifest_name_set)) + "\n").encode()).hexdigest()
payload = {
  "schema":"cellbit-current-tractor-final-quality-v1", "status":"PASS",
  "checkm2_version":"1.0.1", "checkm2_model_mode":"auto", "allmodels":False,
  "minimum_bin_size_filter":None, "total_bins":len(rows), "hq":len(hq),
  "mq_inclusive":len(mq_inc), "mq_only":len(mq_only),
  "lq_or_below":len(rows)-len(hq)-len(mq_only),
  "hq_rule":"Completeness >= 90 and Contamination <= 5",
  "mq_rule":"Completeness >= 50 and Contamination <= 10",
  "quality_report":str(report.resolve()), "quality_report_sha256":hashlib.sha256(report.read_bytes()).hexdigest(),
  "input_view_complete_sha256":hashlib.sha256(view_receipt.read_bytes()).hexdigest(),
  "input_manifest_sha256":manifest_sha256,
  "bin_id_set_sha256":id_set_sha256,
  "exact_input_manifest_checkm2_id_set":True,
  "stage3b_entry_closure_mode":"current_v2_full_entry_closure",
  "stage3b_graph_nodes_equal_assigned_union_unaggregated":True,
  "numeric_fields_finite_and_bounded":True,
}
result_path.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n", encoding="utf-8")
print(json.dumps(payload, sort_keys=True))
PY

end_epoch="$(date +%s)"
cat > "${CONTROL}/RUN_RESULT.json" <<EOF
{
  "schema": "cellbit-attempt003-final-quality-run-v1",
  "status": "PASS",
  "threads": ${THREADS},
  "write_once_control_claim": true,
  "cgroup_relative_path": "${cgroup_rel}",
  "cpu_quota_us": ${cpu_quota},
  "cpu_period_us": ${cpu_period},
  "actual_cpu_quota_cores": ${actual_cpu_quota_cores},
  "checkm2_version_preflight": "${checkm2_version}",
  "start_epoch": ${start_epoch},
  "end_epoch": ${end_epoch},
  "elapsed_seconds": $((end_epoch-start_epoch)),
  "stage3b_complete_sha256": "$(sha256sum "${STAGE3B}/COMPLETE.json" | awk '{print $1}')",
  "view_complete_sha256": "$(sha256sum "${VIEW}/COMPLETE.json" | awk '{print $1}')",
  "quality_result_sha256": "$(sha256sum "${EVAL}/RESULT.json" | awk '{print $1}')"
}
EOF
sha256sum "${STAGE3B}/COMPLETE.json" "${VIEW}/COMPLETE.json" \
  "${VIEW}/INPUT_MANIFEST.tsv" "${EVAL}/CHECKM2_VERSION.txt" \
  "${EVAL}/checkm2_raw/quality_report.tsv" "${EVAL}/RESULT.json" \
  "${CONTROL}/PREFLIGHT.tsv" "${CONTROL}/RUN_RESULT.json" \
  > "${CONTROL}/FINAL_SHA256SUMS.txt"
