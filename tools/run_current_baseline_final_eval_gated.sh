#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1
STAGE3B_UNIT=cellbit-lake-stage3b-exact-baseline-v1-20260904.service
STAGE3B="$ROOT/04_stage3b_exact_baseline_attempt_001"
CONTROL="$ROOT/00_stage3b_control_v1"
VIEW="$ROOT/05_current_baseline_checkm2_input_v1"
EVAL="$ROOT/06_current_baseline_checkm2_attempt_001"
CHECKM2_ENV=/home/data/fyc/biosoft/miniconda3/envs/checkm2_env
CHECKM2_DB=/home/data/shared/software/CheckM2_database/uniref100.KO.1.dmnd

while systemctl --user is-active --quiet "$STAGE3B_UNIT"; do
  sleep 30
done

stage3b_result=$(systemctl --user show "$STAGE3B_UNIT" -p Result --value)
stage3b_status=$(systemctl --user show "$STAGE3B_UNIT" -p ExecMainStatus --value)
if [[ "$stage3b_result" != success || "$stage3b_status" != 0 ]]; then
  echo "Stage3B unit did not finish successfully: Result=$stage3b_result ExecMainStatus=$stage3b_status" >&2
  exit 20
fi

python3 - "$STAGE3B/COMPLETE.json" <<'PY'
import json, pathlib, sys
p = pathlib.Path(sys.argv[1])
if not p.is_file():
    raise SystemExit(f"missing Stage3B receipt: {p}")
d = json.loads(p.read_text(encoding="utf-8"))
required = {
    "schema": "sag-stage3b-complete-v2",
    "status": "PASS",
    "reference_search": "disabled",
    "reference_positive_exclusion": False,
    "triangle_mode": "exact_all_pairs",
    "candidate_or_sketch_prefilter": False,
    "assigned_plus_unaggregated_equals_quality_pass": True,
}
bad = {key: (d.get(key), value) for key, value in required.items() if d.get(key) != value}
if bad:
    raise SystemExit(f"Stage3B current-contract receipt mismatch: {bad}")
PY

test ! -e "$VIEW"
test ! -e "$EVAL"

STAGE3B_CLUSTERS="$STAGE3B/06_subassemble/clusters.tsv"

python3 "$CONTROL/prepare_current_baseline_bin_evaluation.py" \
  --stage3a-groups "$ROOT/01_stage1_3a/03A_subassemble/groups.tsv" \
  --stage3b-clusters "$STAGE3B_CLUSTERS" \
  --out-dir "$VIEW"

mkdir "$EVAL"
export PATH="$CHECKM2_ENV/bin:$PATH"
"$CHECKM2_ENV/bin/checkm2" --version > "$EVAL/CHECKM2_VERSION.txt"
if [[ "$(tr -d '\r\n' < "$EVAL/CHECKM2_VERSION.txt")" != "1.0.1" ]]; then
  echo "unexpected CheckM2 version" >&2
  exit 30
fi

/usr/bin/time -v -o "$EVAL/CHECKM2_TIME.txt" \
  "$CHECKM2_ENV/bin/checkm2" predict \
    --database_path "$CHECKM2_DB" \
    --input "$VIEW/input_fna" \
    --output-directory "$EVAL/checkm2_raw" \
    -x fna \
    --threads 96 \
    >"$EVAL/CHECKM2_STDOUT.log" \
    2>"$EVAL/CHECKM2_STDERR.log"

python3 - "$VIEW/COMPLETE.json" "$EVAL/checkm2_raw/quality_report.tsv" "$EVAL/RESULT.json" <<'PY'
import csv, hashlib, json, pathlib, sys
view_receipt, report, result_path = map(pathlib.Path, sys.argv[1:])
view = json.loads(view_receipt.read_text(encoding="utf-8"))
if view.get("status") != "PASS":
    raise SystemExit("final-bin input view receipt is not PASS")
if view.get("stage3b_entry_closure_mode") != "current_v2_full_entry_closure":
    raise SystemExit("final-bin input view did not close the current Stage3B v2 entry")
if not view.get("checks", {}).get(
    "current_stage3b_graph_nodes_equal_assigned_union_unaggregated"
):
    raise SystemExit("Stage3B graph nodes != assigned union unaggregated")
with report.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
expected = int(view["bin_count"])
if len(rows) != expected:
    raise SystemExit(f"CheckM2 row-count mismatch: expected {expected}, observed {len(rows)}")
names = [r.get("Name", "") for r in rows]
if len(names) != len(set(names)) or any(not x for x in names):
    raise SystemExit("CheckM2 report contains empty or duplicate bin names")
def number(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError) as exc:
        raise SystemExit(f"invalid {key!r} in CheckM2 report") from exc
hq = [r for r in rows if number(r, "Completeness") >= 90 and number(r, "Contamination") <= 5]
mq_inclusive = [r for r in rows if number(r, "Completeness") >= 50 and number(r, "Contamination") <= 10]
hq_names = {r["Name"] for r in hq}
mq_only = [r for r in mq_inclusive if r["Name"] not in hq_names]
payload = {
    "schema": "cellbit-current-tractor-final-quality-v1",
    "status": "PASS",
    "checkm2_version": "1.0.1",
    "checkm2_model_mode": "auto",
    "allmodels": False,
    "minimum_bin_size_filter": None,
    "total_bins": len(rows),
    "hq": len(hq),
    "mq_inclusive": len(mq_inclusive),
    "mq_only": len(mq_only),
    "lq_or_below": len(rows) - len(hq) - len(mq_only),
    "hq_rule": "Completeness >= 90 and Contamination <= 5",
    "mq_rule": "Completeness >= 50 and Contamination <= 10",
    "quality_report": str(report.resolve()),
    "quality_report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
    "input_view_complete_sha256": hashlib.sha256(view_receipt.read_bytes()).hexdigest(),
}
result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(payload, sort_keys=True))
PY

sha256sum \
  "$VIEW/COMPLETE.json" \
  "$VIEW/INPUT_MANIFEST.tsv" \
  "$EVAL/CHECKM2_VERSION.txt" \
  "$EVAL/checkm2_raw/quality_report.tsv" \
  "$EVAL/RESULT.json" \
  > "$EVAL/FINAL_SHA256SUMS.txt"
