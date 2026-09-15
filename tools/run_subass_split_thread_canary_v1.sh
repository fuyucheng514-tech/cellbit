#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/subass_split_threads_v1
SRC=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1/01_stage1_3a/03A_subassemble/Polynucleobacter_sp018882385_/input_subassemblies.fasta
BIN=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/build_speed_v1/subass/cpp-subass
FLYE=/home/data/fyc/biosoft/miniconda3/envs/assemble
EXPECTED_INPUT_SHA=169f9461a0259dd42b27a160902f551a860e97a469488b9970c7ebdbb8bbe41f

if [[ -e "$ROOT" ]]; then
  echo "refusing existing write-once canary root: $ROOT" >&2
  exit 2
fi
[[ "$(sha256sum "$SRC" | awk '{print $1}')" == "$EXPECTED_INPUT_SHA" ]]
mkdir -p "$ROOT"
sha256sum "$SRC" "$BIN" > "$ROOT/AUTHORITY_SHA256.txt"

run_one() {
  local name=$1
  local downstream_threads=$2
  local out="$ROOT/$name"
  mkdir "$out"
  /usr/bin/time -v -o "$ROOT/$name.time" \
    "$BIN" \
      --reads "$SRC" \
      --out-dir "$out" \
      --flye-modules "$FLYE/bin/flye-modules" \
      --minimap2 "$FLYE/bin/flye-minimap2" \
      --samtools "$FLYE/bin/flye-samtools" \
      --package-root "$FLYE/lib/python3.9/site-packages/flye" \
      --config "$FLYE/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg" \
      --threads "$downstream_threads" \
      --assemble-threads 1 \
      --no-overlap-policy passthrough \
      > "$ROOT/$name.stdout" 2> "$ROOT/$name.stderr"
  sha256sum \
    "$out/00-assembly/draft_assembly.fasta" \
    "$out/20-repeat/repeat_graph_edges.fasta" \
    "$out/30-contigger/contigs.fasta" \
    "$out/40-polishing/polished_1.fasta" \
    "$out/assembly.fasta" > "$ROOT/$name.sha256"
}

run_one t1_control 1
run_one t2_repeat_1 2
run_one t2_repeat_2 2
run_one t2_repeat_3 2

python3 - "$ROOT" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
names = ["t1_control", "t2_repeat_1", "t2_repeat_2", "t2_repeat_3"]

def parse_fasta(path):
    records = []
    seq = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    records.append("".join(seq).upper())
                    seq = []
            else:
                seq.append(line)
    if seq:
        records.append("".join(seq).upper())
    payload = "\n".join(sorted(records)).encode()
    return hashlib.sha256(payload).hexdigest(), len(records), sum(map(len, records))

result = {}
for name in names:
    raw = hashlib.sha256((root / name / "assembly.fasta").read_bytes()).hexdigest()
    norm, records, bp = parse_fasta(root / name / "assembly.fasta")
    time_text = (root / f"{name}.time").read_text()
    fields = {}
    for line in time_text.splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", 1)
            fields[key] = value
    result[name] = {
        "downstream_threads": 1 if name == "t1_control" else 2,
        "raw_sha256": raw,
        "normalized_sequence_sha256": norm,
        "records": records,
        "bp": bp,
        "elapsed": fields.get("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "user_seconds": fields.get("User time (seconds)"),
        "max_rss_kb": fields.get("Maximum resident set size (kbytes)"),
    }

t2 = [result[name] for name in names[1:]]
result["checks"] = {
    "t2_raw_identical": len({x["raw_sha256"] for x in t2}) == 1,
    "t2_normalized_identical": len({x["normalized_sequence_sha256"] for x in t2}) == 1,
    "t1_vs_t2_raw_equal": result["t1_control"]["raw_sha256"] == t2[0]["raw_sha256"],
    "t1_vs_t2_normalized_equal": result["t1_control"]["normalized_sequence_sha256"] == t2[0]["normalized_sequence_sha256"],
}
(root / "RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
if not all(result["checks"].values()):
    raise SystemExit("scientific/determinism comparison failed")
(root / "RUN_RESULT.txt").write_text("PASS\n")
PY
