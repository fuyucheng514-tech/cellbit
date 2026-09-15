#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/subass_g0001_split_threads_v1
SRC=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1/04_stage3b_cellbit_negative_exact_attempt_003/06_subassemble/G0001/input_subassemblies.fasta
BIN=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/build_speed_v1/subass/cpp-subass
FLYE=/home/data/fyc/biosoft/miniconda3/envs/assemble
EXPECTED_INPUT_SHA=a49e787f991506cf3abf516c24062a62e7bc2599c63df9b30fe23c16b40a6fb0

if [[ -e "$ROOT" ]]; then
  echo "refusing existing write-once canary root: $ROOT" >&2
  exit 2
fi
[[ "$(sha256sum "$SRC" | awk '{print $1}')" == "$EXPECTED_INPUT_SHA" ]]
mkdir -p "$ROOT"
sha256sum "$SRC" "$BIN" > "$ROOT/AUTHORITY_SHA256.txt"

run_one() {
  local name=$1 threads=$2 out="$ROOT/$1"
  mkdir "$out"
  /usr/bin/time -v -o "$ROOT/$name.time" \
    "$BIN" --reads "$SRC" --out-dir "$out" \
      --flye-modules "$FLYE/bin/flye-modules" \
      --minimap2 "$FLYE/bin/flye-minimap2" \
      --samtools "$FLYE/bin/flye-samtools" \
      --package-root "$FLYE/lib/python3.9/site-packages/flye" \
      --config "$FLYE/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg" \
      --threads "$threads" --assemble-threads 1 --no-overlap-policy passthrough \
      > "$ROOT/$name.stdout" 2> "$ROOT/$name.stderr"
  sha256sum "$out/00-assembly/draft_assembly.fasta" "$out/20-repeat/repeat_graph_edges.fasta" \
    "$out/30-contigger/contigs.fasta" "$out/40-polishing/polished_1.fasta" \
    "$out/assembly.fasta" > "$ROOT/$name.sha256"
}

run_one t1_control 1
run_one t2_repeat_1 2
run_one t2_repeat_2 2

python3 - "$ROOT" <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
names = ["t1_control", "t2_repeat_1", "t2_repeat_2"]

def norm(path):
    seqs, seq = [], []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line: continue
        if line.startswith(">"):
            if seq: seqs.append("".join(seq).upper()); seq=[]
        else: seq.append(line)
    if seq: seqs.append("".join(seq).upper())
    return hashlib.sha256("\n".join(sorted(seqs)).encode()).hexdigest(), len(seqs), sum(map(len,seqs))

result = {}
for name in names:
    fasta = root/name/"assembly.fasta"
    nsha, records, bp = norm(fasta)
    fields = {}
    for line in (root/f"{name}.time").read_text().splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", 1); fields[key]=value
    result[name] = {
        "downstream_threads": 1 if name == "t1_control" else 2,
        "raw_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        "normalized_sequence_sha256": nsha, "records": records, "bp": bp,
        "elapsed": fields.get("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "user_seconds": fields.get("User time (seconds)"),
        "max_rss_kb": fields.get("Maximum resident set size (kbytes)"),
    }
t2=[result[n] for n in names[1:]]
result["checks"]={
  "t2_raw_identical":len({x["raw_sha256"] for x in t2})==1,
  "t2_normalized_identical":len({x["normalized_sequence_sha256"] for x in t2})==1,
  "t1_vs_t2_raw_equal":result["t1_control"]["raw_sha256"]==t2[0]["raw_sha256"],
  "t1_vs_t2_normalized_equal":result["t1_control"]["normalized_sequence_sha256"]==t2[0]["normalized_sequence_sha256"],
}
(root/"RESULT.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
if not all(result["checks"].values()): raise SystemExit("scientific/determinism comparison failed")
(root/"RUN_RESULT.txt").write_text("PASS\n")
PY
