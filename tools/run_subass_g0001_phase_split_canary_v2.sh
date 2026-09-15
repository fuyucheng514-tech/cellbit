#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/subass_g0001_phase_split_v2
SRC=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1/04_stage3b_cellbit_negative_exact_attempt_003/06_subassemble/G0001/input_subassemblies.fasta
GOLD=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/subass_g0001_split_threads_v1/t2_repeat_1
BIN=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/build_speed_v1/subass/cpp-subass
FLYE=/home/data/fyc/biosoft/miniconda3/envs/assemble
[[ ! -e "$ROOT" ]]
mkdir -p "$ROOT/run"
sha256sum "$SRC" "$BIN" > "$ROOT/AUTHORITY_SHA256.txt"
common=(--reads "$SRC" --out-dir "$ROOT/run" --flye-modules "$FLYE/bin/flye-modules"
  --minimap2 "$FLYE/bin/flye-minimap2" --samtools "$FLYE/bin/flye-samtools"
  --package-root "$FLYE/lib/python3.9/site-packages/flye"
  --config "$FLYE/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg"
  --assemble-threads 1 --no-overlap-policy passthrough)
/usr/bin/time -v -o "$ROOT/assemble.time" "$BIN" --threads 1 --phase assemble "${common[@]}" \
  >"$ROOT/assemble.stdout" 2>"$ROOT/assemble.stderr"
/usr/bin/time -v -o "$ROOT/finish.time" "$BIN" --threads 2 --phase finish "${common[@]}" \
  >"$ROOT/finish.stdout" 2>"$ROOT/finish.stderr"
for relative in 00-assembly/draft_assembly.fasta 20-repeat/repeat_graph_edges.fasta \
  30-contigger/contigs.fasta 40-polishing/polished_1.fasta assembly.fasta; do
  cmp "$ROOT/run/$relative" "$GOLD/$relative"
  sha256sum "$ROOT/run/$relative" >> "$ROOT/SCIENTIFIC_SHA256.txt"
done
test -s "$ROOT/run/00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE"
test -s "$ROOT/run/40-polishing/CPP_FINISH_PHASE_COMPLETE"
test -s "$ROOT/run/CPP_FULL_PIPELINE_PASS"
printf 'PASS\n' > "$ROOT/RUN_RESULT.txt"
