#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/18_speed_canary_20260905/subass_phase_split_v1
SRC=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/11_lake_full_current_tractor_baseline_v1/01_stage1_3a/03A_subassemble/Polynucleobacter_sp018882385_/input_subassemblies.fasta
BIN=/home/data/fyc/the_oftware_oiiaioiiiai/c++subass/17_cellbit_sag_pipeline_cpp_speed_v1_20260905/build_speed_v1/subass/cpp-subass
FLYE=/home/data/fyc/biosoft/miniconda3/envs/assemble
[[ ! -e "$ROOT" ]]
mkdir -p "$ROOT"
sha256sum "$SRC" "$BIN" > "$ROOT/AUTHORITY_SHA256.txt"

common=(--reads "$SRC" --flye-modules "$FLYE/bin/flye-modules"
  --minimap2 "$FLYE/bin/flye-minimap2" --samtools "$FLYE/bin/flye-samtools"
  --package-root "$FLYE/lib/python3.9/site-packages/flye"
  --config "$FLYE/lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg"
  --assemble-threads 1 --no-overlap-policy passthrough)

mkdir "$ROOT/all_control"
/usr/bin/time -v -o "$ROOT/all_control.time" "$BIN" --out-dir "$ROOT/all_control" \
  --threads 2 --phase all "${common[@]}" >"$ROOT/all_control.stdout" 2>"$ROOT/all_control.stderr"

for n in split_1 split_2; do
  mkdir "$ROOT/$n"
  /usr/bin/time -v -o "$ROOT/$n.assemble.time" "$BIN" --out-dir "$ROOT/$n" \
    --threads 1 --phase assemble "${common[@]}" >"$ROOT/$n.assemble.stdout" 2>"$ROOT/$n.assemble.stderr"
  /usr/bin/time -v -o "$ROOT/$n.finish.time" "$BIN" --out-dir "$ROOT/$n" \
    --threads 2 --phase finish "${common[@]}" >"$ROOT/$n.finish.stdout" 2>"$ROOT/$n.finish.stderr"
done

for n in all_control split_1 split_2; do
  test -s "$ROOT/$n/00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE"
  test -s "$ROOT/$n/40-polishing/CPP_FINISH_PHASE_COMPLETE"
  test -s "$ROOT/$n/CPP_FULL_PIPELINE_PASS"
  sha256sum "$ROOT/$n/00-assembly/draft_assembly.fasta" \
    "$ROOT/$n/20-repeat/repeat_graph_edges.fasta" "$ROOT/$n/30-contigger/contigs.fasta" \
    "$ROOT/$n/40-polishing/polished_1.fasta" "$ROOT/$n/assembly.fasta" \
    | cut -d' ' -f1 > "$ROOT/$n.scientific_sha256"
done
cmp "$ROOT/all_control.scientific_sha256" "$ROOT/split_1.scientific_sha256"
cmp "$ROOT/split_1.scientific_sha256" "$ROOT/split_2.scientific_sha256"
{
  printf 'mode\tphase\telapsed\n'
  for file in "$ROOT"/*.time; do
    printf '%s\t%s\t%s\n' "$(basename "$file" .time)" time \
      "$(sed -n 's/^\s*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$file")"
  done
} > "$ROOT/TIMING.tsv"
printf 'PASS\n' > "$ROOT/RUN_RESULT.txt"
