#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 INPUT_VIEW_DIR OUTPUT_ROOT" >&2
  exit 2
fi

input_view=$(realpath "$1")
output_root=$(realpath -m "$2")
input_dir="$input_view/input_fna"
input_receipt="$input_view/COMPLETE.json"
checkm2_out="$output_root/checkm2_attempt_001"
gtdbtk_out="$output_root/gtdbtk_identify_attempt_001"

checkm2_bin=/home/data/fyc/biosoft/miniconda3/envs/checkm2_env/bin/checkm2
checkm2_db=/home/data/shared/software/CheckM2_database/uniref100.KO.1.dmnd
gtdbtk_bin=/home/data/fyc/biosoft/miniconda3/envs/gtdbtk/bin/gtdbtk

# Invoking CheckM2 by absolute path does not activate its Conda environment.
# Bind its companion executables (DIAMOND, Prodigal, etc.) explicitly.
export PATH=/home/data/fyc/biosoft/miniconda3/envs/checkm2_env/bin:"$PATH"

[[ -s "$input_receipt" && -d "$input_dir" ]] || {
  echo "input view is not complete" >&2
  exit 2
}
[[ ! -e "$output_root" ]] || {
  echo "write-once output already exists: $output_root" >&2
  exit 2
}
[[ -x "$checkm2_bin" && -s "$checkm2_db" && -x "$gtdbtk_bin" ]] || {
  echo "required CheckM2/GTDB-Tk executable or database is missing" >&2
  exit 2
}

mkdir -p "$output_root"

/usr/bin/time -v -o "$output_root/CHECKM2_TIME.txt" \
  "$checkm2_bin" predict \
    --input "$input_dir" \
    --output-directory "$checkm2_out" \
    --database_path "$checkm2_db" \
    --extension fna \
    --threads 96 \
    --force \
    >"$output_root/CHECKM2_STDOUT.log" \
    2>"$output_root/CHECKM2_STDERR.log"

[[ -s "$checkm2_out/quality_report.tsv" ]] || {
  echo "CheckM2 completed without quality_report.tsv" >&2
  exit 2
}

export GTDBTK_DATA_PATH=/home/data/temp/release232
export PATH=/home/data/fyc/biosoft/miniconda3/envs/gtdbtk/bin:/home/data/fyc/biosoft/miniconda3/envs/eggnog/bin:"$PATH"
/usr/bin/time -v -o "$output_root/GTDBTK_TIME.txt" \
  "$gtdbtk_bin" identify \
    --genome_dir "$input_dir" \
    --out_dir "$gtdbtk_out" \
    -x fna \
    --cpus 96 \
    --force \
    --write_single_copy_genes \
    >"$output_root/GTDBTK_STDOUT.log" \
    2>"$output_root/GTDBTK_STDERR.log"

[[ -s "$gtdbtk_out/identify/gtdbtk.bac120.markers_summary.tsv" ]] || {
  echo "GTDB-Tk identify completed without bac120 marker summary" >&2
  exit 2
}

printf 'PASS\ncheckm2_threads=96\ngtdbtk_cpus=96\n' >"$output_root/UPSTREAM.PASS"
