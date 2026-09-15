#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 QUALITY_FILTERED_SAGS_EFFECTIVE.tsv NEW_OUTPUT_DIR" >&2
  exit 2
fi

source_tsv=$1
out_dir=$2
expected_rows=8785

[[ -f "$source_tsv" ]] || { echo "missing source TSV: $source_tsv" >&2; exit 2; }
if [[ -e "$out_dir" ]]; then
  echo "write-once output already exists: $out_dir" >&2
  exit 2
fi
mkdir -p "$out_dir"

manifest_tmp="$out_dir/LAKE_CONTIGS_8785.tsv.incomplete"
manifest="$out_dir/LAKE_CONTIGS_8785.tsv"
awk -F '\t' 'BEGIN{OFS="\t"}
  NR==1 {
    if ($1!="SAG_id" || $2!="assembly_path") {
      print "unexpected source header" > "/dev/stderr"; exit 2
    }
    print "sag_id", "assembly_fasta"; next
  }
  {print $1, $2}
' "$source_tsv" > "$manifest_tmp"

declare -A seen=()
rows=0
duplicates=0
missing=0
while IFS=$'\t' read -r sag_id assembly_fasta extra; do
  [[ "$sag_id" == "sag_id" ]] && continue
  ((rows+=1))
  if [[ -n "${seen[$sag_id]:-}" ]]; then
    ((duplicates+=1))
  fi
  seen[$sag_id]=1
  [[ -f "$assembly_fasta" ]] || ((missing+=1))
  [[ -z "${extra:-}" ]] || { echo "unexpected third manifest field at $sag_id" >&2; exit 2; }
done < "$manifest_tmp"

if [[ $rows -ne $expected_rows || $duplicates -ne 0 || $missing -ne 0 ]]; then
  echo "input closure failed: rows=$rows duplicates=$duplicates missing=$missing" >&2
  exit 2
fi
mv "$manifest_tmp" "$manifest"

source_sha=$(sha256sum "$source_tsv" | awk '{print $1}')
manifest_sha=$(sha256sum "$manifest" | awk '{print $1}')
source_rows=$(awk 'END{print NR-1}' "$source_tsv")
total_bp=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="total_len")c=i; next} {s+=$c} END{printf "%.0f",s}' "$source_tsv")

{
  printf 'source\tsha256\trows\n'
  printf '%s\t%s\t%s\n' "$source_tsv" "$source_sha" "$source_rows"
  printf '%s\t%s\t%s\n' "$manifest" "$manifest_sha" "$rows"
} > "$out_dir/INPUT_AUTHORITY.tsv"

{
  printf 'status=PASS\n'
  printf 'rows=%s\n' "$rows"
  printf 'duplicates=%s\n' "$duplicates"
  printf 'missing=%s\n' "$missing"
  printf 'total_bp=%s\n' "$total_bp"
  printf 'manifest_sha256=%s\n' "$manifest_sha"
} > "$out_dir/INPUT_PREFLIGHT.txt"

printf 'PASS Lake contig manifest: rows=%s total_bp=%s sha256=%s\n' "$rows" "$total_bp" "$manifest_sha"
