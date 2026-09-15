#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
bin=${1:-"$root/build/gtdb-ani-af"}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" <<'PY'
import pathlib, random, sys
d = pathlib.Path(sys.argv[1])
rng = random.Random(20260904)
a = ''.join(rng.choices('ACGT', k=16000))
b = ''.join(rng.choices('ACGT', k=16000))
(d/'RS_GCF_000000001.1_genomic.fna').write_text('>a\n'+a+'\n')
(d/'GB_GCA_000000002.1_genomic.fna').write_text('>b\n'+b+'\n')
(d/'query.fna').write_text('>q\n'+a[1000:13000]+'\n')
(d/'refs.tsv').write_text(
    'Ref_file\nRS_GCF_000000001.1_genomic.fna\nGB_GCA_000000002.1_genomic.fna\n')
(d/'queries.list').write_text(str(d/'query.fna')+'\n')
(d/'triangle.list').write_text(
    str(d/'RS_GCF_000000001.1_genomic.fna')+'\n'+
    str(d/'GB_GCA_000000002.1_genomic.fna')+'\n')
(d/'old_taxonomy.csv').write_text(
    'RS_GCF_000000001.1,d__Bacteria;p__Old;s__Old_A\n'
    'GB_GCA_000000002.1,d__Bacteria;p__Old;s__Old_B\n')
(d/'new_taxonomy.csv').write_text(
    'accession,gtdb_taxonomy\n'
    'GCF_000000001.1,d__Bacteria;p__New;s__New_A\n'
    'GCA_000000002.1,d__Bacteria;p__New;s__New_B\n')
(d/'bad_taxonomy.csv').write_text(
    'GCF_000000001.1,d__Bacteria;p__New;s__New_A\n'
    'GCA_999999999.1,d__Bacteria;p__New;s__Extra\n')
(d/'duplicate_taxonomy.csv').write_text(
    'GCF_000000001.1,d__Bacteria;p__New;s__New_A\n'
    'RS_GCF_000000001.1,d__Bacteria;p__New;s__Duplicate\n'
    'GCA_000000002.1,d__Bacteria;p__New;s__New_B\n')
(d/'radii.tsv').write_text('representative\tradius\nGCF_000000001.1\t95\n')
(d/'release.txt').write_text('synthetic-release=old-taxonomy-binding\n')
PY

"$bin" index --manifest "$tmp/refs.tsv" --taxonomy "$tmp/old_taxonomy.csv" \
  --radii "$tmp/radii.tsv" --release-metadata "$tmp/release.txt" \
  --out-dir "$tmp/source" --index-format compact-v4 --threads 2 \
  --k 13 --sketch-scale 4 --query-sketch-size 100 --max-posting 10
"$bin" search --index "$tmp/source" --queries "$tmp/queries.list" \
  --out "$tmp/search.before.tsv" --threads 2 --top 2 --deep-verify
"$bin" triangle --list "$tmp/triangle.list" --out "$tmp/triangle.before.tsv" \
  --stats-out "$tmp/triangle.before.stats.json" \
  --threads 2 --profile medium --edge-mode sparse --min-af 0

"$bin" rebind-taxonomy --index "$tmp/source" \
  --taxonomy "$tmp/new_taxonomy.csv" --out-dir "$tmp/rebound"
"$bin" index-info --index "$tmp/rebound" > "$tmp/rebound.info.tsv"
"$bin" search --index "$tmp/rebound" --queries "$tmp/queries.list" \
  --out "$tmp/search.after.tsv" --threads 2 --top 2 --deep-verify
"$bin" triangle --list "$tmp/triangle.list" --out "$tmp/triangle.after.tsv" \
  --stats-out "$tmp/triangle.after.stats.json" \
  --threads 2 --profile medium --edge-mode sparse --min-af 0
cmp "$tmp/search.before.tsv" "$tmp/search.after.tsv"
cmp "$tmp/triangle.before.tsv" "$tmp/triangle.after.tsv"

if "$bin" rebind-taxonomy --index "$tmp/source" \
    --taxonomy "$tmp/bad_taxonomy.csv" --out-dir "$tmp/should_not_exist_bad"; then
  echo 'mismatched taxonomy unexpectedly accepted' >&2; exit 1
fi
test ! -e "$tmp/should_not_exist_bad"
if "$bin" rebind-taxonomy --index "$tmp/source" \
    --taxonomy "$tmp/duplicate_taxonomy.csv" --out-dir "$tmp/should_not_exist_duplicate"; then
  echo 'duplicate taxonomy unexpectedly accepted' >&2; exit 1
fi
test ! -e "$tmp/should_not_exist_duplicate"

python3 - "$tmp" <<'PY'
import csv, hashlib, json, pathlib, re, sys
d = pathlib.Path(sys.argv[1])
src = json.loads((d/'source'/'COMPLETE.json').read_text())
dst = json.loads((d/'rebound'/'COMPLETE.json').read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before_stats = json.loads((d/'triangle.before.stats.json').read_text())
after_stats = json.loads((d/'triangle.after.stats.json').read_text())
for stats, output in ((before_stats, d/'triangle.before.tsv'),
                      (after_stats, d/'triangle.after.tsv')):
    assert stats['schema'] == 'gtdb-ani-af-triangle-stats-v1'
    assert stats['status'] == 'PASS' and stats['triangle_mode'] == 'exact'
    assert stats['nodes'] == 2 and stats['pairs_expected'] == 1
    assert stats['pairs_evaluated'] == 1 and stats['rows_emitted'] == 1
    assert stats['input_list_sha256'] == sha(d/'triangle.list')
    assert stats['output_sha256'] == sha(output)
assert dst['schema'] == src['schema'] == 'gtdb-ani-af-index-v4'
assert dst['status'] == 'PASS' and dst['taxonomy_rebound_from_schema'] == src['schema']
assert dst['source_complete_sha256'] == sha(d/'source'/'COMPLETE.json')
assert sha(d/'rebound'/'SOURCE_COMPLETE.json') == dst['source_complete_sha256']
assert dst['source_taxonomy_sha256'] == src['taxonomy_sha256']
assert dst['taxonomy_sha256'] == sha(d/'new_taxonomy.csv') == sha(d/'rebound'/'TAXONOMY.tsv')
assert dst['taxonomy_sha256'] != src['taxonomy_sha256']
assert dst['taxonomy_entry_count'] == dst['reference_count'] == 2
for name, key in [('REFS.tsv','refs_sha256'), ('SKETCHES.bin','sketches_sha256'),
                  ('POSTINGS.bin','postings_sha256'),
                  ('POSTINGS_LOOKUP.bin','postings_lookup_sha256')]:
    assert sha(d/'source'/name) == sha(d/'rebound'/name) == src[key] == dst[key]
with (d/'search.after.tsv').open() as f:
    hit = next(csv.DictReader(f, delimiter='\t'))
acc = re.search(r'GC[AF]_\d+\.\d+', pathlib.Path(hit['Ref_file']).name).group(0)
def taxonomy(path):
    rows = list(csv.reader(path.open()))
    if rows[0][0] == 'accession': rows = rows[1:]
    return {re.sub(r'^(RS_|GB_)','',r[0]): r[1] for r in rows}
old, new = taxonomy(d/'old_taxonomy.csv'), taxonomy(d/'new_taxonomy.csv')
assert old[acc] != new[acc] and ';p__Old;' in old[acc] and ';p__New;' in new[acc]
assert (d/'rebound.info.tsv').read_text().startswith('schema\treferences\t')
print('PASS rebind: search/triangle byte-exact; core SHA-exact; taxonomy mapping/provenance changed')
PY

printf '\nGCA_999999999.1,d__Bacteria;s__Tampered\n' >> "$tmp/rebound/TAXONOMY.tsv"
if "$bin" index-info --index "$tmp/rebound"; then
  echo 'tampered bound taxonomy unexpectedly accepted by load_index' >&2; exit 1
fi

echo 'PASS gtdb-ani-af taxonomy rebind synthetic regression'
