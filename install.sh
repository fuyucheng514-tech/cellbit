#!/usr/bin/env bash
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PREFIX=${PREFIX:-"$HOME/.local/dna2bit-sag-tractor"}
JOBS=${JOBS:-8}
command -v cmake >/dev/null 2>&1 || { echo "cmake is required" >&2; exit 2; }
command -v c++ >/dev/null 2>&1 || { echo "a C++ compiler is required" >&2; exit 2; }
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$ROOT/build" -j "$JOBS"
test -x "$ROOT/build/dna2bit-sag-pipeline" || { echo "dna2bit-sag-pipeline was not produced" >&2; exit 3; }
mkdir -p "$PREFIX/bin" "$PREFIX/share/dna2bit-sag-tractor"
for name in dna2bit-sag-pipeline dna2bit-embedded-sketch dna2bit-packed-search dna2bit-pack-builder gtdb-ani-af sag-stage3b-tractor; do
  if [ -x "$ROOT/build/$name" ]; then install -m 0755 "$ROOT/build/$name" "$PREFIX/bin/$name"; fi
done
cp -a "$ROOT/python" "$ROOT/tools" "$ROOT/docs" "$PREFIX/share/dna2bit-sag-tractor/"
cp "$ROOT/config/paths.env.example" "$PREFIX/share/dna2bit-sag-tractor/paths.env.example"
printf 'Installed to %s\n' "$PREFIX"
printf 'Main executable: %s/bin/dna2bit-sag-pipeline\n' "$PREFIX"
printf 'Set PATH with: export PATH="%s/bin:$PATH"\n' "$PREFIX"
