#!/usr/bin/env bash
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PREFIX=${PREFIX:-"$HOME/.local"}
JOBS=${JOBS:-8}
command -v cmake >/dev/null 2>&1 || { echo "cmake is required" >&2; exit 2; }
command -v c++ >/dev/null 2>&1 || { echo "a C++ compiler is required" >&2; exit 2; }
cmake_args=(-S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX")
if [ -n "${CONDA_PREFIX:-}" ]; then cmake_args+=("-DHTS_ROOT=$CONDA_PREFIX"); fi
cmake "${cmake_args[@]}"
cmake --build "$ROOT/build" -j "$JOBS"
test -x "$ROOT/build/dna2bit-sag-pipeline" || { echo "dna2bit-sag-pipeline was not produced" >&2; exit 3; }
cmake --install "$ROOT/build"
printf 'Installed to %s\n' "$PREFIX"
printf 'Main executable: %s/bin/microsags\n' "$PREFIX"
printf 'Set PATH with: export PATH="%s/bin:$PATH"\n' "$PREFIX"
