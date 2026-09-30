#!/usr/bin/env bash
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PREFIX=${PREFIX:-"$HOME/.local"}
JOBS=${JOBS:-8}
command -v cmake >/dev/null 2>&1 || { echo "cmake is required" >&2; exit 2; }
command -v c++ >/dev/null 2>&1 || { echo "a C++ compiler is required" >&2; exit 2; }
cmake_args=(-S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX")
if [ -n "${CONDA_PREFIX:-}" ]; then cmake_args+=("-DHTS_ROOT=$CONDA_PREFIX"); fi
if [ "${MICROSAGS_PORTABLE_BUILD:-0}" = 1 ]; then
  cmake_args+=("-DCMAKE_SKIP_RPATH=TRUE")
fi
cmake "${cmake_args[@]}"
cmake --build "$ROOT/build" -j "$JOBS"
test -x "$ROOT/build/dna2bit-sag-pipeline" || { echo "dna2bit-sag-pipeline was not produced" >&2; exit 3; }
cmake --install "$ROOT/build"

# Record the dependency runtime used for this installation. The launcher
# reads it before starting child programs, so users do not need to activate
# the build environment again just to run Microsags.
runtime_file="$PREFIX/share/microsags/runtime.env"
mkdir -p "$(dirname "$runtime_file")"
if [ -n "${CONDA_PREFIX:-}" ]; then
  if [ "${MICROSAGS_PORTABLE_BUILD:-0}" = 1 ]; then
    printf 'RUNTIME_BIN=@PREFIX@/bin\nRUNTIME_LIB=@PREFIX@/lib\n' > "$runtime_file"
  else
    {
      printf 'RUNTIME_BIN=%s/bin\n' "$CONDA_PREFIX"
      printf 'RUNTIME_LIB=%s/lib\n' "$CONDA_PREFIX"
    } > "$runtime_file"
  fi
else
  : > "$runtime_file"
fi
printf 'Installed to %s\n' "$PREFIX"
printf 'Main executable: %s/bin/microsags\n' "$PREFIX"
printf 'Set PATH with: export PATH="%s/bin:$PATH"\n' "$PREFIX"
