#!/usr/bin/env bash

# Build and run the test with the given name for the specified configuration. Defaults to profile.
# The name of a test is the name if its directory under ./tests, assuming it is listed in ./tests/CMakeLists.txt.
# Since this script is minimal, it is quite strict about its arguments, basically the most recent overrides and the test name must be the last argument.

set -euo pipefail

CONFIG="profile"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    *) break ;;
  esac
done

CONFIG="${CONFIG,,}"

case "$CONFIG" in
  debug|profile|release) ;;
  *) echo "Error: unknown config '${CONFIG}'." >&2; exit 1 ;;
esac

if [ $# -eq 0 ]; then
    echo "Usage: $0 [--config <debug|profile|release>] <name>"
    echo "Example: $0 hello-world"
    echo "Example: $0 --config debug hello-world"
    exit 1
fi

cmake --build "./build/cmake-build-${CONFIG}" --target "rtk-test-$1" --parallel "$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"

EXECUTABLE="./build/cmake-build-${CONFIG}/tests/$1/rtk-test-$1"

if [ ! -x "${EXECUTABLE}" ]; then
    echo "Error: Executable not found or not executable:" >&2
    echo "       ${EXECUTABLE}" >&2
    echo "Make sure the target ${TARGET} builds correctly." >&2
    exit 1
fi

"${EXECUTABLE}"
