#!/usr/bin/env bash

# Builds the CMake project for the specified configuration. Defaults to profile.
# Since this script is minimal, it is quite strict about its arguments, basically the most recent overrides.

set -euo pipefail

CONFIG="profile"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

CONFIG="${CONFIG,,}"

case "$CONFIG" in
  debug|profile|release) ;;
  *) echo "Error: unknown config '${CONFIG}'." >&2; exit 1 ;;
esac

cmake --build "./build/cmake-build-${CONFIG}" --target rtk --parallel "$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
