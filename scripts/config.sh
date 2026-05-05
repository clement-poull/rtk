#!/usr/bin/env bash

# Configures the CMake project for the specified build configuration.
# If no configuration is specified, all build configurations are configured.
# Since this script is minimal, it is quite strict about its arguments, basically the most recent overrides.

set -euo pipefail

CONFIG=""
PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

CONFIGS=("debug" "profile" "release")

if [[ -n "$CONFIG" ]]; then
  CONFIGS=("$CONFIG")
fi

for CONFIG in "${CONFIGS[@]}"; do
    CONFIG="${CONFIG,,}"

    case "${CONFIG,,}" in
    debug)   CMAKE_BUILD_TYPE="Debug" ;;
    profile) CMAKE_BUILD_TYPE="RelWithDebInfo" ;;
    release) CMAKE_BUILD_TYPE="Release" ;;
    *) echo "Error: unknown config '${CONFIG}'." >&2; exit 1 ;;
    esac

    if [[ -n "${PREFIX}" ]]; then
        cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_PREFIX_PATH="${PREFIX}" -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S "." -B "./build/cmake-build-${CONFIG}"
    else
        cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S "." -B "./build/cmake-build-${CONFIG}"
    fi
done
