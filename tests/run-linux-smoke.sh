#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

skip_build=0
if [[ "${1:-}" == "--skip-build" ]]; then
    skip_build=1
fi

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "Smoke test must run on Linux, got $(uname -s)" >&2
    exit 1
fi

case "$(uname -m)" in
    aarch64|arm64) ;;
    *)
        echo "Smoke test must run on arm64/aarch64, got $(uname -m)" >&2
        exit 1
        ;;
esac

if [[ "$skip_build" -eq 0 ]]; then
    make cleanOBJ cleanBIN
    make Release
fi

run_dir="$(mktemp -d)"
trap 'rm -rf "$run_dir"' EXIT

cd "$run_dir"
"$repo_root/bin/Release/DiscontinousGalerkin" "$repo_root/tests/linux-smoke-config.yaml" </dev/null

frame_count="$(find . -maxdepth 1 -name '*.png' -print | wc -l | tr -d ' ')"
if [[ "$frame_count" -lt 2 ]]; then
    echo "Expected at least 2 PNG frames, got $frame_count" >&2
    exit 1
fi

if [[ ! -s output.mp4 ]]; then
    echo "Expected non-empty output.mp4" >&2
    exit 1
fi

file "$repo_root/bin/Release/DiscontinousGalerkin"
file output.mp4
echo "Linux arm64 smoke test passed with $frame_count frames."
