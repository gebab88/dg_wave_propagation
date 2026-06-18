#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

image="${IMAGE:-dg-linux-arm64-test}"
platform="${PLATFORM:-linux/arm64}"

docker build \
    --platform "$platform" \
    -f Dockerfile.linux-arm64 \
    -t "$image" \
    .

echo "Docker Linux arm64 test image built successfully: $image"
