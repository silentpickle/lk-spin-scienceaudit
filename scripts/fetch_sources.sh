#!/usr/bin/env bash
# Re-download the pinned arXiv LaTeX sources and check them against docs/SHA256SUMS.
# Each tarball is extracted into its own src/ folder; nothing from it is executed.
set -euo pipefail

cd "$(dirname "$0")/.."

VERSIONS=(1409.3247v1 1411.0739v1 1510.08918v2 1607.03937v1 1610.02626v2)
PAUSE_SECONDS=3   # arXiv asks automated clients to space out requests

for v in "${VERSIONS[@]}"; do
  dir="sources/arxiv/$v"
  mkdir -p "$dir/src"
  curl -sS -L -o "$dir/eprint.bin" "https://arxiv.org/e-print/$v"
  tar -xzf "$dir/eprint.bin" -C "$dir/src"
  sleep "$PAUSE_SECONDS"
done

rg 'sources/arxiv/' docs/SHA256SUMS | shasum -a 256 -c
