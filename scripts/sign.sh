#!/usr/bin/env bash
# Regenerate docs/SHA256SUMS: every git-tracked file (except the manifest itself) plus the pinned arXiv sources.
# Run after `git add` and before committing, so the manifest covers exactly what is committed.
# Check later with:  shasum -a 256 -c docs/SHA256SUMS   (run scripts/fetch_sources.sh first on a fresh clone)
set -euo pipefail

cd "$(dirname "$0")/.."

MANIFEST=docs/SHA256SUMS

{
  git ls-files --cached | rg -v "^${MANIFEST}$" | sort | while IFS= read -r f; do shasum -a 256 "$f"; done
  rg --files --no-ignore sources/arxiv -g 'eprint.bin' | sort | while IFS= read -r f; do shasum -a 256 "$f"; done
} > "$MANIFEST"

shasum -a 256 -c "$MANIFEST" | rg -v ': OK$' || echo "sign: all $(wc -l < "$MANIFEST" | tr -d ' ') files OK"
