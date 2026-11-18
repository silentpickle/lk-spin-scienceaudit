#!/usr/bin/env bash
# Download the raw data archives from the GitHub release, verify them against docs/SHA256SUMS-data.txt, and unpack
# them in place (each archive restores its files under results/). See docs/DATA.md.
# Needs curl, shasum, zstd and tar. Pass archive names to fetch only some, e.g.:  scripts/fetch_data.sh lkspin-data-step5.tar.zst
set -euo pipefail

cd "$(dirname "$0")/.."

RELEASE_URL="https://github.com/silentpickle/lk-spin-scienceaudit/releases/download/v0.1-first-results"
SUMS="docs/SHA256SUMS-data.txt"
DL_DIR="data_download"

mkdir -p "$DL_DIR"
names=("$@")
if [ "$#" -eq 0 ]; then
  while read -r _sum name; do names+=("$name"); done < "$SUMS"   # bash 3.2 (macOS) has no mapfile
fi

PART_LIMIT=100   # the largest archives are uploaded as name.part00, part01, ... (GitHub uploads of ~700 MB kept failing)

fetch_parts() {
  local name=$1 k part
  : > "$DL_DIR/$name.tmp"
  for ((k = 0; k < PART_LIMIT; k++)); do
    part=$(printf '%s.part%02d' "$name" "$k")
    curl -sS -L --fail -o "$DL_DIR/$part" "$RELEASE_URL/$part" 2>/dev/null || break
    cat "$DL_DIR/$part" >> "$DL_DIR/$name.tmp"
    rm "$DL_DIR/$part"
  done
  mv "$DL_DIR/$name.tmp" "$DL_DIR/$name"
}

for name in "${names[@]}"; do
  if [ ! -f "$DL_DIR/$name" ]; then
    curl -sS -L --fail -o "$DL_DIR/$name" "$RELEASE_URL/$name" 2>/dev/null || fetch_parts "$name"
  fi
  rg -F "  $name" "$SUMS" | (cd "$DL_DIR" && shasum -a 256 -c -)
  zstd -d -q -c "$DL_DIR/$name" | tar -xf -
  echo "unpacked $name"
done
