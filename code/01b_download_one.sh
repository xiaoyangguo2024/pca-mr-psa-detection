#!/bin/bash
# Resumable download of one GWAS Catalog harmonised file set (used to fetch the two outcome files in parallel).
# usage: bash code/01b_download_one.sh <GCST range dir> <GCST id>
set -u

# Portable file size and MD5 (macOS stat -f / md5 -r vs GNU stat -c / md5sum)
file_size() { stat -f %z "$1" 2>/dev/null || stat -c %s "$1"; }
md5_of() { if command -v md5 >/dev/null 2>&1; then md5 -r "$@"; else md5sum "$@"; fi; }
cd "$(dirname "$0")/../data/raw/gwas"
rng="$1"; id="$2"
base="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/$rng/$id/harmonised"
for f in "$id.h.tsv.gz" "$id.h.tsv.gz.tbi" md5sum.txt; do
  out="$f"; [ "$f" = md5sum.txt ] && out="$id.md5sum.txt"
  for try in 1 2 3 4 5 6 7 8; do
    curl -s -S -L -C - --retry 3 -o "$out" "$base/$f" && break
    echo "retry $try $out"; sleep 5
  done
  echo "$(date '+%T') done $out $(file_size "$out")"
done
md5_of "$id.h.tsv.gz"; cat "$id.md5sum.txt"
