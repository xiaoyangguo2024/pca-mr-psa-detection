#!/bin/bash
# Downloads the two main outcome GWAS (GWAS Catalog harmonised, GRCh38, with tabix index). Resumable.
set -u

# Portable file size and MD5 (macOS stat -f / md5 -r vs GNU stat -c / md5sum)
file_size() { stat -f %z "$1" 2>/dev/null || stat -c %s "$1"; }
md5_of() { if command -v md5 >/dev/null 2>&1; then md5 -r "$@"; else md5sum "$@"; fi; }
cd "$(dirname "$0")/../data/raw/gwas"
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
for f in \
  GCST90274001-GCST90275000/GCST90274714/harmonised/GCST90274714.h.tsv.gz \
  GCST90274001-GCST90275000/GCST90274714/harmonised/GCST90274714.h.tsv.gz.tbi \
  GCST90274001-GCST90275000/GCST90274714/harmonised/GCST90274714.h.tsv.gz-meta.yaml \
  GCST90461001-GCST90462000/GCST90461907/harmonised/GCST90461907.h.tsv.gz \
  GCST90461001-GCST90462000/GCST90461907/harmonised/GCST90461907.h.tsv.gz.tbi \
  GCST90461001-GCST90462000/GCST90461907/harmonised/md5sum.txt ; do
  out=$(basename "$f"); [ "$out" = md5sum.txt ] && out=GCST90461907.md5sum.txt
  for try in 1 2 3 4 5; do
    curl -s -S -L -C - --retry 3 -o "$out" "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/$f" && break
    echo "retry $try $out"; sleep 5
  done
  echo "$(date '+%T') done $out $(file_size "$out")"
done
curl -s -o GCST90274714.md5sum.txt "$B/GCST90274001-GCST90275000/GCST90274714/harmonised/md5sum.txt"
md5_of GCST90274714.h.tsv.gz GCST90461907.h.tsv.gz
cat GCST90274714.md5sum.txt GCST90461907.md5sum.txt
