#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== all *.out.csv on disk (records, bytes) ==="
for f in *pchunk_*.out.csv; do
  [ -f "$f" ] || continue
  echo "$f -> $(grep -vc '^CHROM' "$f") recs, $(stat -c %s "$f") bytes"
done
echo "=== scored_so_far.csv ==="
[ -f scored_so_far.csv ] && echo "$(grep -vc '^CHROM' scored_so_far.csv) records, $(stat -c %s scored_so_far.csv) bytes"
echo "=== remaining.csv ==="
[ -f remaining.csv ] && grep -vc '^CHROM' remaining.csv
echo "=== other csv files ==="
ls -la *.csv 2>/dev/null | awk '{print $5, $NF}'
