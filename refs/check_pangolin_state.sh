#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== pangolin processes alive ==="
pgrep -fa pangolin | grep -v pgrep || echo "NONE running"
echo "=== records done per chunk (.out.csv) vs chunk size ==="
for i in 00 01 02 03 04; do
  done=$( [ -f pchunk_$i.out.csv ] && wc -l < pchunk_$i.out.csv || echo 0 )
  total=$( wc -l < pchunk_$i 2>/dev/null || echo 0 )
  echo "chunk_$i: done=$done total=$total"
done
echo "=== free mem ==="
free -m | sed -n 2p
