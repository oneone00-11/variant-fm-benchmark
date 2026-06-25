#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
cd /mnt/d/variant-fm-benchmark
OUT=data/raw/scores/pangolin
echo "=== chunk files + any outputs ==="
ls -la "$OUT"/pchunk_* 2>/dev/null | tail -20
echo "=== OOM kills in dmesg? ==="
dmesg 2>/dev/null | grep -iE "out of memory|killed process|oom" | tail -10 || echo "(dmesg unavailable)"
echo "=== single-process test: 200 variants, foreground ==="
head -201 "$OUT/all.csv" > /tmp/p200.csv
/usr/bin/time -v pangolin -d 50 /tmp/p200.csv refs/grch38_subset.fa refs/grch38_subset.gtf.db /tmp/p200_out 2>/tmp/p200.log || { echo "FAILED rc=$?"; tail -25 /tmp/p200.log; }
echo "=== peak RSS ==="
grep -E "Maximum resident" /tmp/p200.log 2>/dev/null
echo "=== output? ==="
head -3 /tmp/p200_out.csv 2>/dev/null
