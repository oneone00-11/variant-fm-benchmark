#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== alive pangolin procs (pid etime chunk) ==="
ps -eo pid,etime,args | grep '[p]angolin -d' | sed 's#/root/miniforge3/envs/pangolin/bin/##g; s#/mnt/d/variant-fm-benchmark/##g'
echo "=== per-chunk done (of 2329) ==="
for i in 00 01 02 03 04 05 06 07; do
  n=$( [ -f cur_pchunk_$i.out.csv ] && grep -vc '^CHROM' cur_pchunk_$i.out.csv || echo NA )
  echo "cur_pchunk_$i: $n"
done
echo "=== tails of chunk logs (look for errors) ==="
for i in 00 01 02 03 04 05 06 07; do
  echo "--- log $i ---"; tail -2 cur_pchunk_$i.log 2>/dev/null
done
echo "=== OOM? ==="; dmesg 2>/dev/null | grep -iE "killed process|out of memory" | tail -5 || echo "(no dmesg)"
echo "=== mem ==="; free -m | sed -n 2p
echo "=== parent job script alive? ==="
pgrep -fa run_pangolin_resume || echo "parent script NOT running"
