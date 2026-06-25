#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== procs alive ==="
ps -eo pid,etime,args | grep '[p]angolin -d' | sed 's#/root/miniforge3/envs/pangolin/bin/##g' | awk '{print $1, $2}'
echo "n_alive=$(pgrep -fc 'pangolin -d')"
echo "=== current done (this pass) ==="
cat cur_pchunk_0[0-7].out.csv 2>/dev/null | grep -vc '^CHROM'
echo "=== last chunk-log lines ==="
for l in cur_pchunk_0[0-7].log; do echo "-- $l --"; tail -1 "$l" 2>/dev/null; done
