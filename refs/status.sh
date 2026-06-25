#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== loop attempts ==="
grep -E 'attempt|ALL DONE|STALLED|FINAL' /mnt/d/variant-fm-benchmark/refs/loop_detached.log 2>/dev/null | tail -8
echo "=== status ==="
echo "banked:   $(grep -vc '^CHROM' scored_so_far.csv 2>/dev/null)"
echo "thispass: $(cat cur_pchunk_0[0-9].out.csv 2>/dev/null | grep -vc '^CHROM')"
echo "procs:    $(pgrep -fc 'pangolin -d')"
if [ -f pangolin_out.csv ]; then echo "FINAL pangolin_out.csv: $(grep -vc '^CHROM' pangolin_out.csv) records"; else echo "FINAL: not yet"; fi
