#!/bin/bash
C=/mnt/d/variant-fm-benchmark/data/raw/scores/alphagenome/ag_cache.tsv
echo "procs: $(pgrep -fc 91_score_alphagenome)"
echo "cache records: $(( $(wc -l < "$C") - 1 ))"
echo "log tail:"; tail -3 /mnt/d/variant-fm-benchmark/refs/ag_detached.log 2>/dev/null
