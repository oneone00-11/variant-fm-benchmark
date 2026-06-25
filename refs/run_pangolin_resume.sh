#!/bin/bash
# Resumable Pangolin: scores data/raw/scores/pangolin/remaining.csv into
# cur_pchunk_*.out.csv. If interrupted (power loss), just re-run
# pangolin_recover.sh (folds ALL *pchunk_*.out.csv into scored_so_far, rebuilds
# remaining) then this script again. No set -e.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
DB=refs/grch38_subset.gtf.db
OUT=data/raw/scores/pangolin
cd "$OUT"
rm -f cur_pchunk_*
tail -n +2 remaining.csv > rem_body.csv
N=8
total=$(wc -l < rem_body.csv)
[ "$total" -eq 0 ] && { echo "nothing remaining"; exit 0; }
chunk=$(( (total + N - 1) / N ))
split -l "$chunk" -d rem_body.csv cur_pchunk_
echo "resume: $total remaining -> $N chunks of $chunk at $(date +%H:%M:%S)"
pids=""
for cf in cur_pchunk_[0-9][0-9]; do
  { echo "CHROM,POS,REF,ALT"; cat "$cf"; } > "$cf.csv"
  pangolin -d 50 "$cf.csv" "/mnt/d/variant-fm-benchmark/$REF" "/mnt/d/variant-fm-benchmark/$DB" "$cf.out" > "$cf.log" 2>&1 &
  pids="$pids $!"
done
echo "launched:$pids ; waiting..."
wait
echo "resume chunks returned at $(date +%H:%M:%S)"
done_now=$(cat cur_pchunk_[0-9][0-9].out.csv 2>/dev/null | grep -vc '^CHROM')
echo "scored this pass: $done_now"
