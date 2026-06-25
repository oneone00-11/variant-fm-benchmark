#!/bin/bash
# Loss-proof resume: each pass writes to its OWN pass_N/ dir, never deletes
# anything. Safe to kill/sleep at any time; just run consolidate.sh then this
# again. Absolute paths so cwd doesn't matter.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
BASE=/mnt/d/variant-fm-benchmark
cd "$BASE/data/raw/scores/pangolin"
RUN=$(( $(ls -1d pass_* 2>/dev/null | wc -l) + 1 ))
mkdir -p "pass_$RUN"
tail -n +2 remaining.csv > "pass_$RUN/body.csv"
N=8
total=$(wc -l < "pass_$RUN/body.csv")
[ "$total" -eq 0 ] && { echo "nothing remaining"; exit 0; }
chunk=$(( (total + N - 1) / N ))
split -l "$chunk" -d "pass_$RUN/body.csv" "pass_$RUN/c"
echo "pass $RUN: $total remaining -> 8 chunks of $chunk at $(date +%H:%M:%S)"
for cf in "pass_$RUN"/c[0-9][0-9]; do
  { echo "CHROM,POS,REF,ALT"; cat "$cf"; } > "$cf.csv"
  pangolin -d 50 "$cf.csv" "$BASE/refs/grch38_subset.fa" "$BASE/refs/grch38_subset.gtf.db" "$cf.out" > "$cf.log" 2>&1 &
done
echo "launched 8 chunks; waiting..."
wait
echo "pass $RUN done at $(date +%H:%M:%S); scored this pass $(cat pass_$RUN/c[0-9][0-9].out.csv 2>/dev/null | grep -vc '^CHROM')"
