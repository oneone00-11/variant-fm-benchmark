#!/bin/bash
# Robust Pangolin run: NO set -e (a chunk error must not abort the others),
# explicit per-chunk logs, parent blocks on wait so the job stays alive.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
DB=refs/grch38_subset.gtf.db
OUT=data/raw/scores/pangolin
mkdir -p "$OUT"
rm -f "$OUT"/pchunk_* "$OUT"/pangolin_out.csv
echo "CHROM,POS,REF,ALT" > "$OUT/all.csv"
grep -v '^#' data/processed/variants_grch38.vcf | awk 'BEGIN{OFS=","}{print $1,$2,$4,$5}' >> "$OUT/all.csv"
tail -n +2 "$OUT/all.csv" > "$OUT/body.csv"
N=8
total=$(wc -l < "$OUT/body.csv")
chunk=$(( (total + N - 1) / N ))
split -l "$chunk" -d "$OUT/body.csv" "$OUT/pchunk_"
echo "split $total into $N chunks of $chunk at $(date +%H:%M:%S)"
pids=""
for cf in "$OUT"/pchunk_[0-9][0-9]; do
  { echo "CHROM,POS,REF,ALT"; cat "$cf"; } > "$cf.csv"
  pangolin -d 50 "$cf.csv" "$REF" "$DB" "$cf.out" > "$cf.log" 2>&1 &
  pids="$pids $!"
  echo "launched $cf -> pid $!"
done
echo "waiting on:$pids"
wait
echo "all chunks returned at $(date +%H:%M:%S)"
# merge whatever completed
hdr_done=0
for of in "$OUT"/pchunk_[0-9][0-9].out.csv; do
  [ -f "$of" ] || continue
  if [ "$hdr_done" = 0 ]; then head -1 "$of" > "$OUT/pangolin_out.csv"; hdr_done=1; fi
  tail -n +2 "$of" >> "$OUT/pangolin_out.csv"
done
echo "MERGED $(( $(wc -l < "$OUT/pangolin_out.csv") - 1 )) records -> $OUT/pangolin_out.csv"
