#!/bin/bash
# Robust, cumulative, self-resuming Pangolin runner.
#  - scored_so_far.csv is the PERSISTENT accumulator (never deleted).
#  - each pass: bank any *pchunk_*.out.csv into scored_so_far (dedup),
#    recompute remaining, score remaining with 5 STAGGERED processes.
#  - loops until remaining==0 or no progress for 2 passes.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
export OMP_NUM_THREADS=3 MKL_NUM_THREADS=3
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
REF=/mnt/d/variant-fm-benchmark/refs/grch38_subset.fa
DB=/mnt/d/variant-fm-benchmark/refs/grch38_subset.gtf.db
[ -f scored_so_far.csv ] || echo "CHROM,POS,REF,ALT,Pangolin" > scored_so_far.csv
prev_rem=-1; stall=0
for attempt in $(seq 1 60); do
  # ---- cumulative bank (dedup by chrom|pos|ref|alt) ----
  : > sf.tmp
  cat *pchunk_*.out.csv 2>/dev/null | grep -v '^CHROM' | awk -F, 'NF>=5 && /Warnings/' >> sf.tmp
  grep -v '^CHROM' scored_so_far.csv >> sf.tmp
  echo "CHROM,POS,REF,ALT,Pangolin" > scored_so_far.csv
  awk -F, '!seen[$1"|"$2"|"$3"|"$4]++' sf.tmp >> scored_so_far.csv
  rm -f sf.tmp
  awk -F, 'NR==FNR{if(FNR>1)s[$1"|"$2"|"$3"|"$4]=1;next} FNR==1{print;next} !s[$1"|"$2"|"$3"|"$4]' \
      scored_so_far.csv all.csv > remaining.csv
  banked=$(grep -vc '^CHROM' scored_so_far.csv)
  nrem=$(grep -vc '^CHROM' remaining.csv)
  echo "[attempt $attempt $(date +%H:%M:%S)] banked=$banked remaining=$nrem"
  [ "$nrem" -le 0 ] && { echo "ALL DONE"; break; }
  if [ "$nrem" -eq "$prev_rem" ]; then stall=$((stall+1)); else stall=0; fi
  [ "$stall" -ge 2 ] && { echo "STALLED (no progress 2 passes) - stopping"; break; }
  prev_rem=$nrem
  # ---- score remaining: 5 staggered processes ----
  rm -f cur_pchunk_*
  tail -n +2 remaining.csv > rb.csv
  N=5; chunk=$(( (nrem + N - 1) / N ))
  split -l "$chunk" -d rb.csv cur_pchunk_
  for cf in cur_pchunk_[0-9][0-9]; do
    { echo "CHROM,POS,REF,ALT"; cat "$cf"; } > "$cf.csv"
    pangolin -d 50 "$cf.csv" "$REF" "$DB" "$cf.out" > "$cf.log" 2>&1 &
    sleep 25   # stagger startup to avoid simultaneous model-load memory spike
  done
  wait
done
cp scored_so_far.csv pangolin_out.csv
echo "FINAL pangolin_out.csv records=$(grep -vc '^CHROM' pangolin_out.csv)"
