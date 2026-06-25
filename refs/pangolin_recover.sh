#!/bin/bash
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "=== pangolin procs before kill ==="
ps -eo pid,etime,args | grep '[p]angolin -d' | sed 's#/root/miniforge3/envs/pangolin/bin/##g'
echo "=== killing all pangolin processes ==="
pkill -9 -f 'pangolin -d' 2>/dev/null
sleep 2
echo "remaining procs: $(pgrep -fc 'pangolin -d' || echo 0)"
# Collect already-scored, well-formed records into one file (dedup by chrom,pos,ref,alt)
echo "CHROM,POS,REF,ALT,Pangolin" > scored_so_far.csv
for f in *pchunk_*.out.csv; do
  [ -f "$f" ] || continue
  # keep only lines with 5 comma-fields and a Warnings: tail (complete records)
  grep -v '^CHROM' "$f" | awk -F, 'NF>=5 && /Warnings/' >> scored_so_far.csv.tmp 2>/dev/null
done
sort -u scored_so_far.csv.tmp >> scored_so_far.csv 2>/dev/null
rm -f scored_so_far.csv.tmp
ndone=$(grep -vc '^CHROM' scored_so_far.csv)
echo "=== unique well-formed scored: $ndone ==="
# Build remaining.csv = all.csv minus scored keys
awk -F, 'NR==FNR{if(FNR>1)seen[$1"|"$2"|"$3"|"$4]=1; next} FNR==1{print; next} !seen[$1"|"$2"|"$3"|"$4]' \
    scored_so_far.csv all.csv > remaining.csv
nrem=$(grep -vc '^CHROM' remaining.csv)
echo "=== remaining to score: $nrem (of 21410) ==="
