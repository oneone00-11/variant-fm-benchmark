#!/bin/bash
# Loss-proof: rebuild scored_so_far ONLY from retained output files (original
# pchunk_*.out.csv + every pass_N/c*.out.csv). Nothing is ever deleted.
cd /mnt/d/variant-fm-benchmark/data/raw/scores/pangolin
echo "CHROM,POS,REF,ALT,Pangolin" > scored_so_far.csv
{
  for f in pchunk_[0-9][0-9].out.csv pass_*/c[0-9][0-9].out.csv; do
    [ -f "$f" ] && grep -v '^CHROM' "$f"
  done
} | awk -F, 'NF>=5 && /Warnings/ && !seen[$1","$2","$3","$4]++' >> scored_so_far.csv
awk -F, 'NR==FNR{if(FNR>1)s[$1"|"$2"|"$3"|"$4]=1;next} FNR==1{print;next} !s[$1"|"$2"|"$3"|"$4]' \
    scored_so_far.csv all.csv > remaining.csv
echo "scored=$(grep -vc '^CHROM' scored_so_far.csv) remaining=$(grep -vc '^CHROM' remaining.csv)"
