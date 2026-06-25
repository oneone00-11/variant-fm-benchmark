#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
export TF_CPP_MIN_LOG_LEVEL=3
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
VCF=data/processed/variants_grch38.vcf
OUT=data/raw/scores/spliceai
mkdir -p "$OUT"
rm -f "$OUT"/chunk_* "$OUT"/spliceai_out.vcf
grep '^#' "$VCF" > "$OUT/header.txt"
grep -v '^#' "$VCF" > "$OUT/body.txt"
N=6
total=$(wc -l < "$OUT/body.txt")
chunk=$(( (total + N - 1) / N ))
split -l "$chunk" -d "$OUT/body.txt" "$OUT/chunk_"
echo "split $total records into chunks of $chunk"
for cf in "$OUT"/chunk_*; do
  [[ "$cf" == *.vcf || "$cf" == *.out.vcf ]] && continue
  cat "$OUT/header.txt" "$cf" > "$cf.vcf"
  ( spliceai -I "$cf.vcf" -O "$cf.out.vcf" -R "$REF" -A grch38 -D 50 >/dev/null 2>&1; echo "finished $cf" ) &
done
wait
cat "$OUT/header.txt" > "$OUT/spliceai_out.vcf"
for of in "$OUT"/chunk_*.out.vcf; do grep -v '^#' "$of" >> "$OUT/spliceai_out.vcf"; done
echo "MERGED $(grep -vc '^#' "$OUT/spliceai_out.vcf") records -> $OUT/spliceai_out.vcf"
