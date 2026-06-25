#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
DB=refs/grch38_subset.gtf.db
echo "CHROM,POS,REF,ALT" > /tmp/pbench.csv
grep -v '^#' data/processed/variants_grch38.vcf | head -100 | awk 'BEGIN{OFS=","}{print $1,$2,$4,$5}' >> /tmp/pbench.csv
echo "=== running pangolin (CSV mode) on 100 variants ==="
start=$(date +%s)
pangolin -d 50 /tmp/pbench.csv "$REF" "$DB" /tmp/pbench_out 2>/tmp/perr.log || { echo "FAILED:"; tail -25 /tmp/perr.log; exit 1; }
end=$(date +%s)
echo "ELAPSED_100=$((end-start)) seconds"
echo "=== output head ==="
head -6 /tmp/pbench_out.csv
