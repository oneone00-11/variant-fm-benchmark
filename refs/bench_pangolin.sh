#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
DB=refs/grch38_subset.gtf.db
VCF=data/processed/variants_grch38.vcf
mkdir -p data/raw/scores/pangolin
grep '^#' "$VCF" > /tmp/pbench.vcf
grep -v '^#' "$VCF" | head -100 >> /tmp/pbench.vcf
echo "pangolin --help (signature):"
pangolin --help 2>&1 | head -25
echo "=== running pangolin on 100 variants ==="
start=$(date +%s)
pangolin -d 50 /tmp/pbench.vcf "$REF" "$DB" /tmp/pbench_out 2>/tmp/pangolin_err.log || tail -20 /tmp/pangolin_err.log
end=$(date +%s)
echo "ELAPSED_100=$((end-start)) seconds"
echo "=== output files ==="; ls -l /tmp/pbench_out* 2>/dev/null
echo "=== sample output ==="; grep -v '^#' /tmp/pbench_out.vcf 2>/dev/null | head -5
