#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
export TF_CPP_MIN_LOG_LEVEL=3
cd /mnt/d/variant-fm-benchmark
REF=refs/grch38_subset.fa
VCF=data/processed/variants_grch38.vcf
grep '^#' "$VCF" > /tmp/bench.vcf
grep -v '^#' "$VCF" | head -100 >> /tmp/bench.vcf
echo "nproc=$(nproc)"
start=$(date +%s)
spliceai -I /tmp/bench.vcf -O /tmp/bench_out.vcf -R "$REF" -A grch38 -D 50 2>/tmp/spliceai_err.log
end=$(date +%s)
echo "ELAPSED_100=$((end-start)) seconds"
echo "=== sample output (SpliceAI INFO field) ==="
grep -v '^#' /tmp/bench_out.vcf | head -5
echo "=== any errors? ==="
tail -3 /tmp/spliceai_err.log
