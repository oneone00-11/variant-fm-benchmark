#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
cd /mnt/d/variant-fm-benchmark/refs
if [ ! -f grch38_subset.fa ]; then
  echo "concatenating chromosomes 2,3,13,16,17 ..."
  zcat Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz \
       Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz \
       Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz \
       Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz \
       Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz > grch38_subset.fa
fi
samtools faidx grch38_subset.fa
echo "=== faidx (sequence names + lengths) ==="
cat grch38_subset.fa.fai
echo "=== VCF check ==="
VCF=/mnt/d/variant-fm-benchmark/data/processed/variants_grch38.vcf
ls -l "$VCF"
echo "first records:"; grep -v '^#' "$VCF" | head -3
echo "n variants:"; grep -vc '^#' "$VCF"
