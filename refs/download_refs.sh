#!/bin/bash
set -e
cd /mnt/d/variant-fm-benchmark/refs
base=https://ftp.ensembl.org/pub/release-112
for c in 2 3 13 16 17; do
  echo "downloading chr$c ..."
  wget -q -c "$base/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$c.fa.gz"
done
echo "downloading GTF ..."
wget -q -c "$base/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
echo "--- done, listing ---"
ls -lh
