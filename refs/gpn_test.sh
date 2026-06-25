#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai   # has htslib/tabix
URL=https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz
echo "=== which tabix ==="; which tabix
echo "=== remote tabix test: BRCA1 region 17:43044295-43044330 ==="
export HTS_S3_HOST=  # ensure plain https
tabix "$URL" 17:43044295-43044330 2>&1 | head -12
echo "exit=$?"
