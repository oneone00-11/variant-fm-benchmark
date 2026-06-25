#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
echo "=== spliceai version ==="
spliceai --help 2>&1 | head -3 || echo "spliceai help failed"
python -c "import spliceai, pysam, pyfaidx, tensorflow as tf; print('spliceai ok; tf', tf.__version__)"
echo "=== bundled annotation files ==="
ann=$(python -c "import os,spliceai; print(os.path.join(os.path.dirname(spliceai.__file__),'annotations'))")
echo "ann dir: $ann"
ls -l "$ann"
echo "=== grch38 annotation head (chrom style?) ==="
head -3 "$ann/grch38.txt"
echo "=== distinct chrom names in annotation ==="
cut -f3 "$ann/grch38.txt" | sort -u | head -30
echo "=== ref download dir ==="
ls -lh /mnt/d/variant-fm-benchmark/refs/
