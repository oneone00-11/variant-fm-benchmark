#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
pip install -q "setuptools<80" 2>&1 | tail -1
echo "=== spliceai import ==="
python - <<'PY'
import os, spliceai
import tensorflow as tf
print("spliceai imported; tf", tf.__version__)
ann = os.path.join(os.path.dirname(spliceai.__file__), "annotations")
print("ann dir:", ann)
print("files:", os.listdir(ann))
p = os.path.join(ann, "grch38.txt")
with open(p) as f:
    head = [next(f) for _ in range(3)]
print("--- grch38.txt head ---")
print("".join(head))
import csv
chroms = set()
with open(p) as f:
    r = csv.reader(f, delimiter="\t")
    hdr = next(r)
    print("columns:", hdr)
    for row in r:
        chroms.add(row[2])
print("distinct chrom names:", sorted(chroms)[:30])
PY
echo "=== refs present ==="
ls -lh /mnt/d/variant-fm-benchmark/refs/*.gz
