#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
echo "=== install pyfastx ==="
pip install -q pyfastx 2>&1 | tail -1
echo "=== pangolin help ==="
pangolin --help 2>&1 | head -25
echo "=== build annotation db ==="
bash /mnt/d/variant-fm-benchmark/refs/build_pangolin_db.sh
echo "=== benchmark pangolin (100 variants) ==="
bash /mnt/d/variant-fm-benchmark/refs/bench_pangolin.sh
