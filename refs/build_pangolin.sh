#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
echo "creating pangolin env..."
mamba create -y -n pangolin -c conda-forge -c bioconda \
  python=3.10 pytorch gffutils pysam pandas numpy biopython pyfaidx 2>&1 | tail -8
conda activate pangolin
echo "pip installing Pangolin from GitHub..."
pip install -q "setuptools<80"
pip install -q git+https://github.com/tkzeng/Pangolin.git 2>&1 | tail -5
echo "=== verify ==="
python -c "import pangolin, torch; print('pangolin import ok; torch', torch.__version__)" 2>&1 | tail -3
which pangolin 2>&1 || echo "pangolin CLI not on PATH"
pangolin --help 2>&1 | head -15 || echo "pangolin --help failed"
