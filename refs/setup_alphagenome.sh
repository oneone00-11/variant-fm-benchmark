#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
mamba create -y -n alphagenome -c conda-forge python=3.11 pip 2>&1 | tail -3
conda activate alphagenome
echo "=== pip install alphagenome ==="
pip install -q alphagenome 2>&1 | tail -6
echo "=== import + version ==="
python - <<'PY'
import alphagenome
print("alphagenome version:", getattr(alphagenome, "__version__", "?"))
from alphagenome.models import dna_client
from alphagenome.data import genome
# list variant-scoring-related entry points
print("dna_client attrs:", [a for a in dir(dna_client) if not a.startswith('_')][:30])
PY
