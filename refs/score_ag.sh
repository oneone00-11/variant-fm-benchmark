#!/bin/bash
REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alphagenome
python "$REPO/scripts/91_score_alphagenome.py"
