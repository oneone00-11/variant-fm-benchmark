#!/bin/bash
REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
# resolve fresh direct CDN URL (signed URLs can expire)
BGZ=$(curl -s -o /dev/null -w '%{url_effective}' -I -L 'https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz')
echo "resolved bgz url (len ${#BGZ})"
python $REPO/scripts/90_score_gpn.py "$BGZ"
