#!/bin/bash
echo "=== AlphaGenome key visible in WSL? (value never printed) ==="
if [ -n "$ALPHAGENOME_API_KEY" ]; then echo "SET (len=${#ALPHAGENOME_API_KEY})"; else echo "NOT SET"; fi
echo "=== network reachability ==="
for u in https://huggingface.co https://hf-mirror.com https://www.googleapis.com https://github.com; do
  printf '%s -> ' "$u"
  curl -s -o /dev/null -w '%{http_code}\n' --max-time 15 "$u" || echo timeout
done
echo "=== conda envs ==="
source ~/miniforge3/etc/profile.d/conda.sh
conda env list | grep -E 'spliceai|pangolin|alphagenome|base'
