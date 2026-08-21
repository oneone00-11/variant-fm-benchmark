#!/bin/bash
REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
# Detach the AlphaGenome scorer (setsid) so it survives the session closing.
# Resumable via ag_cache.tsv (already-scored variants skipped). Key inherited
# from this process's env (forwarded via WSLENV); never written anywhere.
pkill -9 -f 91_score_alphagenome 2>/dev/null
sleep 2
setsid bash -c 'source ~/miniforge3/etc/profile.d/conda.sh; conda activate alphagenome; python "$REPO/scripts/91_score_alphagenome.py" > /mnt/d/variant-fm-benchmark/refs/ag_detached.log 2>&1' < /dev/null &
disown 2>/dev/null
sleep 2
echo "AG detached. log: refs/ag_detached.log"
echo "key present in env: $([ -n "$ALPHAGENOME_API_KEY" ] && echo yes || echo NO)"
