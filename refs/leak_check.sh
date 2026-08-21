#!/bin/bash
REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
# Search for the actual key value across all project files/logs/caches.
# Prints only file names that contain it (should be none). Never echoes the key.
if [ -z "$ALPHAGENOME_API_KEY" ]; then echo "key not in env to check against"; exit 0; fi
HITS=$(grep -rl -F "$ALPHAGENOME_API_KEY" \
  /mnt/d/variant-fm-benchmark/data \
  /mnt/d/variant-fm-benchmark/refs \
  "$REPO/scripts" \
  "$REPO/docs" 2>/dev/null)
if [ -z "$HITS" ]; then echo "LEAK CHECK: clean (key found in 0 files)"; else echo "LEAK FOUND in:"; echo "$HITS"; fi
echo "ag_cache header: $(head -1 /mnt/d/variant-fm-benchmark/data/raw/scores/alphagenome/ag_cache.tsv)"
