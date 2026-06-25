#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spliceai
cd /mnt/d/variant-fm-benchmark/refs
echo "=== download .tbi index (2.65MB) ==="
curl -sL -o gpn_scores.tsv.bgz.tbi 'https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz.tbi'
ls -l gpn_scores.tsv.bgz.tbi
echo "=== resolve direct CDN URL for bgz (HEAD, no body) ==="
BGZ=$(curl -s -o /dev/null -w '%{url_effective}' -I -L 'https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz')
echo "resolved url length: ${#BGZ}"
echo "$BGZ" > /tmp/gpn_bgz_url.txt
echo "=== pysam fetch test: 17:43044295-43044330 ==="
python - "$BGZ" <<'PY'
import sys, pysam
url = sys.argv[1]
tbx = pysam.TabixFile(url, index="/mnt/d/variant-fm-benchmark/refs/gpn_scores.tsv.bgz.tbi")
n = 0
for row in tbx.fetch("17", 43044295, 43044330):
    print(row); n += 1
    if n >= 8: break
print("rows fetched:", n)
PY
