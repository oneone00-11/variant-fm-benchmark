#!/bin/bash
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pangolin
cd /mnt/d/variant-fm-benchmark/refs
# subset GTF to our chroms (2,3,13,16,17) for a fast gffutils db build
if [ ! -f grch38_subset.gtf ]; then
  echo "subsetting GTF..."
  zcat Homo_sapiens.GRCh38.112.gtf.gz | awk 'BEGIN{FS="\t"} /^#/{print;next} $1=="2"||$1=="3"||$1=="13"||$1=="16"||$1=="17"' > grch38_subset.gtf
fi
wc -l grch38_subset.gtf
python - <<'PY'
import gffutils
print("building gffutils db (subset)...")
gffutils.create_db("grch38_subset.gtf", "grch38_subset.gtf.db",
                   force=True, keep_order=True, merge_strategy="create_unique",
                   disable_infer_genes=True, disable_infer_transcripts=True)
print("db created")
PY
ls -lh grch38_subset.gtf.db
