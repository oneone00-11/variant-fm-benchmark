#!/bin/bash
# Fold all completed *pchunk_*.out.csv into scored_so_far, rebuild remaining,
# then resume scoring the rest. Safe to re-run after any interruption.
bash /mnt/d/variant-fm-benchmark/refs/pangolin_recover.sh
bash /mnt/d/variant-fm-benchmark/refs/run_pangolin_resume.sh
