#!/bin/bash
# Stop the currently-tracked loop + its pangolin procs (progress is banked to
# disk), then relaunch the robust loop FULLY DETACHED (setsid) so it survives
# the Claude Code session / wsl relay closing. The loop's first action re-banks
# any partial cur_pchunk output, so nothing is lost.
pkill -f run_pangolin_robust_loop 2>/dev/null
pkill -9 -f 'pangolin -d' 2>/dev/null
sleep 3
cd /mnt/d/variant-fm-benchmark
setsid bash refs/run_pangolin_robust_loop.sh > refs/loop_detached.log 2>&1 < /dev/null &
disown 2>/dev/null
sleep 2
echo "relaunched detached. log: refs/loop_detached.log"
echo "loop procs: $(pgrep -fc run_pangolin_robust_loop)"
