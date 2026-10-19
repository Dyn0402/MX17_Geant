#!/bin/bash
# Yoga driver (2026-10-19): tagged ±5 dumps for the neighbour-asymmetry study. det3 weekend run, then other detectors
# (det2 6-22 FEU 6/8, det4 6-24 FEU 6/8, det6/det7 6-26) for the board-vs-electronics question. Progress in run_dumps.log.
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; S=$(cd $(dirname $0) && pwd)/dump_tagged.py; T=/media/dylan/data/claude_tmp/tmp/nb/tagged
cd ~
for k in g_det3_wknd g_det2 g_det4 g_det7_long g_det6_long; do
  echo "$(date +%T) start $k"; $P $S $k $T/$k.npz --jobs 6 > $T/$k.log 2>&1; echo "$(date +%T) end $k rc $?"
done
echo "$(date +%T) ALL DONE"
