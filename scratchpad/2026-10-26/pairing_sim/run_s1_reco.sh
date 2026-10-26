#!/bin/bash
# NOTE: superseded by condor (pairsim_reco.sh + s1_reco.sub): locally the det3 phys bundle fits at ~5 s/event and
# the run also lacked the production env (WFT_COL_GUARD=0.01 WFT_PIVOT_START=1). Kept for reference only.
# S1 (single cosmic muons, det3 v4c conditions) through the production wft reconstruction with the det3 golden-row
# phys bundle (t0 prior dropped: the sim's trigger timing is its own). One reco at a time (flock), --jobs 2.
# usage: run_s1_reco.sh BUNDLE TAG [angles...]
NT=~/PycharmProjects/nTof_x17; PY=$NT/.venv/bin/python
SIM=/media/dylan/data/claude_tmp/v3wet; OUT=~/x17/cosmic_bench/satstudy/pairing_sim/sim
B=$1; TAG=$2; shift 2
for ang in "$@"; do
  lay=$OUT/layout/v4c_$ang/$TAG; mkdir -p $lay
  [ -e $lay/decoded_root ] || ln -s $SIM/v4c_$ang $lay/decoded_root
  ev=$OUT/events_v4c_${ang}_$TAG.parquet
  [ -e $ev ] && continue
  flock /home/dylan/.claude/jobs/d934b9ce/tmp/heavy.lock $PY $NT/mx17_sim_wft/t13_reco.py --decoded-dir $lay \
     --bundle $B --out $ev --jobs 2 --limit ${LIMIT:-1500} > $OUT/reco_v4c_${ang}_$TAG.log 2>&1
  echo "$ang done $(date +%T)" >> $OUT/run_s1_$TAG.progress
done
echo ALL DONE >> $OUT/run_s1_$TAG.progress
