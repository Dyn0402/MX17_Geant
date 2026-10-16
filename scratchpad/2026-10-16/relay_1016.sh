#!/bin/bash
# DESKTOP relay (2026-10-16). yoga's own ssh to lxplus has no AFS token ("Permission denied" on
# work/), so the desktop pulls each point from AFS (one at a time: the desktop disk has ~4 GB free),
# pushes it to yoga, deletes its copy, and runs T13/T14 on yoga. Then the M3 ledger, the all-event
# shapes and the Y profiles on yoga. Skips points whose t14_summary.json already exists on yoga.
# Start:  setsid nohup bash relay_1016.sh > ~/x17/relay_1016.log 2>&1 < /dev/null &
# Points: condor 4472242 (r5verdict1, rho5M), 4472241 (ionK180/205/225), 4472245 (e17b45/50/55).
AFS=/afs/cern.ch/work/d/dneff/mx17_s1/stageBC
SB=/home/dylan/x17/response_sim/stageB_w2          # on yoga
TMP=~/x17/relay_tmp; mkdir -p $TMP
ANGS="th00 thx10 thx20 thy10 thy20"
run_point() {   # kind tag rho
  local KIND=$1 TAG=$2 RHO=$3
  for ang in $ANGS; do
    name=${KIND}${TAG}_${ang}; dst=$SB/w2_${RHO}_${KIND}${TAG}/${ang}/decoded_root
    if ssh -o BatchMode=yes yoga "test -s $SB/t14_${name}/t14_summary.json"; then echo "skip $name"; continue; fi
    rm -rf $TMP/*
    if ! rsync -a "lxplus:$AFS/${name}/sim_decoded_w2_${RHO}_${name}*.root" $TMP/ 2>/dev/null; then echo "$(date) PULL FAILED $name"; continue; fi
    ssh -o BatchMode=yes yoga "mkdir -p $dst" && rsync -a $TMP/*.root yoga:$dst/ && rm -rf $TMP/*
    echo "$(date) relayed $name"
    for try in 1 2 3; do
      r=$(ssh -o BatchMode=yes yoga "RHO=$RHO bash ~/run_ionspeed_t14.sh $TAG $ang $KIND > ~/t14_${name}.log 2>&1; tail -1 ~/t14_${name}.log")
      echo "  $name try $try: $r"; case "$r" in DONE*) break;; esac; sleep 60
    done
  done
}
run_point r5verdict 1 rho5M
for t in 205 180 225; do run_point ionK $t rho2M; done
for t in 45 50 55; do run_point e17b $t rho2M; done
echo "$(date) all points relayed; analysis on yoga"
ssh -o BatchMode=yes yoga 'cd ~/PycharmProjects/nTof_x17
P=./.venv/bin/python; SB=~/x17/response_sim/stageB_w2
$P ~/ledger_m3.py "w2_rho2M_verdict1/th*"   > ~/ledger_1016_v1.log 2>&1
$P ~/ledger_m3.py "w2_rho5M_r5verdict1/th*" > ~/ledger_1016_r5.log 2>&1
$P ~/ledger_m3.py "w2_rho2M_ionK*/th*"      > ~/ledger_1016_ionK.log 2>&1
$P ~/ledger_m3.py "w2_rho2M_e17b*/th*"      > ~/ledger_1016_e17.log 2>&1
for g in "w2_rho2M_verdict1/th*" "w2_rho5M_r5verdict1/th*" "w2_rho2M_ionK*/th*" "w2_rho2M_e17b*/th*"; do $P ~/shape_all.py "$g"; done > ~/shape_1016.log 2>&1
for p in rho2M_verdict1 rho5M_r5verdict1; do for ang in th00 thy10; do
  $P ~/profile_tails.py $SB/t14_${p#*_}_$ang $SB/w2_$p/$ang > ~/profile_1016_${p}_$ang.log 2>&1
done; done'
echo "$(date) RELAY DONE"
