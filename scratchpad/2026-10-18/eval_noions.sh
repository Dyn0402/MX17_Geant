#!/bin/bash
# Yoga driver (2026-10-18): condor 4480264 (v3 --no-ions diagnostic). Waits for the jobs AND for eval_v4c.sh to finish
# (one CPU-heavy job at a time), pulls, then wft truth-pinned calibration + closure + table.
set -u
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-18; RT=$H/reco_truth
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v3wet; SA=/media/dylan/data/claude_tmp/stageA_e17; C=/media/dylan/data/claude_tmp/closure
while :; do
  n=$(ssh -o BatchMode=yes lxplus 'condor_q 4480264 -af ClusterId 2>/dev/null | wc -l' 2>/dev/null || echo err)
  echo "$(date +%T) jobs left: $n"; [ "$n" = "0" ] && grep -q "EVAL DONE" $H/eval_v4c.log && break; sleep 300
done
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v3ni_*.out' > $H/eval_noions_failed.txt 2>&1
for a in th00 thx05 thx10 thx15 thx20 thy05 thy10 thy15 thy20 txy10; do
  mkdir -p $D/v3ni_$a
  rsync -a "lxplus:work/mx17_s1/stageBC/v3ni_$a/*.root" "lxplus:work/mx17_s1/stageBC/v3ni_$a/*truth*" "lxplus:work/mx17_s1/stageBC/v3ni_$a/stageBC_rho5M.log" $D/v3ni_$a/ || echo "PULL FAILED v3ni_$a"
done
cd $RT
$P sim_calib.py $C/bundle_v3ni $C/calib_work_v3ni $D v3ni $SA --jobs 8 > sim_calib_v3ni.log 2>&1
$P closure.py $C $C/bundle_v3ni $D v3ni $SA --tag simcal --jobs 8 > closure_v3ni_simcal.log 2>&1
$P closure_table.py $C v3ni simcal > closure_table_v3ni_simcal.md 2>&1
echo "$(date +%T) EVAL DONE"
