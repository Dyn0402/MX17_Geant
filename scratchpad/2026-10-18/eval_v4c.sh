#!/bin/bash
# Yoga driver (2026-10-18): wait for condor 4480263 (v4 candidate, stageBC_w4_cand.sub), pull, then ledger / extent /
# neighbour shape / ladder / wft truth-pinned calibration on all 12 angles / closure with that bundle.
set -u
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-18; H17=$H/../2026-10-17; RT=$H/reco_truth
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v3wet; SB=~/x17/response_sim/stageB_w2; SA=/media/dylan/data/claude_tmp/stageA_e17
C=/media/dylan/data/claude_tmp/closure
ANGS="th00 thx05 thx10 thx15 thx20 thx30 thy05 thy10 thy15 thy20 thy30 txy10"
while :; do
  n=$(ssh -o BatchMode=yes lxplus 'condor_q 4480263 -af ClusterId 2>/dev/null | wc -l' 2>/dev/null || echo err)
  echo "$(date +%T) jobs left: $n"; [ "$n" = "0" ] && break; sleep 300
done
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v4c_*.out' > $H/eval_v4c_failed.txt 2>&1
for a in $ANGS; do
  mkdir -p $D/v4c_$a $SB/w2_rho5M_v4c/$a
  rsync -a "lxplus:work/mx17_s1/stageBC/v4c_$a/*.root" "lxplus:work/mx17_s1/stageBC/v4c_$a/*truth*" "lxplus:work/mx17_s1/stageBC/v4c_$a/stageBC_rho5M.log" $D/v4c_$a/ || echo "PULL FAILED v4c_$a"
  ln -sfn $D/v4c_$a $SB/w2_rho5M_v4c/$a/decoded_root
done
cd ~
$P $H/../2026-10-14/ledger_m3.py "w2_rho5M_v4c/th*" > $H/ledger_v4c.log 2>&1
PTS="v3rw v3ws420 v4c" $P $H17/xy_extent_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_extent_v4c.log
mkdir -p /media/dylan/data/claude_tmp/tmp/nb/v4c
$P $H17/dump_pm3_m3.py $D/v4c_th00 /media/dylan/data/claude_tmp/tmp/nb/v4c/waveforms_pm3.npz > /dev/null 2>&1
$P $H/nb_compare.py /media/dylan/data/claude_tmp/tmp/nb/{v3rw,v3ws420,v4c}/waveforms_pm3.npz > $H/nb_compare_v4c.log 2>&1
cd $RT
V_SIM=37.18 $P ladder.py out_ladder_v4c $D v4c > ladder_v4c.log 2>&1
$P sim_calib.py $C/bundle_v4c $C/calib_work_v4c $D v4c $SA --jobs 8 --v 37.18 > sim_calib_v4c.log 2>&1
$P closure.py $C $C/bundle_v4c $D v4c $SA --tag simcal --jobs 8 > closure_v4c_simcal.log 2>&1
$P closure_table.py $C v4c simcal > closure_table_v4c_simcal.md 2>&1
echo "$(date +%T) EVAL DONE"
