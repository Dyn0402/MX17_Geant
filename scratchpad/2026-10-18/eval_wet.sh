#!/bin/bash
# Yoga driver (2026-10-18): wait for condor 4480262 (wet gas + prompt anode spread, stageBC_w3_wetspread.sub), pull,
# then extent / ledger / neighbour shape / ladder / wft truth-pinned calibration for each arm, against v3 (rewired).
# The new points are written on det3's wiring (--wiring connector-inverted), so no rewire step.
set -u
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-18; H17=$H/../2026-10-17; RT=$H/reco_truth
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v3wet; SB=~/x17/response_sim/stageB_w2; SA=/media/dylan/data/claude_tmp/stageA_e17
KS="v3w v3ws300 v3ws420"
mkdir -p $D
while :; do
  n=$(ssh -o BatchMode=yes lxplus 'condor_q 4480262 -af ClusterId 2>/dev/null | wc -l' 2>/dev/null || echo err)
  echo "$(date +%T) jobs left: $n"; [ "$n" = "0" ] && break; sleep 300
done
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v3w*.out' > $H/eval_wet_failed.txt 2>&1
for k in $KS; do for a in th00 thx10 thx20 thy10 thy20; do
  mkdir -p $D/${k}_$a $SB/w2_rho5M_$k/$a
  rsync -a "lxplus:work/mx17_s1/stageBC/${k}_$a/*.root" "lxplus:work/mx17_s1/stageBC/${k}_$a/*truth*" "lxplus:work/mx17_s1/stageBC/${k}_$a/stageBC_rho5M.log" $D/${k}_$a/ || echo "PULL FAILED ${k}_$a"
  ln -sfn $D/${k}_$a $SB/w2_rho5M_$k/$a/decoded_root
done; done
for a in th00 thx10 thx20 thy10 thy20; do ln -sfn /media/dylan/data/claude_tmp/v3rw/v3k265b45_$a $D/v3rw_$a; done
grep -h "gas  \|avalanche\|wiring" $D/v3w*_th00/stageBC_rho5M.log > $H/wet_stageB_heads.log
cd ~
$P $H/../2026-10-14/ledger_m3.py "w2_rho5M_v3w*/th*" > $H/ledger_wet.log 2>&1
PTS="v3rw $KS" $P $H17/xy_extent_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_extent_wet.log
PTS="v3rw $KS" $P $H17/xy_balance_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_balance_wet.log
for k in v3rw $KS; do
  mkdir -p /media/dylan/data/claude_tmp/tmp/nb/$k
  $P $H17/dump_pm3_m3.py $D/${k}_th00 /media/dylan/data/claude_tmp/tmp/nb/$k/waveforms_pm3.npz > /dev/null 2>&1
done
$P $H/nb_compare.py /media/dylan/data/claude_tmp/tmp/nb/{v3rw,v3w,v3ws300,v3ws420}/waveforms_pm3.npz > $H/nb_compare_wet.log 2>&1
cd $RT
V_SIM=37.18 $P ladder.py out_ladder_wet $D $KS > ladder_wet.log 2>&1
for k in v3w v3ws420; do
  $P sim_calib.py /media/dylan/data/claude_tmp/closure/bundle_$k /media/dylan/data/claude_tmp/closure/calib_work_$k $D $k $SA --jobs 8 --v 37.18 > sim_calib_$k.log 2>&1
done
echo "$(date +%T) EVAL DONE"
