#!/bin/bash
# Yoga driver (2026-10-18): wait for condor 4480075 (rho_s 2M x coupling), pull every point from AFS to
# /media/dylan/data/claude_tmp/v3, lay it out as stageB_w2/w2_rho5M_<kind>/<ang>/decoded_root, then ledger / shape /
# X-Y balance / extent. Logs next to this script. Reference points (v2r5k, v2r5kg260) are already on yoga under ../v2.
set -u
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-18; H17=$H/../2026-10-17
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v3; D2=/media/dylan/data/claude_tmp/v2; SB=~/x17/response_sim/stageB_w2
KS="v3r2 v3r2c5t15"
mkdir -p $D
while :; do
  n=$(ssh -o BatchMode=yes lxplus 'condor_q 4480075 -af ClusterId 2>/dev/null | wc -l' 2>/dev/null || echo err)
  echo "$(date +%T) jobs left: $n"; [ "$n" = "0" ] && break; sleep 600
done
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v3r2*.out' > $H/eval_v3r2_failed.txt 2>&1
for k in $KS; do
  for a in th00 thx10 thx20 thy10 thy20; do
    mkdir -p $D/${k}_$a $SB/w2_rho2M_$k/$a
    rsync -a "lxplus:work/mx17_s1/stageBC/${k}_$a/*.root" $D/${k}_$a/ || echo "PULL FAILED ${k}_$a"
    ln -sfn $D/${k}_$a $SB/w2_rho2M_$k/$a/decoded_root
  done
done
cd ~
$P $H17/../2026-10-14/ledger_m3.py "w2_rho2M_v3r2*/th*" > $H/ledger_v3r2.log 2>&1
$P $H17/shape_all_desktop.py "w2_rho2M_v3r2*/th*" 2>&1 | grep -v "events so far" > $H/shape_v3r2.log
PTS="v3k265b45 v3c5t15 $KS" $P $H17/xy_balance_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_balance_v3r2.log
PTS="v3k265b45 v3c5t15 $KS" $P $H17/xy_extent_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_extent_v3r2.log
for k in $KS; do
  mkdir -p /media/dylan/data/claude_tmp/tmp/nb/$k
  $P $H17/dump_pm3_m3.py $D/${k}_th00 /media/dylan/data/claude_tmp/tmp/nb/$k/waveforms_pm3.npz > /dev/null 2>&1
  (cd /media/dylan/data/claude_tmp/tmp/nb/$k && $P $H17/../2026-10-13/xtalk_fit.py) > $H/xtalk_fit_${k}_m3.log 2>&1
done
$P $H/nb_compare.py /media/dylan/data/claude_tmp/tmp/nb/{v3k265b45,v3c5t15,v3r2,v3r2c5t15}/waveforms_pm3.npz > $H/nb_compare_v3r2.log 2>&1
echo "$(date +%T) EVAL DONE"
