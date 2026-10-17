#!/bin/bash
# Yoga copy of eval_v2.sh (2026-10-17): condor 4477507 + 4477508 already finished, so no wait loop; pull every point
# from AFS to /media/dylan/data/claude_tmp (yoga ~/x17 is on that disk too), lay it out as stageB_w2/w2_<rho>_<kind>/<ang>/decoded_root,
# then the desktop ledger / shape / X-Y balance / extent / neighbour-delay evaluation. Logs next to this script.
set -u
CLAUDE_JOB_DIR=/media/dylan/data/claude_tmp   # npz scratch for the neighbour-delay fit
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-17
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v2; SB=~/x17/response_sim/stageB_w2
mkdir -p $D
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v2r*.out' > $H/eval_v2_failed.txt 2>&1
for k in v2r2s v2r2k v2r5s v2r5k v2r5sg279 v2r5kg279 v2r5sg260 v2r5kg260; do
  rho=rho${k:3:1}M
  for a in th00 thx10 thx20 thy10 thy20; do
    mkdir -p $D/${k}_$a $SB/w2_${rho}_$k/$a
    rsync -a "lxplus:work/mx17_s1/stageBC/${k}_$a/*.root" $D/${k}_$a/ || echo "PULL FAILED ${k}_$a"
    ln -sfn $D/${k}_$a $SB/w2_${rho}_$k/$a/decoded_root
  done
done
cd ~
$P $H/../2026-10-14/ledger_m3.py "w2_rho*_v2r*/th*" > $H/ledger_v2.log 2>&1
$P $H/shape_all_desktop.py "w2_rho*_v2r*/th*" 2>&1 | grep -v "events so far" > $H/shape_v2.log
ln -sfn $SB/w2_rho2M_verdict1/th00/decoded_root $D/verdict1_th00; ln -sfn $SB/w2_rho5M_r5verdict1/th00/decoded_root $D/r5verdict1_th00
PTS="v2r2s v2r2k v2r5s v2r5k v2r5sg279 v2r5kg279 v2r5sg260 v2r5kg260" $P $H/xy_balance_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_balance_v2.log
PTS="v2r5s v2r5sg279 v2r5sg260 v2r5k v2r5kg260" $P $H/xy_extent_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_extent_v2.log
for k in v2r5s v2r5sg279 v2r5sg260 v2r5kg260; do
  mkdir -p $CLAUDE_JOB_DIR/tmp/nb/$k
  $P $H/dump_pm3_m3.py $D/${k}_th00 $CLAUDE_JOB_DIR/tmp/nb/$k/waveforms_pm3.npz > /dev/null 2>&1
  (cd $CLAUDE_JOB_DIR/tmp/nb/$k && $P $H/../2026-10-13/xtalk_fit.py) > $H/xtalk_fit_${k}_m3.log 2>&1
done
echo "$(date +%T) EVAL DONE"
