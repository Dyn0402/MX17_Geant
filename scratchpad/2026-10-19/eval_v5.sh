#!/bin/bash
# Yoga driver (2026-10-19): wait for condor 4481160 (v5a/b/c displaced photon feedback, stageBC_w5_fbspread.sub), pull,
# then ledger / extent (5 angles) and, on th00, the prompt profile (nb_prompt.py), per-side timing (nb_side.py) and
# widths (nb_width_clean.py) against the weekend det3 data. Logs next to this script; progress in eval_v5.log.
set -u
H=/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-19; H17=$H/../2026-10-17; H14=$H/../2026-10-14
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python
D=/media/dylan/data/claude_tmp/v3wet; SB=~/x17/response_sim/stageB_w2; T=/media/dylan/data/claude_tmp/tmp/nb/tagged
ANGS="th00 thx10 thx20 thy10 thy20"; TAGS="v5a v5b v5c"
while :; do
  n=$(ssh -o BatchMode=yes lxplus 'condor_q 4481160 -af ClusterId 2>/dev/null | wc -l' 2>/dev/null || echo err)
  echo "$(date +%T) jobs left: $n"; [ "$n" = "0" ] && break; sleep 300
done
ssh -o BatchMode=yes lxplus 'grep -L "=== done" ~/work/mx17_s1/logs/stageBC_v5*.out' > $H/eval_v5_failed.txt 2>&1
for t in $TAGS; do for a in $ANGS; do
  mkdir -p $D/${t}_$a $SB/w2_rho5M_$t/$a
  rsync -a "lxplus:work/mx17_s1/stageBC/${t}_$a/*.root" "lxplus:work/mx17_s1/stageBC/${t}_$a/*truth*" "lxplus:work/mx17_s1/stageBC/${t}_$a/stageBC_rho5M.log" $D/${t}_$a/ || echo "PULL FAILED ${t}_$a"
  ln -sfn $D/${t}_$a $SB/w2_rho5M_$t/$a/decoded_root
done; done
cd ~
for t in $TAGS; do $P $H14/ledger_m3.py "w2_rho5M_$t/th*" > $H/ledger_$t.log 2>&1; done
PTS="v3rw v4c $TAGS" $P $H17/xy_extent_m3.py $D 2>&1 | grep -v "Warn\|it/s\]\|Cutting\|Pre-cut\|Recipe" > $H/xy_extent_v5.log
for t in $TAGS; do $P $H/dump_tagged.py sat_det3 $T/sim_$t.npz --sim $D/${t}_th00 --jobs 2 > /dev/null 2>&1; done
B=100 $P $H/nb_prompt.py $T/g_det3_wknd.npz $T/sim_v3w.npz $T/sim_v3ws300.npz $T/sim_v4c.npz $(for t in $TAGS; do echo $T/sim_$t.npz; done) 2>&1 | grep -v "Warn\|^  " > $H/nb_prompt_v5.log
B=100 $P $H/nb_side.py $(for t in $TAGS; do echo $T/sim_$t.npz; done) 2>&1 | grep -v "Warn\|RuntimeW\|amp = \|M = np\|z = np" > $H/nb_side_v5.log
$P $H/nb_width_clean.py $T/g_det3_wknd.npz $T/sim_v3w.npz $T/sim_v3ws300.npz $T/sim_v4c.npz $(for t in $TAGS; do echo $T/sim_$t.npz; done) 2>&1 | grep -v "Warn\|amp = " > $H/nb_width_v5.log
echo "$(date +%T) EVAL DONE"
