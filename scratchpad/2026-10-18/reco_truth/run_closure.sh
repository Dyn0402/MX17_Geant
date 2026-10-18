#!/bin/bash
# Closure runs, in sequence (one CPU-heavy job at a time on yoga).
cd ~/CLionProjects/MX17_Geant/scratchpad/2026-10-18/reco_truth
P=~/PycharmProjects/nTof_x17/.venv/bin/python; C=/media/dylan/data/claude_tmp/closure; SA=/media/dylan/data/claude_tmp/stageA_e17
R06=~/x17/cosmic_bench/Analysis/mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V/mx17_3/wft/calib_bundle_r06
$P closure.py $C $C/bundle_v3t /media/dylan/data/claude_tmp/v3trw v3t $SA --tag simcal --jobs 8 > closure_v3t_simcal.log 2>&1
$P closure.py $C $R06 /media/dylan/data/claude_tmp/v3trw v3t $SA --tag r06 --jobs 8 > closure_v3t_r06.log 2>&1
echo CLOSURE DONE >> closure_v3t_r06.log
