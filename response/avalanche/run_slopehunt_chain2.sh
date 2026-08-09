#!/bin/bash
# run_slopehunt_chain2.sh — trimmed slope-hunt run (144 slices): drops the
# 92.5/7.5 iso-fraction leg per problem1_8-9-26's Townsend-interpolation
# addendum (gain and slope pull in OPPOSITE directions across 1-12% iso, so
# no fraction fixes both -- the ladder's 90/10 point stands as a cheap
# confirmation, not a real candidate). No gas-table dependency, so this runs
# immediately instead of waiting on Ar_iC4H10_92p5_7p5's Magboltz table
# (left running in the background for completeness/future use, not needed
# here).
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "[slopehunt-chain2] $(date): running the trimmed slope-hunt campaign (144 slices)..."
./run_slopehunt.sh 16 /media/ucla/mx17_response_sim/meshfield_ladder \
    mx17_aval_points_slopehunt_trimmed.txt \
    /media/ucla/mx17_response_sim/avalanche/results_slopehunt
echo "[slopehunt-chain2] $(date): campaign finished"

echo "[slopehunt-chain2] $(date): merging..."
( cd ~/CLionProjects/MX17_Geant && python3 -m response.avalanche.collect \
      /media/ucla/mx17_response_sim/avalanche/results_slopehunt \
      --out response/avalanche/aval_calib_slopehunt.json \
      --figdir /media/ucla/mx17_response_sim/avalanche/figs_slopehunt \
) || echo "[slopehunt-chain2] WARNING: merge failed"

echo "[slopehunt-chain2] $(date): shuttling to EOS..."
klist -s && rsync -av --partial /media/ucla/mx17_response_sim/avalanche/results_slopehunt/ \
    lxplus:/eos/experiment/ntof/data/x17/response_sim/avalanche/raw_slopehunt_20260809/ \
  && rsync -av ~/CLionProjects/MX17_Geant/response/avalanche/aval_calib_slopehunt.json \
    lxplus:/eos/experiment/ntof/data/x17/response_sim/avalanche/aval_calib_slopehunt.json \
  || echo "[slopehunt-chain2] WARNING: shuttle failed (check Kerberos ticket)"

echo "[slopehunt-chain2] $(date): ALL DONE"
