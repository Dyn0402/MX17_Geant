#!/bin/bash
# overnight_chain.sh — orchestrates the unattended overnight pipeline:
# wait for the already-running voltage_ladder.sh -> shuttle the ladder to
# EOS -> run the Ar/iso 95/5 HV scan against it -> merge -> shuttle the HV
# scan + merged calib to EOS. Meant to survive the launching session going
# idle; run it detached (nohup ... & disown).
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "[chain] $(date): waiting for voltage_ladder.sh to finish..."
while pgrep -f voltage_ladder.sh >/dev/null; do sleep 60; done
echo "[chain] $(date): ladder process exited"

LADDER_LOG=/media/ucla/mx17_response_sim/meshfield_ladder_run.log
if ! grep -q "all done" "$LADDER_LOG" 2>/dev/null; then
  echo "[chain] WARNING: $LADDER_LOG does not show 'all done' -- the ladder" \
       "may have crashed. Continuing anyway; check gate failures before" \
       "trusting the HV scan that follows."
fi

echo "[chain] $(date): shuttling the ladder to EOS..."
./shuttle_to_eos.sh ladder || echo "[chain] WARNING: ladder shuttle failed"

echo "[chain] $(date): launching the Ar/iso 95/5 HV scan..."
./run_meshfield_hvscan.sh 8
echo "[chain] $(date): HV scan finished"

echo "[chain] $(date): merging HV scan slices..."
( cd ~/CLionProjects/MX17_Geant && python3 -m response.avalanche.collect \
      /media/ucla/mx17_response_sim/avalanche/results_meshfield_hvscan \
      --out response/avalanche/aval_calib_meshfield_hvscan.json \
      --figdir /media/ucla/mx17_response_sim/avalanche/figs_hvscan \
) || echo "[chain] WARNING: merge failed"

echo "[chain] $(date): shuttling the HV scan + merged calib to EOS..."
./shuttle_to_eos.sh hvscan || echo "[chain] WARNING: hvscan shuttle failed"

echo "[chain] $(date): ALL DONE"
