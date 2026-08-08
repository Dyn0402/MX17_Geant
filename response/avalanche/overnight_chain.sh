#!/bin/bash
# overnight_chain.sh — orchestrates the unattended overnight pipeline:
# wait for the already-running voltage_ladder.sh -> shuttle the ladder to
# EOS -> run the Ar/iso 95/5 HV scan against it -> merge -> shuttle the HV
# scan + merged calib to EOS. Meant to survive the launching session going
# idle; run it detached (nohup ... & disown).
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

LADDER_DIR=/media/ucla/mx17_response_sim/meshfield_ladder
echo "[chain] $(date): waiting for $LADDER_DIR/.done ..."
# A pgrep-on-process-name wait is fragile: the launching "nohup ... &
# disown" wrapper shell itself carries "voltage_ladder.sh" in its cmdline
# and can outlive the real script (hit for real 2026-08-08 -- a stale
# orphaned wrapper kept pgrep matching after the ladder had exited, caught
# and manually killed by a parallel audit session before it could hang this
# chain forever). A marker file the script itself touches on completion
# has no such ambiguity.
while [ ! -f "$LADDER_DIR/.done" ]; do sleep 60; done
echo "[chain] $(date): ladder finished"

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
