#!/bin/bash
# run_slopehunt_chain.sh — wait for the Ar/iC4H10 92.5/7.5 gas table (needed
# by discriminator 1), then run the full slope-hunt campaign, merge, and
# shuttle to EOS. Waits on the OUTPUT FILE existing, not a process-name
# match (pgrep -f is fragile -- see voltage_ladder.sh's history: a launcher
# wrapper's own cmdline can keep matching after the real job exits).
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

GASFILE=~/PycharmProjects/nTof_x17/garfield_sim/gas_tables/Ar_iC4H10_92p5_7p5_Saclay_160m.gas
echo "[slopehunt-chain] $(date): waiting for $GASFILE ..."
while [ ! -f "$GASFILE" ]; do sleep 30; done
echo "[slopehunt-chain] $(date): gas table ready"

echo "[slopehunt-chain] $(date): running the slope-hunt campaign (168 slices)..."
./run_slopehunt.sh 16
echo "[slopehunt-chain] $(date): campaign finished"

echo "[slopehunt-chain] $(date): merging..."
( cd ~/CLionProjects/MX17_Geant && python3 -m response.avalanche.collect \
      /media/ucla/mx17_response_sim/avalanche/results_slopehunt \
      --out response/avalanche/aval_calib_slopehunt.json \
      --figdir /media/ucla/mx17_response_sim/avalanche/figs_slopehunt \
) || echo "[slopehunt-chain] WARNING: merge failed"

echo "[slopehunt-chain] $(date): shuttling to EOS..."
klist -s && rsync -av --partial /media/ucla/mx17_response_sim/avalanche/results_slopehunt/ \
    lxplus:/eos/experiment/ntof/data/x17/response_sim/avalanche/raw_slopehunt_20260809/ \
  && rsync -av ~/CLionProjects/MX17_Geant/response/avalanche/aval_calib_slopehunt.json \
    lxplus:/eos/experiment/ntof/data/x17/response_sim/avalanche/aval_calib_slopehunt.json \
  || echo "[slopehunt-chain] WARNING: shuttle failed (check Kerberos ticket)"

echo "[slopehunt-chain] $(date): ALL DONE"
