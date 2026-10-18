#!/bin/bash
# Rebuild and deploy the living note from design/report/RECO_TRUTH_2026-10-18.md.
set -e
cd ~/CLionProjects/MX17_Geant
T=${CLAUDE_JOB_DIR:-/tmp}/mx17-sim-reco-truth.html
python3 design/report/deck/md_note.py design/report/RECO_TRUTH_2026-10-18.md --out $T --title "Sim truth for track reco" \
  --description "What the det3 response simulation, with its truth, says about the cosmic and beam track reconstruction: hit-time ladder, angle closure, depth residuals. Living record."
python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py $T --slug mx17-sim-reco-truth \
  --tags "X17, simulation, micromegas, reconstruction, waveforms" \
  --summary "Response-sim truth vs the cosmic/beam track reco: hit-time ladder, wiring bug, angle closure, depth residuals (living)" --force --deploy | tail -1
