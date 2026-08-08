#!/bin/bash
# voltage_ladder.sh — gas-agnostic field-map voltage ladder.
#
# The FEM field map is pure electrostatics: gas never enters the solve, so
# ONE map per (V_mesh, E_drift) serves every gas mixture that was actually
# run on the bench/beam. Solving the ladder once here means later per-gas
# avalanche campaigns (Garfield, where the gas DOES matter) just pick the
# right map off the shelf instead of re-solving.
#
# Default range 300-700 V in 10 V steps (41 points) brackets the operating
# voltages actually used across gas configurations; narrow it with
# LO/HI/STEP if the real bench/beam HV records pin it down tighter.
# E_drift is held at the bench 333 V/cm for all of them (a separate map
# axis via --e-drift; not scanned tonight — see FIELD_MAP_RUNBOOK.md and
# the transparency_curve.csv follow-up for its effect on transparency).
#
# Usage: ./voltage_ladder.sh [outdir] [lo] [hi] [step] [jobs]
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

OUTDIR="${1:-/media/ucla/mx17_response_sim/meshfield_ladder}"
LO="${2:-300}"
HI="${3:-700}"
STEP="${4:-10}"
JOBS="${5:-14}"
FIELDMAP_PY=~/CLionProjects/MX17_Geant/.venv-fieldmap/bin/python

mkdir -p "$OUTDIR/logs"
echo "[ladder] V_mesh $LO..$HI step $STEP V, E_drift=333 (bench, held fixed), -> $OUTDIR"

for V in $(seq "$LO" "$STEP" "$HI"); do
  tag=$(printf "vmesh%04d" "$V")
  out_txt="$OUTDIR/meshfield_${tag}.txt"
  if [ -f "$out_txt" ]; then
    echo "[ladder] skip $tag (already done)"
    continue
  fi
  t0=$(date +%s)
  echo "[ladder] V_mesh=$V ..."
  "$FIELDMAP_PY" solve_fieldmap.py --production --v-mesh "$V" --tag "$tag" \
      --jobs "$JOBS" > "$OUTDIR/logs/${tag}.log" 2>&1
  if ! grep -q "ALL SOLVER GATES PASS" "$OUTDIR/logs/${tag}.log"; then
    echo "[ladder] *** GATE FAILURE at V_mesh=$V — see $OUTDIR/logs/${tag}.log ***"
  fi
  mv "meshfield_${tag}.txt" "meshfield_${tag}.json" "$OUTDIR/"
  echo "[ladder] $tag done in $(($(date +%s) - t0))s"
done
echo "[ladder] all done"
