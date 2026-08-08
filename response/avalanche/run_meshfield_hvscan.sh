#!/bin/bash
# run_meshfield_hvscan.sh — T7: Ar/iso 95/5 avalanche campaign against the
# per-voltage field-map ladder (voltage_ladder.sh), one real field map per
# point instead of run_meshfield_campaign.sh's single fixed map. That
# script's campaign accidentally measured the SAME 490V physics at every
# "voltage" label (see MESHFIELD_QUARANTINE_README.md) because it only had
# one map; this one looks up meshfield_vmesh<V>.txt per line so --voltage
# is finally backed by a real per-point field.
#
# Points file default: mx17_aval_points_meshfield.txt (det3's own bench HV
# scan range, 460-530V in 10V steps -- design/RESPONSE_SIM_PLAN.md "their
# HV scan runs 460-530V" -- x8 seeds = 64 slices). Voltages in the points
# file MUST land exactly on the ladder's grid (multiples of 10 within
# [300,700] by default).
#
# Usage: ./run_meshfield_hvscan.sh [jobs] [ladder-dir] [points-file] [out-dir]
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS="${1:-8}"
LADDER_DIR="${2:-/media/ucla/mx17_response_sim/meshfield_ladder}"
POINTS="${3:-mx17_aval_points_meshfield.txt}"
OUT_DIR="${4:-/media/ucla/mx17_response_sim/avalanche/results_meshfield_hvscan}"
GAS_DIR=/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"
source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh

run_one() {
  gasfile=$1; volt=$2; nev=$3; seed=$4; tag=$5
  out="$OUT_DIR/aval_meshfield_${tag}.json"
  if [ -f "$out" ]; then
    echo "[hvscan] skip $tag (already done)"
    return 0
  fi
  vtag=$(printf "vmesh%04d" "$volt")
  field_map="$LADDER_DIR/meshfield_${vtag}.txt"
  if [ ! -f "$field_map" ]; then
    echo "[hvscan] MISSING map for ${volt}V: $field_map -- skipping $tag" \
        | tee "$LOG_DIR/${tag}.log"
    return 1
  fi
  python3 mx17_aval_calib.py \
      --gas-file "$GAS_DIR/$gasfile" \
      --voltage "$volt" --nev "$nev" --seed "$seed" \
      --ion-subsample 50 --field-map "$field_map" \
      --tmax-ns 500 --nbins 2500 \
      --out "$out" > "$LOG_DIR/${tag}.log" 2>&1
  echo "[hvscan] done $tag ($?)"
}
export -f run_one
export GAS_DIR OUT_DIR LOG_DIR LADDER_DIR

echo "[hvscan] ladder: $LADDER_DIR, points: $POINTS, $JOBS-way parallel, $(wc -l < "$POINTS") slices"
awk -F'[, ]+' '{print $1, $2, $3, $4, $5}' "$POINTS" \
  | xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _
echo "[hvscan] all slices submitted; check $LOG_DIR for failures"
