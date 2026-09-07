#!/bin/bash
# run_slopehunt.sh — T14 HV-slope-hunt discriminators (DIAGNOSIS-grid style,
# unconstrained search): the sim's dry-95/5 gain-vs-voltage slope
# (d lnG/dV = 0.296/10V) is ~1.5x too shallow vs det3 data (0.449/10V) --
# not just an offset, the wrong SHAPE. Three axes, one generalized points
# file/runner because they all reduce to "vary one knob, hold the rest of
# mx17_aval_calib.py's invocation fixed":
#
#   1. iso-fraction ladder   -- gas composition is flowmeter-set, never
#                                assayed; more quencher steepens the slope
#   2. Penning rP A/B        -- rP grows with field, so it can move slope
#                                not just offset; literature (Bhattacharya
#                                et al. 2013/2016) puts comparable devices
#                                at rP 40-80%, not just the 30-40% band this
#                                project's dry campaigns have used
#   3. field-map shape A/B   -- does the T6 mesh map's shape (vs a plain
#                                uniform field) itself contribute slope, or
#                                is it exonerated?
#
# Points file columns: GASFILE VOLT NEV SEED TAG PENNING FIELD
#   PENNING: "auto" or a decimal rP (manual mode, gas=ar)
#   FIELD:   "mesh" (per-voltage ladder lookup) or "uniform" (no --field-map)
#
# Usage: ./run_slopehunt.sh [jobs] [ladder-dir] [points-file] [out-dir]
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS="${1:-16}"
LADDER_DIR="${2:-/media/ucla/mx17_response_sim/meshfield_ladder}"
POINTS="${3:-mx17_aval_points_slopehunt.txt}"
OUT_DIR="${4:-/media/ucla/mx17_response_sim/avalanche/results_slopehunt}"
GAS_DIR=/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables
LOG_DIR="$OUT_DIR/logs"
LABEL="DIAGNOSIS-GRID / unconstrained-slope-hunt / T14-HV-slope"

mkdir -p "$OUT_DIR" "$LOG_DIR"
source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh

run_one() {
  gasfile=$1; volt=$2; nev=$3; seed=$4; tag=$5; penning=$6; field=$7
  out="$OUT_DIR/aval_slopehunt_${tag}.json"
  if [ -f "$out" ]; then
    echo "[slopehunt] skip $tag (already done)"
    return 0
  fi

  penning_args=(--penning auto)
  if [ "$penning" != "auto" ]; then
    penning_args=(--penning manual --penning-rp "$penning" --penning-gas ar)
  fi

  field_args=()
  if [ "$field" = "mesh" ]; then
    vtag=$(printf "vmesh%04d" "$volt")
    field_map="$LADDER_DIR/meshfield_${vtag}.txt"
    if [ ! -f "$field_map" ]; then
      echo "[slopehunt] MISSING map for ${volt}V: $field_map -- skipping $tag" \
          | tee "$LOG_DIR/${tag}.log"
      return 1
    fi
    field_args=(--field-map "$field_map" --tmax-ns 500 --nbins 2500)
  fi

  python3 mx17_aval_calib.py \
      --gas-file "$GAS_DIR/$gasfile" \
      --voltage "$volt" --nev "$nev" --seed "$seed" \
      "${penning_args[@]}" "${field_args[@]}" \
      --ion-subsample 50 \
      --campaign-label "$LABEL" \
      --out "$out" > "$LOG_DIR/${tag}.log" 2>&1
  echo "[slopehunt] done $tag ($?)"
}
export -f run_one
export GAS_DIR OUT_DIR LOG_DIR LADDER_DIR LABEL

echo "[slopehunt] points: $POINTS, $JOBS-way parallel, $(wc -l < "$POINTS") slices"
awk -F'[, ]+' '{print $1, $2, $3, $4, $5, $6, $7}' "$POINTS" \
  | xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _
echo "[slopehunt] all slices submitted; check $LOG_DIR for failures"
