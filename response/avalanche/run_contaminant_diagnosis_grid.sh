#!/bin/bash
# run_contaminant_diagnosis_grid.sh — DIAGNOSIS-GRID / unconstrained-contaminant-
# search. NOT a gas assay: no humidity was ever measured on the det3 bench, so
# every water figure here (including 1%, the existing best fit to a slow
# measured v_drift) is a Magboltz fit, not a measurement. This produces gain/
# sigma0/t_arrival for each contaminant candidate at the T14 bench voltage
# (490V) plus a 3-point mini-scan (480/490/500V) for the leading (1% H2O)
# variant, against the existing gas-agnostic field-map ladder. The point is to
# have this ready to interpolate FAST if/when the dry T14 default disagrees
# with data -- it does not pick or bless a winner.
#
# All points use manual Penning at rP=0.40 (the Ar/iC4H10 built-in value,
# treated as the upper bracket -- H2O/N2 IPs both sit above the Ar
# metastables, so neither opens a NEW Penning channel; see
# nTof_x17/garfield_sim/mm_config.py's Ar_iC4H10_H2O_94_5_1 entry for the full
# reasoning). Auto mode would silently run these at rP=0 against a dry
# Ar/iC4H10 95/5 reference that runs at rP=0.40 -- the exact trap that entry
# already documents.
#
# Usage: ./run_contaminant_diagnosis_grid.sh [jobs] [ladder-dir] [points-file] [out-dir]
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS="${1:-8}"
LADDER_DIR="${2:-/media/ucla/mx17_response_sim/meshfield_ladder}"
POINTS="${3:-mx17_aval_points_diagnosis_grid.txt}"
OUT_DIR="${4:-/media/ucla/mx17_response_sim/avalanche/results_diagnosis_grid}"
GAS_DIR=/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables
LOG_DIR="$OUT_DIR/logs"
LABEL="DIAGNOSIS-GRID / unconstrained-contaminant-search"

mkdir -p "$OUT_DIR" "$LOG_DIR"
source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh

run_one() {
  gasfile=$1; volt=$2; nev=$3; seed=$4; tag=$5
  out="$OUT_DIR/aval_diagnosis_grid_${tag}.json"
  if [ -f "$out" ]; then
    echo "[diagnosis-grid] skip $tag (already done)"
    return 0
  fi
  vtag=$(printf "vmesh%04d" "$volt")
  field_map="$LADDER_DIR/meshfield_${vtag}.txt"
  if [ ! -f "$field_map" ]; then
    echo "[diagnosis-grid] MISSING map for ${volt}V: $field_map -- skipping $tag" \
        | tee "$LOG_DIR/${tag}.log"
    return 1
  fi
  python3 mx17_aval_calib.py \
      --gas-file "$GAS_DIR/$gasfile" \
      --voltage "$volt" --nev "$nev" --seed "$seed" \
      --penning manual --penning-rp 0.40 --penning-gas ar \
      --ion-subsample 50 --field-map "$field_map" \
      --tmax-ns 500 --nbins 2500 \
      --campaign-label "$LABEL" \
      --out "$out" > "$LOG_DIR/${tag}.log" 2>&1
  echo "[diagnosis-grid] done $tag ($?)"
}
export -f run_one
export GAS_DIR OUT_DIR LOG_DIR LADDER_DIR LABEL

echo "[diagnosis-grid] ladder: $LADDER_DIR, points: $POINTS, $JOBS-way parallel, $(wc -l < "$POINTS") slices"
awk -F'[, ]+' '{print $1, $2, $3, $4, $5}' "$POINTS" \
  | xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _
echo "[diagnosis-grid] all slices submitted; check $LOG_DIR for failures"
