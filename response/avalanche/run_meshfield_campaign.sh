#!/bin/bash
# run_meshfield_campaign.sh — T7: the S3 avalanche grid against the T6 field
# map (FIELD_MAP_RUNBOOK.md), run directly on the desktop instead of through
# HTCondor: the map itself is 180 MB and does not travel as a condor
# transfer_input_file the way the uniform-field campaign's tiny .gas table
# did, so this is a local GNU-parallel-style fan-out over mx17_aval_points.txt
# instead of mx17_aval.sub.
#
# Usage: ./run_meshfield_campaign.sh [jobs] [field-map]
#   jobs        concurrent workers (default 8; desktop has 16 cores)
#   field-map   path to the map (default ../meshcell/meshfield_production.txt)
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS="${1:-8}"
FIELD_MAP="${2:-../meshcell/meshfield_production.txt}"
GAS_DIR=/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables
OUT_DIR=results_meshfield
LOG_DIR=results_meshfield/logs

mkdir -p "$OUT_DIR" "$LOG_DIR"
source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh

run_one() {
  gasfile=$1; volt=$2; nev=$3; seed=$4; tag=$5
  out="$OUT_DIR/aval_meshfield_${tag}.json"
  if [ -f "$out" ]; then
    echo "[campaign] skip $tag (already done)"
    return 0
  fi
  python3 mx17_aval_calib.py \
      --gas-file "$GAS_DIR/$gasfile" \
      --voltage "$volt" --nev "$nev" --seed "$seed" \
      --ion-subsample 50 --field-map "$FIELD_MAP" \
      --out "$out" > "$LOG_DIR/${tag}.log" 2>&1
  echo "[campaign] done $tag ($?)"
}
export -f run_one
export GAS_DIR OUT_DIR LOG_DIR FIELD_MAP

echo "[campaign] field map: $FIELD_MAP, $JOBS-way parallel, $(wc -l < mx17_aval_points.txt) slices"
awk -F'[, ]+' '{print $1, $2, $3, $4, $5}' mx17_aval_points.txt \
  | xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _
echo "[campaign] all slices submitted; check $LOG_DIR for failures"
