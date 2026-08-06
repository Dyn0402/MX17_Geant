#!/bin/bash
# mx17_aval_wrap.sh — HTCondor wrapper for one S3 avalanche-calibration point.
# Usage: mx17_aval_wrap.sh <GASFILE> <VOLTAGE> <NEV> <SEED> <TAG>
set -e

echo "[wrap] host=$(hostname) gas=$1 V=$2 nev=$3 seed=$4 tag=$5 start=$(date)"

source "$(dirname "${BASH_SOURCE[0]}")/setup_garfield.sh"

python3 -u mx17_aval_calib.py \
    --gas-file "$1" \
    --voltage "$2" \
    --nev "$3" \
    --seed "$4" \
    --ion-subsample 50 \
    --out "aval_$5.json"

echo "[wrap] end=$(date)"
