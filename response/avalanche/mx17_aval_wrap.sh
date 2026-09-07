#!/bin/bash
# mx17_aval_wrap.sh — HTCondor wrapper for one S3 avalanche-calibration point.
# Usage: mx17_aval_wrap.sh <GASFILE> <VOLTAGE> <NEV> <SEED> <TAG> [GAP_UM] [LABEL]
#
# GAP_UM and LABEL are optional and default to mx17_aval_calib.py's own
# defaults, so the original 5-argument call (mx17_aval.sub) is unchanged.
# GAP_UM was added 2026-08-11 for the 135 vs 150 µm gap scan
# (design/report/GAP_SCAN_PREREG_2026-08-11.md); it reaches the config block of
# the output JSON, which is what lets collect.py keep the two arms apart.
set -e

GAP_UM="${6:-}"
LABEL="${7:-}"

echo "[wrap] host=$(hostname) gas=$1 V=$2 nev=$3 seed=$4 tag=$5 gap=${GAP_UM:-default} start=$(date)"

source "$(dirname "${BASH_SOURCE[0]}")/setup_garfield.sh"

extra=()
[ -n "$GAP_UM" ] && extra+=(--gap-um "$GAP_UM")
[ -n "$LABEL" ] && extra+=(--campaign-label "$LABEL")

python3 -u mx17_aval_calib.py \
    --gas-file "$1" \
    --voltage "$2" \
    --nev "$3" \
    --seed "$4" \
    --ion-subsample 50 \
    "${extra[@]}" \
    --out "aval_$5.json"

# Fail loudly if the gap did not reach the record. A slice whose gap silently
# fell back to the 150 µm default looks identical to a real 150 µm slice, and
# the gap scan's whole point is that those two must never merge.
if [ -n "$GAP_UM" ]; then
    python3 - "$GAP_UM" "aval_$5.json" <<'PY'
import json, sys
want, path = float(sys.argv[1]), sys.argv[2]
got = json.load(open(path))["config"].get("gap_um")
if got is None or abs(float(got) - want) > 1e-9:
    sys.exit(f"[wrap] FATAL: asked for gap_um={want}, JSON records {got!r}")
print(f"[wrap] gap_um={got} confirmed in the output record")
PY
fi

echo "[wrap] end=$(date)"
