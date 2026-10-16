#!/bin/bash
# mx17_aval_nc1_wrap.sh — HTCondor wrapper for one T7 meshfield slice with null-collision steps.
# Usage: mx17_aval_nc1_wrap.sh <GASFILE> <VOLTAGE> <NEV> <SEED> <TAG> <EOS_OUT_SUBDIR>
#
# Moved here from the desktop (run_meshfield_campaign.sh NULLCOLL=1) on 2026-10-16, when /media/ucla
# began failing. The 170 MB field map is read from EOS inside the job rather than shipped, and the
# ~55 MB result JSON is written straight to EOS. Everything large is deleted before exit, so condor
# brings back only the done marker and the logs (AFS quota, see ~/.claude/CLAUDE.md).
set -e

EOS=root://eosuser.cern.ch//eos/experiment/ntof/data/x17/response_sim/avalanche
OUT="$EOS/$6"
json="aval_meshfield_$5.json"

echo "[wrap] host=$(hostname) gas=$1 V=$2 nev=$3 seed=$4 tag=$5 out=$OUT start=$(date)"

if xrdfs eosuser.cern.ch stat "/eos/experiment/ntof/data/x17/response_sim/avalanche/$6/$json" >/dev/null 2>&1; then
    echo "[wrap] $json already on EOS, skipping"
    echo skipped > "done_$5.txt"
    exit 0
fi

source "$(dirname "${BASH_SOURCE[0]}")/setup_garfield.sh"

xrdcp -f "$EOS/fieldmaps/meshfield_production.txt" meshfield_production.txt
xrdcp -f "$EOS/fieldmaps/meshfield_production.json" meshfield_production.json
# The unpatched Garfield silently treats the whole map as gas (no wires): stop rather than simulate that.
# (No pipe here: under set -e a failed grep exited the job silently, 2026-10-16.)
rc=0; python3 -u check_grid_active.py meshfield_production.txt > grid_check.log 2>&1 || rc=$?
grep check_grid_active grid_check.log || tail -5 grid_check.log
[ "$rc" = 0 ] || { echo "[wrap] FATAL: active-region check failed (rc=$rc)" >&2; tail -20 grid_check.log >&2; exit 4; }
rm -f grid_check.log

python3 -u mx17_aval_calib.py \
    --gas-file "$1" \
    --voltage "$2" \
    --nev "$3" \
    --seed "$4" \
    --ion-subsample 50 \
    --field-map meshfield_production.txt \
    --null-collision \
    --out "$json"

xrdcp -f "$json" "$OUT/$json"
# Confirm the EOS copy is complete before deleting the local one.
want=$(stat -c %s "$json")
got=$(xrdfs eosuser.cern.ch stat "/eos/experiment/ntof/data/x17/response_sim/avalanche/$6/$json" | awk '/Size:/{print $2}')
if [ "$want" != "$got" ]; then
    echo "[wrap] FATAL: EOS size $got != local $want" >&2
    exit 3
fi
rm -f "$json" meshfield_production.txt meshfield_production.json
rm -rf garfield garfield-*.tar.gz lcg109-* *.gas
echo "$want" > "done_$5.txt"
echo "[wrap] end=$(date)"
