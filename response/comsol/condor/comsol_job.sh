#!/bin/bash
# One COMSOL model run on an HTCondor node. Usage: comsol_job.sh <Model.java> <np> <tag> [VAR=value ...]
# The CSV goes to EOS from inside the job; nothing large is left in the scratch dir
# (condor would copy it back to AFS).
set -e
MODEL=$1; NP=$2; TAG=$3; shift 3
for kv in "$@"; do export "$kv"; done
EOSDIR=/eos/experiment/ntof/data/x17/response_sim/comsol/condor
S=${_CONDOR_SCRATCH_DIR:-$PWD}
export COMSOL_WD=$S/wd OUT=$S/$TAG.csv
./run_lxplus.sh "$MODEL" "$NP" | grep -v '^CSV,' | grep -v '^Memory\|Current Progress' | tail -40 > $TAG.tail.log || true
mkdir -p $EOSDIR
cp "$OUT" "$EOSDIR/$TAG.csv"
cp $TAG.tail.log "$EOSDIR/$TAG.tail.log"
rm -rf "$S/wd" "$OUT"
echo "done $TAG -> $EOSDIR/$TAG.csv"
