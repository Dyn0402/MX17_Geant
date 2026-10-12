#!/bin/bash
# Compile and run one COMSOL Java model on lxplus (or a condor node), keeping
# COMSOL's prefs/tmp/recovery off AFS.  Usage: run_lxplus.sh <Model.java> [np]
# Model parameters come from the environment (see the .java header).
set -e
CS=/cvmfs/projects.cern.ch/engtools/comsol/comsol64/multiphysics/bin/comsol
SRC=$(readlink -f "$1"); NP=${2:-4}
WD=${COMSOL_WD:-/tmp/$USER/comsol_$$}; mkdir -p "$WD"/{prefs,tmp,rec}
# Java models read their parameters from the environment, which
# COMSOL's default security preferences forbid.
cat >"$WD/prefs/comsol.prefs" <<PREFS
security.external.enable=on
security.external.runtimepermission=on
security.external.propertypermission=on
PREFS
cp "$SRC" "$WD/"; cd "$WD"
CLS=$(basename "$SRC" .java)
"$CS" compile -prefsdir "$WD/prefs" "$CLS.java" >compile.log 2>&1 || { cat compile.log; exit 1; }
"$CS" batch -inputfile "$CLS.class" -nosave -np "$NP" \
  -prefsdir "$WD/prefs" -tmpdir "$WD/tmp" -recoverydir "$WD/rec" 2>&1 | tee "$WD/batch.log"
# the model prints its table as "CSV,..." lines (COMSOL's sandbox blocks file writes)
if [ -n "$OUT" ]; then
  grep '^CSV,' "$WD/batch.log" | cut -d, -f2- >"$OUT"
  echo "wrote $OUT ($(wc -l <"$OUT") lines)"
fi
