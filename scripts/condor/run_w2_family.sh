#!/bin/bash
# One W2 Bloch-family job (plan next-in-order item 3). Args: BOX FAM
# BOX in {y, x}; FAM = family index (y: 0..63, x: 0..1).
# Slab -> EOS response_sim/s1_w2_ny512/slabs/. Combine runs later, elsewhere.
set -euo pipefail
BOX="$1"; FAM="$2"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
EOSDIR=/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/slabs
WORK="${TMPDIR:-/tmp}/w2_$$"
mkdir -p "$WORK"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
# scipy is required (dsyevd path) and the system python3 has none — LCG does.
# The LCG setup script is not `set -u`-clean (unbound COMPILER at line 18) —
# it killed the first canary pair instantly. Relax nounset around it only.
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "host $(hostname)  box=${BOX} fam=${FAM}  threads=$OMP_NUM_THREADS"
python3 -u -m response.solver.w2_production family \
    --box "$BOX" --fam "$FAM" --outdir "$WORK" 2>&1

f=$(ls "$WORK"/w2slab_*.npz)
echo "produced $(basename "$f") $(stat -c%s "$f") bytes"
eos mkdir -p "$EOSDIR" 2>/dev/null || true
xrdcp -f "$f" "root://eosuser.cern.ch/${EOSDIR}/$(basename "$f")"
echo "uploaded to $EOSDIR"
rm -rf "$WORK"
