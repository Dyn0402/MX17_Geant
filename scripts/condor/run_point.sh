#!/bin/bash
# One (rho_s, d_k) S1 point at ny=1024 -> EOS. Args: RHO_MOHM DK_UM
set -euo pipefail
RHO="$1"; DK="$2"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
EOSDIR=/eos/experiment/ntof/data/x17/response_sim/s1_ny1024
WORK="${TMPDIR:-/tmp}/s1_$$"
mkdir -p "$WORK"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
# The solve is BLAS-heavy; condor gives one core, so stop OpenBLAS from
# spawning threads it will only fight itself for.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "host $(hostname)  rho=${RHO} MOhm/sq  dk=${DK} um  ny=1024  threads=$OMP_NUM_THREADS"
/usr/bin/time -v python3 -u -m response.solver.kernels \
    --rho "$RHO" --dk "$DK" --ny 1024 --outdir "$WORK" 2>&1

f=$(ls "$WORK"/greens_comb_*.npz)
echo "produced $(basename "$f") $(stat -c%s "$f") bytes"
eos mkdir -p "$EOSDIR" 2>/dev/null || true
xrdcp -f "$f" "root://eosuser.cern.ch/${EOSDIR}/$(basename "$f")"
echo "uploaded to $EOSDIR"
rm -rf "$WORK"
