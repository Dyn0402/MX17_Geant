#!/bin/bash
# LUT caching cert ONLY, under the FIXED harness (f1daf7a: compare at the
# last covered source time + misaligned-axes guard). Args: TAG EOSDIR PROD
set -euo pipefail
TAG="$1"; EOSDIR="$2"; PROD="$3"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
OUTD=/afs/cern.ch/work/d/dneff/mx17_s1/w2_certs
WORK="${TMPDIR:-/tmp}/lutfix_$$"
mkdir -p "$WORK" "$OUTD"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
echo "host $(hostname)  tag=$TAG prod=$PROD"
xrdcp -s -f "root://eosuser.cern.ch/${EOSDIR}/${PROD}" "$WORK/$PROD"
echo "pulled $PROD $(stat -c%s "$WORK/$PROD") bytes"
python3 -u -m response.digitizer.test_lut_vs_solver \
    --kernel "$WORK/$PROD" 2>&1 | tee "$OUTD/lutfix_${TAG}.log"
rm -rf "$WORK"
