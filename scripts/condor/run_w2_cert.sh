#!/bin/bash
# W2 fast-path recertification + T10 slow-path verdict, for ONE rho_s point.
# Arg: RHOTAG (e.g. rho1M, rho2M).
#
# Runs at CERN, not on the laptop, for two independent reasons:
#   * the laptop cannot read /eos/experiment at all (eosexperiment.cern.ch does
#     not resolve off-site — see run_w2_combine.sh), and
#   * RESPONSE_SIM_PLAN §"machine roles" forbids LUT builds on the laptop: the
#     build peaks ~3.7 GB against ~8 GB free on a machine several agents share,
#     and its documented failure mode is a SILENT OOM kill that looks like a
#     hang. Condor gives it a requested, uncontended 16 GB instead.
#
# T10 baseline being reproduced: design/report/t10/t10_prod_ny1024.json used
# greens_comb_rho2M_dk50um_g19um.npz (W1, ny=1024), --mesh-v 490 and the
# aval_calib.json whose md5 is ece7ccd5... — shipped here as
# calib/aval_calib_t10baseline.json so ONLY the boundary model changes.
#
# Optional args 2 and 3 override the product location, which is how the same
# script is dress-rehearsed against the W1 product to reproduce the published
# 8.26 % under LCG_105 before the W2 products exist (arg 1 then only names the
# output files).
set -euo pipefail
RHOTAG="$1"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
BASE=/afs/cern.ch/work/d/dneff/mx17_s1
EOSDIR="${2:-/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products}"
CALIB=$BASE/calib/aval_calib_t10baseline.json
OUTD=$BASE/w2_certs
PROD="${3:-greens_comb_w2_${RHOTAG}_dk50um_g19um.npz}"
WORK="${TMPDIR:-/tmp}/w2cert_$$"
mkdir -p "$WORK" "$OUTD"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "host $(hostname)  rho=$RHOTAG"
xrdcp -s -f "root://eosuser.cern.ch/${EOSDIR}/${PROD}" "$WORK/$PROD"
echo "pulled $PROD $(stat -c%s "$WORK/$PROD") bytes"

echo "=== meta"
python3 -u "$BASE/extract_meta.py" "$WORK/$PROD" > "$OUTD/meta_${RHOTAG}.json"
cat "$OUTD/meta_${RHOTAG}.json"

echo "=== T10 fast-path caching cert (test_lut_vs_solver, bar 2%, W1 got 1e-4)"
python3 -u -m response.digitizer.test_lut_vs_solver \
    --kernel "$WORK/$PROD" 2>&1 | tee "$OUTD/lut_vs_solver_${RHOTAG}.log"

echo "=== T10 slow path (bar 2%, W1 rho2M ny1024 got 8.26%)"
python3 -u -m response.validation.t10_slowpath \
    --kernel "$WORK/$PROD" --calib "$CALIB" --mesh-v 490 \
    --out "$OUTD/t10_w2_${RHOTAG}.json" 2>&1 \
    | tee "$OUTD/t10_w2_${RHOTAG}.log"

echo "=== done"
rm -rf "$WORK"
