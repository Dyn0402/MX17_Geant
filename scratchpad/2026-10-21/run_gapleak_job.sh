#!/bin/bash
# Condor job (2026-10-21): Stage B+C, v5a recipe with a HYBRID kernel (production Y + conducting-gap X, build_hybrid.py)
# and digitizer --gap-funnel F (prompt
# funnelling of gap-landing avalanches onto the strip edge; det3 X asymmetry seventh round); then
# raw_asym.py on X (FEU 07) and Y (FEU 08) inside the job. Decoded ROOT -> EOS only; the small npz -> AFS.
# Args: TAG KERNEL(name under /eos/user/d/dneff/tmp/gx/) FUNNEL ANGLE SEED
set -euo pipefail
TAG="$1"; KER="$2"; FUN="$3"; ANG="$4"; SEED="$5"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src; W2BASE=/afs/cern.ch/work/d/dneff/mx17_s1
EOSPROD=/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products
EOSOUT=/eos/experiment/ntof/data/x17/response_sim/stageB_w2/gapleak/${TAG}
NPZ=$W2BASE/simgapleak_npz; PROD=$KER
CL=mx17_muons_3k_ang_${ANG}_e17_t0.root
WORK="${TMPDIR:-/tmp}/gapleak_$$"; mkdir -p "$WORK/dec" "$NPZ"; cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh; set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-4}
echo "host $(hostname) tag $TAG kernel $KER funnel $FUN ang $ANG seed $SEED"
xrdcp -s -f "root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/${PROD}" "$WORK/$PROD"
cp "$W2BASE/clusters_ang_e17/$CL" "$WORK/$CL"
python3 -u -m response.digitizer.run "$WORK/$CL" \
    --run-config "$W2BASE/calib/run_config_sat_det3.json" --detector mx17_3 --kernel "$WORK/$PROD" \
    --calib "$W2BASE/calib/aval_calib_meshfield_nc1_pooled.json" \
    --noise "$W2BASE/calib/noise_det3_satscan_feu07_cmq.json" \
    --decoded-out "$WORK/dec/sim_${TAG}" --decoded-tag "${TAG}_000" \
    --ion-model reconstructed --ion-mobility K0=2.05 --drift-gap-mm 26.5 --feedback-bg 0.45 --feedback-spread-um 150 \
    --shaper-tp-ns 166 --wiring connector-inverted --seed "$SEED" --gap-funnel "$FUN" ${MAXEV:+--max-events $MAXEV} \
    --out "$WORK/stageBC_${TAG}.json" > "$WORK/stageBC_${TAG}.log" 2>&1
python3 $W2BASE/raw_asym.py "$WORK/dec" 7 "$NPZ/${TAG}_x.npz" --files 1 > "$NPZ/${TAG}_x.log" 2>&1
python3 $W2BASE/raw_asym.py "$WORK/dec" 8 "$NPZ/${TAG}_y.npz" --files 1 > "$NPZ/${TAG}_y.log" 2>&1
eos mkdir -p "$EOSOUT" 2>/dev/null || true
for f in "$WORK"/dec/*.root "$WORK"/stageBC_${TAG}.log; do xrdcp -f "$f" "root://eosuser.cern.ch/${EOSOUT}/$(basename "$f")"; done
tail -2 "$NPZ/${TAG}_x.log"; echo "=== done"
rm -rf "$WORK"
