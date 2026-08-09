#!/bin/bash
# Stage B + C on the W2 kernels (plan §7/§8, handoff §6). Arg: RHOTAG (rho1M, rho2M).
#
# Regenerates sim_decoded ONCE per rho_s point with BOTH new inputs applied
# together — the W2 kernel and the pooled 490 V meshfield calib — so the change
# is one step, not two single-variable regenerations.
#
# Runs at CERN, not the laptop: RESPONSE_SIM_PLAN's machine-roles section
# forbids LUT builds there (~3.7 GB peak, silent-OOM failure mode), and the
# laptop cannot read /eos/experiment at all. The desktop is off limits while
# the T7 chain and its morning collection are running (mx17-geant-6b session).
#
# WHY TWO rho_s POINTS. rho_s is genuinely UNKNOWN for the real detector — that
# is why the plan scans {0.5,1,2,5} MΩ/sq. The T2b spread measurement gives
# 1.5-2.7 MΩ/sq at d_k = 75 µm and 1.1-1.9 at d_k = 50; the production stack's
# EFFECTIVE insulator is 70.5 µm (50 µm kapton + 18.76 µm glue in series), which
# interpolates to ~1.4-2.6 MΩ/sq. So rho2M sits centrally and rho1M sits below
# the data-allowed band, and "rho1M is nominal" dates from the pre-glue d_k =
# 75 µm era. Producing both is a scan point, not indecision; picking one to fit
# T14 later would be tuning, so the choice must rest on T2b, not on the
# comparison.
set -euo pipefail
RHOTAG="$1"
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
W2BASE=/afs/cern.ch/work/d/dneff/mx17_s1     # NOT `BASE`: LCG's setup.sh exports that
EOSPROD=/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products
EOSCLUS=/eos/experiment/ntof/data/x17/response_sim/clusters
EOSOUT=/eos/experiment/ntof/data/x17/response_sim/stageB_w2
CLUSTERS=mx17_muons_3k_spread_t0.root
PROD=greens_comb_w2_${RHOTAG}_dk50um_g19um.npz
OUTD=$W2BASE/stageBC
WORK="${TMPDIR:-/tmp}/stagebc_$$"
mkdir -p "$WORK" "$OUTD"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS

echo "host $(hostname)  rho=$RHOTAG"
xrdcp -s -f "root://eosuser.cern.ch/${EOSPROD}/${PROD}" "$WORK/$PROD"
xrdcp -s -f "root://eosuser.cern.ch/${EOSCLUS}/${CLUSTERS}" "$WORK/$CLUSTERS"
echo "pulled kernel $(stat -c%s "$WORK/$PROD") B, clusters $(stat -c%s "$WORK/$CLUSTERS") B"

# --feu-ids 7 8, NOT the 3/4 default. FEU assignment is per RUN, not per
# detector: mx17_3 sits on 3/4 in the 6-25/6-26 runs but on 7/8 in the 6-27
# saturday scan, which is the T14 target (long_run_resist_490V_drift_1000V)
# per the §0a P1 decision. Verified in that run's own run_config.json
# (detectors[0].dream_feus: x_* -> 7, y_* -> 8) and in qa_config's sat_det3
# (MX17_FEU_X 7, MX17_FEU_Y 8). Wrong ids cost no error — wft finds no files
# and reconstructs zero events.
python3 -u -m response.digitizer.run "$WORK/$CLUSTERS" \
    --feu-ids 7 8 \
    --kernel "$WORK/$PROD" \
    --calib "$SRC/response/avalanche/aval_calib_meshfield_pooled.json" \
    --noise "$W2BASE/calib/noise_det3.json" \
    --decoded-out "$WORK/sim_decoded_w2_${RHOTAG}" \
    --decoded-tag "w2${RHOTAG}_000" \
    --out "$OUTD/stageBC_${RHOTAG}.json" 2>&1 | tee "$OUTD/stageBC_${RHOTAG}.log"

echo "=== outputs"
ls -la "$WORK"/sim_decoded_w2_${RHOTAG}*.root
eos mkdir -p "$EOSOUT" 2>/dev/null || true
for f in "$WORK"/sim_decoded_w2_${RHOTAG}*.root; do
    xrdcp -f "$f" "root://eosuser.cern.ch/${EOSOUT}/$(basename "$f")"
    # AFS copy too: T13 runs wft on the laptop, which cannot reach EOS.
    cp "$f" "$OUTD/" || true
done
echo "=== done"
rm -rf "$WORK"
