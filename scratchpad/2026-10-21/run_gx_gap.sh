#!/bin/bash
# Condor job (2026-10-21, det3 X asymmetry): X-box W2 kernel with conducting ESL gaps (cross-gap leak).
# Both X Bloch families at nx=1560 (20 um) + combine-x, for one gap conductivity fraction GF.
# Products (G_X per rho_s) -> EOS /eos/user/d/dneff/tmp/gx/gap<GF>/; scratch deleted before exit.
# Args: GF
set -euo pipefail
GF="$1"; NX=1560
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
EOSOUT=/eos/user/d/dneff/tmp/gx/gap${GF}_nx${NX}
WORK="${TMPDIR:-/tmp}/gxgap_$$"; mkdir -p "$WORK"; cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh; set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-4} MKL_NUM_THREADS=${OMP_NUM_THREADS:-4}
echo "host $(hostname) gap_frac $GF nx $NX"
for fam in 0 1; do
    python3 -u -m response.solver.w2_production family --box x --fam $fam --nx $NX --gap-frac "$GF" --outdir "$WORK/slabs" 2>&1 | grep -v "propagated"
done
python3 -u -m response.solver.w2_production combine-x --slabdir "$WORK/slabs" --outdir "$WORK/prod" --nx $NX --gap-frac "$GF"
eos mkdir -p "$EOSOUT" 2>/dev/null || true
for f in "$WORK"/prod/*.npz; do xrdcp -f "$f" "root://eosuser.cern.ch/${EOSOUT}/$(basename "$f")"; done
echo "=== done"; rm -rf "$WORK"
