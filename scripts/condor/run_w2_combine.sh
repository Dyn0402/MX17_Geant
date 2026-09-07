#!/bin/bash
# W2 combine stage (plan next-in-order item 3): gather the 66 family slabs from
# EOS, ifft2 them into the 4 production greens_comb_w2_* products, push those
# back to EOS. Runs as ONE condor job rather than on the laptop because
# eosexperiment.cern.ch (the redirect target for /eos/experiment) does not
# resolve off-site — laptop xrdcp cannot read the slabs at all (2026-08-09).
#
# Sizing: _gather() allocates a (61, 512, 3120) complex128 = 1.56 GB out_hat
# plus its ifft2 temporaries, and holds 2 Y kernels (390 MB each, float32) and
# 40 X kernels (12 MB each) live per rho_s -> ~6 GB peak. It re-reads the slab
# set once per (rho_s, drive): 4x2 Y-gathers x 64 slabs + 4x40 X-gathers x 2
# slabs = 832 slab reads, so the slabs go on LOCAL scratch, never EOS fuse.
set -euo pipefail
SRC=/afs/cern.ch/work/d/dneff/mx17_s1/src
EOSBASE=/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512
WORK="${TMPDIR:-/tmp}/w2comb_$$"
SLABS="$WORK/slabs"
OUT="$WORK/products"
mkdir -p "$SLABS" "$OUT"
cd "$SRC"
export MX17_SKIP_HEADER_CHECK=1
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "host $(hostname)  work=$WORK  threads=$OMP_NUM_THREADS"
echo "=== pulling 66 slabs from EOS to local scratch"
n=0
for f in $(xrdfs root://eosuser.cern.ch ls "$EOSBASE/slabs" | xargs -n1 basename); do
    case "$f" in w2slab_*.npz) ;; *) continue ;; esac
    xrdcp -s -f "root://eosuser.cern.ch/${EOSBASE}/slabs/$f" "$SLABS/$f"
    n=$((n+1))
done
echo "pulled $n slabs, $(du -sh "$SLABS" | cut -f1)"
if [ "$n" -ne 66 ]; then echo "FATAL: expected 66 slabs, got $n"; rm -rf "$WORK"; exit 1; fi

echo "=== combine"
python3 -u -m response.solver.w2_production combine \
    --slabdir "$SLABS" --outdir "$OUT" 2>&1

echo "=== uploading products"
eos mkdir -p "$EOSBASE/products" 2>/dev/null || true
for p in "$OUT"/greens_comb_w2_*.npz; do
    echo "upload $(basename "$p") $(stat -c%s "$p") bytes"
    xrdcp -f "$p" "root://eosuser.cern.ch/${EOSBASE}/products/$(basename "$p")"
done
# Keep the products on AFS too: the laptop can rsync them over ssh, which is
# the only path that works off-site (see header).
mkdir -p /afs/cern.ch/work/d/dneff/mx17_s1/w2products
cp "$OUT"/greens_comb_w2_rho1M_dk50um_g19um.npz \
   /afs/cern.ch/work/d/dneff/mx17_s1/w2products/ 2>/dev/null || true
echo "=== done"
rm -rf "$WORK"
