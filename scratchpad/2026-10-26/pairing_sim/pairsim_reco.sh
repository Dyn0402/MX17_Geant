#!/bin/bash
# pairsim_reco.sh <code.tgz> <sim_dir> <tag> <bundle_dir> <limit> <outdir>
# One S1 sim point through the production wft reconstruction (mx17_sim_wft/t13_reco.py), production env
# (WFT_COL_GUARD=0.01, WFT_PIVOT_START=1). Inputs and outputs on EOS; nothing large left in the scratch dir.
CODE=$1; SIMD=$2; TAG=$3; BUNDLE=$4; LIMIT=$5; OUTD=$6
set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc14-opt/setup.sh >/dev/null 2>&1
export WFT_COL_GUARD=0.01 WFT_PIVOT_START=1 MPLBACKEND=Agg
W=${_CONDOR_SCRATCH_DIR:-/tmp/$USER/pairsim_$$}; mkdir -p $W/code $W/sim/run/$TAG/decoded_root $W/out
export X17_BENCH_ROOT=$W/bench X17_ANALYSIS_ROOT=$W/ana
tar xzf $CODE -C $W/code
cp $SIMD/*_0[78].root $W/sim/run/$TAG/decoded_root/
cp -r $BUNDLE $W/bundle
mkdir -p $W/bench/det3/mx17_det3_saturday_scan_6-27-26 && cp $(dirname $CODE)/run_config_sat_det3.json $W/bench/det3/mx17_det3_saturday_scan_6-27-26/run_config.json
NAME=$(basename $SIMD)_$TAG
env -u PYTHONPATH -u PYTHONHOME /eos/user/d/dneff/venvs/ntof-lock-cvmfs/bin/python $W/code/mx17_sim_wft/t13_reco.py \
   --decoded-dir $W/sim/run/$TAG --bundle $W/bundle --out $W/out/events_$NAME.parquet \
   --jobs ${CPUS:-4} --limit $LIMIT > $W/out/reco_$NAME.log 2>&1
rc=$?
mkdir -p $OUTD; cp -r $W/out/* $OUTD/
rm -rf $W/code $W/sim $W/out $W/bundle $W/bench $W/ana
exit $rc
