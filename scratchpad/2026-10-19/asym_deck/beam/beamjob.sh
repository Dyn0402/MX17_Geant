#!/bin/bash
# one (run, sub-run): raw_asym --cm on every MX17 X/Y FEU, 2 files each -> EOS user asym_beam/<run>/<sub>__<det>_<view>.npz
RUN=$1; SUB=$2; B=/eos/experiment/ntof/data/x17/july_beam/runs/$RUN; W=/afs/cern.ch/work/d/dneff/asym929
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
OUT=root://eosuser.cern.ch//eos/user/d/dneff/asym_beam/$RUN
for dv in $(python3 -c "
import json; c=json.load(open('$B/run_config.json'))
for d in c['detectors']:
    if d['name'].startswith('mx17'):
        for v in 'xy': print(f\"{d['name']}:{v}:{d['dream_feus'][v+'_1'][0]}\")"); do
  IFS=: read det v feu <<< "$dv"
  f=${SUB}__${det}_${v}.npz
  python3 $W/raw_asym.py $B/$SUB/decoded_root $feu $f --files 20 --cm > ${f%.npz}.log 2>&1
  [ -s $f ] && xrdcp -f -s $f $OUT/$f && xrdcp -f -s ${f%.npz}.log $OUT/${f%.npz}.log
  rm -f $f ${f%.npz}.log
done
