#!/bin/bash
# Wait for each condor output on AFS, then run the laptop T13/T14 pass for it.
JOBS="ionspeed:146 ionspeed:175 ionspeed:245 fullphys:146 fullphys:175"
for kt in $JOBS; do for ang in th00 thx10 thx20; do
  kind=${kt%%:*}; tag=${kt##*:}
  until ssh -o BatchMode=yes lxplus "ls /afs/cern.ch/work/d/dneff/mx17_s1/stageBC/${kind}${tag}_${ang}/sim_decoded_w2_rho2M_${kind}${tag}_${ang}*_08.root" >/dev/null 2>&1; do sleep 120; done
  echo "$(date) ready $kind $tag $ang"
  ssh -o BatchMode=yes yoga "bash ~/run_ionspeed_t14.sh $tag $ang $kind > ~/t14_${kind}${tag}_${ang}.log 2>&1; tail -1 ~/t14_${kind}${tag}_${ang}.log"
done; done
echo ALL DONE
