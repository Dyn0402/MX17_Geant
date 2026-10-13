#!/bin/bash
# yoga: pull every finished bgtp/noisecm point (job left queue + files present) into the stageB_w2 layout
SB=~/x17/response_sim/stageB_w2; AFS=/afs/cern.ch/work/d/dneff/mx17_s1/stageBC
Q=$(ssh -o BatchMode=yes lxplus "condor_q 4451431 4451432 4451433 4451436 -af Arguments 2>/dev/null" | awk '{print $2}')
for d in $(ssh -o BatchMode=yes lxplus "cd $AFS && ls -d bgtp* noisecm* verdict* 2>/dev/null"); do
  echo "$Q" | grep -qx "$d" && continue
  kind_tag=${d%_*}; ang=${d##*_}
  t=$SB/w2_rho2M_${kind_tag}/${ang}/decoded_root
  ls $t/*_08.root >/dev/null 2>&1 && continue
  mkdir -p $t && rsync -a "lxplus:$AFS/$d/sim_decoded_*.root" $t/ && echo "pulled $d"
done
