#!/bin/bash
# Wait for each condor job to LEAVE THE QUEUE (its files are then complete on AFS), then run
# the laptop T13/T14 pass; skip points whose t14_summary.json already exists.
declare -A CL=([ionspeed]=4450898 [fullphys]=4450899)
JOBS="ionspeed:146 ionspeed:175 ionspeed:245 fullphys:146 fullphys:175"
for kt in $JOBS; do for ang in th00 thx10 thx20; do
  kind=${kt%%:*}; tag=${kt##*:}; name=${kind}${tag}_${ang}
  if ssh -o BatchMode=yes yoga "test -s ~/x17/response_sim/stageB_w2/t14_${name}/t14_summary.json"; then echo "skip $name"; continue; fi
  while ssh -o BatchMode=yes lxplus "condor_q ${CL[$kind]} -af Arguments 2>/dev/null" | grep -q " ${name} "; do sleep 120; done
  until ssh -o BatchMode=yes lxplus "ls /afs/cern.ch/work/d/dneff/mx17_s1/stageBC/${name}/*_08.root" >/dev/null 2>&1; do sleep 120; done
  echo "$(date) ready $name"
  st=/home/dylan/.claude/jobs/913ff27a/tmp/stage/$name; mkdir -p $st
  rsync -a "lxplus:/afs/cern.ch/work/d/dneff/mx17_s1/stageBC/${name}/sim_decoded_w2_rho2M_${name}*.root" $st/ \
    && ssh -o BatchMode=yes yoga "mkdir -p ~/x17/response_sim/stageB_w2/w2_rho2M_${kind}${tag}/${ang}/decoded_root" \
    && rsync -a $st/ yoga:x17/response_sim/stageB_w2/w2_rho2M_${kind}${tag}/${ang}/decoded_root/ && rm -rf $st
  for try in 1 2 3; do
    r=$(ssh -o BatchMode=yes yoga "bash ~/run_ionspeed_t14.sh $tag $ang $kind > ~/t14_${name}.log 2>&1; tail -1 ~/t14_${name}.log")
    echo "  $name try $try: $r"
    case "$r" in DONE*) break;; esac
    sleep 60
  done
done; done
echo ALL DONE
