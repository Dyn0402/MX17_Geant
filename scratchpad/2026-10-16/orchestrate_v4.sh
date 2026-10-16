#!/bin/bash
# Runs ON THE LAPTOP (yoga): it can rsync from lxplus directly. For each point of cluster
# 4451431 (stageBC_w2_bgtp.sub): wait until the job has left the queue and its decoded files
# are on AFS, pull them, then run_ionspeed_t14.sh (T13 + T14, t14_angW arguments).
# Env CL/ANGS/TAGS/KIND override the cluster and angle list (thy add-on: CL=4451432 ANGS="thy10 thy20").
# v4 (2026-10-16): env RHO (default rho2M) for the kernel tag in paths/filenames, e.g. RHO=rho5M.
# Skips points whose t14_summary.json exists. Usage: nohup bash orchestrate_bgtp.sh > ~/orch_bgtp.log &
CL=${CL:-4451431}; SB=~/x17/response_sim/stageB_w2; AFS=/afs/cern.ch/work/d/dneff/mx17_s1/stageBC
KIND=${KIND:-bgtp}; RHO=${RHO:-rho2M}
for tag in ${TAGS:-146b25 146b35 146b45 175b25 175b35 175b45}; do for ang in ${ANGS:-th00 thx10 thx20}; do
  name=${KIND}${tag}_${ang}
  if test -s $SB/t14_${name}/t14_summary.json; then echo "skip $name"; continue; fi
  while ssh -o BatchMode=yes lxplus "condor_q $CL -af Arguments 2>/dev/null" | grep -q " ${name} "; do sleep 180; done
  if ! ssh -o BatchMode=yes lxplus "ls $AFS/${name}/*_08.root" >/dev/null 2>&1; then echo "$(date) MISSING $name (job left queue without output)"; continue; fi
  mkdir -p $SB/w2_${RHO}_${KIND}${tag}/${ang}/decoded_root
  rsync -a "lxplus:$AFS/${name}/sim_decoded_w2_${RHO}_${name}*.root" $SB/w2_${RHO}_${KIND}${tag}/${ang}/decoded_root/
  echo "$(date) pulled $name"
  for try in 1 2 3; do
    RHO=$RHO bash ~/run_ionspeed_t14.sh $tag $ang $KIND > ~/t14_${name}.log 2>&1; r=$(tail -1 ~/t14_${name}.log)
    echo "  $name try $try: $r"; case "$r" in DONE*) break;; esac; sleep 60
  done
done; done
echo ALL DONE
