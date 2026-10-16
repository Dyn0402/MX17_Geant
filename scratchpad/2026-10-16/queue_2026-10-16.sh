#!/bin/bash
# yoga queue (2026-10-16): pull + T13/T14 for condor 4472242 (rho5M verdict 1) and 4472241 (ion K0 scan),
# then the M3 ledger and the all-event shapes for both. Log: ~/queue_1016.log
ANGS="th00 thx10 thx20 thy10 thy20"
CL=4472242 KIND=r5verdict TAGS=1 RHO=rho5M ANGS="$ANGS" bash ~/orchestrate_v4.sh
CL=4472241 KIND=ionK TAGS="205 180 225" ANGS="$ANGS" bash ~/orchestrate_v4.sh
cd ~/PycharmProjects/nTof_x17
./.venv/bin/python ~/ledger_m3.py 'w2_rho2M_verdict1/th*' > ~/ledger_1016_v1.log 2>&1
./.venv/bin/python ~/ledger_m3.py 'w2_rho5M_r5verdict1/th*' > ~/ledger_1016_r5.log 2>&1
./.venv/bin/python ~/ledger_m3.py 'w2_rho2M_ionK*/th*' > ~/ledger_1016_ionK.log 2>&1
for g in 'w2_rho2M_verdict1/th*' 'w2_rho5M_r5verdict1/th*' 'w2_rho2M_ionK*/th*'; do ./.venv/bin/python ~/shape_all.py "$g"; done > ~/shape_1016.log 2>&1
SB=~/x17/response_sim/stageB_w2
for p in rho2M_verdict1 rho5M_r5verdict1; do for ang in th00 thy10; do
  ./.venv/bin/python ~/profile_tails.py $SB/t14_${p#*_}_$ang $SB/w2_$p/$ang > ~/profile_1016_${p}_$ang.log 2>&1
done; done
echo QUEUE DONE
