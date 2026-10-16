#!/bin/bash
# yoga queue b (2026-10-16): after queue_2026-10-16.sh, pull + T13/T14 for condor 4472245
# (verdict 1 on the e17c72c Stage A clusters, bG 0.45/0.50/0.55), then the M3 ledger and shapes.
while pgrep -f queue_2026-10-16.sh >/dev/null; do sleep 120; done
CL=4472245 KIND=e17b TAGS="45 50 55" ANGS="th00 thx10 thx20 thy10 thy20" bash ~/orchestrate_v4.sh
cd ~/PycharmProjects/nTof_x17
./.venv/bin/python ~/ledger_m3.py 'w2_rho2M_e17b*/th*' > ~/ledger_1016_e17.log 2>&1
./.venv/bin/python ~/shape_all.py 'w2_rho2M_e17b*/th*' > ~/shape_1016_e17.log 2>&1
echo QUEUE B DONE
