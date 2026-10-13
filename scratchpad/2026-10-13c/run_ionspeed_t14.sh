#!/bin/bash
# Laptop (yoga) driver: T13 reco + T14 compare for the ion-speed scan (condor 4450898),
# same arguments as the angled ladder t14_angW_* (ANGLED_LADDER_2026-08-09.md §6), so
# the s = 1 baseline is t14_angW_<ang>. Usage: run_ionspeed_t14.sh <tag> <ang>
set -euo pipefail
TAG=$1; ANG=$2
SB=~/x17/response_sim/stageB_w2
NT=~/PycharmProjects/nTof_x17
D=$SB/w2_rho2M_ionspeed$TAG/$ANG
mkdir -p $D/decoded_root
rsync -a "lxplus:/afs/cern.ch/work/d/dneff/mx17_s1/stageBC/ionspeed${TAG}_${ANG}/sim_decoded_w2_rho2M_ionspeed${TAG}_${ANG}*.root" $D/decoded_root/
cd $NT
./.venv/bin/python mx17_sim_wft/t13_reco.py --decoded-dir $D \
    --bundle $SB/data_leg_satdet3/calib_bundle --out $D/events_ionspeed${TAG}_${ANG}.parquet --jobs 6
./.venv/bin/python mx17_sim_wft/t14_compare.py --sim-dir $D --data-dir $SB/data_leg_satdet3 \
    --sim-events $D/events_ionspeed${TAG}_${ANG}.parquet --data-events $SB/data_leg_satdet3/events_data.parquet \
    --data-n-total 47452 --theta-halfwidth 3.0 --out-dir $SB/t14_ionspeed${TAG}_${ANG}
echo DONE $TAG $ANG
