#!/bin/bash
# 490 V meshfield slices 49001-49004 (nev 140 each), default vs null-collision steps, same seeds.
cd /home/dylan/CLionProjects/MX17_Geant/response/avalanche
source ~/Software/root_6_40/bin/thisroot.sh >/dev/null 2>&1; source ~/Software/garfield/install/share/Garfield/setupGarfield.sh; export MX17_GARFIELD_PIN=local-install-root640
G=/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables/Ar_iC4H10_95_5_Saclay_160m.gas
for s in 49001 49002 49003 49004; do for nc in 0 1; do
  ( flag=""; [ $nc = 1 ] && flag="--null-collision"
    nice -n 19 python3 mx17_aval_calib.py --gas-file $G --voltage 490 --nev 140 --seed $s --ion-subsample 50       --field-map ../meshcell/meshfield_production.txt $flag --out /home/dylan/x17/aval_nullcoll_2026-10-16/aval_490_s${s}_nc${nc}.json > /home/dylan/x17/aval_nullcoll_2026-10-16/aval_490_s${s}_nc${nc}.log 2>&1 ) &
done; done
wait; echo AVAL DONE
