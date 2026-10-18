#!/bin/bash
# New Stage A angle points (2026-10-18), same recipe as the 10-16 e17 set
# (mm_sim -g ArIso -p muon -e 4000 -n 3000 --beam-spread 31.2 -t 1), built on yoga
# build-r6404-g4114 (e17c72c+ source, Geant4 11.4.3; the 10-16 set used 11.3.2).
source ~/Software/root_6_40_04/bin/thisroot.sh
source ~/Software/geant4-11.4.3/bin/geant4.sh
B=~/CLionProjects/MX17_Geant/build-r6404-g4114/mm_sim
run(){ tag=$1; seed=$2; shift 2
  $B -g ArIso -p muon -e 4000 -n 3000 --beam-spread 31.2 "$@" -s $seed -t 1 -o mx17_muons_3k_ang_${tag}_e17 > stageA_${tag}.log 2>&1; echo "done $tag $?"; }
run thx05 101606 --theta-x 5 &
run thx15 101607 --theta-x 15 &
run thy05 101608 --theta-y 5 &
run thy15 101609 --theta-y 15 &
run thx30 101610 --theta-x 30 &
run thy30 101611 --theta-y 30 &
run txy10 101612 --theta-x 10 --theta-y 10 &
wait; echo ALL DONE
