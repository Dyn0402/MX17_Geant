#!/bin/bash
# InduceAxi batch on lxplus: two solves at a time, CSVs in ~/MX17_comsol.
# Usage: axi_batch.sh "<RHO_S> <CASE> <VE_UM_NS> <tag>" ...
# Extra env (H_COL_UM, H_MAX_UM, RTOL, DT_MAX_NS) passes through to every run.
cd ~/MX17_comsol
n=0
for spec in "$@"; do
  set -- $spec
  rm -rf /tmp/dneff_axi_$4
  RHO_S=$1 CASE=$2 VE_UM_NS=$3 OUT=$HOME/MX17_comsol/axi_$4.csv COMSOL_WD=/tmp/dneff_axi_$4 \
    ./run_lxplus.sh InduceAxi.java 4 >/dev/null 2>&1 &
  n=$((n + 1)); [ $((n % 2)) -eq 0 ] && wait
done
wait
echo done > axi_batch.done
