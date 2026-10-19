#!/bin/bash
# Desktop driver (2026-10-19): bus_wpot.py scan on the 10 um grid (every feature on cell boundaries). Two at a time.
cd $(dirname $0); P=~/PycharmProjects/nTof_x17/.venv/bin/python; mkdir -p bus_scan
run(){ n=$1; shift; [ -s bus_scan/$n.json ] || $P bus_wpot.py --nrow 4 --dxy 10 --out bus_scan/$n.json "$@" > bus_scan/$n.log 2>&1; }
run d10_nominal & run d10_nobus --no-bus; wait
run d10_L6 --layer 6 & run d10_h132 --h1 132; wait
run d10_h66 --h1 66 & run d10_h33 --h1 33; wait
run d10_gap150 --gap 150 & run d10_gap200 --gap 200; wait
run d10_h15 --h1 15 & run d10_gap250 --gap 250; wait
echo BUS SCAN DONE
