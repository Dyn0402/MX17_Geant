#!/bin/bash
# Mixed-response grid (2026-10-20): X = (1-w)*X(D=0) + w*X(D) per event, from the condor 4505159 arms (th00 + thy*,
# seeds 5/6, FEU 07 decoded files copied to $D). Output $O/m<D>_w<w>_<ang>_s<seed>_x.npz (Y: the xo0 *_y.npz).
D=/media/dylan/data/claude_tmp/tmp/nb/xmix; O=$D/npz; mkdir -p $O
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; S=$(cd $(dirname $0)/../../2026-10-19/asym_deck && pwd)/raw_asym.py
for off in ${OFFS:-100 200 300}; do for w in ${WS:-0.2 0.35 0.5 0.65 0.8}; do
  for ang in th00 thy05 thy10 thy15 thy20 thy30; do for s in 5 6; do
    f=$O/m${off}_w${w}_${ang}_s${s}_x.npz
    [ -s $f ] || $P $S $D/xo0_${ang}_s$s 7 $f --files 1 --mix-with $D/xom${off}_${ang}_s$s --mix-w $w > ${f%.npz}.log 2>&1 &
    while [ $(jobs -r | wc -l) -ge ${NJ:-6} ]; do sleep 2; done
  done; done
done; done
wait; echo MIX DONE
