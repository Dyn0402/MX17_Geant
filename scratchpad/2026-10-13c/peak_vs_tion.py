#!/usr/bin/env python3
"""How much does the ion transit change the PEAK amplitude per unit charge? (FRESH_LOOK §13)
toy_pulse with the T14 geometry (det3, dry 95/5 transport, 166/180 ns shaper), same total
charge, f_ion 0.908: ratio of median peak-strip maximum for T_ion = 140 / 200 ns against the
chain's T = 340 ns (Ar+), per track angle."""
import sys, numpy as np
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
import toy_pulse as T
_ORIG = T.shaper
for tp in (166, 180):
    T.shaper = lambda fam, t_peak=180.0, tp=tp: _ORIG(fam, tp)
    for th in (0.0, 10.7, 20.0):
        pk = {}
        for Ti in (340., 200., 140.):
            rng = np.random.default_rng(7); h = T.shaper('manual_1')
            fine = T.simulate({**T.BASE, 'shaper': 'manual_1', 'f_ion': 0.908, 't_ion': Ti, 'theta': th}, 600, rng, h)
            pk[Ti] = np.median([w.max() for w in fine])
        print(f'shaper {tp} ns, theta {th:4.1f}: peak(T140)/peak(T340) = {pk[140]/pk[340]:.3f},  peak(T200)/peak(T340) = {pk[200]/pk[340]:.3f}')
