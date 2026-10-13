#!/usr/bin/env python3
"""Rise (10-90 %, peak strip, T14 metric) of a point-like or short-column deposit through the
DREAM shaper, vs ion transit T and ion fraction f: what det0's short-column pulses should
show (FRESH_LOOK §13). Uses toy_pulse's shaper and electron_kernel (shared with §7.5/§10.2).
Column: uniform charge over L mm drifting at v (arrival spread L/v), diffusion ignored.
60 ns sampling with random phase, white noise 11 ADC, peak ~1000 ADC (det0 rise sample:
300-3500 ADC)."""
import sys, numpy as np
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
import toy_pulse as T
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-13c')
from rise_all import cross
rng = np.random.default_rng(3)
def rise(tp, f, Tion, spread_ns, n=600, amp=1000.):
    h = T.shaper('manual_1', tp)
    k, _ = T.electron_kernel(h, f, Tion, 0.0, 230.)
    box = np.ones(max(int(spread_ns), 1)) / max(int(spread_ns), 1)
    wf = np.convolve(k, box)[:2000]; wf /= wf.max()
    out = []
    for _ in range(n):
        ph = rng.uniform(0, 60); t = ph + 60 * np.arange(32) - 300
        w = amp * np.interp(t, np.arange(len(wf)), wf, left=0, right=0) + rng.normal(0, 11, 32)
        l10, l90 = cross(w, .1, -1), cross(w, .9, -1)
        out.append((l90 - l10) * 60)
    return np.nanmedian(out)
if __name__ == '__main__':
    for tp in (166, 180):
        for spread in (0, 100):
            print(f'shaper {tp} ns, column arrival spread {spread} ns:')
            for f in (0.0, 0.5, 0.9):
                print('   f=%.1f ' % f + '  '.join(f'T{Ti:.0f}:{rise(tp, f, Ti, spread):.0f}' for Ti in (20., 40., 70., 100., 140., 200., 340.)))
