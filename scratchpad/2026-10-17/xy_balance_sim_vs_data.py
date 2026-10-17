#!/usr/bin/env python3
"""X/Y charge balance, sim vs det3 data, one observable for both (2026-10-17).
Both sides read decoded waveforms with wft.io.FeuReader (pedestal + CNS as in T14). Per event and
view: strips with amp >= 5 sigma and >= 100 ADC; q = sum of strip peak amplitudes, pk = max.
"Unsat" = no strip >= 3500 ADC in either view (the ledger's rail proxy). Data: every trigger in
the files given (all cosmic angles, no M3); sim: each fixed-angle point.
Usage: xy_balance_sim_vs_data.py LABEL:DIR [LABEL:DIR ...]   (DIR holds *_07.root / *_08.root)"""
import os, sys, glob
import numpy as np
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17')]
from wft.io import FeuReader

def per_view(fn, max_files=None):
    q, pk = {}, {}
    r = FeuReader(fn); noise = np.where(r.noise > 0, r.noise, np.inf)
    for e, _t, w in r.iter_events(set(int(x) for x in r.event_ids)):
        amp = w.max(axis=1); over = (amp / noise >= 5) & (amp >= 100)
        if over.any(): q[int(e)] = amp[over].sum(); pk[int(e)] = amp[over].max()
    return q, pk

def main():
    print(f'{"point":22s} {"n":>6s} {"f med":>6s} {"ln qX/qY":>8s} {"pkX/pkY":>7s} | {"n unsat":>7s} {"f med":>6s} {"ln qX/qY":>8s} {"pkX/pkY":>7s}')
    for arg in sys.argv[1:]:
        lab, d = arg.split(':', 1)
        QX, PX, QY, PY = {}, {}, {}, {}
        for fx in sorted(glob.glob(f'{d}/*_07.root')):
            fy = fx[:-8] + '_08.root'
            if not os.path.exists(fy): continue
            qx, px = per_view(fx)
            qy, py = per_view(fy)
            for k in set(qx) & set(qy):
                QX[(fx, k)] = qx[k]; PX[(fx, k)] = px[k]; QY[(fx, k)] = qy[k]; PY[(fx, k)] = py[k]
        keys = list(QX)
        qx = np.array([QX[k] for k in keys]); qy = np.array([QY[k] for k in keys])
        px = np.array([PX[k] for k in keys]); py = np.array([PY[k] for k in keys])
        uns = (px < 3500) & (py < 3500)
        def s(m):
            return (f'{m.sum():6d} {np.median(qx[m] / (qx[m] + qy[m])):6.3f} {np.median(np.log(qx[m] / qy[m])):+8.3f} '
                    f'{np.median(px[m]) / np.median(py[m]):7.3f}')
        print(f'{lab:22s} {s(np.ones(len(qx), bool))} | {s(uns)}', flush=True)


if __name__ == '__main__':
    main()
