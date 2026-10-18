#!/usr/bin/env python3
"""The hit-time ladder with NO readout: Geant4 electrons + drift + diffusion only (2026-10-18, item 3).

Same fit as ladder.py (amplitude-weighted least squares of strip time on strip position, >= 4 strips), but the
strip "amplitude" is the number of electrons that land on it and the strip "time" is the median arrival time of
those electrons. No avalanche, no induction, no resistive sharing, no shaper, no sampling. Whatever compression
is left here comes from the geometry of the measurement itself: each strip integrates a slice of depth, diffusion
smears electrons across strips with no time correlation (errors-in-variables dilution), and the threshold trims
the ends. The difference between this and ladder.py's C on the full sim is the readout's share.

Drift: the digitizer's own DriftGas table at E = 1000 V / 26.5 mm (v ~ 38.0 um/ns), z <= 26.5 mm, sigma_T and
sigma_t from the table. Strips: 0.78 mm pitch on the tilted view's coordinate. Threshold: a strip is kept if it has
>= FRAC of the event's brightest strip (default 0.07, about what 100 ADC is against a typical unsaturated peak).

Usage: ideal_ladder.py OUT.json STAGEA_DIR [FRAC]
"""
import sys, os, glob, json
import numpy as np
sys.path.insert(0, os.path.expanduser('~/CLionProjects/MX17_Geant'))
from response.digitizer.clusters import ClusterFile
from response.digitizer.gas import DriftGas
from response.digitizer import digitize as D

GAP, PITCH = float(os.environ.get('GAP', 26.5)), 0.78
ANGS = {'thx05': ('x', 5), 'thx10': ('x', 10), 'thx15': ('x', 15), 'thx20': ('x', 20), 'thx30': ('x', 30),
        'thy05': ('y', 5), 'thy10': ('y', 10), 'thy15': ('y', 15), 'thy20': ('y', 20), 'thy30': ('y', 30)}


def main():
    out, sdir = sys.argv[1], sys.argv[2]
    frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.07
    gas = DriftGas()
    E = 1000.0 / (GAP * 0.1)
    v = float(np.ravel(gas.v_drift_um_ns(E))[0])
    print(f'E {E:.1f} V/cm  v {v:.2f} um/ns')
    rng = np.random.default_rng(7)
    res = {'v_um_ns': v, 'frac': frac}
    for ang, (view, a) in ANGS.items():
        fn = f'{sdir}/mx17_muons_3k_ang_{ang}_e17_t0.root'
        if not os.path.exists(fn):
            continue
        cf = ClusterFile(fn, volumes=('DriftGas',))
        cf.keep((cf.z > 0) & (cf.z <= GAP))
        tan = np.tan(np.radians(a))
        C = {'med': [], 'mean': [], 'med_xt': [], 'mean_xt': []}
        for _e, g in cf.events():
            n = cf.n_e[g]
            z = np.repeat(cf.z[g], n)
            u = np.repeat((cf.x if view == 'x' else cf.y)[g], n)
            sT = gas.sigma_T_um(E, z) * 1e-3
            st = gas.sigma_t_ns(E, z)
            u = u + rng.normal(0, 1, len(u)) * sT
            t = z * 1e3 / v + rng.normal(0, 1, len(z)) * st
            s = np.floor(u / PITCH).astype(int)
            us = np.unique(s)
            q = np.array([np.sum(s == k) for k in us], float)
            keep = q >= frac * q.max()
            if keep.sum() < 4:
                continue
            us, q = us[keep], q[keep]
            for nm, fn_t in (('med', np.median), ('mean', np.mean)):
                tt = np.array([fn_t(t[s == k]) for k in us])
                x = (us + 0.5) * PITCH
                xm = (q * x).sum() / q.sum(); ym = (q * tt).sum() / q.sum()
                sl = (q * (x - xm) * (tt - ym)).sum() / (q * (x - xm) ** 2).sum()
                C[nm].append(abs(sl) * v * 1e-3 * tan)
                b = (q * (x - xm) * (tt - ym)).sum() / (q * (tt - ym) ** 2).sum()
                C[nm + '_xt'].append(abs(1.0 / b) * v * 1e-3 * tan)
        res[ang] = {nm: {'C': float(np.median(c)), 'q25': float(np.percentile(c, 25)),
                         'q75': float(np.percentile(c, 75)), 'n': len(c)} for nm, c in C.items()}
        print(ang, json.dumps(res[ang]), flush=True)
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
