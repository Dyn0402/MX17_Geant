#!/usr/bin/env python3
"""scenarios.py [--quick] -- bow of the MX17 readout board over the support-plate aperture, scenario scan.

Every scenario is swept in gas overpressure 0-10 mbar (continuation) and reports, at each pressure:
  fit_depth_mm   the bench observable: quadratic surface fitted over the efficient area (399 x 360 mm), depth at
                 150 mm = 0.5 kbar 1.5^2, as nTof_x17 fb_bow_maps/summary.py defines it
  dip150_mm      the true on-axis dip w(150 mm) - w(0)
  centre_mm      centre deflection (away from the cathode)
  slope profile  mean d(height)/d(own axis) in the 9 bins (30-370 mm) of the FB w0 profiles, for the shape test
Bench orientation (cathode up) adds the board's weight to the gas load; at n_TOF the boards stand vertical (no
weight term).  Writes results/scenarios.json and results/maps.npz (deflection maps of selected scenarios).

Materials and the stack (MX17_Geant README stack table, nTof_x17 docs/detector_spec):
  PCB laminate  t (structural: FR4 + 5 Cu layers + 50 um coverlay; the CAD's 1.70 mm also counts the mesh and the
                150 um amplification gap) -- scanned 1.2-2.4 mm, nominal 1.6; FR4 E = 20 GPa (17-24), Cu 117 GPa
                smeared with the in-window coverage of each layer (0, 0.76, 0.458, 0.457, 0.53)
  Rohacell 51   5 mm, E 70 MPa, G 19 MPa, 52 kg/m^3 (grade not known: 31 IG = 36 MPa as a variant)
  back mylar    25 um, E 4.5 GPa, on the Rohacell's far face
  support       8 mm Al plate, 402 mm square aperture; 30 mm Al gas frame (410 aperture, 15 mm wall) on the drift side
"""
import json
import os
import pickle
import sys
import time
from multiprocessing import Pool

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate import Plate, layer_props, fit_depth

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')
G = 9.81
E_CU, E_KAP, E_RESIN = 117e9, 2.5e9, 3.5e9
CU_FILL = (0.0, 0.76, 0.458, 0.457, 0.53)       # L3 guard ring (none inside the window), L4 pads, L5, L6, L7
T_CU, T_COV = 26e-6, 50e-6
ROH = {'51': dict(E=70e6, G=19e6, rho=52.0), '31': dict(E=36e6, G=13e6, rho=32.0)}
E_MY, T_MY, RHO_MY = 4.5e9, 25e-6, 1390.0
K_RIGID = 1e12                                  # N/m^3: a rigid bond (validate.py: within 1 % of the composite)
P_MBAR = [0, 1, 2, 3, 4, 5, 7, 10]
BINS = np.linspace(30, 370, 10)                 # FB profile bins (mm, strip frame; centre at 200)


def pcb_layers(t, E_fr4=20e9):
    g = (t - T_COV - 5 * T_CU) / 5
    L = [(E_KAP, T_COV)]
    for f in CU_FILL:
        L += [(E_CU * f + E_RESIN * (1 - f), T_CU), (E_fr4, g)]
    mass = 1420 * T_COV + sum(8960 * f * T_CU for f in CU_FILL) + 1850 * 5 * g + 7900 * 38e-6 * 0.223  # + mesh
    return L, mass


def backing(kind, grade='51'):
    r = ROH[grade]
    if kind == 'none':
        return None, 0.0
    L = [(r['E'], 5e-3)] + ([(E_MY, T_MY)] if kind != 'bonded_roh' else [])
    mass = r['rho'] * 5e-3 + RHO_MY * T_MY
    return L, mass


def kmap(spec, a):
    """interface stiffness map. spec: ('uniform', k) | ('band', k, width_m) | ('debond_centre', k, R_m) |
    ('patch_centre', k, R_m) | ('strips', k, width_m, pitch_m)."""
    kind = spec[0]
    if kind == 'uniform':
        return spec[1]
    k = spec[1]
    if kind == 'band':
        b = spec[2]; return lambda x, y: np.where(np.maximum(np.abs(x), np.abs(y)) > a / 2 - b, k, 0.0)
    if kind == 'debond_centre':
        R = spec[2]; return lambda x, y: np.where(np.hypot(x, y) > R, k, 0.0)
    if kind == 'patch_centre':
        R = spec[2]; return lambda x, y: np.where(np.hypot(x, y) < R, k, 0.0)
    if kind == 'strips':
        wd, pitch = spec[2], spec[3]
        return lambda x, y: np.where(np.abs(((x + pitch / 2) % pitch) - pitch / 2) < wd / 2, k, 0.0)
    raise ValueError(kind)


def slope_profiles(p, w, other=180.0):
    """mean d(height)/d(own axis) in the FB bins, for the x and y axes [mm/mm]."""
    dwdx = np.gradient(w, p.xs, axis=0); dwdy = np.gradient(w, p.xs, axis=1)
    oth = np.linspace(-other, other, 37) * 1e-3
    prof = {}
    for ax, F in (('x', dwdx), ('y', dwdy)):
        vals = []
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            pos = (np.linspace(lo, hi, 7) - 200.0) * 1e-3
            P, O = np.meshgrid(pos, oth, indexing='ij')
            vals.append(float(np.mean(p.interp(F, P, O) if ax == 'x' else p.interp(F, O, P))))
        prof[ax] = vals
    return prof


def run(sc):
    cache = os.path.join(RES, 'cache', sc['name'] + '.pkl')
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    t0 = time.time()
    a = sc.get('span', 0.402)
    L1l, m1 = pcb_layers(sc['t'], sc.get('E_fr4', 20e9))
    L2l, m2 = backing(sc['backing'], sc.get('grade', '51'))
    L1 = layer_props(L1l)
    L2 = layer_props(L2l) if L2l else None
    k = 0.0
    if L2 is not None and sc['backing'] != 'loose':
        k = kmap(sc.get('bond', ('uniform', K_RIGID)), a)
    edge = sc.get('edge', 'spring')
    K = sc.get('K_rot')
    if edge == 'spring' and K is None:          # foam pinched between the 15 mm frame wall and the plate
        b, Er = 15e-3, ROH[sc.get('grade', '51')]['E']
        K = Er * b ** 3 / (12 * 5e-3) if L2l else 1e7     # no foam -> board directly on the plate: clamped
    p = Plate(a, sc.get('N', 41), L1, L2, k=k, edge=edge, K_rot=K, inplane=sc.get('inplane', 'immovable'),
              theta0=sc.get('theta0', 0.0))
    qg = (m1 + m2) * G if sc.get('gravity', True) else 0.0
    pm = sc.get('p_mbar', P_MBAR)
    ws = p.sweep([-(100.0 * x + qg) for x in pm], m=sc.get('m', 0.0), substeps=1)
    rows = []
    for x, w in zip(pm, ws):
        r = fit_depth(p, w)
        r.update(p_mbar=x, slope=slope_profiles(p, w))
        rows.append(r)
    out = dict(sc, rows=rows, D_pcb=L1['D'], D_backing=(L2['D'] if L2 else 0.0), A_pcb=L1['A'],
               D_composite=layer_props(L1l + (L2l or []))['D'], q_gravity_Pa=qg, K_rot=K, d_m=p.d,
               secs=time.time() - t0)
    if 'bond' in out:
        out['bond'] = list(out['bond'])
    if sc.get('keep_map'):
        out['_maps'] = {f'{x:g}': w for x, w in zip(pm, ws)}
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    pickle.dump(out, open(cache, 'wb'))
    print(f"{sc['name']:38s} {time.time() - t0:6.1f} s   depth@3mbar "
          f"{rows[pm.index(3)]['fit_depth_mm'] if 3 in pm else float('nan'):.2f} mm", flush=True)
    return out


def scenario_list():
    S = []
    # A. main grid: thickness x backing, foam-spring edge, bench
    for t in (1.3, 1.6, 1.7, 2.0):
        for bk in ('none', 'loose', 'bonded_roh', 'bonded_sandwich'):
            S.append(dict(name=f'A_t{t}_{bk}', group='A', t=t * 1e-3, backing=bk, keep_map=(t == 1.6)))
    # A'. finer thickness scan for the two extremes of the backing
    for t in (1.2, 1.4, 1.5, 1.8, 2.4):
        for bk in ('loose', 'bonded_sandwich'):
            S.append(dict(name=f'At_t{t}_{bk}', group='At', t=t * 1e-3, backing=bk))
    # B. partial bonding (t = 1.6, sandwich backing)
    for nm, bond in (('psa_1e8', ('uniform', 1e8)), ('psa_1e9', ('uniform', 1e9)), ('psa_1e10', ('uniform', 1e10)),
                     ('band25', ('band', K_RIGID, 0.025)), ('band60', ('band', K_RIGID, 0.060)),
                     ('debond_R100', ('debond_centre', K_RIGID, 0.100)),
                     ('debond_R150', ('debond_centre', K_RIGID, 0.150)),
                     ('patch_R100', ('patch_centre', K_RIGID, 0.100)),
                     ('strips50_p100', ('strips', K_RIGID, 0.050, 0.100)),
                     ('strips25_p130', ('strips', K_RIGID, 0.025, 0.134))):
        S.append(dict(name=f'B_{nm}', group='B', t=1.6e-3, backing='bonded_sandwich', bond=bond, keep_map=True))
    # C. edges and in-plane restraint (t = 1.6; loose and bonded)
    for bk in ('loose', 'bonded_sandwich'):
        for nm, kw in (('clamped', dict(edge='clamped')), ('K400', dict(edge='spring', K_rot=400.0)),
                       ('K100', dict(edge='spring', K_rot=100.0)), ('ss', dict(edge='ss')),
                       ('free_inplane', dict(inplane='free')), ('span410_ss', dict(edge='ss', span=0.410))):
            S.append(dict(name=f'C_{bk}_{nm}', group='C', t=1.6e-3, backing=bk, keep_map=(nm in ('clamped', 'ss')),
                          **kw))
    # D. materials
    for bk in ('loose', 'bonded_sandwich'):
        S.append(dict(name=f'D_{bk}_roh31', group='D', t=1.6e-3, backing=bk, grade='31'))
        for E in (17e9, 24e9):
            S.append(dict(name=f'D_{bk}_E{E / 1e9:.0f}', group='D', t=1.6e-3, backing=bk, E_fr4=E))
    # E. n_TOF: vertical, no weight term
    for bk in ('loose', 'bonded_sandwich'):
        S.append(dict(name=f'E_{bk}_vertical', group='E', t=1.6e-3, backing=bk, gravity=False))
    # F. non-pressure shapes for the shape test: mesh-tension moment on an SS / soft-spring edge, imposed edge rotation
    for nm, kw in (('mesh10Ncm_ss', dict(edge='ss', m=1000.0 * 0.8e-3)),
                   ('mesh10Ncm_K400', dict(edge='spring', K_rot=400.0, m=1000.0 * 0.8e-3)),
                   ('edge_rot_m5mrad', dict(edge='clamped', theta0=-5e-3)),
                   ('edge_rot_m10mrad', dict(edge='clamped', theta0=-10e-3))):
        S.append(dict(name=f'F_{nm}', group='F', t=1.6e-3, backing='loose', keep_map=True, **kw))
    return S


def main():
    S = scenario_list()
    if '--quick' in sys.argv:
        S = [s for s in S if s['name'] in ('A_t1.6_loose', 'A_t1.6_bonded_sandwich')]
    os.makedirs(RES, exist_ok=True)
    with Pool(int(os.environ.get('NPROC', 6))) as pool:
        res = pool.map(run, S, chunksize=1)
    maps = {}
    for r in res:
        for kx, w in r.pop('_maps', {}).items():
            maps[f"{r['name']}|{kx}"] = w
    json.dump(res, open(f'{RES}/scenarios.json', 'w'), indent=1)
    np.savez_compressed(f'{RES}/maps.npz', **maps)
    print(f'{len(res)} scenarios -> {RES}/scenarios.json')


if __name__ == '__main__':
    main()
