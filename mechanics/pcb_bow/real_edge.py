#!/usr/bin/env python3
"""real_edge.py -- the 2-D scan re-run with the mounting's own edge stiffness (edge1d.py) and with mesh tension.

Edge: rotational spring K_eff from the 1-D cross-section of the STEP geometry (frame bottomed on the board, Rohacell
foundation, O-ring seal at 208.4 mm): 700 / 1110 / 1830 N m/m for a 1.3 / 1.6 / 2.0 mm PCB.  The softer bound is the
frame riding on the O-ring at ~10 mbar (K ~ 500 N m/m at 1.6 mm).

Mesh tension: the mesh is laminated on the board (bulk pillars), so its tension T is balanced inside the board by an
equal compression; the pair carries no net in-plane force (no stress stiffening or softening of the bend), only the
moment m = T e of the mesh about the PCB's neutral surface, e = (depth of the neutral surface below the board top)
+ 0.15 mm (amplification gap).  T = 10 N/m (the stated 5-10 N/m) and 1000 N/m (10 N/cm, typical stretched-mesh
tensions are quoted in N/cm), the moment acting over the whole span, with pressure.  Writes results/real_edge.json.
"""
import json
import os
import sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as S
from plate import layer_props

KEFF = {1.3: 701.0, 1.6: 1113.0, 2.0: 1825.0}
P = [0, 1, 2, 3, 4, 5, 7, 10]


def e_mesh(t):
    L, _ = S.pcb_layers(t)
    return layer_props(L)['zn'] + 0.15e-3


def scen():
    out = []
    for t, K in KEFF.items():
        for bk in ('loose', 'bonded_sandwich'):
            out.append(dict(name=f'G_t{t}_{bk}_Keff', group='G', t=t * 1e-3, backing=bk, edge='spring', K_rot=K,
                            p_mbar=P, keep_map=(t == 1.6)))
    for bk in ('loose', 'bonded_sandwich'):
        out.append(dict(name=f'G_t1.6_{bk}_Koring', group='G', t=1.6e-3, backing=bk, edge='spring', K_rot=500.0,
                        p_mbar=P))
        for T in (10.0, 1000.0):
            out.append(dict(name=f'M_t1.6_{bk}_T{T:g}', group='M', t=1.6e-3, backing=bk, edge='spring',
                            K_rot=KEFF[1.6], m=T * e_mesh(1.6e-3), T_mesh=T, p_mbar=P))
    out.append(dict(name='M_t1.6_loose_T1000_ss', group='M', t=1.6e-3, backing='loose', edge='ss',
                    m=1000.0 * e_mesh(1.6e-3), T_mesh=1000.0, p_mbar=P))
    return out


if __name__ == '__main__':
    with Pool(int(os.environ.get('NPROC', 3))) as pool:
        res = pool.map(S.run, scen(), chunksize=1)
    for r in res:
        r.pop('_maps', None)
    json.dump(res, open(os.path.join(S.RES, 'real_edge.json'), 'w'), indent=1)
    for r in res:
        print(f"{r['name']:32s} " + ' '.join(f"{x['p_mbar']:g}:{x['fit_depth_mm']:.2f}" for x in r['rows']))
