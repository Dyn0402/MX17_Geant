"""Connector-edge test (2026-10-20) on raw_asym.py --edge outputs (yoga edge/*.npz, copied to ~/x17/nb_asym_raw/edge/).
The copy goes into the physically lower column (DAQ n+1 under connector inversion). For a peak at the LOWEST column of a
connector (col % 64 == 0) that column is on the other connector / DREAM chip / cable; the DAQ n+1 is then the far end
of the next chip (127 columns away). Detector-side coupling: excess on physical col-1 survives, DAQ n+1 gets none.
Readout-side coupling: the reverse. Peaks at col % 64 == 63 are the mirror control (col-1 on the same connector).
Uses asym_lib.good + hot veto logic on the column-order window. Usage: edge_test.py [dir]"""
import os, sys, glob, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import asym_lib as L
src = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/x17/nb_asym_raw/edge')
COL = np.array([64 * (c // 64) + 63 - c % 64 for c in range(512)])
def alphas(W):
    c = W[:, 4]; n2 = np.einsum('ij,ij->i', c, c)
    return {k: np.einsum('ij,ij->i', W[:, 4 + k], c) / n2 for k in (-2, -1, 1, 2)}
OUT = {}
for f in sorted(glob.glob(f'{src}/e_*.npz')):
    z = np.load(f); W, Wd, p = z['wf'], z['wf_daq'], z['p'].astype(int)
    g = L.good(W); hv = L.hot(p[g]); g &= ~np.isin(p, hv); W, Wd, p = W[g], Wd[g], p[g]
    pc = COL[p]; lc = pc % 64
    aP, aD = alphas(W), alphas(Wd)      # physical: index -1 = col-1 ; DAQ: +1 = DAQ n+1
    tag = os.path.basename(f)[:-4]; o = {}
    for nm, m in (('col%64==0 (col-1 on OTHER connector)', lc == 0), ('col%64==1..3', (lc >= 1) & (lc <= 3)),
                  ('col%64==63 (col+1 on other connector)', lc == 63), ('col%64==60..62', (lc >= 60) & (lc <= 62))):
        if m.sum() < 15: continue
        phys = aP[-1][m] - aP[1][m]; daq = aD[1][m] - aD[-1][m]
        o[nm] = dict(n=int(m.sum()), phys=[float(phys.mean()), float(phys.std() / np.sqrt(m.sum()))],
                     daq=[float(daq.mean()), float(daq.std() / np.sqrt(m.sum()))],
                     a_colm1=float(aP[-1][m].mean()), a_colp1=float(aP[1][m].mean()), a_daqp1=float(aD[1][m].mean()), a_daqm1=float(aD[-1][m].mean()))
    OUT[tag] = o
    print(f'== {tag}  ({len(p)} edge events)')
    for nm, r in o.items():
        print(f'   {nm:38s} n={r["n"]:4d} | physical α(col−1)−α(col+1) {r["phys"][0]:+.3f}±{r["phys"][1]:.3f}  [{r["a_colm1"]:.2f}/{r["a_colp1"]:.2f}]'
              f' | DAQ α(n+1)−α(n−1) {r["daq"][0]:+.3f}±{r["daq"][1]:.3f}  [{r["a_daqp1"]:.2f}/{r["a_daqm1"]:.2f}]')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'edge_test.json'), 'w'))
