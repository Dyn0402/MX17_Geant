"""r(f) in absolute central-FWHM bins (2026-10-22). The 2026-10-21 near-vertical cut keeps FWHM >= the set's own 25th
percentile, but the data's FWHM is short (det3 May X |dt_nb| < 30 ns: median 423 ns; sim th00: 590 ns): most data events
in the cut are inclined, the central pad seeing only part of the drift column. Here both are binned in absolute FWHM
(|dt_nb| < DTN_MAX), s1 0.2-0.45 and 0.45-0.7. Usage: [DTN_MAX=80 ARMS=cdt0s4,...] fwhm_bins.py"""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import prompt_share as PS, asym_lib as L
RAW = os.path.expanduser(os.environ.get('NB_RAW', '~/x17/nb_asym_raw')); HW = 4; DT = 60.
DTN_MAX = float(os.environ.get('DTN_MAX', 80.)); FWB = (0, 350, 450, 550, 650, 2000)

def load(files):
    W, s1 = [], []
    for f in files:
        try: d = L.load(f)
        except Exception: continue
        p = d['p'] % 512; m = (p % 64 >= 4) & (p % 64 <= 59); W.append(d['W'][m]); s1.append(d['al'][1][m] + d['al'][-1][m])
    W = np.concatenate(W).astype(float); s1 = np.concatenate(s1); c = W[:, HW]; pk = c.argmax(1)
    ts = np.arange(W.shape[2]); pos = lambda w: np.clip(w, 0, None); tc = lambda w: (pos(w) * ts).sum(1) / np.maximum(pos(w).sum(1), 1e-9)
    dtn = (tc(W[:, HW + 1]) - tc(W[:, HW - 1])) * DT; fw = (PS.cross(c, pk, .5, False) - PS.cross(c, pk, .5)) * DT
    return W, s1, dtn, fw, pk

def table(name, files, out=None):
    W, s1, dtn, fw, pk = load(files); c = W[:, HW]; S = W[:, HW - 1] + W[:, HW + 1]
    pts = dict(r05=PS.cross(c, pk, .05), r25=PS.cross(c, pk, .25), peak=pk.astype(float), f50=PS.cross(c, pk, .5, False))
    print(f'\n{name}  (n {len(s1)}, |dt_nb| < {DTN_MAX:.0f})')
    print('   FWHM ns     s1 bin       n    r05          r25          peak         f50          early/pk late/pk')
    rows = []
    for lo, hi in zip(FWB[:-1], FWB[1:]):
        for slo, shi in ((.2, .45), (.45, .7)):
            m = (np.abs(dtn) < DTN_MAX) & (fw >= lo) & (fw < hi) & (s1 >= slo) & (s1 < shi)
            if m.sum() < 25: continue
            v = {}
            for k, t in pts.items():
                nu, de = PS.at(S, t), PS.at(c, t); g = m & np.isfinite(nu) & np.isfinite(de); ii = np.flatnonzero(g)
                bs = PS.RNG.integers(0, len(ii), (200, len(ii)))
                v[k] = (nu[g].sum() / de[g].sum(), np.std(nu[ii][bs].sum(1) / de[ii][bs].sum(1)))
            rows.append(dict(fw=(lo, hi), s1=(slo, shi), n=int(m.sum()), **{k: list(x) for k, x in v.items()}))
            print(f'   {lo:4d}-{hi:<5d}  {slo:.2f}-{shi:.2f}  {m.sum():5d}  ' + '  '.join(f'{v[k][0]:.3f}±{v[k][1]:.3f}' for k in pts)
                  + f'   {v["r05"][0] / v["peak"][0]:.2f}  {v["f50"][0] / v["peak"][0]:.2f}')
    return rows

if __name__ == '__main__':
    import json; R = {}
    for view in ('x', 'y'):
        for n, t in (('det3 May', 'd3may'), ('det3 June-15', 'd3jun'), ('det2', 'd2')):
            R[f'{n} {view.upper()}'] = table(f'{n} {view.upper()}', [f'{RAW}/{t}_{view}.npz'])
        pat = 'glg0_th*_x.npz' if view == 'x' else 'glg0_th00_*_y.npz'
        R[f'sim glg0 {view.upper()}'] = table(f'sim glg0 {view.upper()}', sorted(glob.glob(f'{RAW}/simgapleak/{pat}')))
        for tag in filter(None, os.environ.get('ARMS', '').split(',')):          # cloud arms (simcloud/c<tag>_*)
            R[f'sim {tag} {view.upper()}'] = table(f'sim {tag} {view.upper()}', sorted(glob.glob(f'{RAW}/simcloud/{pat.replace("glg0", tag)}')))
    json.dump(R, open(os.path.join(HERE, f'fwhm_bins_dtn{DTN_MAX:.0f}.json'), 'w'), indent=1)
