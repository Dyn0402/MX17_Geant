"""Late across-strip growth, X vs Y view, per session (2026-10-22). G = r(f50 falling) / r(peak) for each session's own
long pulses (central FWHM >= 0.8 x the session/view p90, |dt_nb| < 80 ns, 0.1 < s1 < 1), bootstrap errors; the ratio
G_X / G_Y within a session cancels most of the angle-mix and pulse-length differences between setups. Listed with the
session's X asymmetry amplitude a (survey) and drift voltage. Sim glg0 (vertical, matched FWHM) for reference.
Usage: late_xy.py [out.json]"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt')); sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import prompt_share as PS, r5lib
from fwhm_bins import load, HW
RNG = np.random.default_rng(11)

def growth(files, fw_cut=None):
    W, s1, dtn, fw, pk = load(files)
    p90 = np.nanpercentile(fw, 90); m = (fw >= (fw_cut if fw_cut else .8 * p90)) & (np.abs(dtn) < 80) & (s1 > .1) & (s1 < 1.)
    c = W[:, HW]; S = W[:, HW - 1] + W[:, HW + 1]; tl = PS.cross(c, pk, .5, False); tp = pk.astype(float)
    nl, dl, npk, dp = PS.at(S, tl), PS.at(c, tl), PS.at(S, tp), PS.at(c, tp)
    g = m & np.isfinite(nl) & np.isfinite(dl); ii = np.flatnonzero(g)
    if len(ii) < 25: return None
    G = lambda j: (nl[j].sum() / dl[j].sum()) / (npk[j].sum() / dp[j].sum())
    bs = [G(ii[RNG.integers(0, len(ii), len(ii))]) for _ in range(300)]
    return dict(G=float(G(ii)), err=float(np.std(bs)), n=int(len(ii)), s1=float(s1[ii].mean()), p90=float(p90))

if __name__ == '__main__':
    out = {}; rows = {}
    for r in r5lib.rows():
        if r['build'] not in ('det3 June-15', 'det3 May', 'det2'): continue
        k = (r['build'], r['run'][-14:], r['sub'][:24], r['mesh'], r['drift'])
        try: rows.setdefault(k, {})[r['view']] = (growth([r['f']]), r.get('a'))
        except Exception: pass
    print('build        | session                                  | mesh drift |  G_X            G_Y           G_X/G_Y        a_X')
    for k, v in rows.items():
        if 'X' not in v or 'Y' not in v or v['X'][0] is None or v['Y'][0] is None: continue
        gx, gy = v['X'][0], v['Y'][0]; q = gx['G'] / gy['G']; qe = q * np.hypot(gx['err'] / gx['G'], gy['err'] / gy['G'])
        out[' | '.join(map(str, k))] = dict(X=gx, Y=gy, ratio=[q, qe], a=v['X'][1])
        print(f'{k[0]:12s} | {k[1]} {k[2]:24s} | {k[3]} {k[4]:>5} | {gx["G"]:.2f}±{gx["err"]:.2f} ({gx["n"]:3d})  {gy["G"]:.2f}±{gy["err"]:.2f} ({gy["n"]:3d})  '
              f'{q:.2f}±{qe:.2f}   {v["X"][1]:+.3f}')
    RAW = os.path.expanduser('~/x17/nb_asym_raw')
    for fwc in (450, 550):
        gx = growth(sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_x.npz')), fwc); gy = growth(sorted(glob.glob(f'{RAW}/simgapleak/glg0_th00_*_y.npz')), fwc)
        q = gx['G'] / gy['G']; out[f'sim glg0 FWHM>={fwc}'] = dict(X=gx, Y=gy, ratio=[q, q * np.hypot(gx['err'] / gx['G'], gy['err'] / gy['G'])])
        print(f'sim glg0 FWHM >= {fwc} (vertical)                               | {gx["G"]:.2f}±{gx["err"]:.2f} ({gx["n"]:4d}) {gy["G"]:.2f}±{gy["err"]:.2f} ({gy["n"]:3d})  {q:.2f}')
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'late_xy.json'), 'w'), indent=1)
