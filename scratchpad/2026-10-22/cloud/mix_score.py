"""Does the data's inclination mix make its 'flat early sharing'? (2026-10-22) The sim's near-vertical X set has zero
x-tilt; the data's |dt_nb| < 80 ns, FWHM >= p25 cut keeps short (inclined) pulses. Sim X at x-tilt 0 (glg0 th00+thy),
1/2/3/5 deg (ctx*, sheared th00) and Stage A thx05/thx10 (cg0_thx*): non-negative weights per angle fitted so the sim's
central-FWHM histogram (|dt_nb| < 80) matches the data's, then the 2026-10-21 cut (FWHM >= the weighted p25) and the
r(f) profile with the weights. Per-angle profiles too. Usage: mix_score.py [data tag, default d3may] [out.json]"""
import os, sys, glob, json, numpy as np
from scipy.optimize import nnls
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import prompt_share as PS
from fwhm_bins import load, HW
RAW = os.path.expanduser(os.environ.get('NB_RAW', '~/x17/nb_asym_raw'))
ANG = {0: sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_x.npz')), 1: 'ctx1_th00', 2: 'ctx2_th00', 3: 'ctx3_th00',
       5: 'ctx5_th00', '5g4': 'cg0_thx05', 10: 'cg0_thx10'}
FR = dict(r05=.05, r10=.1, r25=.25); EDGES = np.arange(100, 1001, 50); BINS = ((.2, .45), (.45, .7))

def prof(W, s1, pk, w, m):
    c = W[:, HW]; S = W[:, HW - 1] + W[:, HW + 1]
    pts = {k: PS.cross(c, pk, f) for k, f in FR.items()}; pts['peak'] = pk.astype(float); pts['f50'] = PS.cross(c, pk, .5, False)
    out = {}
    for lo, hi in BINS:
        b = m & (s1 >= lo) & (s1 < hi); row = {'n': int(b.sum()), 'n_eff': float(w[b].sum() ** 2 / max((w[b] ** 2).sum(), 1e-12))}
        for k, t in pts.items():
            nu, de = PS.at(S, t), PS.at(c, t); g = b & np.isfinite(nu) & np.isfinite(de)
            if g.sum() > 20: row[k] = float((w[g] * nu[g]).sum() / (w[g] * de[g]).sum())
        out[f'{lo}-{hi}'] = row
    return out

def line(name, pr):
    for b, r in pr.items():
        print(f'   {name:22s} s1 {b:9s} n {r["n"]:5d} (eff {r["n_eff"]:6.0f})  ' + '  '.join(f'{k} {r[k]:.3f}' for k in ('r05', 'r10', 'r25', 'peak', 'f50') if k in r)
              + (f'   early/pk {r["r05"] / r["peak"]:.2f} late/pk {r["f50"] / r["peak"]:.2f}' if 'r05' in r and 'peak' in r and 'f50' in r else ''))

if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'd3may'
    Wd, sd, dd, fd, pd = load([f'{RAW}/{tag}_x.npz']); md = np.abs(dd) < 80
    hd = np.histogram(fd[md], EDGES)[0].astype(float)
    S = {}
    for a, f in ANG.items():
        files = f if isinstance(f, list) else sorted(glob.glob(f'{RAW}/simcloud/{f}_*_x.npz'))
        if files: S[a] = load(files)
    keys = list(S); H = np.array([np.histogram(S[a][3][np.abs(S[a][2]) < 80], EDGES)[0] for a in keys], float).T
    wts, res = nnls(H, hd); print(f'{tag} X: FWHM-histogram weights per sim event ' + ', '.join(f'{a}: {x:.3g}' for a, x in zip(keys, wts))
                                     + f'  (chi2-ish {res ** 2 / max(hd.sum(), 1):.3f})')
    print('   FWHM bins       data  ' + ' '.join(f'{e:4.0f}' for e in EDGES[:-1]))
    print('   data counts           ' + ' '.join(f'{v:4.0f}' for v in hd)); print('   fitted mixture        ' + ' '.join(f'{v:4.0f}' for v in H @ wts))
    # data: the 2026-10-21 cut and the matched-FWHM cut
    R = {'weights': dict(zip(map(str, keys), map(float, wts)))}
    m21 = md & (fd >= np.nanpercentile(fd, 25)); R['data'] = prof(Wd, sd, pd, np.ones(len(sd)), m21); line(f'data {tag} (10-21 cut)', R['data'])
    W = np.concatenate([S[a][0] for a in keys]); s1 = np.concatenate([S[a][1] for a in keys]); dt = np.concatenate([S[a][2] for a in keys])
    fw = np.concatenate([S[a][3] for a in keys]); pk = np.concatenate([S[a][4] for a in keys])
    w = np.concatenate([np.full(len(S[a][1]), x) for a, x in zip(keys, wts)])
    ok = np.isfinite(fw); o = np.argsort(fw[ok]); cw = np.cumsum(w[ok][o]); p25 = fw[ok][o][np.searchsorted(cw, .25 * cw[-1])]
    R['sim mix'] = prof(W, s1, pk, w, (np.abs(dt) < 80) & (fw >= p25)); line('sim mix (10-21 cut)', R['sim mix'])
    R['sim tilt0'] = prof(*S[0][:1], S[0][1], S[0][4], np.ones(len(S[0][1])), (np.abs(S[0][2]) < 80) & (S[0][3] >= np.nanpercentile(S[0][3], 25)))
    line('sim tilt 0 (10-21 cut)', R['sim tilt0'])
    for a in keys:
        Wa, sa, da, fa, pa = S[a]; R[f'tilt {a}'] = prof(Wa, sa, pa, np.ones(len(sa)), np.abs(da) < 80)
        print(f'   tilt {a}: FWHM p50 {np.nanmedian(fa[np.abs(da) < 80]):.0f} ns, |dt_nb| p50 {np.nanmedian(np.abs(da)):.0f}'); line(f'tilt {a} (|dt|<80)', R[f'tilt {a}'])
    json.dump(R, open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, f'mix_score_{tag}.json'), 'w'), indent=1)
