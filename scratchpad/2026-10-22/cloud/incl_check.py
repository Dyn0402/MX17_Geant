"""Is the data's flat early sharing residual x-inclination? (2026-10-22) The sim's near-vertical X sample has exactly
zero x-inclination (th00, thy*), the data's |dt_nb| < 80 ns sample has some. Inclined: the shallow (early) electrons sit
on one side of the peak pad and the deep (late) ones on the other, so the early neighbour excess is ONE-SIDED and
anti-correlated with the late one. Per event, d(t) = (W+1 - W-1)/W_c at r25 and at f50; split by |dt_nb|.
Usage: incl_check.py   (NB_RAW as coupling_checks)"""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NB_RAW', os.path.expanduser('~/x17/nb_asym_raw'))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import coupling_checks as CC, prompt_share as PS, asym_lib as L
RAW = os.environ['NB_RAW']; HW = 4; DT = 60.

def load(files):
    W, s1 = [], []
    for f in files:
        d = L.load(f); p = d['p'] % 512; m = (p % 64 >= 4) & (p % 64 <= 59)
        W.append(d['W'][m]); s1.append(d['al'][1][m] + d['al'][-1][m])
    W = np.concatenate(W).astype(float); s1 = np.concatenate(s1); c = W[:, HW]; pk = c.argmax(1)
    ts = np.arange(W.shape[2]); pos = lambda w: np.clip(w, 0, None)
    tc = lambda w: (pos(w) * ts).sum(1) / np.maximum(pos(w).sum(1), 1e-9)
    dtn = (tc(W[:, HW + 1]) - tc(W[:, HW - 1])) * DT
    fwhm = (PS.cross(c, pk, .5, False) - PS.cross(c, pk, .5)) * DT
    keep = fwhm >= np.nanpercentile(fwhm, CC.FW_PCT)
    return W[keep], s1[keep], dtn[keep]

def run(name, files):
    W, s1, dtn = load(files); c = W[:, HW]; pk = c.argmax(1); D = W[:, HW + 1] - W[:, HW - 1]; S = W[:, HW + 1] + W[:, HW - 1]
    te, tl, t5 = PS.cross(c, pk, .25), PS.cross(c, pk, .5, False), PS.cross(c, pk, .05)
    de, dl = PS.at(D, te) / PS.at(c, te), PS.at(D, tl) / PS.at(c, tl)
    print(f'\n{name}: n {len(s1)}')
    print('   |dt_nb| ns    n    corr(d_early,d_late)   r05    r25    peak   f50     (s1 0.2-0.45)')
    for lo, hi in ((0, 30), (30, 80), (80, 200), (200, 1e4)):
        m = (np.abs(dtn) >= lo) & (np.abs(dtn) < hi) & np.isfinite(de) & np.isfinite(dl) & (np.abs(de) < 3) & (np.abs(dl) < 3)
        b = m & (s1 >= .2) & (s1 < .45)
        rr = []
        for t in (t5, te, pk.astype(float), tl):
            g = b & np.isfinite(t); nu, dd = PS.at(S, t)[g], PS.at(c, t)[g]; rr.append(nu.sum() / dd.sum() if g.sum() > 20 else np.nan)
        cc = np.corrcoef(de[m], dl[m])[0, 1] if m.sum() > 20 else np.nan
        print(f'   {lo:4.0f}-{hi:<5.0f} {m.sum():6d}    {cc:+.2f}                 ' + '  '.join(f'{v:5.3f}' for v in rr) + f'   [{b.sum()}]')

if __name__ == '__main__':
    for n, t in (('det3 May X', 'd3may_x'), ('det3 June-15 X', 'd3jun_x'), ('det3 May Y', 'd3may_y'), ('det3 June-15 Y', 'd3jun_y')):
        run(n, [f'{RAW}/{t}.npz'])
    run('sim glg0 X (th00+thy)', sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_x.npz')))
    run('sim glg0 Y (th00)', sorted(glob.glob(f'{RAW}/simgapleak/glg0_th00_*_y.npz')))
    run('sim glg0 Y (thy, inclined in Y)', sorted(glob.glob(f'{RAW}/simgapleak/glg0_thy*_y.npz')))
