"""Score models on the x-position fingerprint (2026-10-20): Y-pad events, E[alpha(n+1) | alpha(n-1)] in Y-share x u bins
(xpos_probe.py binning). Effect = det3 June-15 (or a session) minus the copy-free reference (det2 + det3 May);
model effect = model minus the sim D0 arm. Also prints each model's effect row for Ysh 0.4-0.85 pooled.
Usage: eval_xpos.py"""
import os, glob, re, json, numpy as np
import r5lib as R, xpos_probe as XP
from eval_mix import arm, SIM, MIX, SESS
def ev_from(a):
    Wx, Wy, px, py = [np.asarray(v) for v in a]; Wx = Wx.astype(float); Wy = Wy.astype(float)
    cx, cy = Wx[:, 4], Wy[:, 4]; nx = np.einsum('ij,ij->i', cx, cx); ny = np.einsum('ij,ij->i', cy, cy)
    alx = np.stack([np.einsum('ij,ij->i', Wx[:, 4 + k], cx) / nx for k in range(-4, 5)], 1)
    ysh = (np.einsum('ij,ij->i', Wy[:, 3], cy) + np.einsum('ij,ij->i', Wy[:, 5], cy)) / ny; par = (px + py) % 2
    return dict(alx=alx, ysh=ysh, ypad=par == int(ysh[par == 1].mean() < ysh[par == 0].mean()))
def pooled(E, ylo=.4, yhi=.85):
    u = E['alx'][:, 3]; v = E['alx'][:, 5]; m = E['ypad'] & (E['ysh'] >= ylo) & (E['ysh'] < yhi); r = []
    for a, b in zip(XP.UB[:-1], XP.UB[1:]):
        k = m & (u >= a) & (u < b); r.append((v[k].mean(), v[k].std() / np.sqrt(k.sum())) if k.sum() >= 30 else (np.nan, np.nan))
    return np.array(r)
if __name__ == '__main__':
    ref = pooled({k: np.concatenate([R.ev(R.BUILDS[b])[k] for b in ('det2', 'det3 May')]) for k in ('alx', 'ysh', 'ypad')})
    data = {k: pooled(R.ev(SESS[k])) for k in ('6-22/23 1000V', '6-23/24 600V', '6-27/28 1000V', 'all June-15')}
    eff = {k: np.c_[v[:, 0] - ref[:, 0], np.hypot(v[:, 1], ref[:, 1])] for k, v in data.items()}
    f = lambda r: ' '.join('  -   ' if np.isnan(x) else f'{x:+.3f}' for x in r[:, 0])
    print('u bins:', XP.UB.tolist(), ' (Y-pad, Ysh 0.4-0.85; effect vs copy-free)')
    for k, e in eff.items(): print(f'DATA {k:16s} {f(e)}')
    pats = {'D0': f'{SIM}/xo0_{{ang}}_s{{s}}_x.npz'}
    for off in (100, 200, 300, 450, 600, 780): pats[f'D-{off} w1'] = f'{SIM}/xom{off}_{{ang}}_s{{s}}_x.npz'
    for fx in sorted(glob.glob(f'{MIX}/m*_th00_s5_x.npz')):
        m = re.search(r'm(\d+)_w([\d.]+)_th00', fx); pats[f'D-{m.group(1)} w{m.group(2)}'] = f'{MIX}/m{m.group(1)}_w{m.group(2)}_{{ang}}_s{{s}}_x.npz'
    M = {}
    for nm, p in pats.items():
        a = arm(p)
        if a is not None: M[nm] = pooled(ev_from(a))
    a0 = arm(pats['D0']); E0 = ev_from(a0)
    for c in (.15, .25):
        Wx = a0[0].astype(float).copy(); Wx[:, 1:] += c * a0[0][:, :-1]; M[f'copy {c}'] = pooled(ev_from((Wx, *a0[1:])))
    OUT = {}
    for nm, r in M.items():
        e = np.c_[r[:, 0] - M['D0'][:, 0], np.hypot(r[:, 1], M['D0'][:, 1])]
        chi = {k: float(np.nansum((e[:, 0] - d[:, 0]) ** 2 / (e[:, 1] ** 2 + d[:, 1] ** 2 + .02 ** 2))) for k, d in eff.items()}
        OUT[nm] = dict(eff=e.tolist(), chi=chi)
        print(f'MODEL {nm:14s} {f(e)} | chi2 ' + ' '.join(f'{k[:5]} {v:5.0f}' for k, v in chi.items()))
    json.dump(dict(UB=XP.UB.tolist(), data={k: v.tolist() for k, v in eff.items()}, models=OUT),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eval_xpos.json'), 'w'), indent=1)
