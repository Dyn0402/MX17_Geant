"""(C2) Charge bookkeeping at matched width (2026-10-20). Per Y-share bin and pad class: X(n-1)/SY, X(n+1)/SY,
Xc/SY, SX/SY (integrals of positive samples, 9-channel windows), for
 data builds; det2 / det3 May with a plain copy x(n+1) += a x(n) injected (a = 0.2, 0.3), peak re-found;
 sim arms xo0 and xom100/200/300 (th00 + thy*) with and without an added copy.
A copy adds to n+1 and leaves n-1; an offset moves charge from n-1 to n+1. Usage: probe_charge_bins.py"""
import os, re, glob, json, numpy as np
import r5lib as R
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
def classify(Wx, Wy, px, py):
    """peak re-found inside the window (after an injection), parity, Y share, Y-pad class read from the data."""
    c = Wx.max(2).argmax(1); sh = c - 4
    ok = np.abs(sh) <= 1
    idx = np.clip(np.arange(9)[None, :] + sh[:, None], 0, 8)
    Wx2 = np.take_along_axis(Wx, idx[:, :, None], 1); Wx2[sh == 1, 8] = 0; Wx2[sh == -1, 0] = 0
    px2 = px + sh
    cy = Wy[:, 4]; n2 = np.einsum('ij,ij->i', cy, cy)
    ysh = (np.einsum('ij,ij->i', Wy[:, 3], cy) + np.einsum('ij,ij->i', Wy[:, 5], cy)) / n2
    par = (px2 + py) % 2; yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean())
    return Wx2[ok], Wy[ok], ysh[ok], yp[ok]
def table(Wx, Wy, ysh, yp):
    qx = np.clip(Wx, 0, None).sum(2); qy = np.clip(Wy, 0, None).sum(2); SY = qy.sum(1); b = np.digitize(ysh, ED) - 1
    out = {}
    for cls, cm in (('Y', yp), ('X', ~yp)):
        rows = []
        for k in range(len(ED) - 1):
            m = cm & (b == k) & (SY > 0)
            if m.sum() < 25: rows.append(None); continue
            f = lambda v: (float(np.mean(v[m] / SY[m])), float(np.std(v[m] / SY[m]) / np.sqrt(m.sum())))
            rows.append(dict(n=int(m.sum()), m1=f(qx[:, 3]), p1=f(qx[:, 5]), c=f(qx[:, 4]), sx=f(qx.sum(1))))
        out[cls] = rows
    return out
def show(name, t):
    for cls in ('Y', 'X'):
        s = ' | '.join('       -        ' if r is None else f'{r["m1"][0]:.3f} {r["p1"][0]:.3f} {r["c"][0]:.2f}' for r in t[cls])
        print(f'{name:22s} {cls}pad | {s}')
def main():
    OUT = {}
    print('per Y-share bin', ED.tolist(), ': X(n-1)/SY  X(n+1)/SY  Xc/SY')
    for b in ('det3 June-15', 'det3 May', 'det2', 'det1'):
        E = R.ev(R.BUILDS[b]); args = (E['Wx'], E['Wy'], E['px'], E['py'])
        t = table(*classify(*args)); OUT[b] = t; show(b, t)
        if b in ('det3 May', 'det2'):
            for a in (.2, .3):
                Wx = E['Wx'].copy(); Wx[:, 1:] += a * E['Wx'][:, :-1]
                t = table(*classify(Wx, E['Wy'], E['px'], E['py'])); OUT[f'{b} +copy {a}'] = t; show(f'{b} +copy {a}', t)
    # sim arms
    SRC = os.path.expanduser('~/x17/nb_asym_raw/simxoff')
    def sim(arm):
        Wx, Wy, px, py = [], [], [], []
        for fx in sorted(glob.glob(f'{SRC}/{arm}_*_x.npz')):
            if not re.search('_th00_|_thy', fx): continue
            fy = fx[:-6] + '_y.npz'
            zx, zy = np.load(fx), np.load(fy)
            if zx['wf'].ndim != 3 or zy['wf'].ndim != 3: continue
            dx, dy = R.L.load(fx), R.L.load(fy)
            kx = dx['fi'].astype(np.int64) * 10**9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10**9 + dy['ev']
            _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
            Wx.append(dx['W'][ix]); Wy.append(dy['W'][iy]); px.append(dx['p'][ix]); py.append(dy['p'][iy])
        return [np.concatenate(v).astype(np.float32) for v in (Wx, Wy, px, py)]
    for arm in ('xo0', 'xom100', 'xom200', 'xom300'):
        Wx0, Wy, px, py = sim(arm)
        for a in (0, .2):
            Wx = Wx0.copy()
            if a: Wx[:, 1:] += a * Wx0[:, :-1]
            t = table(*classify(Wx, Wy, px, py)); OUT[f'sim {arm} +copy {a}'] = t; show(f'sim {arm} +copy {a}', t)
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'probe_charge_bins.json'), 'w'))


if __name__ == '__main__':
    main()
