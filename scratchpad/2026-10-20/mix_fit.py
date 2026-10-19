"""Can ANY fixed linear channel mixing reproduce det3 June-15? (2026-10-20)
A readout-side effect (M1 card, connector, cable, FEU) is channel-level: x'(n) = x(n) + sum_k m_k x(n+k), the same for
every event whatever the charge does on the pads (linearity is measured). Inject that into copy-free boards' raw
windows, run det3's analysis (raw peak, X/Y parity, Y-share bins), fit m_{-1} (the observed n -> n+1 copy), m_{+1}
(the reverse, may be negative: an antisymmetric, derivative-like mixing) and m_{+-2}. Observables per Y-share bin:
class asymmetry, class +-1 sum, class fraction. Usage: mix_fit.py"""
import os, sys, json, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uncopy_parity as UP
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
def analyse(W, px0, py, ysh, wts):
    A_ = W.max(2); sh = A_[:, 3:6].argmax(1) - 1; c0 = 4 + sh; ii = np.arange(len(W))
    c = W[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
    ap = np.einsum('ij,ij->i', W[ii, c0 + 1], c) / n2; am = np.einsum('ij,ij->i', W[ii, c0 - 1], c) / n2
    px = px0 + sh; par = (px + py) % 2; ok = (px % 64 >= 4) & (px % 64 <= 59)
    yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean()); b = np.digitize(ysh, ED) - 1; r = []
    for k in range(len(ED) - 1):
        mb = ok & (b == k)
        for cl in (True, False):
            m = mb & (yp == cl); w = wts[m]; sw = w.sum()
            for v in (ap - am, ap + am):
                mu = (w * v[m]).sum() / sw if sw > 0 else np.nan
                r.append((mu, np.sqrt((w * (v[m] - mu) ** 2).sum() / sw / max(m.sum(), 1)) if sw > 0 else np.nan))
        fy = wts[mb & yp].sum() / wts[mb].sum(); r.append((fy, np.sqrt(fy * (1 - fy) / max(mb.sum(), 1))))
    return np.array(r)
def mix(W, m):
    out = W.copy()
    for k, v in m.items():
        if not v: continue
        if k > 0: out[:, :-k] += v * W[:, k:]
        else: out[:, -k:] += v * W[:, :k]
    return out
if __name__ == '__main__':
    Wd, pxd, pyd, yshd = UP.events('det3 June-15')
    data = analyse(Wd.astype(float), pxd, pyd, yshd, np.ones(len(Wd)))
    hd = np.histogram(yshd, ED)[0] / len(yshd)
    OUT = dict(data=data.tolist(), fits={})
    def chi(t, sel=None):
        ok = np.isfinite(t[:, 0]) & np.isfinite(data[:, 0])
        return float((((t[ok, 0] - data[ok, 0]) ** 2) / (t[ok, 1] ** 2 + data[ok, 1] ** 2)).sum()), int(ok.sum())
    for tb in ('det2', 'det3 May', 'det1'):
        W, px, py, ysh = UP.events(tb); W = W.astype(float)
        ht = np.histogram(ysh, ED)[0] / len(ysh); wts = (hd / np.maximum(ht, 1e-9))[np.clip(np.digitize(ysh, ED) - 1, 0, len(ED) - 2)]
        t0 = analyse(W, px, py, ysh, wts); c0, nd = chi(t0)
        res = []
        for mm1, mp1 in itertools.product(np.arange(0, .61, .05), np.arange(-.3, .31, .05)):
            t = analyse(mix(W, {-1: mm1, 1: mp1}), px, py, ysh, wts); res.append((chi(t)[0], mm1, mp1, 0., 0., t))
        res.sort(key=lambda v: v[0]); b = res[0]
        # refine with +-2 around the best
        for mm2, mp2 in itertools.product(np.arange(-.1, .11, .05), np.arange(-.1, .11, .05)):
            t = analyse(mix(W, {-1: b[1], 1: b[2], -2: mm2, 2: mp2}), px, py, ysh, wts); res.append((chi(t)[0], b[1], b[2], mm2, mp2, t))
        res.sort(key=lambda v: v[0]); b = res[0]
        uni = min((r for r in res if abs(r[2]) < 1e-9 and r[3] == 0 and r[4] == 0), key=lambda v: v[0])
        print(f'template {tb:9s}: no mixing chi2 {c0:.0f}/{nd} | plain copy best m-1={uni[1]:.2f}: {uni[0]:.0f} | '
              f'best general m-1={b[1]:.2f} m+1={b[2]:+.2f} m-2={b[3]:+.2f} m+2={b[4]:+.2f}: chi2 {b[0]:.0f}/{nd}')
        t = b[5]
        for k in range(len(ED) - 1):
            o = 5 * k
            print(f'   Ysh {ED[k]:.2f}-{ED[k+1]:.2f}: Ypad asym {data[o,0]:+.2f}/{t[o,0]:+.2f} sum {data[o+1,0]:.2f}/{t[o+1,0]:.2f} | '
                  f'Xpad asym {data[o+2,0]:+.2f}/{t[o+2,0]:+.2f} sum {data[o+3,0]:.2f}/{t[o+3,0]:.2f} | fY {data[o+4,0]:.2f}/{t[o+4,0]:.2f}')
        OUT['fits'][tb] = dict(none=c0, ndf=nd, plain=dict(m=uni[1], chi2=uni[0]),
                               best=dict(m_m1=b[1], m_p1=b[2], m_m2=b[3], m_p2=b[4], chi2=b[0], model=t.tolist()))
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mix_fit.json'), 'w'), default=float)
