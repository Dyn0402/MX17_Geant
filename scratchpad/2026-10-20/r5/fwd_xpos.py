"""Fast forward model (static solved weighting fields, fwd_variants.py machinery) scored on all fingerprints at once
(2026-10-20): Y-pad / X-pad class asymmetry and Y-pad fraction in Y-share QUANTILE bins (the static model's sharing
scale differs from data), and the x-position conditional E[v | u] for Y-pad events (u = clean-side share, v = n+1 side
share) in u-quantile bins, as an effect against the symmetric model (data: det3 June-15 minus det2 + det3 May).
Mechanisms ('dn' = line k-1 = data's DAQ n+1):
  copy a      S[k-1] += a S[k]
  shift D     X lines read the cloud at x + D (D < 0 toward k-1)
  mix D,w     (1-w) S(0) + w S(D)
  skip t      line k hands a fraction t of its signal to line k-2... i.e. S'[k-1] = S[k-1] + t S[k+1], S'[k+1] = (1-t) S[k+1]
              applied to every line: x'(j) = (1-t) x(j) + t x(j+2)  (channel transfer two lines over, toward k-1)
Usage: fwd_xpos.py [--n 60000]"""
import os, sys, json, argparse, numpy as np
from scipy.interpolate import RegularGridInterpolator
ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=60000); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, '..', '..', '2026-10-19', 'bus_par')
P = 780.; dg = 20.; EX = np.load(f'{B}/d20_nobus_esl.npy').astype(float)
nx, ny = EX.shape; nrow = int(round(ny * dg / P / 2))
xc = (np.arange(nx) + .5) * dg - P / 2; yc = (np.arange(ny) + .5) * dg - nrow * P
Ex = np.concatenate([EX[-10:], EX, EX[:10]]); xx = np.concatenate([xc[-10:] - 2 * P, xc, xc[:10] + 2 * P])
_f = RegularGridInterpolator((xx, yc), Ex, bounds_error=False, fill_value=0.)
fW = lambda x, y: _f(np.stack([((x + P / 2) % (2 * P)) - P / 2, y], -1))
g = np.linspace(-3, 3, 13); GX, GY = np.meshgrid(g, g, indexing='ij'); W = np.exp(-.5 * (GX ** 2 + GY ** 2)); W /= W.sum()
K = np.arange(-5, 6)
rng = np.random.default_rng(11); n = A.n
cx = rng.uniform(-P, P, n); cy = rng.uniform(-P / 2, P / 2, n); s = rng.uniform(250, 650, n)
px = cx[:, None, None] + s[:, None, None] * GX[None]; py = cy[:, None, None] + s[:, None, None] * GY[None]
SY = np.stack([(fW(py - P * c + P / 2, px + P / 2 - P * c) * W).sum((1, 2)) for c in K], 1)
_c = {}
def SX(D=0.):
    if D not in _c: _c[D] = np.stack([(fW(px - P * k, py - P * k + D) * W).sum((1, 2)) for k in K], 1)
    return _c[D]
NQ = 6; NU = 6
def analyse(S):
    ii = np.arange(n); Sx = S * rng.normal(1, .03, S.shape); Sy = SY * rng.normal(1, .03, SY.shape)
    kx = Sx[:, 2:-2].argmax(1) + 2; ky = Sy[:, 2:-2].argmax(1) + 2
    c = Sx[ii, kx]; v = Sx[ii, kx - 1] / c; u = Sx[ii, kx + 1] / c; dd = v - u
    ysh = (Sy[ii, ky - 1] + Sy[ii, ky + 1]) / Sy[ii, ky]; par = (K[kx] + K[ky]) % 2
    yp = par == (0 if ysh[par == 0].mean() < ysh[par == 1].mean() else 1)
    return dict(dd=dd, u=u, v=v, ysh=ysh, yp=yp)
def fingerprint(e, qy=None, qu=None):
    qy = np.quantile(e['ysh'], np.linspace(0, 1, NQ + 1)) if qy is None else qy
    b = np.clip(np.digitize(e['ysh'], qy[1:-1]), 0, NQ - 1); o = dict(AY=[], AX=[], fY=[])
    for k in range(NQ):
        m = b == k
        o['AY'].append(e['dd'][m & e['yp']].mean()); o['AX'].append(e['dd'][m & ~e['yp']].mean()); o['fY'].append(e['yp'][m].mean())
    # x-position conditional: Y-pad events in the middle Y-share quantiles 1..4, u quantiles
    m = e['yp'] & (b >= 1) & (b <= 4)
    qu = np.quantile(e['u'][m], np.linspace(0, 1, NU + 1)) if qu is None else qu
    bu = np.clip(np.digitize(e['u'], qu[1:-1]), 0, NU - 1)
    o['vu'] = [e['v'][m & (bu == k)].mean() for k in range(NU)]; o['uq'] = qu.tolist()
    o['a1'] = float(e['dd'].mean())
    return {k: (np.array(v) if k not in ('a1', 'uq') else v) for k, v in o.items()}
# ---- data in the same estimator (amplitude ratios at the peak sample are close to the in-phase alpha) ----
sys.path.insert(0, HERE)
import r5lib as R
def data_fp(Es):
    e = dict(dd=np.concatenate([E['alx'][:, 5] - E['alx'][:, 3] for E in Es]), u=np.concatenate([E['alx'][:, 3] for E in Es]),
             v=np.concatenate([E['alx'][:, 5] for E in Es]), ysh=np.concatenate([E['ysh'] for E in Es]), yp=np.concatenate([E['ypad'] for E in Es]))
    return e
D3 = data_fp([R.ev(R.BUILDS['det3 June-15'])]); REF = data_fp([R.ev(R.BUILDS['det2']), R.ev(R.BUILDS['det3 May'])])
fd = fingerprint(D3); fr = fingerprint(REF)
# for the data effect in vu, use det3's own u quantiles on both (same x gauge)
fr_u = fingerprint(REF, qu=np.array(fd['uq']))
DAT = dict(AY=fd['AY'], AX=fd['AX'], dfY=fd['fY'] - fr['fY'], dvu=fd['vu'] - fr_u['vu'])
f = lambda a: ' '.join(f'{x:+.2f}' for x in a)
print(f'{"DATA det3 June-15":18s} a1 {fd["a1"]:+.3f} | AY [{f(DAT["AY"])}] AX [{f(DAT["AX"])}] dfY [{f(DAT["dfY"])}] dvu [{f(DAT["dvu"])}]')
def main():
    base = analyse(SX(0.)); fb = fingerprint(base)
    def score(fm, fmu):
        r = dict(AY=fm['AY'], AX=fm['AX'], dfY=fm['fY'] - fb['fY'], dvu=fmu['vu'] - fingerprint(base, qu=np.array(fmu['uq']))['vu'])
        chi = {k: float(np.sum((r[k] - DAT[k]) ** 2) / .03 ** 2) for k in DAT}
        return r, chi
    MODELS = [('none', SX(0.))]
    for a in (.1, .2, .3):
        S = SX(0.).copy(); S[:, :-1] += a * SX(0.)[:, 1:]; MODELS.append((f'copy {a}', S))
    for D in (-100., -200., -300.): MODELS.append((f'shift {D:.0f}', SX(D)))
    for D in (-200., -390., -585., -780.):
        for w in (.15, .3, .5):
            MODELS.append((f'mix {D:.0f},{w}', (1 - w) * SX(0.) + w * SX(D)))
    for t in (.1, .2, .3, .4):
        S0 = SX(0.); S = (1 - t) * S0.copy(); S[:, :-2] += t * S0[:, 2:]; MODELS.append((f'skip {t}', S))
    for D, w, a in ((-200., 1., .1), (-780., .3, .1), (-390., .5, .1)):
        S0 = (1 - w) * SX(0.) + w * SX(D); S = S0.copy(); S[:, :-1] += a * S0[:, 1:]; MODELS.append((f'mix {D:.0f},{w}+copy {a}', S))
    OUT = dict(data={k: v.tolist() for k, v in DAT.items()}, models={})
    for nm, S in MODELS:
        e = analyse(S); fm = fingerprint(e); fmu = fm   # u quantiles of the model itself
        r, chi = score(fm, fmu); tot = sum(chi.values())
        OUT['models'][nm] = dict(a1=fm['a1'], chi=chi, **{k: v.tolist() for k, v in r.items()})
        print(f'{nm:22s} a1 {fm["a1"]:+.3f} | AY [{f(r["AY"])}] AX [{f(r["AX"])}] dfY [{f(r["dfY"])}] dvu [{f(r["dvu"])}] | chi ' +
              ' '.join(f'{k} {v:4.0f}' for k, v in chi.items()) + f' tot {tot:5.0f}')
    json.dump(OUT, open(os.path.join(HERE, 'fwd_xpos.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
