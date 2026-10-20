"""(A) Where across the Y pad (in x) does the n+1 excess come from? (2026-10-20)
Y-pad events in Y-share bins; x gauge u = alpha(n-1) (clean side: small u = charge toward the n+1 side).
excess(u) = E[alpha(n+1) | u] on det3 June-15 minus the same on copy-free boards (det2 + det3 May pooled) in the same
Y-share bin. A gap-local coupling peaks at small u (charge near the n+1-side gap); a pad-wide one is flat in u.
Also: E[alpha(n+1) | u] for det2 + injected plain copy 0.25 (the copy shape for reference).
Usage: xpos_probe.py"""
import json, os, numpy as np
import r5lib as R
YB = [(.4, .55), (.55, .7), (.7, .85), (.85, 1.0)]; UB = np.array([-.05, .05, .1, .15, .2, .25, .3, .4, .6])
def cond(E, mask):
    u = E['alx'][:, 3]; v = E['alx'][:, 5]; out = []
    for lo, hi in YB:
        m = mask & (E['ysh'] >= lo) & (E['ysh'] < hi); row = []
        for a, b in zip(UB[:-1], UB[1:]):
            k = m & (u >= a) & (u < b)
            row.append((float(v[k].mean()), float(v[k].std() / np.sqrt(k.sum())), int(k.sum())) if k.sum() >= 30 else (np.nan, np.nan, int(k.sum())))
        out.append(row)
    return np.array(out)   # (ybin, ubin, 3)
def main():
    D3 = R.ev(R.BUILDS['det3 June-15']); C2 = R.ev(R.BUILDS['det2']); CM = R.ev(R.BUILDS['det3 May'])
    ref = {k: np.concatenate([C2[k], CM[k]]) for k in ('alx', 'ysh', 'ypad')}
    r3 = cond(D3, D3['ypad']); rr = cond(ref, ref['ypad'])
    # plain copy on det2: alpha(+1) gets + a * (central's in-phase = 1), alpha(-1) gets + a * alpha(-2)
    a = .25; cp = {k: C2[k].copy() for k in ('alx', 'ysh', 'ypad')}
    cp['alx'][:, 5] = C2['alx'][:, 5] + a; cp['alx'][:, 3] = C2['alx'][:, 3] + a * C2['alx'][:, 2]
    rc = cond(cp, cp['ypad']); r2 = cond(C2, C2['ypad'])
    print('u bins (alpha(n-1)):', UB.tolist())
    for j, (lo, hi) in enumerate(YB):
        print(f'Ysh {lo}-{hi}')
        for nm, r in (('det3 June-15', r3), ('ref det2+May', rr), ('det2 +copy .25', rc)):
            print(f'  {nm:14s} E[a+1|u] ' + ' '.join('  -   ' if np.isnan(x[0]) else f'{x[0]:.3f}' for x in r[j]) + '   n ' + ' '.join(f'{int(x[2]):4d}' for x in r[j]))
        print(f'  {"excess det3":14s}          ' + ' '.join('  -   ' if np.isnan(x) else f'{x:+.3f}' for x in r3[j, :, 0] - rr[j, :, 0]))
        print(f'  {"excess copy":14s}          ' + ' '.join('  -   ' if np.isnan(x) else f'{x:+.3f}' for x in rc[j, :, 0] - r2[j, :, 0]))
    json.dump(dict(UB=UB.tolist(), YB=YB, det3=r3.tolist(), ref=rr.tolist(), copy=rc.tolist(), det2=r2.tolist()),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'xpos_probe.json'), 'w'))


if __name__ == '__main__':
    main()
