"""Paired comparison of sim arms on the strip-registration odd amplitude (2026-10-21).
For every (angle, seed) present in the baseline arm AND each variant arm: X-only events, delta per peak column as in
sim_delta.py ('col' orientation), odd amplitude O = <sgn(delta) a1> over |delta| >= 50 um, all events and the narrow
half (s1 below the baseline's median). Same seeds share clusters, diffusion, gain and noise draws (the funnel has its
own rng stream), so the difference is the variant's effect; errors are from the per-(angle, seed) spread of
the paired differences.
Usage: paired_odd.py <base prefix> <base dir> <variant prefix> <variant dir> [<variant prefix> <variant dir> ...]"""
import os, sys, glob, re, json, numpy as np
import r5lib as R
P = 780.; ORG = -(512 - 1) / 2 * P
def arm(pre, d):
    out = {}
    for f in sorted(glob.glob(os.path.expanduser(f'{d}/{pre}_*_x.npz'))):
        m = re.search(rf'{re.escape(pre)}_(th00|thy\d+)_s(\d+)_x\.npz$', f)
        if not m or np.load(f)['wf'].ndim != 3: continue
        out[(m.group(1), int(m.group(2)))] = f
    return out
def odd(f, med=None):
    d = R.L.load(f); p = d['p'] % 512; ok = (p % 64 >= 4) & (p % 64 <= 59)
    col = 64 * (p // 64) + 63 - p % 64; delta = ((ORG + col * P + 400.) % 800.) - 400.
    s1 = d['al'][1] + d['al'][-1]; a1 = d['a1']; m = ok & (np.abs(delta) >= 50)
    med = np.median(s1[m]) if med is None else med
    o = np.sign(delta[m]) * a1[m]; nar = s1[m] < med
    return dict(all=float(o.mean()), nar=float(o[nar].mean()), n=int(m.sum()), s1=float(s1[m].mean())), med
B = arm(sys.argv[1], sys.argv[2]); OUT = {}
for k in range(3, len(sys.argv), 2):
    pre, d = sys.argv[k], sys.argv[k + 1]; V = arm(pre, d); keys = sorted(set(B) & set(V)); rows = []
    for key in keys:
        b, med = odd(B[key]); v, _ = odd(V[key], med); rows.append((key, b, v))
    if not rows: print(f'{pre}: no common (angle, seed)'); continue
    for q in ('all', 'nar'):
        bb = np.array([r[1][q] for r in rows]); vv = np.array([r[2][q] for r in rows]); w = np.array([r[1]['n'] for r in rows], float)
        mb, mv = np.average(bb, weights=w), np.average(vv, weights=w); dd = vv - bb
        e = np.std(dd) / np.sqrt(len(dd)) if len(dd) > 1 else np.nan
        print(f'{pre:6s} vs {sys.argv[1]} [{q}]  {len(rows)} pairs  base O {mb:+.4f}  variant O {mv:+.4f}  '
              f'ratio {mv / mb:.2f}  diff {np.average(dd, weights=w):+.4f} ± {e:.4f}')
        OUT.setdefault(pre, {})[q] = dict(pairs=len(rows), base=mb, var=mv, ratio=mv / mb, diff=float(np.average(dd, weights=w)), err=float(e))
    sb = np.mean([r[1]['s1'] for r in rows]); sv = np.mean([r[2]['s1'] for r in rows])
    print(f'{"":6s}    X sharing s1: base {sb:.3f}  variant {sv:.3f}'); OUT[pre]['s1'] = dict(base=sb, var=sv)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paired_odd.json'), 'w'), indent=1)
