"""Even term of every sim arm (2026-10-21, audit next step 1). Strip registration makes a1 odd in the strip offset
delta and the X +-1 sharing s1 = alpha(+1) + alpha(-1) even in it. Per arm (X-only events, th00 + thy*, edge cut, delta
per peak column as in paired_odd.py): a1 and s1 in 100 um delta bins, the odd amplitude O = <sgn(delta) a1> over
|delta| >= 50 um, the even contrast E = s1(150 <= |delta| <= 350) - s1(|delta| < 100), and the slope k of
s1 = c + k |a1| across the delta bins (the data's observable, since data have no delta: det2 cells give k = -0.19 +- 0.06,
|a1| up to 0.2-0.3 at s1X - s1Y flat to +-0.04).
Usage: sim_even.py DIR [DIR ...]"""
import os, sys, glob, re, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L
P = 780.; ORG = -(512 - 1) / 2 * P; DB = np.arange(-400, 401, 100)
def tc(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
arms = {}
for D in sys.argv[1:]:
    for f in sorted(glob.glob(f'{D}/*_x.npz')):
        m = re.search(r'/([^/]+)_(th00|thy\d+)_s(\d+)_x\.npz$', f)
        if m: arms.setdefault(m.group(1), []).append(f)
OUT = {}
print(f'{"arm":14s} {"n":>6s}  {"O":>7s}  {"E":>7s}  {"E/O":>5s}  {"k(s1 vs |a1|)":>14s}  s1 per delta bin (-400..400 by 100)')
for pre, fs in sorted(arms.items()):
    A1, S1, DL = [], [], []
    for f in fs:
        if np.load(f)['wf'].ndim != 3: continue
        d = L.load(f); p = d['p'] % 512; ok = (p % 64 >= 4) & (p % 64 <= 59)
        DL.append((((ORG + tc(p) * P + 400.) % 800.) - 400.)[ok]); A1.append(d['a1'][ok]); S1.append((d['al'][1] + d['al'][-1])[ok])
    if not A1: continue
    a1, s1, dl = map(np.concatenate, (A1, S1, DL)); se = lambda v: v.std() / np.sqrt(len(v))
    m = np.abs(dl) >= 50; O = (np.sign(dl[m]) * a1[m]).mean(); Oe = se(np.sign(dl[m]) * a1[m])
    far = (np.abs(dl) >= 150) & (np.abs(dl) <= 350); near = np.abs(dl) < 100
    E = s1[far].mean() - s1[near].mean(); Ee = np.hypot(se(s1[far]), se(s1[near]))
    j = np.digitize(dl, DB) - 1; B = []
    for k in range(len(DB) - 1):
        q = j == k
        if q.sum() >= 30: B.append((DB[k] + 50, a1[q].mean(), se(a1[q]), s1[q].mean(), se(s1[q]), int(q.sum())))
    B = np.array(B); w = 1 / B[:, 4] ** 2; X = np.stack([np.ones(len(B)), np.abs(B[:, 1])], 1)
    c, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], B[:, 3] * np.sqrt(w), rcond=None)
    cov = np.linalg.inv((X * w[:, None]).T @ X); r = B[:, 3] - X @ c; s2 = max((w * r ** 2).sum() / max(len(B) - 2, 1), 1)
    ke = np.sqrt(cov[1, 1] * s2)
    OUT[pre] = dict(n=int(len(a1)), O=float(O), O_e=float(Oe), E=float(E), E_e=float(Ee), k=float(c[1]), k_e=float(ke),
                    bins=[dict(d=float(b[0]), a1=float(b[1]), a1_e=float(b[2]), s1=float(b[3]), s1_e=float(b[4]), n=int(b[5])) for b in B])
    print(f'{pre:14s} {len(a1):6d}  {O:+.4f}  {E:+.4f}  {E / O:5.1f}  {c[1]:+6.2f}±{ke:4.2f}    ' + ' '.join(f'{b[3]:.2f}' for b in B))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_even.json'), 'w'), indent=1)
