"""Test 1 follow-up (2026-10-23): could a pillar-driven a1 of the data's size hide from pillar.py's a1 search? For pillars to
make a field of |a1| ~ 0.1-0.2 they must act on a large share of events, not the 1.2 % on top of them. Inject into the real
survey events a1 += A for events within r of a synthetic lattice (s 1.002, t 0.2 deg), with A's sign following the
local pillar/pad phase (sign flips between the 4 x 4 patches in a checkerboard, as a moire would), and smear the event
position used for the injection by 0.53 mm (the centroid-to-track scatter, test 2) so the search sees what the data would.
Usage: pillar_inject.py -> pillar_inject.json"""
import os, json, numpy as np, pillar as PL
OUT = {}; rng = np.random.default_rng(5)
for b in ('det3 June-15', 'det2'):
    E = PL.R.ev(PL.R.BUILDS[b]); pc = lambda p: 64 * (p // 64) + 63 - p % 64
    x, y = 0.78 * pc(E['px']), 0.78 * pc(E['py']); a1 = E['alx'][:, 5] - E['alx'][:, 3]; s1 = E['alx'][:, 5] + E['alx'][:, 3]
    m = (np.abs(a1) < 3) & (np.abs(s1) < 3); x, y, a1 = x[m], y[m], a1[m]
    null = max(PL.zgrid(x, y, a1 - a1.mean(), PL.PIT * f, 4).max() for f in PL.NULLP)
    xt, yt = x + rng.normal(0, 0.53, len(x)), y + rng.normal(0, 0.53, len(x))
    sgn = np.where(((x // (PL.L / 4)) + (y // (PL.L / 4))) % 2 == 0, 1, -1); r_ = {}
    for r in (0.6, 1.0, 1.4):
        near = PL.pillars_near(xt, yt, r=r)
        for A in (0.1, 0.2, 0.5):
            a = a1 + A * sgn * near; z = PL.zgrid(x, y, a - a.mean(), PL.PIT, 4).max()
            r_[f'r{r}_A{A}'] = dict(frac=float(near.mean()), field=float(A * near.mean()), z=float(z))
            print(f'{b:13s} r {r} mm ({100 * near.mean():4.1f} % of events) A {A}: patch field {A * near.mean():.3f}, Zmax {z:6.1f} (off-lattice null max {null:.1f})', flush=True)
    OUT[b] = dict(null_max=float(null), inject=r_)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pillar_inject.json'), 'w'), indent=1)
