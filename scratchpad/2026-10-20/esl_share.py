"""Do the ESL stripes show up AT ALL? Shares (X +-1, Y +-1, A_x/A_y) vs true column phase mod 40, per build/class.
If the shares show a 40-column modulation the stripe/gap position is resolved and the copy test has power."""
import os, numpy as np
z = np.load(os.path.expanduser('~/x17/nb_asym_raw/esl_cache.npz'), allow_pickle=True)
col = lambda ch: 64 * (ch % 512 // 64) + 63 - ch % 64
def sine(y, ph):
    mu = np.array([y[ph == i].mean() for i in range(40)]); n = np.array([(ph == i).sum() for i in range(40)])
    se = np.array([y[ph == i].std() for i in range(40)]) / np.sqrt(np.maximum(n, 1)); ok = n > 5
    th = 2 * np.pi * np.arange(40) / 40; w = 1 / se[ok] ** 2
    X = np.c_[np.ones(40), np.cos(th), np.sin(th)][ok]; A = X.T @ (w[:, None] * X); c = np.linalg.solve(A, X.T @ (w * mu[ok]))
    cov = np.linalg.inv(A); amp = np.hypot(c[1], c[2]); ea = np.sqrt((c[1]**2*cov[1,1] + c[2]**2*cov[2,2])) / amp
    ph0 = np.degrees(np.arctan2(c[2], c[1])) % 360
    return c[0], amp, ea, ph0, (((mu[ok] - c[0]) / se[ok]) ** 2).sum() / (ok.sum() - 1)
import sys
blds = sys.argv[1:] or ['det3 June-15', 'det3 May', 'det2', 'det6', 'det1']
for bu in blds:
    for cls in ('ypad', 'xpad'):
        m = (z['build'] == bu) & (z['ypad'] if cls == 'ypad' else ~z['ypad']) & (z['px'] % 64 >= 4) & (z['px'] % 64 <= 59)
        ph = col(z['px'][m]) % 40
        obs = dict(xsh=z['ax'][m][:, 3] + z['ax'][m][:, 5], ysh=z['ay'][m][:, 3] + z['ay'][m][:, 5],
                   lr=np.log(z['Ax'][m] / z['Ay'][m]), copy=z['ax'][m][:, 5] - z['ax'][m][:, 3], x2=z['ax'][m][:, 2] + z['ax'][m][:, 6])
        s = f'{bu:13s} {cls} n={m.sum():5d}'
        for k, y in obs.items():
            c0, amp, ea, p0, chi = sine(y, ph); s += f' | {k} {c0:.3f} A={amp:.3f}±{ea:.3f} @{p0:3.0f} χ²={chi:.1f}'
        print(s)
