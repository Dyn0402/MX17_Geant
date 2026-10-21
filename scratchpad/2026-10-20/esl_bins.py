"""ESL strip phase, binned for the deck (2026-10-20): Y-pad / X-pad class copy and the Y +-1 share vs true pad column
mod 40 (connector-reversed map), in 10 bins of 4 columns = 80 um steps of the strip offset under the pad.
Plus the sine-fit amplitude over all 40 phases. Writes esl_bins.json."""
import os, json, numpy as np
z = np.load(os.path.expanduser('~/x17/nb_asym_raw/esl_cache.npz'), allow_pickle=True)
col = lambda ch: 64 * (ch % 512 // 64) + 63 - ch % 64
def sine(y, ph):
    th = 2 * np.pi * ph / 40; X = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    c, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ c; cov = np.linalg.inv(X.T @ X) * r.var()
    amp = np.hypot(c[1], c[2]); return float(amp), float(np.sqrt((c[1] ** 2 * cov[1, 1] + c[2] ** 2 * cov[2, 2])) / amp)
OUT = {}
for bu in ('det3 June-15', 'det3 May', 'det2'):
    o = {}
    for cls, m0 in (('ypad', z['ypad']), ('xpad', ~z['ypad'])):
        m = (z['build'] == bu) & m0 & (z['px'] % 64 >= 4) & (z['px'] % 64 <= 59)
        ph = col(z['px'][m]) % 40; cp = z['ax'][m][:, 5] - z['ax'][m][:, 3]; ys = z['ay'][m][:, 3] + z['ay'][m][:, 5]
        b = ph // 4
        o[cls] = dict(n=int(m.sum()), copy=[[float(cp[b == k].mean()), float(cp[b == k].std() / np.sqrt((b == k).sum()))] for k in range(10)],
                      ysh=[[float(ys[b == k].mean()), float(ys[b == k].std() / np.sqrt((b == k).sum()))] for k in range(10)],
                      sine_copy=sine(cp, ph), sine_ysh=sine(ys, ph), mean=float(cp.mean()))
    OUT[bu] = o
    print(bu, {c: (o[c]['n'], round(o[c]['mean'], 3), np.round(o[c]['sine_copy'], 3), np.round(o[c]['sine_ysh'], 3)) for c in o})
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'esl_bins.json'), 'w'))
