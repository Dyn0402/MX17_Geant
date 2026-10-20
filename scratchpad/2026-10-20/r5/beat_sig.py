"""Single-frequency significance of the stripe beat (2026-10-20), from periodogram.py's inputs.
For o with std s over N events, a Fourier amplitude 2|sum o e^{i phi}|/N is Rayleigh with scale s*sqrt(2/N) under no
signal; z = A / (s sqrt(2/N)). Reports z and phase at f = 1/40 (P = 800 um) and the best f in 0.015..0.040
(P ~ 792..812 um), for s1 (sharing) and d1 (asymmetry), per build, X view (and Y view for a control: the Y view should
not carry the X-column beat). Both observables are detrended in column first (detrend(), sigma 30 columns) so the smooth gain/width
variation over the board does not leak into f ~ 0.02. Usage: beat_sig.py"""
import json, os, numpy as np
import r5lib as R
F = np.arange(0.015, 0.0401, 0.00025); OUT = {}
def detrend(col, o, sig=30.):
    # subtract a smooth trend in column (Gaussian, sigma 30 columns: transfer exp(-2 pi^2 sig^2 f^2) ~ 1e-5 at f = 1/40)
    n = np.bincount(col, minlength=512).astype(float); s = np.bincount(col, o, minlength=512)
    k = np.exp(-0.5 * (np.arange(-120, 121) / sig) ** 2)
    tr = np.convolve(s, k, 'same') / np.maximum(np.convolve(n, k, 'same'), 1e-9)
    return o - tr[col]
def amp(col, o, f):
    o = o - o.mean(); c = (o * np.exp(2j * np.pi * f * col)).sum() * 2 / len(o); return np.abs(c), np.angle(c)
for b in ('det3 June-15', 'det3 May', 'det2', 'det1', 'det4', 'det6', 'det7'):
    for view in ('X', 'Y'):
        p, s1, d1 = [], [], []
        for r in R.rows():
            if r['build'] != b or r['view'] != view: continue
            d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
            p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
        if not p: continue
        p, s1, d1 = map(np.concatenate, (p, s1, d1)); col = 64 * (p // 64) + 63 - p % 64
        m = (np.abs(s1) < 3) & (np.abs(d1) < 3); col, s1, d1 = col[m], s1[m], d1[m]; N = len(col)
        line = f'{b:13s} {view} n={N:6d}'
        for nm, v in (('s1', detrend(col, s1)), ('d1', detrend(col, d1))):
            sc = v.std() * np.sqrt(2 / N); A0, ph0 = amp(col, v, 1 / 40)
            As = np.array([amp(col, v, f)[0] for f in F]); i = As.argmax()
            OUT[f'{b}|{view}|{nm}'] = dict(z40=float(A0 / sc), A40=float(A0), ph40=float(ph0), fbest=float(F[i]), zbest=float(As[i] / sc), scale=float(sc))
            line += f' | {nm}: 1/40 A={A0:.4f} z={A0/sc:4.1f} ph={np.degrees(ph0):+4.0f} ; best f={F[i]:.4f} (P={780/(1-F[i]):.1f}um) z={As[i]/sc:4.1f}'
        print(line)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'beat_sig.json'), 'w'), indent=1)
