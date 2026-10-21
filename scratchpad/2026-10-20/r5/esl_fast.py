"""Fast strip-registration model (2026-10-21): can a cross-gap leak or gap funnelling raise the sim's a1(delta) from
~1/3 of the data to the data's size?  Static X-line weighting at the ESL plane (bus_par/d20_nobus_esl.npy, the
fwd_xpos.py field). Deterministic: gaussian clouds of width sig on a grid of centres (cx over one column, cy over the
2P comb period). Locally pitch-locked strips (780 um pitch, the 550/800 fill) at offset s from the column centres.
Charge relocation at peak time:
  strip charge   spread uniformly across its strip's width (fast across-width equalisation)
  gap charge     stays where it lands (the sim: sigma_gap = 0) ... or a fraction f is funnelled onto the nearer strip
  leak eps       a fraction eps of every strip's charge moves to its two neighbour strips (eps/2 each)
  sy             optional along-strip (y) gaussian smear of strip charge [um]
a1 = mean(alpha(k-1) - alpha(k+1)) over clouds whose X peak is line 0 (k-1 = data's DAQ n+1, the fwd_xpos.py
convention). s > 0 puts the strip on the k+1 side of the column centre.
Usage: esl_fast.py"""
import os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, '..', '..', '2026-10-19', 'bus_par')
P = 780.; dg = 20.; EX = np.load(f'{B}/d20_nobus_esl.npy').astype(float)
nx, ny = EX.shape; nrow = int(round(ny * dg / P / 2))
X0 = P / 2; WS = 550. * P / 800.
# strip-averaged map: unit charge spread over a strip centred at x (boxcar of WS along x, 20 um cells)
nb = int(round(WS / dg)); EXs = np.zeros_like(EX); pad = np.concatenate([np.zeros((nb, ny)), EX, np.zeros((nb, ny))])
for t in range(nb): EXs += pad[nb - nb // 2 + t: nb - nb // 2 + t + nx]
EXs /= nb
def ysmear(M, sy):
    if sy <= 0: return M
    k = np.arange(ny); k = np.minimum(k, ny - k) * dg; w = np.exp(-.5 * (k / sy) ** 2); w /= w.sum()
    return np.real(np.fft.ifft(np.fft.fft(M, axis=1) * np.fft.fft(w)[None], axis=1))
def lut(M):
    Mp = np.concatenate([M, np.zeros((1, ny))])
    def f(x, y):
        i = np.floor((x + P / 2) / dg).astype(int); i = np.where((i >= 0) & (i < nx), i, nx)
        return Mp[i, (np.floor((y + nrow * P) / dg).astype(int)) % ny]
    return f
fW = lut(EX)
K = np.arange(-4, 5); NK = len(K); i0 = 4
g = np.linspace(-3, 3, 25); GX, GY = np.meshgrid(g, g, indexing='ij'); Wg = np.exp(-.5 * (GX ** 2 + GY ** 2)); Wg /= Wg.sum()
CX = np.arange(-P / 2, P / 2, 20.) + 10; CY = np.arange(-P, P, 60.) + 30
SIG = np.array([250., 350., 450., 550., 650.])
cx, cy, sg = [a.ravel() for a in np.meshgrid(CX, CY, SIG, indexing='ij')]; n = len(cx)
px = (X0 + cx[:, None, None] + sg[:, None, None] * GX[None]).reshape(n, -1)
py = (cy[:, None, None] + sg[:, None, None] * GY[None]).reshape(n, -1); wq = Wg.ravel()
def signals(s, f=0., eps=0., sy=0., prompt=False):
    if prompt: return np.stack([(fW(px - k * P, py - k * P) * wq).sum(1) for k in K], 1)
    fS = lut(ysmear(EXs, sy))
    u = px - X0 - s; j = np.round(u / P); on = np.abs(u - j * P) < WS / 2; xs = X0 + s + j * P
    w_gap = np.where(on, 0., wq * (1 - f)); w_str = np.where(on, wq, wq * f)
    S = np.zeros((n, NK))
    for ik, k in enumerate(K):
        st = fS(xs - k * P, py - k * P)
        if eps: st = (1 - eps) * st + eps / 2 * (fS(xs - P - k * P, py - k * P) + fS(xs + P - k * P, py - k * P))
        S[:, ik] = (w_gap * fW(px - k * P, py - k * P) + w_str * st).sum(1)
    return S
def score(S):
    pk = S.argmax(1) == i0; c = S[:, i0]
    dd = (S[:, i0 - 1] - S[:, i0 + 1]) / c; sh = (S[:, i0 - 1] + S[:, i0 + 1]) / c
    m = pk; nar = m & (sg <= 350); wid = m & (sg >= 550)
    return float(dd[m].mean()), float(dd[nar].mean()), float(dd[wid].mean()), float(sh[m].mean())
SS = [-325., -225., -150., -75., 0., 75., 150., 225., 325.]
MODELS = [('prompt (no ESL move)', dict(prompt=True)), ('sim-like', dict()), ('sim-like sy400', dict(sy=400.)),
          ('sim-like sy1500', dict(sy=1500.)), ('funnel 0.5', dict(f=.5)), ('funnel 1.0', dict(f=1.)),
          ('leak 0.1', dict(eps=.1)), ('leak 0.3', dict(eps=.3)),
          ('funnel 1 + leak 0.1', dict(f=1., eps=.1)), ('funnel 1 + leak 0.3', dict(f=1., eps=.3)),
          ('funnel 1 sy1500', dict(f=1., sy=1500.))]
if __name__ == '__main__':
    OUT = {}
    print(f'{n} clouds; columns: a1 / narrow (sig<=350) / wide (sig>=550) / sharing, per strip offset s (um)')
    print(f'{"":22s}     ' + ' '.join(f'{s:+6.0f}' for s in SS))
    for nm, kw in MODELS:
        rows = [score(signals(s, **kw)) for s in SS]
        OUT[nm] = dict(s=SS, a1=[r[0] for r in rows], narrow=[r[1] for r in rows], wide=[r[2] for r in rows], share=[r[3] for r in rows])
        for lab, q in (('a1', 0), ('nar', 1), ('wide', 2), ('sh', 3)):
            print(f'{nm if q == 0 else "":22s} {lab:4s} ' + ' '.join(f'{r[q]:+.3f}' for r in rows), flush=True)
    json.dump(OUT, open(os.path.join(HERE, 'esl_fast.json'), 'w'), indent=1)
