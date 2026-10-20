"""(T) Time profile of the n+1 excess for Y-pad events (2026-10-20), data vs models.
Aligned on the Y central's peak sample, normalised per event by the Y central amplitude; mean waveforms of Y central,
X central, X(n-1) (clean side), and the excess e = X(n+1) - X(n-1). Reports peak time (parabolic interpolation, ns,
relative to the Y central peak), 50 % rise time, and the fit e = a Xc + b Yc + c X(n-1) (shape decomposition).
Usage: imported: prof(Wx, Wy, ypad) -> dict."""
import numpy as np
def _align(W, ref_pk, Z=8, T=40):
    out = np.zeros((len(W), W.shape[1], T))
    for i in range(len(W)):
        s = Z - ref_pk[i]; lo, hi = max(0, -s), min(32, T - s)
        out[i, :, lo + s:hi + s] = W[i, :, lo:hi]
    return out
def _pk(w):
    k = int(np.argmax(w)); 
    if 0 < k < len(w) - 1:
        a, b, c = w[k - 1], w[k], w[k + 1]; d = (a - 2 * b + c)
        return (k + (0.5 * (a - c) / d if d != 0 else 0)) * 60.0
    return k * 60.0
def _t50(w):
    m = w.max(); i = int(np.flatnonzero(w >= 0.5 * m)[0])
    if i == 0: return 0.
    return (i - 1 + (0.5 * m - w[i - 1]) / (w[i] - w[i - 1])) * 60.0
def prof(Wx, Wy, ypad, Z=8):
    Wx, Wy = Wx[ypad].astype(float), Wy[ypad].astype(float)
    pk = Wy[:, 4].argmax(1); A = Wy[:, 4].max(1)[:, None, None]
    X = _align(Wx / A, pk, Z); Y = _align(Wy / A, pk, Z)
    yc, xc, xm, xp = Y[:, 4].mean(0), X[:, 4].mean(0), X[:, 3].mean(0), X[:, 5].mean(0); e = xp - xm
    o = dict(n=int(len(Wx)))
    for nm, w in (('Yc', yc), ('Xc', xc), ('Xm1', xm), ('exc', e)):
        o[nm] = dict(peak=_pk(w) - Z * 60, t50=_t50(w) - Z * 60, amp=float(w.max()))
    B = np.stack([xc, yc, xm], 1); c, *_ = np.linalg.lstsq(B, e, rcond=None); r = e - B @ c
    o['fit'] = dict(a_Xc=float(c[0]), b_Yc=float(c[1]), c_Xm1=float(c[2]), resid=float(r.std() / e.std()))
    o['w'] = dict(yc=yc.tolist(), xc=xc.tolist(), xm=xm.tolist(), exc=e.tolist())
    return o
def line(name, o):
    f = lambda k: f'{o[k]["peak"]:+5.0f}/{o[k]["t50"]:+5.0f}'
    F = o['fit']
    return (f'{name:22s} n={o["n"]:6d} | peak/t50 ns: Yc {f("Yc")}  Xc {f("Xc")}  X(n-1) {f("Xm1")}  excess {f("exc")} (amp {o["exc"]["amp"]:.3f})'
            f' | e = {F["a_Xc"]:+.2f} Xc {F["b_Yc"]:+.2f} Yc {F["c_Xm1"]:+.2f} Xn-1 (resid {F["resid"]:.2f})')
