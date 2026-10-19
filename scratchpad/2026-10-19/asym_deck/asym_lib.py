"""Observables for the X-asymmetry deck (2026-10-19), from raw_asym.py outputs (channel-order 9x32 windows, p-4..p+4).
good():  200 < A_c < 3500 ADC (not railed), central peak sample 4..22, flat pre-samples (|mean 0-2| < 0.1 A_c).
hot():   channels that are the peak far more often than their neighbours (> 6x the median count of busy channels).
load(tag) -> dict with W (wf), A (max per channel), p, ev, fi, c; all observables are in DAQ channel numbers:
 L = A[p-1], R = A[p+1]; alpha_k = in-phase content of channel p+k (projection on the central pulse);
 copy a1 = alpha_+1 - alpha_-1."""
import os, numpy as np
RAW = os.path.expanduser(os.environ.get('NB_RAW', '~/x17/nb_asym_raw'))
HW = 4

def good(W):
    c = W[:, HW]; A = c.max(1); pk = c.argmax(1); pre = c[:, :3].mean(1)
    return (A > 200) & (A < 3500) & (pk >= 4) & (pk <= 22) & (np.abs(pre) < 0.1 * A)

def hot(p):
    n = np.bincount(p, minlength=512); busy = n[n > 0]
    return np.flatnonzero(n > 6 * max(np.median(busy), 1)) if len(busy) else np.array([], int)

def load(path, veto=True):
    z = np.load(path if path.endswith('.npz') else f'{RAW}/{path}.npz')
    W = z['wf']; g = good(W); p = z['p'].astype(int)
    hv = hot(p[g]) if veto else np.array([], int)
    g &= ~np.isin(p, hv)
    W = W[g]; A = W.max(2); c = W[:, HW]; n2 = np.einsum('ij,ij->i', c, c)
    al = {k: np.einsum('ij,ij->i', W[:, HW + k], c) / n2 for k in range(-HW, HW + 1)}
    return dict(W=W, A=A, p=p[g], ev=z['ev'][g], fi=z['fi'][g], al=al, a1=al[1] - al[-1], hot=hv, z=z,
                n=int(g.sum()), nall=len(z['p']))

def coin(d):
    R = d['A'][:, HW + 1] > d['A'][:, HW - 1]; f = R.mean(); return f, np.sqrt(f * (1 - f) / max(len(R), 1))

def mean_se(x):
    return float(np.mean(x)), float(np.std(x) / np.sqrt(max(len(x), 1)))

def uncopy(W, a):
    """Undo channel n+1 += a * channel n inside the window (channel order), recursively upward."""
    out = W.astype(float).copy()
    for k in range(1, W.shape[1]):
        out[:, k] = W[:, k] - a * out[:, k - 1]
    return out

def fit_copy(d):
    """a that zeroes the mean alpha_+1 - alpha_-1 after uncopy (bisection)."""
    def asym(a):
        U = uncopy(d['W'], a); c = U[:, HW]; n2 = np.einsum('ij,ij->i', c, c)
        return np.mean(np.einsum('ij,ij->i', U[:, HW + 1] - U[:, HW - 1], c) / n2)
    lo, hi = -0.2, 0.5
    for _ in range(40):
        m = (lo + hi) / 2
        if asym(m) > 0: lo = m
        else: hi = m
    return (lo + hi) / 2

def profile(W):
    """Mean in-phase content alpha_k, k = -4..4, of an array of windows."""
    c = W[:, HW]; n2 = np.einsum('ij,ij->i', c, c)
    return np.array([np.mean(np.einsum('ij,ij->i', W[:, HW + k], c) / n2) for k in range(-HW, HW + 1)])
