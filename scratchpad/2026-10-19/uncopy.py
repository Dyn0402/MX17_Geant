"""Undo det3's one-sided copy (2026-10-19, NB_ASYM.md): channel n+1 carries a * channel n (same chip, same sample),
so x'_{n+1} = x_{n+1} - a x'_n, recursively upward in channel number within each 64-channel chip.
uncopy_feu(W, a): W[512, n_sample] in CHANNEL order -> corrected copy.
uncopy_window(W, ch, a): a position-ordered window (k = -5..5, ch = channel ids) as in the tagged dumps; with
connector-inverted wiring the channel predecessor of position k is position k+1, so it runs from the top down.
Run as a script: closure test on a tagged dump (copy before/after, + side unchanged).  Usage: uncopy.py NPZ A"""
import sys, numpy as np
def uncopy_feu(W, a):
    out = W.astype(float).copy()
    for n in range(1, W.shape[0]):
        if n % 64: out[n] = W[n] - a * out[n - 1]
    return out
def uncopy_window(W, ch, a):
    """W[n_ev, 2HW+1, ns], ch[n_ev, 2HW+1]; position k's channel predecessor is ch - 1, normally at position k + 1."""
    out = W.astype(float).copy(); hw = W.shape[1] // 2
    for k in range(2 * hw - 1, -1, -1):
        prev = (ch[:, k + 1] == ch[:, k] - 1) & (ch[:, k] % 64 != 0)
        out[prev, k] = W[prev, k] - a * out[prev, k + 1]
    return out
if __name__ == '__main__':
    HW = 5; z = np.load(sys.argv[1]); a = float(sys.argv[2])
    for v in 'xy':
        W = np.nan_to_num(z[f'{v}_wf']); ch = z[f'{v}_ch'].astype(int); c = W[:, HW].max(1); m = (c > 200) & (c < 3500)
        for nm, X in (('raw', W[m]), (f'uncopy a={a}', uncopy_window(W[m], ch[m], a))):
            wc = X[:, HW]; n2 = np.einsum('ij,ij->i', wc, wc); al = lambda k: np.einsum('ij,ij->i', X[:, HW + k], wc) / n2
            print(f'{v.upper()} {nm:14s}: alpha-1 {al(-1).mean():.4f} alpha+1 {al(1).mean():.4f} copy {np.mean(al(-1) - al(1)):+.4f}  '
                  f'alpha-2 {al(-2).mean():.4f} alpha+2 {al(2).mean():.4f}')
