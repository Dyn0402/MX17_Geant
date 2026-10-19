"""Offline emulation of laterally displaced photon feedback (2026-10-19).

Hypothesis: a UV photon from an avalanche at height ~0 frees a photoelectron at the mesh (height H above the
resistive layer). Isotropic emission gives P(R > r) = H / sqrt(r^2 + H^2) in the plane; its projection on one strip
axis is a Cauchy (Lorentzian) of half-width H. The secondaries avalanche promptly, so their charge has the primary's
time structure. Stage B places feedback where the primary landed (digitize.add_feedback). Emulation on the sim's mean
aligned waveforms M_k(t) (k = strip offset, central-normalised per event as nb_side.aligned):
    M'_k = (1 - f) M_k + f sum_j w_j M_{k-j},   w_j = P(displacement moves the charge j strips),
with the parent uniform within its strip and f = the feedback fraction of the total charge (bG of the sim point).
Strips beyond ±5 are taken as 0. Then compares the early-window prompt alpha_k, the peak delays and the peak shares
with the data (X + side; Y folded).
Usage: nb_fb_emul.py DATA.npz SIM.npz F H1,H2,...   (H in um)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import load, HW
from nb_side import aligned, grid, DT
P = 780.; win = (grid * DT >= -200) & (grid * DT <= 300)
def wj(H, jmax=10):
    x = (np.arange(400) + .5) / 400 * P - P / 2; w = np.zeros(2 * jmax + 1)
    for xi in x:                                                   # Cauchy CDF between strip edges
        e = (np.arange(-jmax, jmax + 2) - .5) * P - xi
        cdf = .5 + np.arctan(e / H) / np.pi; w += np.diff(cdf)
    return w / len(x)
def means(f, v):
    A = aligned(load(f, v)['W']); cov = np.isfinite(A[:, HW]).mean(0) >= .5
    M = np.nanmean(A, 0); M[:, ~cov] = np.nan; return M
def summary(M, side):
    """side: '+' (X data: clean side only) or 'fold'. Returns per k = 1..5: alpha, peak delay, peak share."""
    c = M[HW]; tc = grid[np.nanargmax(c)]; out = []
    for k in range(1, 6):
        y = M[HW + k] if side == '+' else .5 * (M[HW + k] + M[HW - k])
        a = np.nansum(y[win] * c[win]) / np.nansum(c[win] ** 2)
        out.append((a, (grid[np.nanargmax(y)] - tc) * DT, np.nanmax(y) / np.nanmax(c)))
    return out
data, sim, F = sys.argv[1], sys.argv[2], float(sys.argv[3]); Hs = [float(h) for h in sys.argv[4].split(',')]
for v in 'xy':
    D = summary(means(data, v), '+' if v == 'x' else 'fold'); M = means(sim, v); M0 = np.nan_to_num(M)
    rows = [('data', D), ('sim', summary(M, 'fold'))]
    for H in Hs:
        w = wj(H); jm = (len(w) - 1) // 2; Mp = np.zeros_like(M0)
        for k in range(-HW, HW + 1):
            Mp[k + HW] = (1 - F) * M0[k + HW] + F * sum(w[j + jm] * M0[k - j + HW] for j in range(-jm, jm + 1) if abs(k - j) <= HW)
        Mp[:, np.isnan(M[HW])] = np.nan
        rows.append((f'sim+fb H={H:.0f}', summary(Mp, 'fold')))
    print(f'\n{os.path.basename(sim)} {v.upper()}, f = {F}  (data {os.path.basename(data)}: {"+ side" if v == "x" else "folded"})')
    print('   ' + ' ' * 16 + '  '.join(f'   ±{k}: alpha   delay  share ' for k in range(1, 4)) + '    ±4 alpha  ±5 alpha')
    for nm, S in rows:
        print(f'   {nm:16s}' + '  '.join(f'{S[k][0]:8.4f} {S[k][1]:+6.0f} {S[k][2]:6.3f}' for k in range(3)) + f'   {S[3][0]:8.4f}  {S[4][0]:8.4f}')
