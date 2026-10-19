"""Neighbour delay vs strip distance, X and Y separately, each side separately (-k, +k) and folded (±k).
Same alignment as 2026-10-18/nb_compare.py (t = 0 at the central's leading 50 %, each event normalised to its
central peak, railed centrals dropped), on the ±5 dumps of dump_pmN_m3.py. Per strip distance: peak time of the
mean waveform minus the central's, the mean's peak and minimum (fraction of central, after a pre-signal baseline) and its 50 % leading-edge time minus the
central's. Errors: bootstrap over events (B resamples). Strips past the board edge are NaN and skipped (nanmean).
A peak still rising at the end of the record prints as inf; '(x% late)' is the bootstrap fraction where it was.
The 50 % edge crossing is the last one before the peak (for a far strip, how late its rise reaches half its peak).
Usage: nb_far.py NPZ...   (B env, default 300)"""
import os, sys, numpy as np
DT = 60.; grid = np.arange(-8, 24, .25); B = int(os.environ.get('B', 300))
def cross(w, frac):
    ipk = int(np.nanargmax(w)); lvl = frac * w[ipk]
    b = np.flatnonzero(~(w[:ipk + 1] >= lvl))
    if not len(b): return np.nan
    i = b[-1]; j = i + 1
    if j >= len(w): return np.nan
    return i + (lvl - w[i]) / (w[j] - w[i]) if w[j] != w[i] else float(i)
def aligned(W):
    hw = W.shape[1] // 2; W = W[np.nanmax(W[:, hw], 1) < 3500]; x = np.arange(W.shape[2]); out = []
    for w in W:
        l = cross(w[hw], .5)
        if np.isnan(l): continue
        pk = w[hw].max()
        out.append([np.interp(grid + l, x, s / pk, left=np.nan, right=np.nan) if np.isfinite(s).all() else np.full(len(grid), np.nan) for s in w])
    return np.array(out)                                    # events x strips x grid
def stats(A, hw, idx):
    """A: events x strips x grid; idx: strip rows to average (one side or both). Returns delay, peak, edge delay."""
    cov = np.isfinite(A[:, hw]).mean(0) >= .5                  # grid points inside the record for most events
    c = np.where(cov, np.nanmean(A[:, hw], 0), np.nan); n = np.where(cov, np.nanmean(np.nanmean(A[:, idx], 1), 0), np.nan)
    tc = grid[np.nanargmax(c)]; tn = grid[np.nanargmax(n)]
    if np.nanargmax(n) == np.flatnonzero(cov)[-1]: tn = np.inf    # still rising at the end of the record
    ec = cross(c, .5); en = cross(n, .5)
    return (tn - tc) * DT, np.nanmax(n), np.nanmin(n), (np.interp(en, np.arange(len(grid)), grid) - np.interp(ec, np.arange(len(grid)), grid)) * DT
rng = np.random.default_rng(1)
def boot(A, hw, idx):
    v = np.array(stats(A, hw, idx)); bs = np.array([stats(A[rng.integers(0, len(A), len(A))], hw, idx) for _ in range(B)])
    fr = np.isinf(bs[:, 0]).mean(); bs[~np.isfinite(bs)] = np.nan
    return v, np.nanstd(bs, 0), fr
for f in sys.argv[1:]:
    z = np.load(f); name = f.split('/')[-2]
    for v in 'xy':
        D = aligned(z[f'data_{v}_wf']); S = aligned(z[f'sim_{v}_wf']); hw = D.shape[1] // 2
        print(f'\n{name} {v.upper()}  (data {len(D)} ev, sim {len(S)} ev)   peak delay [ns] | peak/central | min/central | 50% edge delay [ns]')
        for k in range(1, hw + 1):
            for lab, idx in ((f'-{k}', [hw - k]), (f'+{k}', [hw + k]), (f'±{k}', [hw - k, hw + k])):
                (dd, dp, dm, de), (sd, sp, sm, se), fd = boot(D, hw, idx)
                (md, mp, mm, me), (smd, smp, smm, sme), fs = boot(S, hw, idx)
                print(f'  {lab:>3s}  delay data {dd:+5.0f}±{sd:3.0f} ({fd:.0%} late) sim {md:+5.0f}±{smd:3.0f} ({fs:.0%} late) | '
                      f'peak data {dp:.3f}±{sp:.3f} sim {mp:.3f}±{smp:.3f} | edge data {de:+5.0f}±{se:3.0f} sim {me:+5.0f}±{sme:3.0f}')
