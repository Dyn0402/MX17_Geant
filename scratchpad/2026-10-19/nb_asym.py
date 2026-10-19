"""Left/right neighbour asymmetry on the tagged dumps (dump_tagged.py), 2026-10-19.
Per event: n_k = max of strip k (position order) / max of the peak strip; A1 = n_-1 - n_+1, A2 = n_-2 - n_+2.
Data: |ax|, |ay| < AMAX (M3); sim: every event. Railed peak strips (>= 3500) and weak ones (< 200 ADC) dropped.
Splits: beat phase (rank mod 40, the 39 ESL : 40 pad superperiod), channel within the 64-channel connector, connector,
rank parity, and the track's position within the peak strip u (M3 or truth, offset fitted per file so that the
mean of u is 0 over events), in [-0.5, 0.5] strip pitches.
Usage: nb_asym.py NPZ...   (AMAX env, default 3)"""
import os, sys, numpy as np
HW = 5; PITCH = 0.78; AMAX = float(os.environ.get('AMAX', 3))

def load(f, v):
    z = np.load(f); W = z[f'{v}_wf']; m3 = z[f'{v}_m3']; amp = np.nanmax(W, 2); c = amp[:, HW]
    sim = np.all(m3[:, 2:] == 0)
    sel = (c < 3500) & (c > 200) & np.isfinite(amp[:, HW - 2:HW + 3]).all(1)
    if not sim: sel &= (np.abs(m3[:, 2]) < AMAX) & (np.abs(m3[:, 3]) < AMAX)
    n = amp[sel] / c[sel, None]
    pos = z[f'{v}_pos'][z[f'{v}_peak_ch'][sel]]
    # the coordinate along the view's strips: data det3 is rotated 90 deg to M3 (X strips ~ ref_y, Y strips ~ -ref_x)
    if sim: t = m3[sel, 0] if v == 'x' else m3[sel, 1]
    else:   t = m3[sel, 1] if v == 'x' else -m3[sel, 0]
    # sign: pick the one that correlates positively with the strip position
    if np.corrcoef(pos, t)[0, 1] < 0: t = -t
    d = pos - t; off = np.median(d); u = -(pos - t - off) / PITCH                      # u > 0: track toward +1
    return dict(n=n, rank=z[f'{v}_rank'][sel].astype(int), ch=z[f'{v}_peak_ch'][sel].astype(int), u=u,
                W=W[sel] / c[sel, None, None], sim=sim, off=off)

def mean_err(x): return (np.nanmean(x), np.nanstd(x) / np.sqrt(max(np.isfinite(x).sum(), 1)))

def table(name, key, A, labels):
    out = []
    for k in labels:
        m = key == k
        if m.sum() < 5: out.append(f'{k}:  –  '); continue
        a, e = mean_err(A[m]); out.append(f'{k}:{a:+.3f}±{e:.3f}')
    print(f'    by {name}: ' + '  '.join(out))

def sine_fit(phase, A, period):
    """A ~ a0 + a1 cos + b1 sin over the phase; returns a0, amplitude, phase of max, chi2 improvement."""
    w = 2 * np.pi * phase / period; M = np.stack([np.ones_like(w), np.cos(w), np.sin(w)], 1)
    c, res, *_ = np.linalg.lstsq(M, A, rcond=None)
    r0 = np.sum((A - A.mean()) ** 2); r1 = np.sum((A - M @ c) ** 2); s2 = r1 / (len(A) - 3)
    cov = s2 * np.linalg.inv(M.T @ M); amp = np.hypot(c[1], c[2])
    eamp = np.sqrt((c[1] ** 2 * cov[1, 1] + c[2] ** 2 * cov[2, 2]) / amp ** 2)
    return c[0], amp, eamp, (np.degrees(np.arctan2(c[2], c[1])) % 360) / 360 * period, (r0 - r1) / s2

if __name__ == '__main__':
    for f in sys.argv[1:]:
        for v in 'xy':
            D = load(f, v); n = D['n']; A1 = n[:, HW - 1] - n[:, HW + 1]; A2 = n[:, HW - 2] - n[:, HW + 2]
            print(f'\n{os.path.basename(f)} {v.upper()}  n={len(n)}  (u offset {D["off"]:+.2f} mm)')
            for k in (1, 2, 3):
                (l, el), (r, er) = mean_err(n[:, HW - k]), mean_err(n[:, HW + k])
                print(f'    ±{k}: -{k} {l:.4f}±{el:.4f}  +{k} {r:.4f}±{er:.4f}  diff {l - r:+.4f}±{np.hypot(el, er):.4f}')
            for nm, A in (('A1', A1), ('A2', A2)):
                print(f'  {nm}:')
                table('rank parity', D['rank'] % 2, A, range(2))
                table('beat phase (rank%40)//5', (D['rank'] % 40) // 5, A, range(8))
                table('channel in connector //8', (D['ch'] % 64) // 8, A, range(8))
                table('connector', D['ch'] // 64, A, range(8))
                table('u bin (-0.5..0.5 in 0.1)', np.clip(((D['u'] + .5) * 10).astype(int), 0, 9), A, range(10))
                a0, amp, eamp, ph, dchi = sine_fit(D['rank'] % 40, A, 40)
                print(f'    beat sine fit: mean {a0:+.4f}, amplitude {amp:.4f}±{eamp:.4f} (max at rank%40 = {ph:.1f}), Δχ² {dchi:.1f} for 2 dof')
            # the u-symmetric part: A(u) + A(-u) cancels the geometry and leaves the asymmetric response
            ub = np.abs(D['u']) < .25
            print(f'  A1 for central tracks |u| < 0.25: {mean_err(A1[ub])[0]:+.4f}±{mean_err(A1[ub])[1]:.4f}  (n {ub.sum()}),'
                  f' mean u {D["u"][ub].mean():+.3f}')
