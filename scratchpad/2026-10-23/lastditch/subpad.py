"""Test 2 (2026-10-23): a software collimated scan. M3 points to the DUT to 0.21 mm in X (nTof_x17 M3_SELF_RESOLUTION),
so the sub-pad phase phi = M3 x offset from the nearest X-channel centre (pads at pos = 0.78 x column), in units of the
0.78 mm pitch, keeps ~25 % contrast at the pad fundamental. Per board: an 8 x 8 cell map of a1x (by M3 position) splits
the events into a strong-field region and a reference region (det3: cells with mean a1 below / above the board median;
det2: negative / positive cells). In each region, a1x, s1x and the centroid residual vs phi (8 bins).
The difference D(phi) = a1_strong(phi) - a1_ref(phi) separates the readings:
  - position-independent transfer (n-1 -> n+1 skip, skewed response): D flat in phi;
  - a displacement of the whole cloud by d: D = d x da1_ref/dphi (the derivative shape, largest where a1 is steep);
  - an attractor / a local feature: D concentrated at some phi, with s1 changing there.
Fit D(phi) = c0 + k * da1_ref/dphi (k in units of pitch); report c0, k (as um), chi2 of flat / derivative / both.
The M3 point is moved to the cluster plane (x += L tan_x, L fitted per board: ~10 mm), and only |tan_x| < TAN (0.03)
is kept: the max channel follows the track's projection, so inclined tracks smear the phase.
Usage: [TAN=0.03] subpad.py -> subpad.json"""
import os, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); P = 0.78; NB = 8; TAN = float(os.environ.get('TAN', 0.03)); OUT = {}
T = {os.path.basename(f)[:-4]: dict(np.load(f)) for f in sorted(glob.glob(f'{HERE}/tables/*.npz'))}

def cellmap(t, m, nc=8):
    ix = np.clip((t['ref_x'] / 399.36 * nc).astype(int), 0, nc - 1); iy = np.clip((t['ref_y'] / 399.36 * nc).astype(int), 0, nc - 1)
    cid = ix * nc + iy; n = np.bincount(cid[m], minlength=nc * nc); s = np.bincount(cid[m], t['a1x'][m], minlength=nc * nc)
    return cid, np.where(n >= 15, s / np.maximum(n, 1), np.nan), n

def prof(phi, v, m):
    b = np.clip(((phi[m] + 0.5) * NB).astype(int), 0, NB - 1); n = np.bincount(b, minlength=NB)
    mu = np.bincount(b, v[m], minlength=NB) / np.maximum(n, 1); sq = np.bincount(b, v[m] ** 2, minlength=NB) / np.maximum(n, 1)
    return mu, np.sqrt(np.maximum(sq - mu ** 2, 0) / np.maximum(n - 1, 1)), n

for board, tags in (('det3 June-15', ['d3_622', 'd3_627s', 'd3_627']), ('det2', ['d2_622', 'd2_622l']), ('det4', ['d4_624'])):
    t = {k: np.concatenate([T[g][k] for g in tags]) for k in T[tags[0]] if k != 'build'}
    sess = np.concatenate([np.full(len(T[g]['ev']), i) for i, g in enumerate(tags)])
    r0 = t['xcen'] - t['ref_x']; m0 = (np.abs(t['tan_x']) < 0.1) & (np.abs(r0) < 4)
    L = np.polyfit(t['tan_x'][m0], r0[m0], 1)[0]; t['ref_x'] = t['ref_x'] + L * t['tan_x']   # M3 plane -> cluster plane
    ok = (np.abs(t['tan_x']) < TAN) & (np.abs(t['a1x']) < 3) & (np.abs(t['s1x']) < 3) & (np.abs(t['xcen'] - t['ref_x']) < 3)
    t['a1n'] = t['a1x'].copy()
    for i in range(len(tags)):                                            # remove each session's board mean (per-session scale)
        mm = ok & (sess == i); t['a1n'][mm] -= t['a1x'][mm].mean()
    cid, M, n = cellmap(t, ok); med = np.nanmedian(M)
    strong = ok & (M[cid] < med) if board != 'det2' else ok & (M[cid] < 0)
    ref = ok & (M[cid] >= med) if board != 'det2' else ok & (M[cid] >= 0)
    phi = ((t['ref_x'] / P + 0.5) % 1) - 0.5
    res = t['xcen'] - t['ref_x']; R = {}
    for nm, v in (('a1', t['a1x']), ('s1', t['s1x']), ('res_um', res * 1e3)):
        a, ea, na = prof(phi, v, strong); b, eb, nb = prof(phi, v, ref)
        R[nm] = dict(strong=a.tolist(), e_strong=ea.tolist(), ref=b.tolist(), e_ref=eb.tolist(), n_strong=na.tolist(), n_ref=nb.tolist())
    a, ea = np.array(R['a1']['strong']), np.array(R['a1']['e_strong']); b, eb = np.array(R['a1']['ref']), np.array(R['a1']['e_ref'])
    D, eD = a - b, np.hypot(ea, eb); ph = (np.arange(NB) + 0.5) / NB - 0.5
    der = (np.roll(b, -1) - np.roll(b, 1)) / (2 / NB)                    # da1_ref/dphi, periodic central difference
    w = 1 / eD ** 2
    def fit(X):
        cf = np.linalg.solve((X.T * w) @ X, (X.T * w) @ D); cov = np.linalg.inv((X.T * w) @ X)
        return cf, np.sqrt(np.diag(cov)), float((((D - X @ cf) / eD) ** 2).sum())
    f0, e0, c0 = fit(np.c_[np.ones(NB)]); f1, e1, c1 = fit(np.c_[der]); f2, e2, c2 = fit(np.c_[np.ones(NB), der])
    # first harmonic of D (a localized feature shows up here beyond the derivative shape)
    h = 2 * np.mean((D - D.mean()) * np.exp(-2j * np.pi * ph)); eh = np.sqrt(2 / NB) * np.sqrt(np.mean(eD ** 2)) * np.sqrt(2)
    OUT[board] = dict(lever_mm=float(L), tan_cut=TAN, res_mad_um=float(1482.6 * np.median(np.abs(res[ok] - np.median(res[ok])))), cells=M.reshape(8, 8).tolist(), median=float(med), n_strong=int(strong.sum()), n_ref=int(ref.sum()),
                      mean_a1_strong=float(t['a1x'][strong].mean()), mean_a1_ref=float(t['a1x'][ref].mean()),
                      phi=ph.tolist(), D=D.tolist(), eD=eD.tolist(), dA1ref_dphi=der.tolist(),
                      flat=dict(c0=float(f0[0]), e=float(e0[0]), chi2=c0, dof=NB - 1),
                      deriv=dict(k_um=float(f1[0] * P * 1e3), e_um=float(e1[0] * P * 1e3), chi2=c1, dof=NB - 1),
                      both=dict(c0=float(f2[0]), e_c0=float(e2[0]), k_um=float(f2[1] * P * 1e3), e_k_um=float(e2[1] * P * 1e3), chi2=c2, dof=NB - 2),
                      harm1=dict(amp=float(abs(h)), phase=float(np.angle(h)), e=float(eh)), prof=R)
    q = OUT[board]
    print(f"{board}: lever {L:.1f} mm, |tan_x| < {TAN}, centroid - M3 MAD {q['res_mad_um']:.0f} um")
    print(f"   strong n {q['n_strong']} <a1> {q['mean_a1_strong']:+.3f} | ref n {q['n_ref']} <a1> {q['mean_a1_ref']:+.3f}")
    print('   phi        ' + ' '.join(f'{x:+6.2f}' for x in ph))
    for nm in ('a1', 's1', 'res_um'):
        for reg in ('strong', 'ref'):
            print(f'   {nm:6s} {reg:6s}' + ' '.join(f'{x:+6.3f}' if nm != 'res_um' else f'{x:+6.0f}' for x in R[nm][reg]))
    print('   D            ' + ' '.join(f'{x:+6.3f}' for x in D) + '   +- ' + ' '.join(f'{x:.3f}' for x in eD))
    print(f"   flat: c0 {q['flat']['c0']:+.3f} +- {q['flat']['e']:.3f}  chi2 {q['flat']['chi2']:.1f}/{NB - 1};  derivative only: k {q['deriv']['k_um']:+.0f} +- {q['deriv']['e_um']:.0f} um  chi2 {q['deriv']['chi2']:.1f}/{NB - 1};"
          f"  both: c0 {q['both']['c0']:+.3f} +- {q['both']['e_c0']:.3f}, k {q['both']['k_um']:+.0f} +- {q['both']['e_k_um']:.0f} um  chi2 {q['both']['chi2']:.1f}/{NB - 2}", flush=True)
json.dump(OUT, open(f'{HERE}/subpad.json', 'w'), indent=1)
