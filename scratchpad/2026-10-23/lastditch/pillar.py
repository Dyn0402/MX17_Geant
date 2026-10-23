"""Test 1 (2026-10-23): is the pillar lattice (Dynamask, dia 0.6 mm on a 4.68 mm = 6-pad grid, 3498A_bulk.gbr) visible in
the June boards' data, and does a1 follow it?  The 10-19 "flat mod 6" check averaged column mod 6 over the whole board,
which cancels a lattice whose phase drifts against the pads (a moire), so it never tested this.
Hypothesis: a pillar film laid out on the pad grid but printed ~0.2 % off scale / a few mrad rotated would make a smooth
two-signed X-only a1 field (beat 0.78 mm / |s - 1| ~ 400 mm for s - 1 = 0.002, the map fit's ~500-column period).
Data:
  production wft events.parquet (all June bench runs, x_ok & y_ok & quality, <= 1 track): positions x_p0, y_p0 (raw strip
    frame, sigma ~0.55 / 0.69 mm), weights: occupancy, log charge (x_q_sum + y_q_sum, centred and clipped);
  raw survey matched X/Y events (r5lib.ev): central pads -> pos = 0.78 x physical column/row, weight a1x (DAQ order).
Statistic: 2-D Fourier power at the two lattice fundamentals k1 = (cos t, sin t) / (4.68 s), k2 = (-sin t, cos t) / (4.68 s),
Hann-windowed over the active area, normalised by the null expectation sum(h^2 w^2) (Z ~ Exp(1) per vector; Z1 + Z2 ~
Gamma(2)); occupancy and charge coherent over the board (the lattice phase is global), a1 incoherent over 4 x 4 patches (its
sign follows the local pillar/pad phase). Search s in [0.975, 1.025], t in [-1.5, 1.5] deg; the look-elsewhere null is the
same search at off-lattice pitches (4.68 x {0.86, 0.9, 1.1, 1.14}, between the 49.92 mm chip harmonics).
Sensitivity: the same search after injecting synthetic pillars into the real events (a lattice at s = 1.002, t = 0.2 deg;
events whose true position, unsmeared estimate = measured, falls within 0.3 mm of a pillar lose a fraction `depth` of
their weight / are dropped with probability `depth`), reported as the depth that reaches Z_max at the null's 99th percentile.
Usage: pillar.py -> pillar.json"""
import os, sys, glob, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); AN = '/media/dylan/data/x17/cosmic_bench/Analysis'
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5')); import r5lib as R
PIT, L = 4.68, 399.36; S = np.arange(0.975, 1.02501, 0.0025); TH = np.radians(np.arange(-1.5, 1.501, 0.1))
NULLP = [0.86, 0.90, 1.10, 1.14]; rng = np.random.default_rng(23); OUT = {}
DETS = {'det3 June-15': 'mx17_3', 'det2': 'mx17_2', 'det4': 'mx17_4', 'det6': 'mx17_6', 'det7': 'mx17_7'}

def hann(x, y): return np.sin(np.pi * np.clip(x / L, 0, 1)) ** 2 * np.sin(np.pi * np.clip(y / L, 0, 1)) ** 2

def zgrid(x, y, w, pitch, npatch=1):
    """Z1 + Z2 over the (s, t) grid; w already centred (occupancy: w = 1)."""
    Z = np.zeros((len(S), len(TH)))
    edges = np.linspace(0, L, npatch + 1)
    for a in range(npatch):
        for b in range(npatch):
            m = (x >= edges[a]) & (x < edges[a + 1]) & (y >= edges[b]) & (y < edges[b + 1])
            if m.sum() < 50: continue
            xs, ys = x[m] - edges[a], y[m] - edges[b]; lp = L / npatch
            h = (np.sin(np.pi * xs / lp) * np.sin(np.pi * ys / lp)) ** 2; hw = h * w[m]; norm = (hw ** 2).sum()
            for i, s in enumerate(S):
                for kx, ky in ((np.cos(TH), np.sin(TH)), (-np.sin(TH), np.cos(TH))):
                    ph = (2j * np.pi / (pitch * s)) * (np.outer(xs, kx) + np.outer(ys, ky))
                    Z[i] += np.abs(hw @ np.exp(ph)) ** 2 / norm
    return Z

def search(x, y, w, npatch=1):
    Z = zgrid(x, y, w, PIT, npatch); i, j = np.unravel_index(Z.argmax(), Z.shape)
    nul = [zgrid(x, y, w, PIT * f, npatch).max() for f in NULLP]
    return dict(zmax=float(Z.max()), s=float(S[i]), t_deg=float(np.degrees(TH[j])), null_max=[float(v) for v in nul],
                z_at_nominal=float(Z[np.argmin(abs(S - 1)), np.argmin(abs(TH))])), Z

def pillars_near(x, y, s=1.002, t=np.radians(0.2), r=0.3, x0=1.1, y0=2.3):
    """distance to the nearest point of a synthetic lattice (pitch PIT * s, rotation t, origin x0, y0)"""
    c, sn = np.cos(t), np.sin(t); u = ((x - x0) * c + (y - y0) * sn) / (PIT * s); v = (-(x - x0) * sn + (y - y0) * c) / (PIT * s)
    du, dv = (u - np.round(u)) * PIT * s, (v - np.round(v)) * PIT * s
    return np.hypot(du, dv) < r

def prod(det):
    fs = sorted(glob.glob(f'{AN}/*/*/{det}/wft/events.parquet')); D = []
    for f in fs:
        P = pd.read_parquet(f, columns=['x_ok', 'y_ok', 'x_quality_ok', 'y_quality_ok', 'n_tracks', 'x_p0', 'y_p0', 'x_q_sum', 'y_q_sum', 'spark'])
        P = P[P.x_ok & P.y_ok & P.x_quality_ok & P.y_quality_ok & (P.n_tracks <= 1) & ~P.spark]
        P = P[(P.x_p0 > 2) & (P.x_p0 < L - 2) & (P.y_p0 > 2) & (P.y_p0 < L - 2) & (P.x_q_sum > 0) & (P.y_q_sum > 0)]
        D.append(P[['x_p0', 'y_p0', 'x_q_sum', 'y_q_sum']].to_numpy(float))
    return np.concatenate(D), len(fs)

for b, det in (DETS.items() if __name__ == '__main__' else ()):   # importable by pillar_inject.py
    d, nf = prod(det); x, y = d[:, 0], d[:, 1]; lq = np.log(d[:, 2] + d[:, 3]); lq = np.clip(lq - np.median(lq), -2, 2)
    res = dict(n=len(x), runs=nf)
    occ, Zo = search(x, y, np.ones(len(x))); chg, Zc = search(x, y, lq - lq.mean())
    res['occupancy'], res['charge'] = occ, chg
    # sensitivity: inject pillars into the real events (occupancy: drop with prob depth; charge: scale by 1 - depth)
    near = pillars_near(x, y); res['inject_frac_near'] = float(near.mean()); inj = {}
    for depth in (0.25, 0.5, 1.0):
        keep = ~(near & (rng.random(len(x)) < depth)); zo = zgrid(x[keep], y[keep], np.ones(keep.sum()), PIT).max()
        lq2 = lq + np.where(near, np.log(max(1 - depth, 0.05)), 0); lq2 = np.clip(lq2, -3, 2); zc = zgrid(x, y, lq2 - lq2.mean(), PIT).max()
        inj[str(depth)] = dict(occ=float(zo), charge=float(zc))
    res['inject'] = inj; OUT[b] = res
    print(f"{b:13s} production n {len(x):6d} ({nf} runs): occupancy Zmax {occ['zmax']:5.1f} at s {occ['s']:.4f} t {occ['t_deg']:+.1f}  "
          f"(nominal {occ['z_at_nominal']:.1f}; off-lattice maxima {', '.join(f'{v:.1f}' for v in occ['null_max'])}); "
          f"charge Zmax {chg['zmax']:5.1f} at s {chg['s']:.4f} t {chg['t_deg']:+.1f} (off-lattice {', '.join(f'{v:.1f}' for v in chg['null_max'])})", flush=True)
    print(f"{'':13s} injected pillars ({100 * res['inject_frac_near']:.1f} % of events within 0.3 mm): " +
          '; '.join(f"depth {k}: occ Z {v['occ']:.1f}, charge Z {v['charge']:.1f}" for k, v in inj.items()), flush=True)

# a1 on the raw survey (matched X/Y, central pads)
for b in (('det3 June-15', 'det2', 'det3 May') if __name__ == '__main__' else ()):
    E = R.ev(R.BUILDS[b]); pc = lambda p: 64 * (p // 64) + 63 - p % 64
    x, y = 0.78 * pc(E['px']), 0.78 * pc(E['py']); a1 = E['alx'][:, 5] - E['alx'][:, 3]; s1 = E['alx'][:, 5] + E['alx'][:, 3]
    m = (np.abs(a1) < 3) & (np.abs(s1) < 3); x, y, a1, s1 = x[m], y[m], a1[m], s1[m]; r = dict(n=int(m.sum()))
    for nm, w, npch in (('a1_patch4', a1 - a1.mean(), 4), ('a1_board', a1 - a1.mean(), 1), ('s1_board', s1 - s1.mean(), 1), ('occ_pads', np.ones(len(x)), 1)):
        q, _ = search(x, y, w, npch); r[nm] = q
    near = pillars_near(x, y); r['inject_frac_near'] = float(near.mean())
    a1i = a1 + np.where(near, 0.5, 0); r['inject_a1_0.5_patch4'] = float(zgrid(x, y, a1i - a1i.mean(), PIT, 4).max())
    OUT[f'survey {b}'] = r
    print(f"survey {b:13s} n {r['n']:6d}: " + '; '.join(f"{nm} Zmax {r[nm]['zmax']:.1f} at s {r[nm]['s']:.4f} t {r[nm]['t_deg']:+.1f} (off-lattice {max(r[nm]['null_max']):.1f})"
                                                   for nm in ('a1_patch4', 'a1_board', 's1_board', 'occ_pads')) +
          f"; injected a1 +0.5 near pillars ({100 * r['inject_frac_near']:.1f} %): Z {r['inject_a1_0.5_patch4']:.1f}", flush=True)
if __name__ == '__main__': json.dump(OUT, open(f'{HERE}/pillar.json', 'w'), indent=1)
