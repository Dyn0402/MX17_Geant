#!/usr/bin/env python3
"""Gain ladders for every bench chamber x gas we have, one recipe (FRESH_LOOK §13).

Datasets (threshold-free peak-strip max, M3-good rays, same population at every V):
  det0  Ar/iso 95/5     4-7-26   extract_bench.py with track window; fiducial ref_x 50..150, ref_y -200..150
  det0  He/C2H6 96.5/3.5 4-2-26  same chamber, same fiducial
  det3o Ar/iso 95/5     5-5-26   old det3 (bulked 5-5-26); fiducial x -200..100, y -200..50
  det3o Ar/CF4 90/10    5-6-26   same chamber
  det3  Ar/iso 95/5     6-27-26  June det3 = T14 (nTof_x17 hv_slope mesh_ladder.csv, as in §9)
Clean-point rule: quantile >= FLOOR_MULT x the view's noise-only peak max (~25 ADC det0,
~40 ADC det3o) and <= 0.70 x rail (2700 ADC). Excluded for time drift: det0 Ar/iso 490 V
(gain +40 % within its own 3 h, first point after the He -> Ar/iso change), det3o Ar/iso 510 V
(+25 %), det0 He 505 V (taken first, ~6 h out of the descending sequence) and final_505 (separate 24 h run); det3o Ar/CF4 above 630 V (every quantile
flattens: a separate saturation regime above the fitted runaway).
Fits per (dataset, view, quantile) with >= 4 clean points, errors floored at 2 %:
  exp    ln A = a + s (V - V0)
  gfb    A = c G(V) / (1 - beta G(V)/G(V0)), G = Garfield uniform-gap shape for that gas
         (nTof_x17 garfield_sim/results/<gas>_Saclay_160m.json, log-linear fit) -> bG(V0), V*
  ffb    same with a free base slope
Writes ladders_all.csv (points) and fits_all.csv.
"""
import glob, json, numpy as np, pandas as pd
from scipy.optimize import least_squares
C = '/home/dylan/.claude/jobs/913ff27a/tmp/c/'
GR = '/home/dylan/PycharmProjects/nTof_x17/garfield_sim/results/'
RAIL = 3872.0
def garf_slope(name):
    d = json.load(open(f'{GR}{name}_Saclay_160m.json'))
    V, G = np.array(d['voltages'], float), np.array(d['gain_mean'], float)
    return np.polyfit(V, np.log(G), 1)[0]

SETS = [  # name, parquet glob, fiducial, gas json, exclude volts, floor, quantiles, V0
    ('det0 Ar/iso', 'det0__*Ar_Iso*', (50, 150, -200, 150), 'Ar_iC4H10_95_5', {490}, 25, (.5, .75, .9), 470),
    ('det0 He/C2H6', 'det0__*He_HV_Scan_4-2*', (50, 150, -200, 150), 'He_C2H6_96p5_3p5', {505}, 25, (.5, .75, .9), 500),
    ('det0 He/C2H6 (4-1 repeat)', 'det0__*He_HV_Scan_4-1*', (50, 150, -200, 150), 'He_C2H6_96p5_3p5', set(), 25, (.5, .75, .9), 490),
    ('det3o Ar/iso', 'det3old__*HV_Scan_5-5*', (-200, 100, -200, 50), 'Ar_iC4H10_95_5', {510}, 40, (.1, .25, .5), 480),
    ('det3o Ar/CF4', 'det3old__*ArCF4*', (-200, 100, -200, 50), 'Ar_CF4_90_10', set(), 40, (.1, .25, .5), 610),
]
pts = []
for name, pat, (x0, x1, y0, y1), gas, excl, floor, qs, V0 in SETS:
    df = pd.concat([pd.read_parquet(f) for f in glob.glob(C + pat + '.parquet')])
    df = df[(df.ref_x > x0) & (df.ref_x < x1) & (df.ref_y > y0) & (df.ref_y < y1) & ~df.subrun.str.startswith('final')]
    rng = np.random.default_rng(1)
    for (view, v), g in df.groupby(['view', 'volt']):
        a = g.peak_amp.values
        for q in qs:
            val = np.quantile(a, q)
            bs = np.quantile(a[rng.integers(0, len(a), (200, len(a)))], q, axis=1)
            ok = (val >= 5 * floor) and (val <= 0.7 * RAIL) and (v not in excl) and not (name.endswith('Ar/CF4') and v > 630)
            pts.append(dict(set=name, gas=gas, view=view, volt=v, q=f'p{int(q*100):02d}', A=val,
                            lnerr=np.std(np.log(np.maximum(bs, 1))), ok=ok, V0=V0, n=len(a)))
J = pd.read_csv('/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-13/mesh_ladder_det3_2026-06-27.csv')
for _, r in J.iterrows():
    for q in ('p02', 'p10', 'p25', 'p50'):
        pts.append(dict(set='det3 (June, T14) Ar/iso', gas='Ar_iC4H10_95_5', view=r['view'], volt=r['volt'], q=q,
                        A=r[q], lnerr=r[q + '_lnerr'], ok=bool(r[q + '_ok']), V0=490, n=r['n'], scan=r['scan']))
P = pd.DataFrame(pts); P.to_csv('ladders_all.csv', index=False)

fits = []
for (name, gas, view, q), g in P[P.ok].groupby(['set', 'gas', 'view', 'q']):
    g = g.groupby('volt').agg(A=('A', 'mean'), lnerr=('lnerr', 'max'), V0=('V0', 'first')).reset_index()
    if len(g) < 4: continue
    V, y = g.volt.values.astype(float), np.log(g.A.values)
    e = np.hypot(g.lnerr.values, 0.02); V0 = g.V0.iloc[0]; sG = garf_slope(gas)
    ex = np.polyfit(V - V0, y, 1, w=1 / e); chi_ex = np.sum(((np.polyval(ex, V - V0) - y) / e) ** 2)
    qd = np.polyfit(V - V0, y, 2, w=1 / e)
    def gfb(p, s=None):
        s = sG if s is None else s
        u = np.exp(p[1] + s * (V - V0))
        return np.where(u < 1, p[0] + s * (V - V0) - np.log1p(-np.minimum(u, 0.999999)), 1e3)
    best = min((least_squares(lambda p: (gfb(p) - y) / e, [y.mean(), lb]) for lb in (-4, -2, -1, -0.3)), key=lambda r: r.cost)
    bestf = min((least_squares(lambda p: (gfb(p[:2], p[2]) - y) / e, [y.mean(), lb, ex[0] * 0.7]) for lb in (-4, -2, -1, -0.3)), key=lambda r: r.cost)
    b, bf = np.exp(best.x[1]), np.exp(bestf.x[1])
    loc = np.diff(y) / np.diff(V) * 10
    fits.append(dict(set=name, view=view, q=q, V=f'{int(V.min())}-{int(V.max())}', n=len(V),
                     slope_lo=loc[:2].mean(), slope_hi=loc[-2:].mean(), garf_slope=sG * 10,
                     data_over_garf=ex[0] / sG, curv=qd[0] * 100,
                     chi_exp=chi_ex, chi_gfb=2 * best.cost, bG_V0=b, Vstar_gfb=V0 - np.log(b) / sG,
                     chi_ffb=2 * bestf.cost, s0_ffb=bestf.x[2] * 10, bG_ffb=bf, Vstar_ffb=V0 - np.log(bf) / bestf.x[2], V0=V0))
F = pd.DataFrame(fits); F.to_csv('fits_all.csv', index=False)
pd.set_option('display.width', 300)
print(F.round(3).to_string(index=False))
