#!/usr/bin/env python3
"""Per-voltage quantile ladder for the May old-det3 Ar/iso and Ar/CF4 scans
(parquets from extract_may.py). Fiducial chosen once from the top-voltage
efficiency maps of both gases (x -200..100, y -200..50 mm, ref frame at z=232),
applied at every voltage. Quantiles used only while below 0.70 of the rail
(the hv_slope rule). Writes ladder_may.csv.

    python ladder_may.py PARQUET_DIR
"""
import sys, glob, numpy as np, pandas as pd
RAIL, SAT_FRAC = 3872.0, 0.70
d = sys.argv[1]
df = pd.concat([pd.read_parquet(f) for f in glob.glob(f'{d}/*.parquet')])
df = df[(df.ref_x > -200) & (df.ref_x < 100) & (df.ref_y > -200) & (df.ref_y < 50)]
df['gas'] = np.where(df.scan.str.contains('CF4'), 'ArCF4_90_10', 'ArIso_95_5')
rng = np.random.default_rng(1)
rows = []
for (gas, view, v), g in df.groupby(['gas', 'view', 'volt']):
    a = g.peak_amp.values
    r = dict(gas=gas, view=view, volt=v, n=len(a), eff=(g.n_over > 0).mean(), frail=(a > 0.88 * RAIL).mean())
    for q in (0.10, 0.25, 0.50):
        val = np.quantile(a, q)
        bs = np.quantile(a[rng.integers(0, len(a), (300, len(a)))], q, axis=1)
        r[f'p{int(q*100):02d}'] = val if val < SAT_FRAC * RAIL else np.nan
        r[f'e{int(q*100):02d}'] = np.std(np.log(bs))
    rows.append(r)
L = pd.DataFrame(rows).sort_values(['gas', 'view', 'volt'])
L.to_csv('ladder_may.csv', index=False)
pd.set_option('display.width', 200)
print(L.round(3).to_string(index=False))
