#!/usr/bin/env python3
"""Local slopes and a free-slope feedback fit on the May old-det3 ladders
(ladder_may.csv). Model: ln A = ln c + s0 (V - V0) - ln(1 - beta exp(s0 (V - V0))),
beta = bG(V0); runaway at V* = V0 - ln(beta)/s0. Compared with a plain
exponential and a free quadratic (same/one-more parameters).
Ranges: Ar/iso <= 500 V (510 V drifted +25 % within its run, first after install);
Ar/CF4 from where the quantile is >= 150 ADC (noise-max floor ~35-40 ADC) to 630 V
(above 640 V every quantile flattens: a separate saturation regime).
"""
import numpy as np, pandas as pd
from scipy.optimize import least_squares
L = pd.read_csv('ladder_may.csv')
FLOOR = 0.02   # 2 % error floor, as in feedback_fit.py
def fb(p, V, V0): c, s, lb = p; return c + s * (V - V0) - np.log1p(-np.exp(lb + s * (V - V0)))
out = []
for (gas, view), g in L.groupby(['gas', 'view']):
    vmax = 500 if gas.startswith('ArIso') else 630
    V0 = 490 if gas.startswith('ArIso') else 610
    for q in ('p10', 'p25', 'p50'):
        h = g[(g.volt <= vmax) & (g[q] >= 150) & g[q].notna()]
        if len(h) < 4: continue
        V, y = h.volt.values.astype(float), np.log(h[q].values)
        e = np.hypot(h['e' + q[1:]].values, FLOOR)
        loc = np.diff(y) / np.diff(V) * 10
        lin = np.polyfit(V - V0, y, 1, w=1 / e); chi_lin = np.sum(((np.polyval(lin, V - V0) - y) / e) ** 2)
        qd = np.polyfit(V - V0, y, 2, w=1 / e); chi_q = np.sum(((np.polyval(qd, V - V0) - y) / e) ** 2)
        best = None
        for lb0 in (-3, -2, -1, -0.5):
            try:
                r = least_squares(lambda p: (fb(p, V, V0) - y) / e, [y.mean(), lin[0] * 0.7, lb0],
                                  bounds=([-np.inf, 0.001, -12], [np.inf, 1, np.log(0.999) - (V.max() - V0) * 0.0]))
            except Exception: continue
            if best is None or r.cost < best.cost: best = r
        c, s, lb = best.x
        beta = np.exp(lb); vstar = V0 - lb / s
        out.append(dict(gas=gas, view=view, q=q, V=f'{int(V.min())}-{int(V.max())}', n=len(V),
                        local_slopes=' '.join(f'{x:.2f}' for x in loc),
                        curv_per_10V2=qd[0] * 100, chi_exp=chi_lin, chi_quad=chi_q, chi_fb=2 * best.cost,
                        s0_per10V=s * 10, bG_V0=beta, gain_factor=1 / (1 - beta), V_runaway=vstar))
O = pd.DataFrame(out)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
print(O.round(3).to_string(index=False))
