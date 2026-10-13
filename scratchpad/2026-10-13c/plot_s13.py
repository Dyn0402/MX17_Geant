#!/usr/bin/env python3
"""FRESH_LOOK §13 figures: (a) gain ladders, local slope per 10 V vs voltage for every
chamber x gas with Garfield's slope for that gas; (b) det0 Ar/iso mean pulse against the
He/C2H6 template folded with fast (fitted) and Ar+ (340 ns) ion signals."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import he_template_fit as H
P = pd.read_csv('ladders_all.csv')
COL = {'Ar_iC4H10_95_5': '#2a6fdb', 'Ar_CF4_90_10': '#d9711c', 'He_C2H6_96p5_3p5': '#2e9e5b'}
GS = {'Ar_iC4H10_95_5': 0.318, 'Ar_CF4_90_10': 0.208, 'He_C2H6_96p5_3p5': 0.190}
MK = {'det3 (June, T14) Ar/iso': 'o', 'det3o Ar/iso': 's', 'det0 Ar/iso': '^', 'det3o Ar/CF4': 'D', 'det0 He/C2H6': 'v'}
QS = {'det3 (June, T14) Ar/iso': 'p10', 'det3o Ar/iso': 'p25', 'det0 Ar/iso': 'p75', 'det3o Ar/CF4': 'p50', 'det0 He/C2H6': 'p75'}
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
for s, q in QS.items():
    g = P[(P.set == s) & (P.q == q) & P.ok & (P.view == 'x')]
    if s.startswith('det3 (June'): g = g[g.scan == 'scan1'] if 'scan' in g else g
    g = g.groupby('volt').A.mean()
    V = g.index.values.astype(float); y = np.log(g.values)
    if len(V) < 3: continue
    gas = P[P.set == s].gas.iloc[0]
    ax.plot(0.5 * (V[1:] + V[:-1]), np.diff(y) / np.diff(V) * 10, marker=MK[s], color=COL[gas], lw=1.5, label=f'{s} (X {q})')
for gas, sl in GS.items():
    ax.axhline(sl, color=COL[gas], ls='--', lw=1)
ax.set_xlabel('mesh (resist) voltage [V]'); ax.set_ylabel('local slope d lnA / dV  [per 10 V]')
ax.set_title('Gain ladders; dashed = Garfield slope for that gas', fontsize=10)
ax.set_ylim(0, 2.5); ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25)
ax = axs[1]
He = H.collect('det0__*He_HV_Scan_4-2*', 495, 515, 0, .05); Ar = H.collect('det0__*Ar_Iso*', 440, 480, 0, .05)
g, th = H.fine_avg(He); _, ta = H.fine_avg(Ar); w = (g >= -200) & (g <= 800)
ax.plot(g[w], th[w], color=COL['He_C2H6_96p5_3p5'], lw=1.5, label='det0 He/C2H6 (template)')
ax.plot(g[w], ta[w], color=COL['Ar_iC4H10_95_5'], lw=2.5, label='det0 Ar/iso data')
ax.plot(g[w], H.model(g, th, 0.96, 150, 1)[w], 'k--', lw=1.2, label='template (x) ion f 0.96, T 150 ns (fit)')
ax.plot(g[w], H.model(g, th, 0.9, 340, 6)[w], color='#b03030', ls=':', lw=1.5, label='template (x) Ar+ ions f 0.9, T 340 ns')
ax.plot(g[w], H.model(g, th, 0.2, 340, 8)[w], color='#888', ls='-.', lw=1.2, label='template (x) suppressed ions f 0.2, T 340 ns')
ax.set_xlabel('t from leading 50 % [ns]'); ax.set_ylabel('mean peak-strip pulse (norm.)')
ax.set_title('det0 vertical tracks: shaper-free ion test', fontsize=10); ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig('../../design/report/figs/s13_ladders_and_ion_2026-10-13.png', dpi=140)
