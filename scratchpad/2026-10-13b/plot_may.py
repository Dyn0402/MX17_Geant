#!/usr/bin/env python3
"""Figure: May old-det3 ladders, Ar/iso 95/5 vs Ar/CF4 90/10 (same chamber), median
peak-strip amplitude vs mesh voltage, with Garfield's uniform-gap slope for each gas
(nTof_x17 garfield_sim/results/*_Saclay_160m.json) anchored at the lowest clean point."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
L = pd.read_csv('ladder_may.csv')
R = '/home/dylan/PycharmProjects/nTof_x17/garfield_sim/results/'
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
C = {'x': '#2a6fdb', 'y': '#d9711c'}
for ax, gas, gf, lo, hi in ((axs[0], 'ArIso_95_5', 'Ar_iC4H10_95_5', 460, 500),
                            (axs[1], 'ArCF4_90_10', 'Ar_CF4_90_10', 590, 630)):
    for view in ('x', 'y'):
        g = L[(L.gas == gas) & (L.view == view)]
        ax.plot(g.volt, g.p50, 'o', ms=4, color=C[view], alpha=0.35)
        h = g[(g.volt >= lo) & (g.volt <= hi) & g.p50.notna()]
        ax.plot(h.volt, h.p50, 'o-', ms=6, lw=2, color=C[view], label=f'data {view.upper()} median (fit range)')
    d = json.load(open(R + gf + '_Saclay_160m.json'))
    V, G = np.array(d['voltages']), np.array(d['gain_mean'])
    s = np.polyfit(V, np.log(G), 1)[0]
    a0 = L[(L.gas == gas) & (L.view == 'x') & (L.volt == lo)].p50.iloc[0]
    vv = np.linspace(lo, hi + 40, 50)
    ax.plot(vv, a0 * np.exp(s * (vv - lo)), '--', color='#555', lw=1.5, label=f'Garfield slope {10*s:.2f}/10 V')
    ax.set_yscale('log'); ax.set_xlabel('mesh (resist) voltage [V]'); ax.grid(alpha=0.25, which='both')
    ax.set_title({'ArIso_95_5': 'Ar/iC4H10 95/5 (5 May)', 'ArCF4_90_10': 'Ar/CF4 90/10 (6-7 May)'}[gas], fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc='upper left')
axs[0].set_ylabel('median peak-strip amplitude [ADC]')
fig.suptitle('Old det3 (bulked 5-5-26), same chamber, two gases: faded = noise-floor, drifting or railing points', fontsize=10)
fig.tight_layout(); fig.savefig('../../design/report/figs/may_two_gas_ladder_2026-10-13.png', dpi=140)
