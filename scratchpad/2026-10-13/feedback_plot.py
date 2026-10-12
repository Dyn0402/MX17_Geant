#!/usr/bin/env python3
"""Figure for FRESH_LOOK §9: the bench gain ladder against Townsend-only Garfield and
Garfield x photon feedback G/(1-bG). Args: mesh_ladder.csv aval_calib_meshfield_hvscan.json out.png"""
import json, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

d = pd.read_csv(sys.argv[1]); pts = json.load(open(sys.argv[2]))["points"]
sv = np.array([v["voltage_V"] for k, v in pts.items() if "95_5" in k])
sg = np.array([v["polya"]["gain_mean"] for k, v in pts.items() if "95_5" in k])
o = np.argsort(sv); sv, sg = sv[o], sg[o]
P = np.polyfit(sv - 490, np.log(sg), 1)
Gs = lambda V: np.exp(np.polyval(P, V - 490))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
cols = {"p50": "#2a6fdb", "p10": "#d9822b", "p02": "#2f9e44"}
Vf = np.linspace(420, 515, 300)
for est in ("p50", "p10", "p02"):
    t = d[(d.view == "x") & d[est + "_ok"]]
    V = t.volt.values.astype(float); y = np.log(t[est].values)
    e = np.maximum(t[est + "_lnerr"].values, 0.02)
    f = least_squares(lambda p: (y - (p[0] + np.log(Gs(V) / Gs(490)) - np.log(1 - p[1] * Gs(V) / Gs(490)))) / e,
                      [y.max(), 0.3], bounds=([-np.inf, 0], [np.inf, 0.999]))
    c, bb = f.x
    fs = least_squares(lambda p: (y - (p[0] + np.log(Gs(V) / Gs(490)))) / e, [y.max()])
    ax[0].errorbar(V, np.exp(y), yerr=np.exp(y) * e, fmt="o", ms=4, color=cols[est], label=f"data X {est}")
    g = Gs(Vf) / Gs(490)
    ok = bb * g < 0.98
    ax[0].plot(Vf[ok], np.exp(c) * g[ok] / (1 - bb * g[ok]), "-", color=cols[est], lw=1.5,
               label=f"Garfield x feedback, bG(490)={bb:.2f}")
    ax[0].plot(Vf, np.exp(fs.x[0]) * g, ":", color=cols[est], lw=1)
    # local slope
    Vm = 0.5 * (V[1:] + V[:-1]); sl = np.diff(y) / np.diff(V) * 10
    ax[1].plot(Vm, sl, "o", ms=4, color=cols[est])
    ax[1].plot(Vf[ok], P[0] * 10 / (1 - bb * g[ok]), "-", color=cols[est], lw=1.5)
ax[0].set_yscale("log"); ax[0].set_xlabel("mesh voltage [V]"); ax[0].set_ylabel("peak-strip amplitude [ADC]")
ax[0].set_title("det3 27 June ladder (rail-safe points)\nsolid: Garfield x G/(1-bG), dotted: Garfield alone (scale fitted)", fontsize=9)
ax[0].legend(fontsize=7)
ax[1].axhline(P[0] * 10, color="k", ls=":", label=f"Garfield (Townsend only) {P[0]*10:.3f}")
ax[1].set_xlabel("mesh voltage [V]"); ax[1].set_ylabel("d ln A / dV  per 10 V")
ax[1].set_title("local slope: data (points, finite differences) vs feedback model (lines)", fontsize=9)
ax[1].set_ylim(0, 1.0); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(sys.argv[3], dpi=130)
print("wrote", sys.argv[3])
