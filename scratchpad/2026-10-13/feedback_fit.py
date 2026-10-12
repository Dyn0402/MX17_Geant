#!/usr/bin/env python3
"""
feedback_fit.py -- does secondary-avalanche feedback, G_eff = G/(1 - b G), explain
the T14 gain slope AND amplitude at once?

Why: d ln G_eff/dV = s0/(1 - bG) and G_eff/G = 1/(1 - bG). One mechanism multiplies the
gain and the slope by the SAME factor. The bench has slope x1.52 (HV_SLOPE_2026-08-09)
and amplitude x1.77-1.82 (T14 sim/data 0.565 / 0.549 at 490 V), and its local slope
rises with voltage. Pure Townsend (the sim) has a flat slope.

Data: det3 saturday scan, mesh_ladder.csv (threshold-free peak-strip max, M3 selected,
rail-safe points only), from nTof_x17 mx17_sim_wft/hv_slope (laptop ~/x17/response_sim/hv_slope).
Sim:  MX17_Geant response/avalanche/aval_calib_meshfield_hvscan.json, dry 95/5, 460-530 V.

The fit uses only the SHAPE of the data ladder (a free scale c), so it is independent of
any ADC/fC calibration. The amplitude ratio it predicts at 490 V is then compared with
T14's measured one.
"""
import json, sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ladder = sys.argv[1]
simj = sys.argv[2]

pts = json.load(open(simj))["points"]
sv, sg = [], []
for k, v in pts.items():
    if "95_5" in k:
        sv.append(v["voltage_V"]); sg.append(v["polya"]["gain_mean"])
sv, sg = np.array(sv), np.array(sg)
o = np.argsort(sv); sv, sg = sv[o], sg[o]
# sim: ln G linear in V (its local slopes scatter 0.26-0.34 per 10 V around a flat mean)
P = np.polyfit(sv - 490, np.log(sg), 1, w=None)
P2 = np.polyfit(sv - 490, np.log(sg), 2)
print(f"sim: ln G = {P[1]:.4f} + {P[0]*10:.4f}/10V (V-490);  G(490) = {np.exp(P[1]):.0f}; "
      f"quadratic term {P2[0]*100:+.4f} per (10V)^2")
Gs = lambda V: np.exp(np.polyval(P, V - 490))

d = pd.read_csv(ladder)
res_all = {}
for view in ("x", "y"):
    for est in ("p50", "p25", "p10", "p02"):
        t = d[(d.view == view) & (d[est + "_ok"])]
        V = t.volt.values.astype(float); A = t[est].values; e = t[est + "_lnerr"].values
        e = np.maximum(e, 0.02)          # floor: gas drift between the two passes is 2-4 %
        y = np.log(A)

        def r_fb(p):   # c, b' = b*G(490)
            c, bb = p
            g = Gs(V) / Gs(490.0)
            return (y - (c + np.log(g) - np.log(1 - bb * g))) / e
        def r_exp(p):
            return (y - (p[0] + p[1] * (V - 490))) / e
        def r_quad(p):
            return (y - (p[0] + p[1] * (V - 490) + p[2] * (V - 490) ** 2)) / e
        def r_sim(p):  # sim shape, free scale only
            return (y - (p[0] + np.log(Gs(V)))) / e
        f1 = least_squares(r_fb, [np.log(A.max()), 0.3], bounds=([-np.inf, 0], [np.inf, 0.999]))
        f2 = least_squares(r_exp, [np.log(A.max()), 0.04])
        f3 = least_squares(r_quad, [np.log(A.max()), 0.04, 0.0])
        f4 = least_squares(r_sim, [np.log(A.max())])
        n = len(V)
        chi = lambda f: float(np.sum(f.fun ** 2))
        bb = f1.x[1]
        # uncertainty on bb from the Jacobian
        J = f1.jac; cov = np.linalg.pinv(J.T @ J); sbb = float(np.sqrt(cov[1, 1]))
        enh = 1 / (1 - bb)
        slope490 = P[0] * 10 / (1 - bb)
        res_all[(view, est)] = dict(n=n, V=(V.min(), V.max()), bb=bb, sbb=sbb, enh=enh,
                                    chi_fb=chi(f1), chi_exp=chi(f2), chi_quad=chi(f3), chi_sim=chi(f4),
                                    slope_exp=f2.x[1] * 10, slope490_fb=slope490, quad=f3.x[2] * 100)
        print(f"{view} {est}: {n:2d} pts {V.min():.0f}-{V.max():.0f} V | feedback bG(490) = {bb:.3f} +- {sbb:.3f}"
              f" -> gain x{enh:.2f}, slope@490 {slope490:.3f}/10V | chi2: feedback {chi(f1):6.1f} "
              f"(ndf {n-2}), free-exp {chi(f2):6.1f} (ndf {n-2}, slope {f2.x[1]*10:.3f}), "
              f"quad {chi(f3):6.1f} (ndf {n-3}, curv {f3.x[2]*100:+.4f}), sim-shape {chi(f4):7.1f} (ndf {n-1})")

print("\nT14 measured amplitude deficit at 490 V (data/sim): X 1/0.565 = 1.77 (published), "
      "1/0.549 = 1.82 (rail-bias corrected); views are identical unselected.")

# ---- control: can ANY Townsend alpha(E) give the measured curvature? ----------------
# parallel-plate Townsend: ln G = K exp(-D/V); with the free ADC scale c that is 3 params
# (one more than the feedback model).  Curvature sign is (D/V - 2): for D/V < 2 it is
# negative, so a rising slope needs D/V > 2 -- check what the fit does.
print("\nTownsend-form control: ln A = c + K exp(-D/V)   [3 params]")
for view in ("x", "y"):
    for est in ("p50", "p10", "p02"):
        t = d[(d.view == view) & (d[est + "_ok"])]
        V = t.volt.values.astype(float); y = np.log(t[est].values)
        e = np.maximum(t[est + "_lnerr"].values, 0.02)
        def r_t(p):
            c, lk, D = p
            return (y - (c + np.exp(lk) * np.exp(-D / V))) / e
        best = None
        for D0 in (300, 600, 900, 1200, 1500, 2000):
            lnG490 = 10.0
            lk0 = np.log(lnG490) + D0 / 490
            f = least_squares(r_t, [np.log(t[est].values.max()) - lnG490, lk0, D0], max_nfev=20000)
            if best is None or np.sum(f.fun**2) < np.sum(best.fun**2): best = f
        c, lk, D = best.x
        lnG490 = np.exp(lk) * np.exp(-D / 490)
        print(f"{view} {est}: chi2 {np.sum(best.fun**2):6.1f} (ndf {len(V)-3}); D/V(490) = {D/490:.2f}; "
              f"implied ln G(490) = {lnG490:.1f} (G = {np.exp(lnG490):.3g}); sim ln G(490) = {P[1]:.2f}")
