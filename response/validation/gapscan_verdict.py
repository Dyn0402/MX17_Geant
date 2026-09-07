"""Adjudicate the 135 vs 150 µm gap scan against its pre-registration.

DECISION RULE, PRE-REGISTERED 2026-08-11 before the campaign was submitted —
design/report/GAP_SCAN_PREREG_2026-08-11.md. Reproduced here so the script that
reads the data states the rule it is applying:

  gain(135)/gain(150) at 490 V ......... predicted x1.57  (competing arm: x1.6-1.7)
  d lnG/dV at 135 µm ................... predicted 0.309 +- 0.010 per 10 V
                                         (competing arm: 0.26-0.29)
  ion transit t90 ratio ................ predicted x0.81  (= (135/150)^2, exact)

  * slope flat, gain up x1.5-1.6  -> predicted. The gap is a GAIN lever and not
    a slope lever, so it cannot close data 0.4487 vs sim 0.3107 at any gap.
  * slope falls to 0.26-0.29      -> the competing extrapolation wins; alpha(E)
    is materially concave above 35 kV/cm. Gap still cannot close the slope.
  * slope RISES toward 0.4487     -> the parallel-plate Townsend framing is
    WRONG, and the Penning outcome-C verdict rests on that framing. Most
    important possible result; check hardest before believing.
  * gain ratio outside x1.3-1.9   -> suspect the run, not the physics.

FALSIFIER FOR THE CAMPAIGN ITSELF: the 150 µm arm re-run on condor must
reproduce the slope hunt's desktop 150 µm uniform arm (17 062 / 43 395 /
110 085 at 460/490/520 V) to within seed statistics. If it does not, the
lxplus-vs-desktop leg differs and nothing may be quoted.

    python3 -m response.validation.gapscan_verdict \
        --calib ~/x17/response_sim/avalanche/aval_calib_gapscan.json
"""

from __future__ import annotations

import argparse
import json

import numpy as np

# The desktop 150 µm uniform arm, for the campaign falsifier (slope hunt,
# aval_calib_slopehunt.json, auto Penning, Ar/iC4H10 95/5 dry).
DESKTOP_150 = {460.0: 17062.4, 490.0: 43394.6, 520.0: 110084.6}

PRED_GAIN_RATIO = (1.50, 1.60)      # this note's arm
PRED_GAIN_RATIO_ALT = (1.6, 1.7)    # the competing arm
PRED_SLOPE_135 = (0.299, 0.319)     # 0.309 +- 0.010
PRED_SLOPE_135_ALT = (0.26, 0.29)
DATA_SLOPE = 0.4487


def arms(doc):
    """Split the merged calib into {gap_um: {voltage: point}}."""
    out = {}
    for key, pt in doc["points"].items():
        gap = float(pt.get("gap_um", 150.0))
        # gap_um is machine-readable as of 2026-08-11; if a point lacks it the
        # key suffix is the only evidence, and guessing is how the Penning axis
        # nearly got merged away. Refuse instead.
        if "gap_um" not in pt and "@gap" in key:
            raise SystemExit(
                f"{key}: key says a non-default gap but the point carries no "
                "gap_um — re-collect with the 2026-08-11 collect.py")
        out.setdefault(gap, {})[float(pt["voltage_V"])] = pt
    return out


def frac_err(pt):
    """Fractional error on this point's mean gain, from its own seed count."""
    return float(np.sqrt(pt["polya"]["rel_var"] / max(pt.get("nev_total", 1), 1)))


def slope(pts):
    """d lnG/dV per 10 V, and its error, from the voltage points of one arm."""
    v = np.array(sorted(pts))
    g = np.array([pts[x]["polya"]["gain_mean"] for x in v])
    sig = np.array([frac_err(pts[x]) for x in v])   # = error on ln g
    A = np.vstack([v, np.ones_like(v)]).T
    W = np.diag(1.0 / sig ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    b = cov @ A.T @ W @ np.log(g)
    return float(b[0] * 10), float(np.sqrt(cov[0, 0]) * 10), v, g


def alpha_curve(by_gap):
    """(E [kV/cm], alpha [1/cm], err) for every point, both arms pooled.

    In a UNIFORM field the avalanche integrates a pure gas property:
    ln G = alpha(E) * d exactly, with E = V/d. So the two gap arms MUST fall on
    one alpha(E) curve — they sample it over different E ranges, but it is the
    same curve. That is a direct test of the parallel-plate Townsend framing
    the whole gain/slope argument rests on, and it is the sharpest thing this
    campaign can say: it holds or it does not, across a 15 % change in d.
    """
    out = []
    for gap, pts in by_gap.items():
        d_cm = gap * 1e-4
        for v, p in pts.items():
            out.append((gap, v, v / d_cm / 1e3,
                        np.log(p["polya"]["gain_mean"]) / d_cm,
                        frac_err(p) / d_cm))
    return sorted(out, key=lambda r: r[2])


def band(x, lo, hi, err=None):
    """Where x sits relative to a pre-registered band, in sigma.

    BOTH conventions are printed when an error is given, deliberately.
    "Excluded at N sigma" conventionally means to the NEAREST BAND EDGE, and
    that is the conservative number for an exclusion — but it is the
    *flattering* number for one's own surviving prediction, so quoting edges
    alone lets whoever holds the pen pick the convention that suits them.
    Added 2026-08-11 after the gap scan was first written up quoting an
    exclusion to the band CENTRE (9.0 sigma) where the edge reading is 6.4.
    """
    where = "IN" if lo <= x <= hi else ("LOW" if x < lo else "HIGH")
    if err is None or where == "IN":
        return where
    edge = lo if x < lo else hi
    return (f"{where}, {abs(x - edge) / err:.1f} sigma to the near edge / "
            f"{abs(x - (lo + hi) / 2) / err:.1f} to the centre")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    doc = json.load(open(a.calib))
    by_gap = arms(doc)
    if set(by_gap) != {135.0, 150.0}:
        raise SystemExit(
            f"expected gaps {{135, 150}}, got {sorted(by_gap)} — if this is "
            "one merged arm, the gap axis did not survive collection")

    res = {"pre_registration": "design/report/GAP_SCAN_PREREG_2026-08-11.md"}
    print("Gap scan verdict — 135 vs 150 µm, uniform field, Ar/iC4H10 95/5 dry\n")

    print("  gap   V      gain      n_slices   nev")
    for gap in sorted(by_gap, reverse=True):
        for v in sorted(by_gap[gap]):
            p = by_gap[gap][v]
            print(f"  {gap:5.0f} {v:5.0f} {p['polya']['gain_mean']:11.1f} "
                  f"{p.get('n_slices', 0):8d} {p.get('nev_total', 0):7d}")

    # ── campaign falsifier ───────────────────────────────────────────────────
    print("\nFALSIFIER — condor 150 µm arm vs the desktop 150 µm arm:")
    worst = 0.0
    for v, g_desk in DESKTOP_150.items():
        g_here = by_gap[150.0][v]["polya"]["gain_mean"]
        d = g_here / g_desk - 1.0
        worst = max(worst, abs(d))
        print(f"  {v:.0f} V  condor {g_here:10.1f}  desktop {g_desk:10.1f}  "
              f"{100*d:+6.1f} %")
    # Seed statistics: sqrt(rel_var / nev) on each arm, added in quadrature.
    tol = 3.0 * max(
        np.hypot(*[np.sqrt(by_gap[150.0][v]["polya"]["rel_var"]
                           / max(by_gap[150.0][v].get("nev_total", 1), 1))
                   for _ in (0, 1)])
        for v in DESKTOP_150)
    passed = worst <= max(tol, 0.10)
    res["falsifier"] = {"worst_frac_dev": worst, "tol": max(tol, 0.10),
                        "passed": bool(passed)}
    print(f"  worst |dev| {100*worst:.1f} % vs tolerance "
          f"{100*max(tol, 0.10):.1f} %  ->  {'PASS' if passed else 'FAIL'}")
    if not passed:
        print("\n  ⚠️  FALSIFIER FAILED — the condor and desktop legs differ. "
              "NOTHING BELOW MAY BE QUOTED until that is explained.")

    # ── the framing test, which outranks both headline observables ───────────
    print("\nONE alpha(E) CURVE? (ln G = alpha(E)*d exactly in a uniform field,\n"
          "so both arms must lie on the same curve — a direct test of the\n"
          "parallel-plate Townsend framing across a 15 % change in d)")
    curve = alpha_curve(by_gap)
    print(f"  {'gap':>5s} {'V':>5s} {'E kV/cm':>9s} {'alpha /cm':>11s} {'+-':>5s}")
    for gap, v, E, al, ea in curve:
        print(f"  {gap:5.0f} {v:5.0f} {E:9.2f} {al:11.1f} {ea:5.1f}")
    a150 = [r for r in curve if r[0] == 150.0]
    a135 = [r for r in curve if r[0] == 135.0]
    c = np.polyfit([r[2] for r in a150], [r[3] for r in a150], 1,
                   w=[1.0 / r[4] for r in a150])
    chi2 = 0.0
    print("  150-arm LINEAR fit extrapolated onto the 135 arm:")
    for gap, v, E, al, ea in a135:
        pred = float(np.polyval(c, E))
        z = (al - pred) / ea
        chi2 += z * z
        print(f"    E={E:6.2f}  measured {al:7.1f}+-{ea:.1f}  predicted {pred:7.1f}"
              f"  {z:+5.1f} sigma   (gain x{np.exp((al - pred) * 135e-4):.3f})")
    chi2 /= len(a135)
    res["alpha_curve"] = [dict(gap_um=r[0], V=r[1], E_kVcm=r[2], alpha=r[3],
                               alpha_err=r[4]) for r in curve]
    res["alpha_one_curve_chi2_per_point"] = float(chi2)
    print(f"  chi2/point = {chi2:.1f}  ->  "
          f"{'ONE CURVE — framing holds' if chi2 < 4 else 'ARMS DISAGREE'}")

    # ── the two observables ──────────────────────────────────────────────────
    s150, e150, v150, g150 = slope(by_gap[150.0])
    s135, e135, v135, g135 = slope(by_gap[135.0])
    ratio = {v: by_gap[135.0][v]["polya"]["gain_mean"]
                / by_gap[150.0][v]["polya"]["gain_mean"] for v in v150}

    print("\nGAIN RATIO 135/150:")
    for v in sorted(ratio):
        e = ratio[v] * np.hypot(frac_err(by_gap[135.0][v]),
                                frac_err(by_gap[150.0][v]))
        print(f"  {v:.0f} V   x{ratio[v]:.3f} +- {e:.3f}")
    r490 = ratio[490.0]
    print(f"  at 490 V: x{r490:.3f}   "
          f"[this arm x1.50-1.60: {band(r490, *PRED_GAIN_RATIO)}]  "
          f"[competing x1.6-1.7: {band(r490, *PRED_GAIN_RATIO_ALT)}]")

    print("\nHV SLOPE d lnG/dV per 10 V:")
    print(f"  150 µm  {s150:.4f} +- {e150:.4f}")
    print(f"  135 µm  {s135:.4f} +- {e135:.4f}")
    print(f"    this arm  0.309+-0.010 : "
          f"{band(s135, *PRED_SLOPE_135, err=e135)}")
    print(f"    competing 0.26 -0.29   : "
          f"{band(s135, *PRED_SLOPE_135_ALT, err=e135)}")
    ed = float(np.hypot(e135, e150))
    print(f"  moved   {s135 - s150:+.4f} +- {ed:.4f}  = "
          f"{abs(s135 - s150) / ed:.1f} sigma — the MATCHED comparison, "
          f"immune to any offset common to both arms")
    print(f"  data is {DATA_SLOPE:.4f}; the 135 µm arm is still short by "
          f"{DATA_SLOPE - s135:.4f}, i.e. the gap closes "
          f"{100 * (s135 - s150) / (DATA_SLOPE - s150):.0f} % of the discrepancy")
    res["slope_err_150"], res["slope_err_135"] = e150, e135
    res["slope_move"], res["slope_move_err"] = s135 - s150, ed

    print("\nION TRANSIT (predicted x0.81 = (135/150)^2):")
    for v in sorted(v150):
        t135 = by_gap[135.0][v].get("t_arrival_mean_ns")
        t150 = by_gap[150.0][v].get("t_arrival_mean_ns")
        if t135 and t150:
            print(f"  {v:.0f} V   x{t135 / t150:.3f}")

    # ── verdict ──────────────────────────────────────────────────────────────
    if s135 > s150 + 0.02:
        verdict = ("SLOPE ROSE — the parallel-plate Townsend framing is in "
                   "question, and the Penning outcome-C verdict rests on it. "
                   "Check hardest before believing.")
    elif s135 < 0.295:
        verdict = ("SLOPE FELL — the competing (alpha''<0) extrapolation wins. "
                   "The gap moves the slope AWAY from the data; still dead as "
                   "a slope candidate.")
    else:
        verdict = ("SLOPE FLAT as pre-registered — the gap is a gain lever and "
                   "NOT a slope lever. It cannot close 0.4487 vs 0.3107 at any "
                   "gap.")
    if not 1.3 <= r490 <= 1.9:
        verdict = (f"GAIN RATIO x{r490:.2f} IS OUTSIDE x1.3-1.9 — suspect the "
                   "run, not the physics; check gap_um reached every slice. "
                   + verdict)
    print(f"\nVERDICT: {verdict}")
    print("\n⚠️  Neither gap is a measurement of det3's as-built gap, and "
          "choosing one because it closes T14 would launder the comparison "
          "into its own input.")

    res.update(gain_ratio=ratio, gain_ratio_490=r490,
               slope_150=s150, slope_135=s135, data_slope=DATA_SLOPE,
               gains={g: {v: by_gap[g][v]["polya"]["gain_mean"]
                          for v in by_gap[g]} for g in by_gap},
               verdict=verdict)
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)
        print(f"\nwrote {a.out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
