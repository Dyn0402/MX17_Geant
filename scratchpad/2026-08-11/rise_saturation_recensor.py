#!/usr/bin/env python3
"""Re-derive OVERNIGHT_2026-08-10 §4 with the saturation-biased rise estimator
censored.

THE DEFECT. `t14_compare.extract_view` measured the 10-90 % rise against the
waveform's own maximum with no saturation handling. On a railed waveform that
maximum IS the rail, so the 90 % level sits at 0.9 x rail — reached early on
the still-rising edge — and the pulse reads FAST. Fine if it hit both legs
alike; it does not. The data legs rail on ~24-26 % of events, the sim legs on
~3-5 %, so it is very nearly a pure data-leg speed-up, in the direction that
makes the sim look slow. Every §4 offset is therefore biased.

THE FIX HERE. Read-only on the FROZEN parquets: drop events with
peak_amp >= 3500 ADC on BOTH legs, identically, and recompute every quantile.
`t14_compare.py` now applies the same censor at source (SAT_ADC), so a re-run
reproduces this without the post-hoc step.

Nothing is fitted and no leg is treated differently. Run:

    ./.venv/bin/python scratchpad/2026-08-11/rise_saturation_recensor.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path.home() / "x17/response_sim/stageB_w2"
OUT = Path(__file__).resolve().parents[2] / "design/report"
QS = [5, 10, 25, 50, 75, 90]
SAT_ADC = 3500.0
VIEW = "x"

ANGLED = [
    ("10 deg", "t14_angW_thx10", "t14_DIAGNOSIS_noions_thx10"),
    ("20 deg", "t14_angW_thx20", "t14_DIAGNOSIS_noions_thx20"),
]
VERTICAL = [
    ("fion030", "t14_DIAGNOSIS_fion030", 0.30),
    ("fion050", "t14_DIAGNOSIS_fion050", 0.50),
    ("fion070", "t14_DIAGNOSIS_fion070", 0.70),
    ("default", "t14_compare", 0.9006),
]
VERT_NOIONS = "t14_DIAGNOSIS_noions"


def rise(dirname, leg, censor):
    """Rise quantiles for one leg, with the censor on or off."""
    t = pd.read_parquet(BASE / dirname / f"wf_{leg}_{VIEW}.parquet")
    n_all = len(t)
    sat = (t.peak_amp >= SAT_ADC)
    keep = ~sat if censor else pd.Series(True, index=t.index)
    r = t.loc[keep, "rise_ns"].to_numpy()
    finite = np.isfinite(r)
    return dict(
        q=np.percentile(r[finite], QS),
        n=int(finite.sum()),
        n_all=n_all,
        sat_frac=float(sat.mean()),
        # NaN rise among the KEPT events — a real shape statement (no 10 % or
        # 90 % crossing in the window), not a cut, so it is reported apart.
        nan_frac=float(1.0 - finite.mean()) if len(r) else float("nan"),
    )


def fmt(v):
    return "  ".join(f"{x:7.1f}" for x in v)


def main():
    doc = {"sat_adc": SAT_ADC, "view": VIEW, "quantiles": QS,
           "note": ("re-derivation of OVERNIGHT_2026-08-10 §4 with railed "
                    "events (peak_amp >= 3500 ADC) dropped identically on "
                    "both legs; read-only on the frozen parquets")}

    print(f"=== ANGLED closure, view {VIEW.upper()}, censor peak_amp < {SAT_ADC:.0f} ===\n")
    doc["angled"] = {}
    for name, d_ions, d_noions in ANGLED:
        blk = {}
        # Each sim leg is compared against the data leg extracted ALONGSIDE it,
        # not against one shared data leg: the per-directory event cap makes
        # the two data legs differ slightly, and the original
        # `noions_angled_closure.py` offsets it this way. Keeping the
        # convention is what makes this a re-derivation rather than a new
        # analysis that happens to disagree.
        for lbl, dirn, leg in (("sim_ions", d_ions, "sim"),
                               ("sim_noions", d_noions, "sim"),
                               ("data", d_ions, "data"),
                               ("data_noions_leg", d_noions, "data")):
            blk[lbl] = {k: rise(dirn, leg, k == "censored")
                        for k in ("raw", "censored")}
        print(f"-- {name}")
        print(f"{'leg':22s} {'p5':>7s} {'p10':>8s} {'p25':>8s} {'p50':>8s} "
              f"{'p75':>8s} {'p90':>8s}   n   sat%   nan%")
        for lbl in ("sim_ions", "sim_noions", "data", "data_noions_leg"):
            for mode in ("raw", "censored"):
                b = blk[lbl][mode]
                print(f"  {lbl:16s} {mode:9s} {fmt(b['q'])}  {b['n']:5d} "
                      f"{100*b['sat_frac']:5.1f} {100*b['nan_frac']:5.1f}")
        for mode in ("raw", "censored"):
            off_i = blk["sim_ions"][mode]["q"] - blk["data"][mode]["q"]
            off_n = (blk["sim_noions"][mode]["q"]
                     - blk["data_noions_leg"][mode]["q"])
            blk[f"offset_ions_{mode}"] = list(off_i)
            blk[f"offset_noions_{mode}"] = list(off_n)
            blk[f"mean_offset_ions_{mode}"] = float(off_i.mean())
            blk[f"mean_offset_noions_{mode}"] = float(off_n.mean())
            print(f"  offset with ions   {mode:9s} {fmt(off_i)}  "
                  f"mean {off_i.mean():+7.1f}  span {off_i.ptp():5.1f}")
            print(f"  offset ions removed{mode:9s} {fmt(off_n)}  "
                  f"mean {off_n.mean():+7.1f}  span {off_n.ptp():5.1f}")
        # f_eff implied by linear interpolation between "ions removed" (f=0)
        # and "with ions" (f=0.9006) on the MEAN offset: the data sits where
        # the offset would be zero.
        for mode in ("raw", "censored"):
            a = blk[f"mean_offset_noions_{mode}"]
            b = blk[f"mean_offset_ions_{mode}"]
            f = 0.9006 * (0.0 - a) / (b - a) if b != a else float("nan")
            blk[f"f_eff_{mode}"] = float(f)
            print(f"  implied f_eff      {mode:9s} {f:6.3f}")
        for k, v in list(blk.items()):
            if isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, dict) and "q" in vv:
                        vv["q"] = list(vv["q"])
        doc["angled"][name] = blk
        print()

    print(f"=== VERTICAL f_ion rigidity (shift vs noions), view {VIEW.upper()} ===\n")
    doc["vertical"] = {}
    base = {m: rise(VERT_NOIONS, "sim", m == "censored")
            for m in ("raw", "censored")}
    print(f"{'leg':22s} {'p5':>7s} {'p10':>8s} {'p25':>8s} {'p50':>8s} "
          f"{'p75':>8s} {'p90':>8s}   span")
    for name, dirn, fion in VERTICAL:
        blk = {"f_ion": fion}
        for mode in ("raw", "censored"):
            b = rise(dirn, "sim", mode == "censored")
            sh = b["q"] - base[mode]["q"]
            blk[mode] = dict(q=list(b["q"]), shift=list(sh),
                             span=float(sh.ptp()), n=b["n"],
                             sat_frac=b["sat_frac"])
            print(f"  {name:12s} {mode:9s} {fmt(sh)}  {sh.ptp():6.1f}")
        doc["vertical"][name] = blk
    doc["vertical"]["_noions_base"] = {
        m: dict(q=list(base[m]["q"]), n=base[m]["n"],
                sat_frac=base[m]["sat_frac"]) for m in base}

    # The vertical DATA leg's own NaN-rise rate, which the peer flagged as a
    # one-sided ~9 % loss. It is a property of the data leg, not of the censor.
    v = pd.read_parquet(BASE / "t14_compare" / f"wf_data_{VIEW}.parquet")
    s = pd.read_parquet(BASE / "t14_compare" / f"wf_sim_{VIEW}.parquet")
    doc["vertical_nan_rise"] = dict(
        data=float(v.rise_ns.isna().mean()), sim=float(s.rise_ns.isna().mean()),
        data_sat=float((v.peak_amp >= SAT_ADC).mean()),
        sim_sat=float((s.peak_amp >= SAT_ADC).mean()))
    print(f"\nvertical NaN-rise:  data {100*doc['vertical_nan_rise']['data']:.1f} %  "
          f"sim {100*doc['vertical_nan_rise']['sim']:.1f} %   "
          f"sat: data {100*doc['vertical_nan_rise']['data_sat']:.1f} %  "
          f"sim {100*doc['vertical_nan_rise']['sim_sat']:.1f} %")

    p = OUT / "rise_saturation_recensor_2026-08-11.json"
    json.dump(doc, open(p, "w"), indent=1)
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
