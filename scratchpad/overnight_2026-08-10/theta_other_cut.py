#!/usr/bin/env python3
"""Is the vertical leg's non-rigid rise offset a one-view-window artifact?

The angled ladder shows the sim's rise distribution is the data's rigidly
displaced by ~65-80 ns at every quantile -- at the two INCLINED points. At the
vertical point alone the offset is not rigid: +90 ns at p5 but -33 ns at p90,
i.e. the data's vertical leg has extra spread in both tails.

Suspected cause (T14_CAMPAIGN open thread #4): t14_compare's _sel_ids windows
the data on ONE view's theta. So the "vertical" X leg is cut to |theta_x|<3 deg
with theta_y left completely free, admitting cosmics steeply inclined in y --
while the sim gun is a true pencil beam, vertical in both views. Inclined
tracks are faster, so the data leg gets a fast contamination the sim cannot
have.

Test: intersect the vertical data leg with the June wft analysis of the SAME
run, which carries per-view theta, and additionally require |theta_other|<3 deg.
If the non-rigidity is the artifact, the offset should flatten toward the ~70 ns
rigid value the inclined points show.

Falsifier: if cutting theta_other leaves the offset just as non-rigid
(still strongly negative at p90), the vertical broadening is real and needs a
physical explanation, not a selection fix.

Caveat carried into the result: the two legs are only partially overlapping
(the June wft pass covers 7093 of the run's 47452 events), so this runs on a
subsample and the control below is what makes it interpretable -- the
comparison that matters is overlap-with-cut vs overlap-WITHOUT-cut, both drawn
from the same 703 events, not against the full 2500-event leg.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
THETA = Path("/media/dylan/data/x17/cosmic_bench/Analysis/"
             "mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V/"
             "mx17_3/wft/events.parquet")
QS = [5, 10, 25, 50, 75, 90]
HALFWIDTH = 3.0


def qs(r):
    r = np.asarray(r, dtype=float)
    r = r[np.isfinite(r)]
    return np.percentile(r, QS), r.size


def main():
    ev = pd.read_parquet(THETA, columns=["event_id", "x_theta_deg",
                                         "y_theta_deg"])
    out = {}
    for view, other in (("x", "y"), ("y", "x")):
        sim = pd.read_parquet(BASE / "t14_angW_th00" / f"wf_sim_{view}.parquet")
        dat = pd.read_parquet(BASE / "t14_angW_th00" / f"wf_data_{view}.parquet")
        m = dat.merge(ev, on="event_id", how="inner")
        keep = m[f"{other}_theta_deg"].abs() < HALFWIDTH

        q_sim, n_sim = qs(sim["rise_ns"])
        q_full, n_full = qs(dat["rise_ns"])          # the leg as published
        q_ovl, n_ovl = qs(m["rise_ns"])              # overlap, no extra cut
        q_cut, n_cut = qs(m.loc[keep, "rise_ns"])    # overlap + |theta_other|<3

        print(f"\n===== vertical point, {view.upper()} view "
              f"(cutting |theta_{other}| < {HALFWIDTH} deg) =====")
        print(f"{'leg':>28} {'n':>6} " + " ".join(f"{'p'+str(p):>7}" for p in QS))
        for lab, q, n in (("sim (pencil beam)", q_sim, n_sim),
                          ("data, published leg", q_full, n_full),
                          ("data, overlap only", q_ovl, n_ovl),
                          (f"data, overlap + theta_{other} cut", q_cut, n_cut)):
            print(f"{lab:>28} {n:>6d} " + " ".join(f"{v:>7.1f}" for v in q))

        print(f"{'':>28} {'':>6} " + " ".join(f"{'-'*7}" for _ in QS))
        for lab, q in (("offset vs published", q_sim - q_full),
                       ("offset vs overlap", q_sim - q_ovl),
                       ("offset vs overlap+cut", q_sim - q_cut)):
            span = q.max() - q.min()
            print(f"{lab:>28} {'':>6} " + " ".join(f"{v:>7.1f}" for v in q)
                  + f"   span {span:6.1f} ns")

        out[view] = dict(
            n=dict(sim=n_sim, published=n_full, overlap=n_ovl, cut=n_cut),
            quantiles=QS,
            offset_published=(q_sim - q_full).tolist(),
            offset_overlap=(q_sim - q_ovl).tolist(),
            offset_cut=(q_sim - q_cut).tolist(),
            span_published=float((q_sim - q_full).max() - (q_sim - q_full).min()),
            span_overlap=float((q_sim - q_ovl).max() - (q_sim - q_ovl).min()),
            span_cut=float((q_sim - q_cut).max() - (q_sim - q_cut).min()),
        )

    p = BASE / "t14_ang_trend" / "theta_other_cut.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {p}")
    print("\nReminder: 'span' is the spread of the sim-data offset across "
          "p5..p90. Small span = rigid displacement. The inclined points sit "
          "at span ~10-13 ns.")


if __name__ == "__main__":
    main()
