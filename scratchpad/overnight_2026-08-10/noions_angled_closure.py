#!/usr/bin/env python3
"""The closing test: at INCLINED incidence, does removing the ion term land the
sim's rise distribution on the data at every quantile?

Chain of reasoning this closes:

  1. Angled ladder -- at 10 and 20 deg the sim-vs-data rise mismatch is a RIGID
     ~65-80 ns offset at every quantile, not a missing fast population.
  2. f_ion rigidity (vertical legs) -- the f_ion dial does NOT act as a delay
     there: 0 -> 0.9006 moves p5 by +115 ns but p90 by -28 ns. It COMPRESSES.
  3. So the obvious objection: if f_ion is not a rigid delay, how can it be the
     ~70 ns rigid offset?

The resolution has to be measured, not argued. At vertical the sim's rise
distribution is wide (p5-p90 spans ~280 ns) and the ion floor only bites on the
fast half. At 20 deg the distribution is already compressed to a ~63 ns span --
essentially everything sits ON the floor -- so there the same dial should act
much more nearly as a uniform shift.

Test: the DIAGNOSIS_noions legs at thx10/thx20 are on disk. Compare, at each
inclined point, sim-with-ions vs sim-without-ions vs data, quantile by
quantile.

Pre-stated outcomes:
  * If noions lands ON the data at every quantile at 10 and 20 deg (offsets
    collapsing from ~70 ns to near zero with small span), the ion term IS the
    whole rise discrepancy at inclined incidence, and the f_ion contradiction
    is confirmed on the full distribution rather than on p5 alone.
  * If noions OVERSHOOTS (sim now faster than data), the ion term is too strong
    but something else is also wrong.
  * If noions still leaves a rigid offset, the offset is NOT the ion term and
    the fixed-delay candidates (shaper group delay, template t0, seeding
    latency) move to the front.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
QS = [5, 10, 25, 50, 75, 90]

POINTS = [
    ("10 deg", "t14_angW_thx10", "t14_DIAGNOSIS_noions_thx10"),
    ("20 deg", "t14_angW_thx20", "t14_DIAGNOSIS_noions_thx20"),
]
VIEW = "x"


def q(d, leg):
    r = pd.read_parquet(BASE / d / f"wf_{leg}_{VIEW}.parquet")["rise_ns"].to_numpy()
    r = r[np.isfinite(r)]
    return np.percentile(r, QS), r.size


out = {}
for name, d_ions, d_noions in POINTS:
    qi, ni = q(d_ions, "sim")
    qn, nn = q(d_noions, "sim")
    qd, nd = q(d_ions, "data")
    qd2, nd2 = q(d_noions, "data")

    print(f"\n===== X view, {name} =====")
    print(f"{'leg':>22} {'n':>6} " + " ".join(f"{'p'+str(p):>7}" for p in QS))
    print(f"{'sim, with ions':>22} {ni:>6d} " + " ".join(f"{v:>7.1f}" for v in qi))
    print(f"{'sim, ions REMOVED':>22} {nn:>6d} " + " ".join(f"{v:>7.1f}" for v in qn))
    print(f"{'data (ions leg)':>22} {nd:>6d} " + " ".join(f"{v:>7.1f}" for v in qd))
    if not np.allclose(qd, qd2, rtol=0, atol=1.0):
        print(f"{'data (noions leg)':>22} {nd2:>6d} "
              + " ".join(f"{v:>7.1f}" for v in qd2)
              + "   <-- data legs differ, offsets below use each leg's own data")

    off_i = qi - qd
    off_n = qn - qd2
    print(f"{'':>22} {'':>6} " + " ".join("-------" for _ in QS))
    print(f"{'offset, with ions':>22} {'':>6} "
          + " ".join(f"{v:>7.1f}" for v in off_i)
          + f"   span {off_i.max()-off_i.min():6.1f}  mean {off_i.mean():6.1f}")
    print(f"{'offset, ions removed':>22} {'':>6} "
          + " ".join(f"{v:>7.1f}" for v in off_n)
          + f"   span {off_n.max()-off_n.min():6.1f}  mean {off_n.mean():6.1f}")
    print(f"{'ion term buys':>22} {'':>6} "
          + " ".join(f"{v:>7.1f}" for v in (qn - qi))
          + f"   span {(qn-qi).max()-(qn-qi).min():6.1f}")

    out[name] = dict(quantiles=QS,
                     sim_ions=qi.tolist(), sim_noions=qn.tolist(),
                     data=qd.tolist(), data_noions_leg=qd2.tolist(),
                     offset_ions=off_i.tolist(), offset_noions=off_n.tolist(),
                     mean_offset_ions=float(off_i.mean()),
                     mean_offset_noions=float(off_n.mean()),
                     span_offset_ions=float(off_i.max() - off_i.min()),
                     span_offset_noions=float(off_n.max() - off_n.min()))

p = BASE / "t14_ang_trend" / "noions_angled_closure.json"
p.write_text(json.dumps(out, indent=1))
print(f"\nwrote {p}")
