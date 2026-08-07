#!/usr/bin/env python3
"""
test_lut_vs_solver.py — T10 fast-path certification, channel by channel.

THE TWO-TIER CLAIM (plan principle 3): a rigorous slow path validates a fast
path that is the production digitizer. Same physics, different caching. Until
now the fast path had never been checked against anything but itself.

THE REFERENCE. `response.solver.kernels.charge_budget_y` / `charge_budget_x`
read the per-channel charge straight out of the S1 product with the indexing
the solver itself uses — the same code path whose closed-form sum rule passes
at 4.8e-07. That is the right reference precisely because hand-rolled indexing
of these arrays is where the mistakes live: the Y kernels are stored with y = 0
at ny//2 (`_to_y0_origin`), rows alternate parity, and X kernels are indexed by
absolute column mod the 40-pad superperiod. Three separate ad-hoc attempts to
reproduce the sum rule on 2026-08-07 all got it wrong.

THE THING UNDER TEST is everything the LUT does to that product on the way to
the digitizer: window it in y, stride it in x, resample a 61-point log time
axis onto a uniform 1 ns grid, differentiate G into a current, and re-index the
X kernels from absolute column to channel offset. Any of those can silently
scramble or lose charge.

WHAT THE RIGHT TOTAL IS, since this cost an afternoon. The closed-form
`prompt_sum_rule` = S(0)/C(0) = 0.875 applies ONLY to the fictitious
pitch-sized (0.78 mm) pads that `check_sum_rule` substitutes so that a closed
form exists. Production kernels use the REAL 0.68 mm pads, which do not tile:
(0.68/0.78)^2 = 0.760 of the plane is pad and the rest is inter-pad gap, whose
image charge lands on no readout channel. The solver already computes the
correct figure and stores it in every product as `channel_capture_prompt`
(0.665 for this geometry). Asserting against 0.875 produces a phantom 23 %
"loss" that is really the inter-pad copper.

    python3 -m response.digitizer.test_lut_vs_solver --kernel <s1.npz>
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from ..common import constants as C
from .kernel_lut import CombKernelLUT
from ..solver import kernels as K

TOL = 0.02          # plan T10: < 2 % residual, fast vs rigorous


class _Grid:
    """The bits of a solver object that charge_budget_* actually reads."""

    def __init__(self, x, y):
        self.x, self.y = x, y
        self.nx, self.ny = len(x), len(y)
        self.lx = float(self.nx * (x[1] - x[0]))
        self.ly = float(self.ny * (y[1] - y[0]))


def certify(path, dmax=3, n_probe=6, seed=3):
    d = np.load(path)
    meta = json.loads(str(d["meta"]))
    t_src = d["t"]
    sy = _Grid(d["x"], d["y_Y"])
    sx = _Grid(d["x"], d["y_X"])
    gy = (d["G_Y_even"], d["G_Y_odd"])
    GX = d["G_X"]
    gx = {c: GX[c] for c in range(C.N_PAD_PER_SUPER)}

    lut = CombKernelLUT(path)
    # Compare at the LAST time the LUT covers, where every approximation it
    # makes has had maximum opportunity to accumulate.
    t_ref = float(lut.t[-1])
    it = int(np.argmin(np.abs(t_src - t_ref)))
    t_ref = float(t_src[it])
    kt = int(np.argmin(np.abs(lut.t - t_ref)))

    print(f"  product   {path.split('/')[-1]}")
    print(f"  compare at t = {t_ref*1e9:.0f} ns "
          f"(source index {it}, LUT sample {kt})")
    print(f"  solver's own channel_capture_prompt = "
          f"{meta['channel_capture_prompt']:.6f}  "
          f"(NOT sum_rule_expect = {meta['sum_rule_expect']:.3f}: that is the "
          f"tiling-pad closed form)\n")

    rng = np.random.default_rng(seed)
    worst = 0.0
    ok = True

    for _ in range(n_probe):
        c0 = int(rng.integers(0, C.N_PAD_PER_SUPER))
        x0 = float(K.pad_x(c0))
        par = 0                                  # deposit on an even row

        ref_y = K.charge_budget_y(sy, gy, x0, dmax=dmax, row0_parity=par)
        ref_x = K.charge_budget_x(sx, gx, x0, dmax=dmax, y0_m=0.0)

        # --- the same quantities out of the LUT ---------------------------
        # LUT holds CURRENT on a uniform grid, so charge = cumulative integral.
        ix = int(np.argmin(np.abs(lut.x - np.mod(x0, lut.x[-1] - lut.x[0]
                                                 + (lut.x[1] - lut.x[0])))))
        got_y, got_x = {}, {}
        for dd in range(-dmax, dmax + 1):
            # Y channel d rows away: its kernel parity is (par+d)%2 and the
            # source sits at y_rel = -d * pad pitch relative to that channel.
            iy = lut.iy_Y(-dd * C.PAD_PITCH_M)
            got_y[dd] = float(lut.I_Y[(par + dd) % 2, iy, ix, :kt + 1].sum()
                              * lut.dt)
            jd = dd + lut.n_side
            if 0 <= jd < lut.I_X.shape[0]:
                jy = lut.iy_X(0.0)
                got_x[dd] = float(lut.I_X[jd, jy, ix, :kt + 1].sum() * lut.dt)

        for tag, ref, got in (("Y", ref_y, got_y), ("X", ref_x, got_x)):
            ds = sorted(got)
            r = np.array([ref[dd][it] for dd in ds])
            g = np.array([got[dd] for dd in ds])
            scale = float(np.abs(r).max())
            if scale <= 0:
                continue
            err = float(np.abs(g - r).max() / scale)
            worst = max(worst, err)
            good = err <= TOL
            ok &= good
            if not good or _ == 0:
                print(f"  col {c0:2d} {tag}  ref " +
                      " ".join(f"{v:7.4f}" for v in r))
                print(f"         {' '*len(tag)}  lut " +
                      " ".join(f"{v:7.4f}" for v in g) +
                      f"   max rel {err:.4f}  {'OK' if good else 'FAIL'}")

    print(f"\n  worst per-channel residual over {n_probe} source positions: "
          f"{worst:.4f}   tol {TOL}")
    print("\n" + ("PASS" if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--n-probe", type=int, default=6)
    a = ap.parse_args()
    return 0 if certify(a.kernel, n_probe=a.n_probe) else 1


if __name__ == "__main__":
    raise SystemExit(main())
