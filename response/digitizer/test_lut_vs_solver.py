#!/usr/bin/env python3
"""
test_lut_vs_solver.py — T10 fast-path certification, channel by channel.

THE TWO-TIER CLAIM (plan principle 3): a rigorous slow path validates a fast
path that is the production digitizer. Same physics, different caching. Until
now the fast path had never been checked against anything but itself.

THE REFERENCE, AND EXACTLY WHAT IT IS WORTH (corrected 2026-08-07, audit A4).
`response.solver.kernels.charge_budget_y` / `charge_budget_x` read the
per-channel charge straight out of the S1 product with the solver's own
conventions: the Y kernels are stored with y = 0 at ny//2 (`_to_y0_origin`),
rows alternate parity, and X kernels are indexed by absolute column mod the
40-pad superperiod. Three separate ad-hoc attempts to reproduce the sum rule on
2026-08-07 all got that wrong, which is why the reference is the solver's code
and not fresh indexing.

This docstring used to claim `charge_budget_*` was "the same code path whose
closed-form sum rule passes at 4.8e-07". IT IS NOT. The sum rule runs through
`sum_over_rows` / `sum_over_columns` (kernels.py), which is separate indexing
from `charge_budget_y` / `charge_budget_x`. So what this test certifies is
narrower, and worth stating precisely — the trust hierarchy is three-tier:

  1. CLOSED FORM. `check_sum_rule` certifies S(0)/C(0) = 0.875 exactly, but
     only for the FICTITIOUS pitch-sized (0.78 mm) pads that tile the plane.
     Production geometry is not this.
  2. REAL-PAD ANCHOR. The production 0.68 mm capture, 0.665023, against the
     semi-analytic (PAD_SIZE/PAD_PITCH)^2 x S(0)/C(0). Asserted in
     `test_charge_audit.py` — it was a comment until the audit.
  3. THIS TEST. Two INDEPENDENTLY WRITTEN indexings of the same product agree
     to 1e-4. That is a strong check on the LUT's transformations and a weak
     one on the product itself: a bug shared by both indexings survives it.

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

    # THE SOURCE POSITION MUST EXIST IN BOTH GRIDS. The LUT keeps every
    # x_stride-th x sample, so pad centres generally fall BETWEEN LUT samples —
    # evaluating the reference at a pad centre and the LUT at the nearest
    # sample compares two different source positions, up to 20 um apart. That
    # is not a small error here: the kernel's dependence on sub-pad position is
    # strong (the d=0 charge ranges 0.22-0.35 across columns), and a first
    # version of this test read 59 % "residual" that was entirely this.
    # lut.x == xfull[::x_stride], so every LUT sample is exactly a full-res
    # sample: pick sources FROM the LUT grid and index the reference at
    # x_stride times that index.
    stride = int(round((lut.x[1] - lut.x[0]) / (sy.x[1] - sy.x[0])))
    assert np.allclose(lut.x, sy.x[::stride][:len(lut.x)]), \
        "LUT x grid is not a subsample of the product's x grid"

    rng = np.random.default_rng(seed)
    worst = 0.0
    ok = True

    for _ in range(n_probe):
        c0 = int(rng.integers(0, C.N_PAD_PER_SUPER))
        # nearest LUT sample to this pad centre, then the identical full-res x
        ix = int(np.argmin(np.abs(lut.x - float(K.pad_x(c0)))))
        x0 = float(lut.x[ix])
        ix_ref = ix * stride
        # the reference's own channel assignment, at that exact x
        c0 = K.nearest_column(x0)
        par = 0                                  # deposit on an even row

        ref_y = K.charge_budget_y(sy, gy, x0, dmax=dmax, row0_parity=par)
        ref_x = K.charge_budget_x(sx, gx, x0, dmax=dmax, y0_m=0.0)
        assert int(np.argmin(np.abs(sy.x - np.mod(x0, sy.lx)))) == ix_ref, \
            "reference did not land on the same x sample as the LUT"

        # --- the same quantities out of the LUT ---------------------------
        # LUT holds CURRENT on a uniform grid, so charge = cumulative integral.
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


def certify_off_grid(path, n_probe=20000):
    """
    C2: channel booking must be self-consistent for sources drawn OFF the LUT
    grid — which is every real avalanche.

    `certify` above deliberately probes ON the LUT grid, because comparing
    per-channel charge at two different source positions is meaningless. That
    blinds it to the booking question, which is not about charge values at all:
    the LUT evaluates the kernel at the nearest stored x SAMPLE, so the column
    it is serving as d = 0 is the one nearest THAT sample. Book the charge
    against the column nearest the TRUE x instead — two roundings, onto a
    40 µm grid and a 780 µm grid whose boundaries do not coincide — and near
    every second pad boundary the two disagree by a whole pad.

    So the invariant, which is exact and needs no tolerance, is that the booked
    column is the one whose kernel the band actually holds:

        lut.col_at(x)  ==  lut.col_of_x[lut.ix(x)]   (mod the superperiod)

    and the pre-fix rule (round the true x onto the pad lattice) is reported
    alongside, to show its disagreement rate is real rather than hypothetical.

    NOT tested against `nearest_column` of the sample: at a pad boundary the
    two pad centres are exactly equidistant, and `argmin` takes the lower index
    while `rint` takes the even one. That is a tie-break difference on a source
    sitting exactly on the boundary — both answers are equally right, and 14 of
    80 boundary probes hit it. What must not differ is which kernel was used.
    """
    lut = CombKernelLUT(path)
    rng = np.random.default_rng(17)

    # Sources at pad boundaries +- 5 um, where the two conventions part, plus a
    # uniform sample for the overall rate.
    edges = (K.pad_x(np.arange(C.N_PAD_PER_SUPER)) + C.PAD_PITCH_M / 2)
    probes = np.concatenate([edges - 5e-6, edges + 5e-6,
                             rng.uniform(0, C.SUPERPERIOD_M, n_probe)])

    ix = lut.ix(probes)
    got = lut.col_at(probes, ix=ix)
    band = lut.col_of_x[ix]                       # the kernel actually served
    old = np.rint((probes - K.PAD_ORIGIN_M) / C.PAD_PITCH_M).astype(int)

    n_bad = int((((got - band) % C.N_PAD_PER_SUPER) != 0).sum())
    n_moved = int((((old - got) % C.N_PAD_PER_SUPER) != 0).sum())
    n_edge = 2 * len(edges)
    n_uni = int((((old - got) % C.N_PAD_PER_SUPER)[n_edge:] != 0).sum())

    print(f"\n  off-grid booking (C2), {len(probes)} probes "
          f"({n_edge} at pad boundaries +-5 um, {n_probe} uniform)")
    print(f"    col_at serves the band's own column: {n_bad} disagreements"
          f"   {'OK' if n_bad == 0 else 'FAIL'}")
    print(f"    the pre-fix rule (round the TRUE x) differs on "
          f"{n_moved}/{len(probes)} — {n_uni}/{n_probe} "
          f"({n_uni/n_probe:.1%}) of the UNIFORM probes, which is the rate at "
          f"which real avalanches were booked a pad off")
    return n_bad == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--n-probe", type=int, default=6)
    a = ap.parse_args()
    ok = certify(a.kernel, n_probe=a.n_probe)
    ok &= certify_off_grid(a.kernel)
    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
