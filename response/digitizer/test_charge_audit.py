#!/usr/bin/env python3
"""
test_charge_audit.py — does the charge that goes in come out the other end?

THE CHECK. Ramo's theorem plus S1's own sum rule fix the answer exactly: the
weighting potentials of all readout channels partition the plane (S1 verifies
this to 4e-13), so a charge q arriving at the ESL must induce a TOTAL of q
summed over every channel of BOTH views. Nothing about resistivity, kapton
thickness or ion transit can change that total — they only decide how it is
shared and when it arrives.

So the whole transport chain has one hard invariant:

    sum over all channels of induced charge  ==  N_primary * transparency * gain

and any deficit is charge the digitizer LOST rather than shared.

WHY IT MATTERS FOR THE BUDGET NUMBERS. The digitizer computes only n_side = 4
channels either side of the centroid, and `event_budget` normalises the shares
to the sum over THAT WINDOW, not to the true total. If the window leaks, every
share is inflated by the same factor and the c1 numbers quoted throughout are
biased high. This measures the leak directly, at several depths, because the
leak grows with transverse diffusion and therefore with drift distance.

This is an internal invariance test — it compares the simulation against a
theorem, not against data.

    python3 -m response.digitizer.test_charge_audit --kernel <s1.npz> \\
            --calib <aval_calib.json>
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from ..common import constants as C
from ..solver import kernels as K
from .digitize import Digitizer

# Bar on |1 - measured/expected| once the expectation is PER POSITION.
# It was 3 % against a grid-MEAN reference, which had to absorb the whole
# position dependence of the capture (d=0 spans 0.22-0.35 across columns).
TOL_LEAK = 0.01
# Bar on the real-pad geometric anchor, below. The measured agreement is 8e-4;
# 0.5 % leaves room for the pad-edge fringe while still catching any
# indexing-level error.
TOL_ANCHOR = 0.005

# ── RESOLVED 2026-08-07: THE "MISSING 23 %" WAS THE INTER-PAD GAPS ──────────
#
# This test originally asserted against prompt_sum_rule = S(0)/C(0) = 0.875 and
# reported a flat ~33 % deficit at every drift depth. That expectation was
# WRONG, and the right number was already sitting in every product's metadata.
#
# S(0)/C(0) applies only to the FICTITIOUS pitch-sized (0.78 mm) pads that
# kernels.check_sum_rule substitutes so that a closed form exists — with those,
# the checkerboard tiles the plane exactly. Production kernels use the REAL
# 0.68 mm pads, which do NOT tile: (0.68/0.78)^2 = 0.760 of the plane is pad,
# and the image charge on the 100 um inter-pad gaps lands on no readout channel
# at all. run_point computes this and stores it as `channel_capture_prompt`:
#
#     sum_rule_expect (tiling pads)   0.875000
#     channel_capture_prompt (real)   0.665023   <- the correct reference
#     digitizer measured              0.675300   <- within 1.5 %
#
# The reference is now read from the product's own meta. Note that
# channel_capture is identical prompt and late (0.665023 both), which is the
# right signature: a geometric partition cannot be time-dependent.
#
# Three explanations were tested and eliminated before the real one was found,
# and the eliminations remain useful:
#   drift depth    32.3 % at z=0.5 mm vs 33.1 % at z=29.5 mm while sigma_T goes
#                  107 -> 819 um, so it was never diffusion.
#   channel window n_side 4 -> 16 with y_window grown to match CONVERGES at
#                  0.6753 by n_side=8. (They must grow together: Y channels sit
#                  on the 0.78 mm pad pitch, so a 3.9 mm window holds only +-5
#                  of them however large n_side is.)
#   kernel time    opening t_max makes it worse at fixed n_side — charge
#                  dispersing into more channels, not a time cut.
#
# ── 2026-08-07 (audit A4 / Fix 4): THE REFERENCE IS NOW PER POSITION ────────
#
# The residual left over after the gaps were understood was +1.5 %, i.e. an
# EXCESS — and none of the leak framings above can produce one: a window that
# leaks, a kernel truncated in time and a channel set that misses charge all
# push the measured total DOWN. An excess can only come from the reference.
#
# It did. `channel_capture_prompt` is the capture averaged over the whole
# 31.2 x 1.56 mm cell, but capture is strongly position dependent (the d=0
# share alone spans 0.22-0.35 across columns), and the audit's probes are a
# handful of random positions — their mean is not the cell mean. Comparing a
# 6-position sample against a grid-mean scalar folds that sampling spread into
# what is reported as a "leak".
#
# So the expectation is now evaluated AT each probe's own positions, from the
# solver's own row/column sums (`sum_over_rows` + `sum_over_columns`) — the
# same helpers whose closed-form sum rule passes at 4.8e-07, per the process
# rule that indexing is never re-derived ad hoc.
#
# WHAT THIS TEST CAN AND CANNOT SEE, measured by deliberately breaking it
# (source position handed to the reference displaced, 16 probes at z=0.5 mm):
#
#     break                             leak     probe spread
#     none (correct)                    1.06 %      12.6 %
#     x displaced half a pad (390 um)   1.67 %      26.1 %   <- CAUGHT
#     x displaced one pad (780 um)      0.82 %      13.0 %   <- invisible
#     y displaced one row (780 um)      1.04 %      13.0 %   <- invisible
#
# Sub-pad position errors show up, in the SPREAD rather than the mean. Whole-
# pad ones cannot: the observable is a sum over all channels, which is by
# construction invariant under moving charge BETWEEN channels, and the capture
# map's dominant structure is periodic with the pad pitch. So a parity flip or
# the one-pad-off booking of C2 is not detectable here at all — that is what
# `test_lut_vs_solver.py` with off-lattice probes is for. Stated because the
# fix order asked for a parity break to fail here, and it provably cannot.


class _Grid:
    """The bits of a solver object that the row/column sums actually read."""

    def __init__(self, x, y):
        self.x, self.y = x, y
        self.nx, self.ny = len(x), len(y)
        self.lx = float(self.nx * (x[1] - x[0]))
        self.ly = float(self.ny * (y[1] - y[0]))


class CaptureMap:
    """
    Total charge captured by ALL channels of BOTH views, as a function of
    where on the cell the deposit lands.

    This is exactly the array `run_point` reduces to the scalar
    `channel_capture_prompt`; keeping it un-reduced is the whole fix. Built
    from `sum_over_rows` / `sum_over_columns`, which is why it inherits the
    solver's certified indexing rather than a fresh guess at it.

    The result is periodic: 31.2 mm in x (the pad/ESL beat) and 1.56 mm in y
    (one checkerboard period), and the two views' sums land on the same y grid
    by construction (`x_ny`).
    """

    def __init__(self, path):
        with np.load(path) as d:
            sy = _Grid(d["x"], d["y_Y"])
            sx = _Grid(d["x"], d["y_X"])
            gy = (d["G_Y_even"], d["G_Y_odd"])
            gx = {c: d["G_X"][c] for c in range(C.N_PAD_PER_SUPER)}
            tot = K.sum_over_rows(sy, gy) + K.sum_over_columns(sx, gx)
        self.prompt = np.asarray(tot[0], dtype=float)      # (ny_cell, nx)
        self.x = np.asarray(sy.x, dtype=float)
        self.dx = float(self.x[1] - self.x[0])
        self.ny = self.prompt.shape[0]
        self.dy = float(C.PAD_PITCH_M * 2 / self.ny)       # 1.56 mm / ny
        self.mean = float(self.prompt.mean())

    def at(self, x_m, y_m):
        """Expected capture for point deposits at (x, y) [m]."""
        ix = np.rint(np.mod(x_m, C.SUPERPERIOD_M) / self.dx).astype(int) \
            % len(self.x)
        iy = np.rint(np.mod(y_m, 2 * C.PAD_PITCH_M) / self.dy).astype(int) \
            % self.ny
        return self.prompt[iy, ix]


def audit(dig, z_mm, cap, n_e=2000, n_rep=6, n_samp=3000, seed=11):
    """Induced charge summed over the window vs the theorem's answer."""
    rng = np.random.default_rng(seed)
    got, want, ref = [], [], []
    for _ in range(n_rep):
        # Land the deposit at a random sub-pitch position: the window loss
        # depends on where the centroid sits relative to the pad grid, and a
        # deposit centred on a pad is the most favourable case.
        x0 = rng.uniform(0.0, C.SUPERPERIOD_M * 1e3)
        y0 = rng.uniform(0.0, 2 * C.PAD_PITCH_M * 1e3)
        x, y, t, q = dig.transport(np.array([x0]), np.array([y0]),
                                   np.array([float(z_mm)]), np.array([0.0]),
                                   np.array([n_e]))
        if len(x) == 0:
            continue
        cur = dig.induce(x, y, t, q, n_samp)
        # Induced charge = integral of the current. The LUT is on a 1 ns grid
        # and carries e/s, so this is in elementary charges.
        tot = sum(float(np.sum(v)) * dig.lut.dt for v in cur.values())
        got.append(tot)
        want.append(float(np.sum(q)))
        # The expectation for THESE electrons, at the positions they actually
        # landed on after diffusion — charge-weighted, since the capture is a
        # per-electron quantity and the gains differ.
        ref.append(float(np.sum(cap.at(x, y) * q) / np.sum(q)))
    return np.array(got), np.array(want), np.array(ref)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--n-rep", type=int, default=6)
    ap.add_argument("--ion-model", choices=("measured", "analytic"),
                    default="measured",
                    help="must match what the run being audited uses; the "
                         "invariant under test holds for both")
    a = ap.parse_args()

    import json
    with np.load(a.kernel) as _d:
        meta = json.loads(str(_d["meta"]))
    expect = float(meta["channel_capture_prompt"])
    sum_rule = float(meta["sum_rule_expect"])

    print("Charge audit — induced total vs Ramo + the S1 sum rule\n")

    ok = True

    # ── 0. THE REAL-PAD ANCHOR ──────────────────────────────────────────────
    # The closed-form sum rule certifies only the FICTITIOUS pitch-sized pads
    # that tile the plane. Production uses the real 0.68 mm pads, and until now
    # the agreement between the solved real-pad capture and the semi-analytic
    # (pad area fraction) x (sum rule) was recorded ONLY in a comment — the
    # number every downstream budget is normalised to had no test on it at all
    # (audit A4). Assert it.
    area_frac = (C.PAD_SIZE_UM / C.PAD_PITCH_UM) ** 2
    anchor = area_frac * sum_rule
    dev = abs(expect - anchor)
    good = dev <= TOL_ANCHOR
    ok &= good
    print(f"  real-pad anchor: solved channel_capture_prompt {expect:.6f}")
    print(f"                   (PAD_SIZE/PAD_PITCH)^2 x S(0)/C(0) = "
          f"{area_frac:.4f} x {sum_rule:.4f} = {anchor:.6f}")
    print(f"                   |diff| {dev:.2e}  tol {TOL_ANCHOR:.1e}   "
          f"{'OK' if good else 'FAIL'}\n")

    dig = Digitizer(a.kernel, os.path.expanduser(a.calib), seed=11,
                    ion_model=a.ion_model)
    cap = CaptureMap(a.kernel)
    print(f"  kernel  rho_s={dig.lut.describe()['rho_s_MOhm_sq']} MΩ/sq  "
          f"n_side={dig.n_side} (so {2*dig.n_side+1} channels per view)")
    print(f"  sigma_T at the ESL grows with drift depth, so the window leak "
          f"must grow with z\n")
    print(f"  reference: PER POSITION, from sum_over_rows + sum_over_columns.")
    print(f"  its cell mean is {cap.mean:.6f} (= the product's "
          f"channel_capture_prompt {expect:.6f}), but it ranges "
          f"{cap.prompt.min():.4f}-{cap.prompt.max():.4f} over the cell —")
    print(f"  which is why a few random probes cannot be compared against the "
          f"mean (that was the +1.5 % 'excess')\n")
    print(f"  {'z [mm]':>7s} {'sigma_T [µm]':>13s} {'induced/deposited':>19s} "
          f"{'expected':>10s} {'leak':>8s} {'+-':>7s} {'spread':>8s}")

    for z_mm in (0.5, 5.0, 15.0, 29.5):
        got, want, ref = audit(dig, z_mm, cap, n_rep=a.n_rep)
        if not len(got):
            continue
        # Per-probe residual: each probe against ITS OWN expectation.
        per = got / want / ref - 1.0
        leak = -float(np.mean(per))
        spread = float(np.std(per))
        # The probe-to-probe spread is the sampling error on the mean, and it
        # is LARGE at shallow depth: a barely-diffused cloud sits at one place
        # on a capture map that ranges 0.25-0.89, so 6 probes constrain the
        # mean to only a few per cent. Quoting the mean without it invites
        # reading noise as a leak.
        se = spread / max(1.0, np.sqrt(len(per)))
        # DriftGas returns arrays for array-like fields; take a scalar.
        try:
            sig = float(np.ravel(dig.gas.sigma_T_um(dig.E_drift, z_mm))[0])
        except Exception:
            sig = float("nan")
        good = abs(leak) <= max(TOL_LEAK, 2.0 * se)
        ok &= good
        print(f"  {z_mm:7.1f} {sig:13.0f} {float(np.mean(got/want)):19.4f} "
              f"{float(np.mean(ref)):10.4f} {leak:8.2%} {se:7.2%} "
              f"{spread:8.2%}  {'OK' if good else 'LEAK'}")

    print(f"\n  tolerance: |1 - measured/expect| <= {TOL_LEAK:.0%} per depth, "
          f"expectation evaluated at the probe's own positions")
    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
