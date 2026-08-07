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
from .digitize import Digitizer

TOL_LEAK = 0.03          # vs the product's own channel_capture prediction

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

def audit(dig, z_mm, n_e=2000, n_rep=6, n_samp=3000, seed=11):
    """Induced charge summed over the window vs the theorem's answer."""
    rng = np.random.default_rng(seed)
    got, want = [], []
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
    return np.array(got), np.array(want)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--n-rep", type=int, default=6)
    a = ap.parse_args()

    import json
    with np.load(a.kernel) as _d:
        meta = json.loads(str(_d["meta"]))
    expect = float(meta["channel_capture_prompt"])
    dig = Digitizer(a.kernel, os.path.expanduser(a.calib), seed=11)
    print("Charge audit — induced total vs Ramo + the S1 sum rule\n")
    print(f"  kernel  rho_s={dig.lut.describe()['rho_s_MOhm_sq']} MΩ/sq  "
          f"n_side={dig.n_side} (so {2*dig.n_side+1} channels per view)")
    print(f"  sigma_T at the ESL grows with drift depth, so the window leak "
          f"must grow with z\n")
    print(f"  reference: this product's channel_capture_prompt = {expect:.6f}")
    print(f"  (the real 0.68 mm pads cover (0.68/0.78)^2 = "
          f"{(0.68/0.78)**2:.3f} of the plane; the rest is inter-pad gap)\n")
    print(f"  {'z [mm]':>7s} {'sigma_T [µm]':>13s} {'induced/deposited':>19s} "
          f"{'vs expect':>10s}")

    ok = True
    for z_mm in (0.5, 5.0, 15.0, 29.5):
        got, want = audit(dig, z_mm, n_rep=a.n_rep)
        if not len(got):
            continue
        ratio = float(np.mean(got / want))
        leak = 1.0 - ratio / expect
        # DriftGas returns arrays for array-like fields; take a scalar.
        try:
            sig = float(np.ravel(dig.gas.sigma_T_um(dig.E_drift, z_mm))[0])
        except Exception:
            sig = float("nan")
        good = abs(leak) <= TOL_LEAK
        ok &= good
        print(f"  {z_mm:7.1f} {sig:13.0f} {ratio:19.4f} {leak:10.2%}  "
              f"{'OK' if good else 'LEAK'}")

    print(f"\n  tolerance: |1 - measured/expect| <= {TOL_LEAK:.0%}")
    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
