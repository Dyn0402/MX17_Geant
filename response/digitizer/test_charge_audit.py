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

TOL_LEAK = 0.02          # >2 % lost outside the channel window is a problem

# ── STATUS 2026-08-07: THIS TEST FAILS, AND THE CAUSE IS NOT YET KNOWN ───────
#
# Measured: the digitizer's total induced charge is 0.675 of the deposited
# charge, against the S1 closed form's 0.875. Three candidate explanations were
# tested and ELIMINATED:
#
#   drift depth       leak is 32.3 % at z = 0.5 mm and 33.1 % at z = 29.5 mm
#                     while sigma_T goes 107 -> 819 um. An 8x change in the
#                     transverse spread moves it by 0.8 %, so it is not
#                     diffusion pushing charge out of the window.
#   channel window    n_side 4 -> 16 with y_window grown to match (the two must
#                     grow together: Y channels sit on the 0.78 mm pad pitch,
#                     so a 3.9 mm window holds only +-5 of them however large
#                     n_side is) CONVERGES at 0.6753 by n_side = 8. More
#                     channels cannot recover it.
#   kernel time cut   opening t_max 1000 -> 4000 ns makes it WORSE at fixed
#                     n_side (0.672 -> 0.625), which is the expected behaviour
#                     of charge dispersing into more channels, not a time cut.
#
# So ~23 % is unaccounted for. NOT yet established as a digitizer bug: the
# solver's own check_sum_rule passes on this product at 4.8e-07, and three
# separate hand-rolled attempts to reproduce that sum from the raw npz were all
# wrong (the verified summation is kernels.sum_over_rows / sum_over_columns,
# which apply a _to_y0_origin re-origining that ad-hoc indexing misses).
#
# NEXT STEP, and it is the proper T10 fast-path certification rather than more
# probing: compare the LUT lookup channel-by-channel against
# kernels.charge_budget_y / sum_over_rows at the same source position, so the
# reference is the code that already passes its own closed-form check.
#
# WHAT IT DOES AND DOES NOT INVALIDATE. Share observables (c1, peak ratios) are
# ratios within the window and are unaffected if the loss is uniform across
# channels -- which the window-convergence above suggests but does not prove.
# ABSOLUTE charge and therefore ADC scale ARE affected, and that is entangled
# with the separate observation that simulated MIPs saturate the 12-bit ADC.


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

    dig = Digitizer(a.kernel, os.path.expanduser(a.calib), seed=11)
    print("Charge audit — induced total vs Ramo + the S1 sum rule\n")
    print(f"  kernel  rho_s={dig.lut.describe()['rho_s_MOhm_sq']} MΩ/sq  "
          f"n_side={dig.n_side} (so {2*dig.n_side+1} channels per view)")
    print(f"  sigma_T at the ESL grows with drift depth, so the window leak "
          f"must grow with z\n")
    print(f"  {'z [mm]':>7s} {'sigma_T [µm]':>13s} {'induced/deposited':>19s} "
          f"{'leak':>8s}")

    ok = True
    for z_mm in (0.5, 5.0, 15.0, 29.5):
        got, want = audit(dig, z_mm, n_rep=a.n_rep)
        if not len(got):
            continue
        ratio = float(np.mean(got / want))
        leak = 1.0 - ratio
        # DriftGas returns arrays for array-like fields; take a scalar.
        try:
            sig = float(np.ravel(dig.gas.sigma_T_um(dig.E_drift, z_mm))[0])
        except Exception:
            sig = float("nan")
        good = leak <= TOL_LEAK
        ok &= good
        print(f"  {z_mm:7.1f} {sig:13.0f} {ratio:19.4f} {leak:8.2%}  "
              f"{'OK' if good else 'LEAK'}")

    print(f"\n  tolerance: leak <= {TOL_LEAK:.0%}")
    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
