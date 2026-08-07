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
