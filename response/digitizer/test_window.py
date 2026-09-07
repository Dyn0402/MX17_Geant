#!/usr/bin/env python3
"""
test_window.py — how much of the late halo did the 1 us LUT window throw away?

Audit A1 / Fix 1. `kernel_lut` defaulted `t_max_ns = 1000` while the DAQ area
window is 32 x 60 ns = 1.92 us. Sheet transport feeds the far neighbours on
exactly that timescale (D ~ 1 m^2/s at rho_s = 2 MOhm/sq, so ~0.6 us to reach
d = +-2 and past 1 us for d = +-3), so a 1 us kernel cannot deliver charge the
real detector delivers between 1.0 and 1.92 us. That biases the d = +-2/+-3
area shares LOW — the same direction as the headline §9 "missing broad slow
halo" tension, which is why it had to be measured before anything is concluded
about the halo.

This quantifies it directly, with no digitizer and no Geant4 in the way: build
the same product's LUT at several t_max, integrate each channel's current over
the window, and report the d-profile.

Two separate effects are reported and must not be conflated:

  * WINDOW  — how much more charge each channel has collected by the later
    time. Every channel grows; the halo grows fastest, which is the point.
  * SHARE   — the normalised d-profile. This is what plan §9 compares against
    data, and it is the number that was biased.

    python3 -m response.digitizer.test_window --kernel <s1.npz>
"""

from __future__ import annotations

import argparse

import numpy as np

from ..common import constants as C
from .kernel_lut import CombKernelLUT
from ..solver import kernels as K

# The windows worth comparing: the pre-fix default, the DAQ area window, the
# new default, and the SPS-config frame.
WINDOWS_NS = (1000.0, 1920.0, 3000.0, 4200.0)


def profile(path, t_max_ns, dmax=3, n_probe=12, seed=5, ion_ns=None):
    """Charge on d = -dmax..dmax, summed over probe positions, at one t_max."""
    lut = CombKernelLUT(path, t_max_ns=t_max_ns)
    if ion_ns:
        # Fold the analytic ion transit so the comparison includes the
        # longitudinal smear that pushes charge further toward the window edge.
        lut.apply_ion_transit(0.0923, ion_ns)
    rng = np.random.default_rng(seed)
    xs = rng.uniform(0.0, C.SUPERPERIOD_M, n_probe)

    accY = np.zeros(2 * dmax + 1)
    accX = np.zeros(2 * dmax + 1)
    for x0 in xs:
        ix = int(lut.ix(x0))
        for j, d in enumerate(range(-dmax, dmax + 1)):
            iy = lut.iy_Y(-d * C.PAD_PITCH_M)
            accY[j] += float(lut.I_Y[d % 2, iy, ix, :].sum()) * lut.dt
            jd = d + lut.n_side
            if 0 <= jd < lut.I_X.shape[0]:
                accX[j] += float(
                    lut.I_X[jd, lut.iy_X(0.0), ix, :].sum()) * lut.dt
    return accX / n_probe, accY / n_probe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--n-probe", type=int, default=12)
    ap.add_argument("--ion-ns", type=float, default=340.0,
                    help="fold the analytic ion transit (0 = bare kernel)")
    a = ap.parse_args()

    dmax = 3
    ds = list(range(-dmax, dmax + 1))
    print("LUT window vs the late halo (audit A1 / Fix 1)\n")
    print(f"  product {a.kernel.split('/')[-1]}, {a.n_probe} source positions, "
          f"ion transit {a.ion_ns:.0f} ns")
    print(f"  DAQ area window is 32 x 60 = 1920 ns; the pre-fix LUT stopped at "
          f"1000\n")

    rows = {}
    for t in WINDOWS_NS:
        rows[t] = profile(a.kernel, t, dmax=dmax, n_probe=a.n_probe,
                          ion_ns=(a.ion_ns or None))
        print(f"  built t_max = {t:6.0f} ns", flush=True)

    for tag, k in (("X", 0), ("Y", 1)):
        print(f"\n  {tag} view — CHARGE per channel (unit deposit)")
        print(f"  {'t_max':>8s}" + "".join(f"{d:>9d}" for d in ds)
              + f"{'total':>10s}")
        for t in WINDOWS_NS:
            q = rows[t][k]
            print(f"  {t:8.0f}" + "".join(f"{v:9.4f}" for v in q)
                  + f"{q.sum():10.4f}")

        print(f"\n  {tag} view — SHARE (what §9 compares), and the shift vs "
              f"the pre-fix 1000 ns")
        print(f"  {'t_max':>8s}" + "".join(f"{d:>9d}" for d in ds))
        base = rows[1000.0][k] / rows[1000.0][k].sum()
        for t in WINDOWS_NS:
            q = rows[t][k]
            s = q / q.sum()
            print(f"  {t:8.0f}" + "".join(f"{v:9.4f}" for v in s))
        s3 = rows[3000.0][k] / rows[3000.0][k].sum()
        d_abs = s3 - base
        print(f"  {'delta':>8s}" + "".join(f"{v:+9.4f}" for v in d_abs))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(np.abs(base) > 1e-9, d_abs / base, np.nan)
        print(f"  {'rel':>8s}" + "".join(
            ("      n/a" if not np.isfinite(v) else f"{v:+8.1%}") for v in rel))

    # The headline: c1 and c2, the two numbers plan §9 quotes.
    print("\n  track-relevant summary (mean of d = +-1 and d = +-2 shares)")
    print(f"  {'t_max':>8s} {'c1_X':>8s} {'c2_X':>8s} {'c1_Y':>8s} {'c2_Y':>8s}")
    for t in WINDOWS_NS:
        qx, qy = rows[t]
        sx, sy = qx / qx.sum(), qy / qy.sum()
        c1x = 0.5 * (sx[dmax - 1] + sx[dmax + 1])
        c2x = 0.5 * (sx[dmax - 2] + sx[dmax + 2])
        c1y = 0.5 * (sy[dmax - 1] + sy[dmax + 1])
        c2y = 0.5 * (sy[dmax - 2] + sy[dmax + 2])
        print(f"  {t:8.0f} {c1x:8.4f} {c2x:8.4f} {c1y:8.4f} {c2y:8.4f}")

    print("\n  This is a PREDICTION REVISION, not a tuning step: no §9 number "
          "was\n  consulted in choosing the window (plan §9 firewall).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
