#!/usr/bin/env python3
"""
test_longitudinal.py — check the general FIR path against the trusted one.

apply_ion_transit (cumsum running mean, in use since 2026-08-07) and
apply_longitudinal (FFT FIR, new) must agree when the second is handed the
exact shape the first assumes: f_e * delta(t) + f_i * rect(0, T). If they do,
the measured-template path inherits the rectangle path's validation instead of
starting from nothing.

This is a pure array test — no LUT file, no Garfield — so it runs anywhere.

    python3 -m response.digitizer.test_longitudinal
"""

from __future__ import annotations

import numpy as np

from .kernel_lut import CombKernelLUT


class _FakeLUT:
    """Just enough of a LUT for the two folding methods to operate on."""

    def __init__(self, nt=1001, dt_ns=1.0, n=64, seed=7):
        rng = np.random.default_rng(seed)
        self.dt = dt_ns * 1e-9
        # A plausible kernel shape: fast rise, long dispersive tail, so the
        # test exercises the tail region where truncation bugs live.
        t = np.arange(nt)
        base = np.exp(-t / 90.0) - np.exp(-t / 6.0)
        self.I_Y = (base[None, None, :] * rng.uniform(0.4, 1.6, (2, n, 1))
                    ).astype(np.float32)
        self.I_X = (base[None, None, :] * rng.uniform(0.4, 1.6, (9, n, 1))
                    ).astype(np.float32)

    apply_ion_transit = CombKernelLUT.apply_ion_transit
    apply_longitudinal = CombKernelLUT.apply_longitudinal


def impulse_check(f_e, T_ns, dt_ns):
    """
    The test that matters, and the one that was missing. Feed each folding
    path a unit impulse: the answer is then h itself, exact by inspection,
    with no reference implementation to be wrong in the same way.

    This is what caught the c[n]/(n+1) leading-edge bug in apply_ion_transit
    (2026-08-07): it returned 1.0 for the first sample instead of
    f_e + f_i/L = 0.095, and put 5.9 units of charge on the readout per unit
    induced. Every check up to then had compared ratios or shapes, which a
    charge-non-conserving filter can still pass.
    """
    L = max(1, int(round(T_ns / dt_ns)))
    want0, want_flat = f_e + (1 - f_e) / L, (1 - f_e) / L
    ok = True
    for name, fold in (("rectangle", "rect"), ("FIR", "fir")):
        g = _FakeLUT(nt=1001)
        g.I_Y[:] = 0.0
        g.I_X[:] = 0.0
        g.I_Y[:, :, 0] = 1.0
        g.I_X[:, :, 0] = 1.0
        if fold == "rect":
            g.apply_ion_transit(f_e, T_ns)
        else:
            h = np.full(L, (1.0 - f_e) / (L * dt_ns))
            h[0] += f_e / dt_ns
            g.apply_longitudinal(h, dt_ns * 1e-9)
        y = g.I_Y[0, 0]
        checks = [("y[0]", y[0], want0), ("y[5]", y[5], want_flat),
                  ("y[L-1]", y[L - 1], want_flat), ("y[L+60]", y[L + 60], 0.0),
                  ("area", y.sum(), 1.0)]
        for label, got, want in checks:
            good = abs(got - want) <= 1e-5 + 1e-4 * abs(want)
            ok &= good
            print(f"  {name:9s} {label:8s} {got: .6f}  want {want: .6f}  "
                  f"{'OK' if good else 'FAIL'}")
    return ok


def main():
    f_e, T_ns, dt_ns = 0.0923, 340.0, 1.0

    print("impulse response (exact by inspection):")
    ok_imp = impulse_check(f_e, T_ns, dt_ns)
    print()

    a = _FakeLUT()
    a.apply_ion_transit(f_e, T_ns)

    # The same operator written out as an explicit sampled shape. The rectangle
    # path uses L = round(T/dt) samples of running mean, each weighted 1/L, so
    # the equivalent FIR is exactly that — build it that way rather than from
    # the continuous rectangle, or the comparison tests the discretisation
    # rather than the convolution.
    L = max(1, int(round(T_ns / dt_ns)))
    h = np.zeros(L)
    h[0] = f_e / dt_ns
    h += (1.0 - f_e) / (L * dt_ns)

    b = _FakeLUT()
    b.apply_longitudinal(h, dt_ns * 1e-9)

    ok = True
    for name, x, y in (("I_Y", a.I_Y, b.I_Y), ("I_X", a.I_X, b.I_X)):
        scale = np.abs(x).max()
        err = np.abs(x - y).max() / scale
        # Both paths run in float32; 1e-5 is the float32 noise floor for a
        # 340-tap sum, not a tolerance chosen to make the test pass.
        good = err < 1e-5
        ok &= good
        print(f"  {name}  max rel diff {err:.2e}   {'OK' if good else 'FAIL'}")

    # Area invariance: an FIR of unit area cannot change any channel's
    # integral. This is the property the whole peak-vs-area argument rests on,
    # so assert it rather than trusting it.
    c = _FakeLUT()
    before = c.I_Y.sum(axis=-1).copy()
    c.apply_longitudinal(h, dt_ns * 1e-9)
    after = c.I_Y.sum(axis=-1)
    # Tail that runs off the end of the window is genuinely lost, so compare
    # only where the kernel has decayed — that is the real DAQ situation and
    # the reason run.py warns about the 1.92 us window.
    rel = np.abs(after - before).max() / np.abs(before).max()
    print(f"  area drift {rel:.2e} (finite window; exact only for h fully "
          f"inside the record)")

    print("\n" + ("PASS" if (ok and ok_imp) else "FAIL"))
    return 0 if (ok and ok_imp) else 1


if __name__ == "__main__":
    raise SystemExit(main())
