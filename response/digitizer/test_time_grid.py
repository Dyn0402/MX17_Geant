#!/usr/bin/env python3
"""
test_time_grid.py — is the S1 log time grid fine enough for the fast path?

THE APPROXIMATION UNDER TEST. S1 products are solved on a 61-point LOG time
axis (`kernels.log_times`, 0.1 ns to 10 us). `CombKernelLUT._resample`
interpolates the weighting response G LINEARLY IN t onto a uniform 1 ns grid,
and `_to_current` then differentiates. Linear interpolation of G means the
induced CURRENT is piecewise CONSTANT between source samples — a staircase.

Measured on the production rho_s = 2 MOhm/sq, d_k = 75 um product: the 1001-
sample current series takes only **79 distinct values**, with steps up to
**114 ns** wide. The DREAM peaking time is 180 ns. A staircase whose step is
comparable to the shaping constant is exactly the regime where "the shaper
smooths it out" stops being obviously true, so it is measured rather than
assumed.

THE TEST. Re-solve one point with a 4x denser time grid (--nt 240) and compare
the fast path built from each, after the DREAM shaper, since the shaped
waveform is what the DAQ samples. The plan's T10 acceptance is a < 2 % waveform
residual, so that is the bar used here.

If this fails, the fix is not a smarter interpolator: it is to re-solve the S1
grid with more time points, because the information is genuinely absent from
the coarse product.

    python3 -m response.digitizer.test_time_grid <coarse.npz> <dense.npz>
"""

from __future__ import annotations

import argparse

import numpy as np

from .kernel_lut import CombKernelLUT
from ..dream.shaper import DreamShaper

TOL_RESIDUAL = 0.02        # plan T10: < 2 % waveform residual
TOL_PEAK = 0.02
TOL_PEAK_TIME_NS = 5.0


def staircase_report(lut, label):
    """How many distinct current levels the 1 ns series actually carries."""
    iy = len(lut.y_Y) // 2
    ix = lut.I_Y.shape[2] // 2
    w = lut.I_Y[0, iy, ix, :]
    nz = np.nonzero(np.abs(np.diff(w)) > 1e-12)[0]
    gaps = np.diff(nz) if len(nz) > 1 else np.array([0])
    print(f"  {label:8s} {len(w):5d} samples, "
          f"{len(np.unique(np.round(w, 12))):4d} distinct levels, "
          f"max step {gaps.max():4d} ns")
    return w


def compare(coarse_path, dense_path, n_probe=24, seed=5):
    a = CombKernelLUT(coarse_path)
    b = CombKernelLUT(dense_path)
    if a.I_Y.shape != b.I_Y.shape or a.I_X.shape != b.I_X.shape:
        raise SystemExit(f"LUT shapes differ ({a.I_Y.shape} vs {b.I_Y.shape}) "
                         "— the two products must share nx/ny, only nt may "
                         "differ, or this compares spatial resolution instead")

    print("  time-grid staircase, central Y channel:")
    staircase_report(a, "coarse")
    staircase_report(b, "dense")
    print()

    sh = DreamShaper()
    rng = np.random.default_rng(seed)
    ok = True
    rows = []

    # Probe a spread of channel offsets and source positions rather than one
    # convenient point: the staircase is worst where the kernel still has
    # structure late, which is the NEIGHBOUR channels, not the central one.
    for arr_name in ("I_Y", "I_X"):
        A, B = getattr(a, arr_name), getattr(b, arr_name)
        n0, n1, n2 = A.shape[0], A.shape[1], A.shape[2]
        for _ in range(n_probe):
            i0 = int(rng.integers(0, n0))
            i1 = int(rng.integers(max(0, n1 // 2 - 6), min(n1, n1 // 2 + 7)))
            i2 = int(rng.integers(0, n2))
            wa = sh.apply(A[i0, i1, i2, :].astype(float), dt_ns=a.dt * 1e9)
            wb = sh.apply(B[i0, i1, i2, :].astype(float), dt_ns=b.dt * 1e9)
            scale = np.abs(wb).max()
            if scale <= 0:
                continue
            res = float(np.sqrt(np.mean((wa - wb) ** 2)) / scale)
            dpk = float(abs(np.abs(wa).max() - scale) / scale)
            dtp = float(abs(np.argmax(wa) - np.argmax(wb)) * a.dt * 1e9)
            rows.append((arr_name, res, dpk, dtp))

    for name in ("I_Y", "I_X"):
        r = [x for x in rows if x[0] == name]
        if not r:
            continue
        res = np.array([x[1] for x in r])
        dpk = np.array([x[2] for x in r])
        dtp = np.array([x[3] for x in r])
        for label, v, tol, unit in (
                ("rms residual", res, TOL_RESIDUAL, ""),
                ("peak error", dpk, TOL_PEAK, ""),
                ("peak-time shift", dtp, TOL_PEAK_TIME_NS, " ns")):
            good = v.max() <= tol
            ok &= good
            print(f"  {name}  {label:16s} median {v.mean():8.4f}{unit}  "
                  f"max {v.max():8.4f}{unit}  tol {tol}{unit}  "
                  f"{'OK' if good else 'FAIL'}")
        print()

    print("PASS" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("coarse")
    ap.add_argument("dense")
    ap.add_argument("--n-probe", type=int, default=24)
    a = ap.parse_args()
    return 0 if compare(a.coarse, a.dense, a.n_probe) else 1


if __name__ == "__main__":
    raise SystemExit(main())
