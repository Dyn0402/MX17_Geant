#!/usr/bin/env python3
"""
test_ny_grid.py — is the S1 transverse (y) grid fine enough? (audit C6)

THE APPROXIMATION UNDER TEST. Production S1 products are solved at ny = 512
over the 49.92 mm Y box, i.e. 97.5 um samples, and a readout pad is 680 um wide
on a 780 um pitch. The PROMPT kernel is the sharpest thing in the product — at
t = 0 the sheet is a plain dielectric and the response is the bare
electrostatic image, so its pad-edge shoulder is a near-step in y. A near-step
sampled at 97.5 um is exactly where a grid is least converged, and nothing had
ever measured it: `test_time_grid` certifies the TIME axis, `test_lut_vs_solver`
compares two indexings of the SAME product and so cannot see a grid error at
all.

THE TEST, mirroring test_time_grid: re-solve the nominal point at ny = 1024
(one doubling, ~2x cost) and difference the prompt kernels.

WHY THE ANSWER IS REPORTED TWICE. The raw pad-edge difference is not what the
simulation delivers. Every avalanche arrives smeared by transverse diffusion —
at least ~150 um sigma even for a cluster born at the mesh, and 300-800 um over
a real drift — and that smear is applied AFTER the kernel lookup. So:

  * RAW      the worst-case per-sample difference, which lands on the pad edge;
  * SMEARED  the same difference after a Gaussian of the given sigma, which is
             what any real avalanche actually sees.

The certificate the plan needs is the smeared one: if it is <= 0.3 %, ny = 512
stays with this run recorded as its evidence; if not, production needs ny = 1024.

    python3 -m response.digitizer.test_ny_grid <ny512.npz> <ny1024.npz>
"""

from __future__ import annotations

import argparse
import json

import numpy as np

TOL_SMEARED = 0.003        # 0.3 %, the bar the plan asks for
SMEAR_UM = (150.0, 300.0, 800.0)


def _prompt(path):
    """Prompt (t = 0) Y kernels and the y axis they live on."""
    with np.load(path) as d:
        meta = json.loads(str(d["meta"]))
        y = d["y_Y"].copy()
        g = np.stack([d["G_Y_even"][0].astype(np.float64),
                      d["G_Y_odd"][0].astype(np.float64)])   # (par, ny, nx)
    return y, g, meta


def _smear(g, y, sigma_m):
    """Gaussian smear along y, normalised, done by FFT on the periodic axis."""
    dy = float(y[1] - y[0])
    n = g.shape[1]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dy)
    H = np.exp(-0.5 * (k * sigma_m) ** 2)
    return np.fft.ifft(np.fft.fft(g, axis=1) * H[None, :, None], axis=1).real


def compare(path_coarse, path_dense):
    y_c, g_c, m_c = _prompt(path_coarse)
    y_d, g_d, m_d = _prompt(path_dense)

    if g_c.shape[2] != g_d.shape[2]:
        raise SystemExit(f"nx differs ({g_c.shape[2]} vs {g_d.shape[2]}); the "
                         "two products must share nx so only ny is compared")
    if g_d.shape[1] % g_c.shape[1]:
        raise SystemExit(f"ny_dense={g_d.shape[1]} is not a multiple of "
                         f"ny_coarse={g_c.shape[1]}")
    step = g_d.shape[1] // g_c.shape[1]

    print("  coarse  ny = %4d  (%.1f um)   %s" %
          (g_c.shape[1], (y_c[1] - y_c[0]) * 1e6, path_coarse.split('/')[-1]))
    print("  dense   ny = %4d  (%.1f um)   %s" %
          (g_d.shape[1], (y_d[1] - y_d[0]) * 1e6, path_dense.split('/')[-1]))
    print("  both at rho_s = %.1f MOhm/sq, d_k = %.0f um\n" %
          (m_c["rho_s_ohm_sq"] / 1e6, m_c["d_kapton_m"] * 1e6))

    # Sub-sample the dense solve onto the coarse grid. Every coarse y sample is
    # exactly a dense sample (both boxes are the same length and ny doubles),
    # so this compares the SAME physical positions and no interpolation enters.
    sub = g_d[:, ::step, :]
    assert np.allclose(y_c, y_d[::step]), "y grids are not nested"

    scale = float(np.abs(g_d).max())
    raw = float(np.abs(sub - g_c).max())
    print(f"  prompt kernel peak value            {scale:.4f}")
    print(f"  RAW max |ny512 - ny1024|            {raw:.4f}  "
          f"({raw/scale:.2%} of peak)")
    print(f"  (this lands on the pad edge, where a near-step is sampled)\n")

    ok = True
    print(f"  {'smear sigma':>12s} {'max |diff|':>12s} {'vs peak':>10s} "
          f"{'tol':>8s}")
    for s_um in SMEAR_UM:
        a = _smear(g_c, y_c, s_um * 1e-6)
        b = _smear(sub, y_c, s_um * 1e-6)
        d = float(np.abs(a - b).max())
        rel = d / scale
        good = rel <= TOL_SMEARED
        # Only the smallest smear is a real requirement: it is the floor any
        # avalanche gets. Larger ones are reported to show the trend.
        if s_um == min(SMEAR_UM):
            ok &= good
        print(f"  {s_um:10.0f} um {d:12.5f} {rel:10.3%} {TOL_SMEARED:8.1%}  "
              f"{'OK' if good else 'FAIL'}")

    print(f"\n  certificate applies to the MINIMUM smear ({min(SMEAR_UM):.0f} "
          f"um), the floor every avalanche gets.")
    print("  " + ("ny = 512 is adequate; keep it, with this run as the evidence."
                  if ok else
                  "ny = 512 is NOT adequate; production needs ny = 1024."))
    print("\n" + ("PASS" if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("coarse")
    ap.add_argument("dense")
    a = ap.parse_args()
    return 0 if compare(a.coarse, a.dense) else 1


if __name__ == "__main__":
    raise SystemExit(main())
