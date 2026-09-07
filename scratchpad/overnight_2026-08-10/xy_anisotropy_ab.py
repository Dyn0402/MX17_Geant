#!/usr/bin/env python3
"""Is the sim's X/Y amplitude asymmetry a Y-path BUG or the strip model itself?

The data says the detector is X/Y symmetric to 0.02 % (unselected amplitudes
2606.5 vs 2607.0 ADC). The sim is not: the frozen T14 default reads Y/X = 0.79
on peak amplitude, and RHO_S_SENSITIVITY_2026-08-09.md found the asymmetry
shrinks with rho_s and vanishes near 5 MOhm/sq -- out of the T2b band -- and
closed with "the Y response MODEL is the more likely culprit than the rho_s
NUMBER". This tests that.

The hypothesis being tested
---------------------------
The asymmetry is not a defect in the Y code path. It is forced by the sheet
geometry, and any correct implementation of that geometry would show it:

  * The ESL is 550 um strips on an 800 um pitch, patterned in x and UNIFORM in
    y (wpot.py: "y is uniform in the model, so its modes decouple completely").
    So charge spreads freely ALONG y and is blocked ACROSS x by the 250 um
    gaps.
  * An X channel is a pad COLUMN -- a comb running along y, parallel to the
    strips. Charge diffusing along the strip stays over the same column, so
    the X channel keeps it.
  * A Y channel is a pad ROW -- running along x, perpendicular to the strips.
    Charge diffusing along the strip walks straight off it onto neighbouring
    rows, so the Y channel loses it.

If that is the mechanism, then replacing the patterned sheet with a UNIFORM
sheet of the same rho_s -- the only change -- must collapse the asymmetry.

Falsifier, stated before running: if the uniform sheet leaves Y/X essentially
unchanged (say it stays below 0.90), the anisotropy is NOT the cause and the
asymmetry really does live somewhere in the Y code path, which is then where
to look next.

Both arms use the same solver, the same box, the same grid, the same pad
patterns and the same shaper; the ONLY difference is sigma_s(x). The uniform
arm uses solve_uniform_analytic (the V1 closed form, exact for a uniform
sheet); the strip arm uses the full Bloch solve. Absolute values at this
reduced grid are not production numbers -- the strip-vs-uniform CONTRAST is
the result, and the discretisation is common to both arms.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np

from response.common import constants as C
from response.solver import kernels as K
from response.solver.wpot import WeightingSolver
from response.dream.shaper import DreamShaper

DT_NS = 2.0
T_MAX_NS = 2400.0
DMAX = 6
RHO_S = 2.0e6
D_K = 50e-6
F_ELECTRON = 0.092
TRANSIT_NS = 340.0


def make_solver(ly, nx, ny, uniform):
    return WeightingSolver(rho_s_ohm_sq=RHO_S, d_kapton_m=D_K, nx=nx, ny=ny,
                           ly_m=ly, phase_m=K.ESL_PHASE_CENTERED_M,
                           tau_drain_s=None, uniform_sheet=uniform)


def run_solve(s, pattern, times_s, uniform):
    return (s.solve_uniform_analytic(pattern, times_s) if uniform
            else s.solve(pattern, times_s))


def ion_fold(q_t, t_s, f_e, transit_s):
    """Induced charge for a realistic arrival: f_e prompt + (1-f_e) ramp.

    The kernel is the response to unit charge present since t=0. The electron
    fraction arrives promptly; the ion fraction is delivered linearly over the
    transit. Superposing shifted copies of the kernel gives the ion arm.
    """
    out = f_e * q_t
    n = 24
    for k in range(n):
        dt = (k + 0.5) / n * transit_s
        out = out + (1.0 - f_e) / n * np.interp(t_s - dt, t_s, q_t, left=0.0)
    return out


def shaped_peaks(budget, t_s, shaper, t_grid_ns, use_ions):
    """Shaped peak amplitude per channel offset."""
    peaks = {}
    for d, q in budget.items():
        q = np.asarray(q, dtype=float)
        if use_ions:
            q = ion_fold(q, t_s, F_ELECTRON, TRANSIT_NS * 1e-9)
        # kernel is on a log time axis -> resample to the shaper's uniform grid
        q_u = np.interp(t_grid_ns, t_s * 1e9, q, left=0.0, right=q[-1])
        # I = dq/dt, and q(0) != 0: the prompt capture is a genuine delta at
        # t=0. A centred gradient smears it across two samples and loses half
        # of it off the front of the array; the backward difference against an
        # implicit q(t<0)=0 puts the whole prompt term in sample 0, which is
        # the correct delta discretisation at this dt against a 180 ns shaper.
        cur = np.diff(np.concatenate([[0.0], q_u])) / DT_NS
        peaks[d] = float(shaper.apply(cur, dt_ns=DT_NS).max())
    return peaks


def profile_stats(peaks):
    ds = np.array(sorted(peaks))
    a = np.array([peaks[d] for d in ds])
    a = np.clip(a, 0.0, None)
    if a.max() <= 0:
        return dict(amp0=0.0, n2=0, rms=float("nan"))
    w = a / a.sum()
    mu = float(np.sum(ds * w))
    return dict(amp0=float(peaks[0]),
                n2=int(np.sum(a / a.max() > 0.02)),
                rms=float(np.sqrt(np.sum(w * (ds - mu) ** 2))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=1560)
    ap.add_argument("--ny-y", type=int, default=256)
    ap.add_argument("--n-xpos", type=int, default=8)
    ap.add_argument("--out", default="xy_anisotropy_ab.json")
    a = ap.parse_args()

    ny_x = K.x_ny(a.ny_y)
    times = K.log_times(n=48, t_min=1e-10, t_max=3e-6)
    t_grid = np.arange(0.0, T_MAX_NS + DT_NS, DT_NS)
    shaper = DreamShaper(dt_ns=DT_NS)

    # deposit x positions: 8 across two pad pitches, as rho_s_sensitivity does
    x0s = 5 * C.PAD_PITCH_M + np.linspace(
        0.0, 2 * C.PAD_PITCH_M, a.n_xpos, endpoint=False)

    results = {}
    for arm, uniform in (("strips", False), ("uniform", True)):
        t0 = time.time()
        sy = make_solver(K.Y_BOX_M, a.nx, a.ny_y, uniform)
        gy = {p: run_solve(sy, K.y_channel_pattern(sy, p), times, uniform)
              for p in (0, 1)}
        sx = make_solver(K.X_BOX_M, a.nx, ny_x, uniform)
        gx = {c: run_solve(sx, K.x_channel_pattern(sx, c), times, uniform)
              for c in range(C.N_PAD_PER_SUPER)}
        print(f"[{arm}] solved in {time.time()-t0:.1f} s", flush=True)

        for tag, use_ions in (("prompt", False), ("ion", True)):
            accx, accy = [], []
            for x0 in x0s:
                by = K.charge_budget_y(sy, gy, x0, dmax=DMAX, row0_parity=0)
                bx = K.charge_budget_x(sx, gx, x0, dmax=DMAX, y0_m=0.0)
                accy.append(profile_stats(
                    shaped_peaks(by, times, shaper, t_grid, use_ions)))
                accx.append(profile_stats(
                    shaped_peaks(bx, times, shaper, t_grid, use_ions)))

            def avg(rows):
                return {k: float(np.mean([r[k] for r in rows]))
                        for k in rows[0]}

            mx, my = avg(accx), avg(accy)
            results[f"{arm}_{tag}"] = dict(X=mx, Y=my,
                                           YX_ratio=my["amp0"] / mx["amp0"])
            print(f"  {tag:6s}  X amp0={mx['amp0']:.6g} rms={mx['rms']:.2f}  |  "
                  f"Y amp0={my['amp0']:.6g} rms={my['rms']:.2f} n>2%={my['n2']:.1f}"
                  f"  |  Y/X = {my['amp0']/mx['amp0']:.3f}", flush=True)
        del sy, gy, sx, gx

    print("\n=== VERDICT ===")
    for tag in ("prompt", "ion"):
        rs = results[f"strips_{tag}"]["YX_ratio"]
        ru = results[f"uniform_{tag}"]["YX_ratio"]
        print(f"{tag:6s}  Y/X strips = {rs:.3f}   uniform = {ru:.3f}   "
              f"asymmetry removed: {(1-abs(1-ru)/abs(1-rs))*100:5.1f} %"
              if abs(1 - rs) > 1e-9 else "")
        print(f"        Y sharing rms strips = "
              f"{results[f'strips_{tag}']['Y']['rms']:.2f} -> uniform "
              f"{results[f'uniform_{tag}']['Y']['rms']:.2f} strips")

    results["_meta"] = dict(nx=a.nx, ny_y=a.ny_y, ny_x=ny_x, rho_s=RHO_S,
                            d_k=D_K, dmax=DMAX, n_xpos=a.n_xpos,
                            f_electron=F_ELECTRON, transit_ns=TRANSIT_NS)
    with open(a.out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
