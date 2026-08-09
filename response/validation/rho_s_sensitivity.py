#!/usr/bin/env python3
"""
rho_s_sensitivity.py — could a wrong ESL surface resistivity explain the T14
discrepancies? (2026-08-09, answer: no — see the report for the full argument.)

THE QUESTION
============
ρ_s has never been measured on our modules (it is a scan axis frozen at
2 MΩ/sq, `design/report/T14_FREEZE_QUEUE_2026-08-09.md`). The natural worry
after T14: higher ρ_s would raise amplitudes, narrow sharing, and speed the
rise — the three directions the sim misses in. This module measures all three
derivatives on the REAL chain components rather than arguing them: unit point
charge on the ESL → `CombKernelLUT` per-channel induced currents →
`DreamShaper` (180 ns peaking, β = 0.75), at every ρ_s the S1 kernel archive
holds, prompt-only and with the analytic ion rectangle folded in.

WHAT IT FOUND (design/report/RHO_S_SENSITIVITY_2026-08-09.md)
=============================================================
* amplitude: ×10 in ρ_s → shaped peak ×1.08–1.15 (X), ×1.61 (Y). Cannot reach
  the ×1.7–1.9 T14 deficit, which is anyway voltage-dependent (HV-slope) —
  a static sheet property is excluded by signature, not just by size.
* rise: SLOWER with higher ρ_s (ion-folded Y 196→228 ns over the ×10) — the
  opposite of what data wants; the rise floor is the ion term (noions).
* sharing: the one observable with the right signature (higher ρ_s narrows
  it), but halving the spread needs ~8 MΩ/sq, outside the T2b band — and T2b
  pins exactly the sharing-relevant product ρ_s·c′.
* nugget: the sim-side X/Y amplitude asymmetry shrinks with ρ_s and vanishes
  at 5 MΩ/sq in this test, in direct tension with the T2b band.

CAVEATS
=======
Point deposit, not a track; the Y-channel parity assignment alternates with
row offset and is averaged over deposit x (second order for cross-ρ_s
ratios); the ρ_s ladder is the W1 dk50 family because W2 exists only at
rho2M — the W1-vs-W2 rho2M pair is reported so the boundary/glue offset
(~12 % on amplitude) can be seen to be ρ_s-independent.

    python3 -m response.validation.rho_s_sensitivity \
        --s1 ~/x17/response_sim/s1 \
        --w2 ~/x17/response_sim/s1_ny1024/greens_comb_rho2M_dk50um_g19um.npz \
        --out design/report/rho_s_sensitivity.json
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re

import numpy as np

from ..common import constants as C
from ..digitizer.kernel_lut import CombKernelLUT
from ..dream.shaper import DreamShaper

DT_NS = 2.0
T_MAX_NS = 2400.0            # covers the 1.92 µs DAQ window + margin
N_SIDE = 8
F_ELECTRON = 0.092           # 1 - f_ion, ions.py defended value
TRANSIT_NS = 340.0


def rise_10_90(t_ns, w):
    ip = int(np.argmax(w))
    p = w[ip]
    if p <= 0:
        return np.nan
    seg = w[: ip + 1]
    t10 = np.interp(0.1 * p, seg, t_ns[: ip + 1])
    t90 = np.interp(0.9 * p, seg, t_ns[: ip + 1])
    return t90 - t10


def metrics(shaper, chan_currents, t_ns):
    shaped = {n: shaper.apply(cur, dt_ns=DT_NS)
              for n, cur in chan_currents.items()}
    ns = np.array(sorted(shaped))
    amps = np.array([shaped[n].max() for n in ns])
    frac = amps / amps.max()
    w = amps / amps.sum()
    mu = np.sum(ns * w)
    return dict(amp0=float(shaped[0].max()),
                rise=float(rise_10_90(t_ns, shaped[0])),
                n5=int(np.sum(frac > 0.05)),
                n2=int(np.sum(frac > 0.02)),
                rms=float(np.sqrt(np.sum(w * (ns - mu) ** 2))))


def channel_currents(lut, ix):
    """X and Y per-channel currents for a deposit at LUT x sample ix, y=0."""
    jy0 = lut.iy_X(0.0)
    xcur = {int(d): lut.I_X[j, jy0, ix, :].astype(float)
            for j, d in enumerate(lut.ds)}
    ycur = {}
    for n in range(-N_SIDE, N_SIDE + 1):
        p = abs(n) % 2
        iy = lut.iy_Y(n * C.PAD_PITCH_M)
        ycur[n] = lut.I_Y[p, iy, ix, :].astype(float)
    return xcur, ycur


def avg_metrics(shaper, lut):
    """Metrics averaged over 8 deposit x positions across two pad pitches."""
    t_ns = lut.t * 1e9
    x0 = 5 * C.PAD_PITCH_M
    offs = np.linspace(0.0, 2 * C.PAD_PITCH_M, 8, endpoint=False)
    accx, accy = [], []
    for o in offs:
        ix = int(lut.ix(x0 + o))
        xcur, ycur = channel_currents(lut, ix)
        accx.append(metrics(shaper, xcur, t_ns))
        accy.append(metrics(shaper, ycur, t_ns))

    def avg(rows):
        return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}

    return avg(accx), avg(accy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", default=os.path.expanduser("~/x17/response_sim/s1"),
                    help="directory of W1 greens_comb_*.npz products")
    ap.add_argument("--dk", default="dk50um",
                    help="which d_k family to ladder over (default dk50um)")
    ap.add_argument("--w2", default=os.path.expanduser(
                        "~/x17/response_sim/s1_ny1024/"
                        "greens_comb_rho2M_dk50um_g19um.npz"),
                    help="production W2 kernel for the boundary/glue "
                         "reference row ('' to skip)")
    ap.add_argument("--out", default="rho_s_sensitivity.json")
    a = ap.parse_args()

    paths = {}
    for p in sorted(glob.glob(os.path.join(a.s1, f"greens_comb_*_{a.dk}.npz"))):
        m = re.search(r"greens_comb_(rho[\d.]+M)_", os.path.basename(p))
        if m:
            paths[f"W1 {m.group(1)} {a.dk}"] = p
    if a.w2 and os.path.exists(a.w2):
        paths["W2 " + os.path.basename(a.w2)] = a.w2
    if not paths:
        raise SystemExit(f"no kernels found under {a.s1}")

    shaper = DreamShaper(dt_ns=DT_NS)
    results = {}
    for name, path in paths.items():
        print(f"=== {name} ===", flush=True)
        lut = CombKernelLUT(path, dt_ns=DT_NS, t_max_ns=T_MAX_NS,
                            n_side=N_SIDE)
        out = {"path": path}
        out["prompt_X"], out["prompt_Y"] = avg_metrics(shaper, lut)
        lut.apply_ion_transit(F_ELECTRON, TRANSIT_NS)
        out["ion_X"], out["ion_Y"] = avg_metrics(shaper, lut)
        del lut
        results[name] = out
        for tag in ("prompt", "ion"):
            for v in ("X", "Y"):
                m = out[f"{tag}_{v}"]
                print(f"  {tag:6s} {v}: amp0={m['amp0']:.4f}  "
                      f"rise={m['rise']:6.1f} ns  n>5%={m['n5']:.1f}  "
                      f"n>2%={m['n2']:.1f}  rms={m['rms']:.2f} strips",
                      flush=True)

    with open(a.out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
