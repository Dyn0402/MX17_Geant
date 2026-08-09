#!/usr/bin/env python3
"""
psi_readout.py — the READOUT weighting potential through the real woven mesh,
and the ion/electron charge split it implies.

Work-list item 2 of design/report/HANDOFF_S3_ION_2026-08-09.md:

    "Re-derive f_ion on the readout electrode with the T6 production field map
     + the local Garfield ComponentGrid patch (mesh screening: in-gap split
     != through-mesh split)."

WHY THIS IS THE QUESTION. Everything downstream of S3 uses the parallel-plate
ansatz psi(z) = 1 - z/gap. `mx17_aval_calib.py` states it plainly and, crucially,
KEEPS it even in the meshfield branch: the DRIFT field is the realistic woven
mesh loaded through ComponentGrid, but the WEIGHTING field is still a
`ComponentConstant` with `SetWeightingField(0, 0, 1/gap)`. So the measured
S3 v2 template's f_ion = 0.9006 is the parallel-plate in-gap split evaluated on
a meshfield-driven avalanche depth profile -- it is NOT a through-mesh readout
split, and it has never been checked against one. That is the gap this closes.

WHAT IS ACTUALLY SOLVED. `solve_fieldmap.py` already computes the answer and
nobody noticed: its unit problem `u_A` has anode = +V_MESH, wires = 0,
top = 0, which IS the readout weighting potential up to the scale V_MESH.
This script reuses that solve verbatim -- same geometry, same mesher, same
element -- and does not re-derive any geometry. psi = u_A / V_MESH.

The physical solution exported as `meshfield_production.txt` is u_A + c*u_B,
so the shipped map cannot be used for this: u_B is exactly the drift-side
admixture that a weighting potential must NOT have.

THE OBSERVABLE. For a pair born at height z above the ESL, the electron falls
to the anode and the ion climbs to the mesh, so with the true psi

    f_electron(z) = psi(anode) - psi(z) = 1 - psi(z)
    f_ion(z)      = psi(z) - psi(end)

and the split is averaged over the MEASURED alpha_z profile from the S3 calib,
not over its mean -- f_ion is linear in psi, so for the parallel-plate case
mean-vs-distribution makes no difference, but for a curved psi it does, and
assuming it away is how the 5 um / 13.84 um error happened the first time.

psi(end) is where the ion stops CONTRIBUTING inside a waveform, not where it
stops existing. `funnel_ion_endpoints.json` (production map, 3000 ions)
measured 94.8 % absorbed on the mesh wires (psi = 0 exactly) and 5.2 % escaping
into the drift bulk. An escaped ion sits in E_drift = 333 V/cm where
v = mu*E ~ 5 nm/ns, so it moves ~5 um in a microsecond and is frozen on the
timescale of any DREAM waveform: its residual psi is charge the readout never
sees. That is a REAL loss channel and it is included below, weighted 5.2 %.

Run (~1 min at smoke resolution, on the laptop):

    ~/PycharmProjects/nTof_x17/.venv/bin/python psi_readout.py --smoke

Resolution is checked, not assumed: away from the weave the transverse average
of psi is EXACTLY linear in z (Laplace + periodic cell kills every harmonic but
the constant, and they decay as exp(-2*pi*d/pitch), i.e. e^-12.8 at the ion
birth height). The script fits that line over the amp bulk and reports the
residual; if the mesh is too coarse to reproduce an exactly-linear function,
the number is not to be used. That gate is what makes a smoke-resolution solve
admissible here, where it is not admissible for S3.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from response.meshcell import solve_fieldmap as S  # noqa: E402


def plane_average_psi(ev, z_um, n_side, pitch):
    """Transverse average of psi over one unit cell at height z.

    Averaged over the GAS only: a point inside a wire is not a place an ion can
    be, and including the wires' psi = 0 would bias the average toward zero in
    exactly the z-range where the answer matters. Returns (mean, gas_fraction);
    a gas fraction < 1 flags that the plane cuts the weave.
    """
    q = (np.arange(n_side) + 0.5) / n_side - 0.5
    X, Y = np.meshgrid(q * pitch, q * pitch, indexing="ij")
    pts = np.vstack([X.ravel(), Y.ravel(), np.full(X.size, float(z_um))])
    gas = ~S.in_wire(pts)
    if not gas.any():
        return np.nan, 0.0
    val, _, ok = ev(pts[:, gas])
    ok = np.asarray(ok, bool)
    if not ok.any():
        return np.nan, 0.0
    return float(np.mean(val[ok])), float(gas.sum()) / gas.size


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--smoke", action="store_true",
                   help="lc_wire=2.0, lc_far=14.0 — ~1 min. Admissible here "
                        "because of the linearity gate; NOT for S3.")
    g.add_argument("--production", action="store_true",
                   help="lc_wire=0.8, lc_far=8.0 — minutes, more memory.")
    g.add_argument("--lc", type=float, nargs=2, metavar=("WIRE", "FAR"),
                   help="explicit (lc_wire, lc_far), for the resolution "
                        "convergence check the quoted number rests on")
    ap.add_argument("--n-side", type=int, default=48,
                    help="transverse samples per side per z-plane")
    ap.add_argument("--nz", type=int, default=241,
                    help="z-planes from the ESL to the drift cut")
    ap.add_argument("--calib", default=os.path.join(
        HERE, "..", "avalanche", "aval_calib_meshfield_pooled.json"),
        help="S3 calib supplying the MEASURED alpha_z profile to average over")
    ap.add_argument("--funnel", default=os.path.join(
        HERE, "funnel_ion_endpoints.json"),
        help="T6 ion-endpoint fractions (absorbed vs escaped)")
    ap.add_argument("--out", default=os.path.join(HERE, "psi_readout.json"))
    a = ap.parse_args()

    if a.lc:
        lc_wire, lc_far = a.lc
        tag = f"lc{lc_wire:g}_{lc_far:g}"
    else:
        lc_wire, lc_far = (2.0, 14.0) if a.smoke else (0.8, 8.0)
        tag = "smoke" if a.smoke else "production"
    t0 = time.time()

    msh = os.path.join(tempfile.gettempdir(), f"psi_readout_{tag}.msh")
    print(f"[psi] meshing (lc_wire={lc_wire}, lc_far={lc_far}) ...")
    S.build_mesh(lc_wire, lc_far, msh)
    m = S.load_mesh(msh)
    print(f"[psi] solving ... ({time.time()-t0:.0f}s)")
    basis, u_A, _u_B = S.solve(m)

    # psi = u_A / V_MESH: u_A is anode=+V_MESH, wires=0, top=0, which is the
    # readout weighting problem exactly. Normalising by V_MESH (not V_ANODE,
    # though they are equal) keeps it correct under --v-mesh.
    ev = S.P2Evaluator(m, basis, u_A / S.V_MESH)

    zs = np.linspace(S.Z_ANODE, S.Z_CATH, a.nz)
    psi, gasfrac = np.zeros(a.nz), np.zeros(a.nz)
    for i, z in enumerate(zs):
        psi[i], gasfrac[i] = plane_average_psi(ev, z, a.n_side, S.PITCH)
    print(f"[psi] sampled {a.nz} planes ({time.time()-t0:.0f}s)")

    # ── Gate: <psi> must be EXACTLY linear in the amp bulk ───────────────────
    # Harmonics decay as exp(-2*pi*d/pitch), so 30 um below the weave underside
    # they are e^-2.8 = 6 %, and 60 um below they are 0.4 %. Fit over the lower
    # two thirds of the gap, where the exact answer is a straight line.
    z_rel = zs - S.Z_ANODE                       # height above the ESL [um]
    bulk = (z_rel >= 2.0) & (z_rel <= 0.66 * S.AMP_GAP)
    cf = np.polyfit(z_rel[bulk], psi[bulk], 1)
    resid = psi[bulk] - np.polyval(cf, z_rel[bulk])
    lin_rms = float(np.sqrt(np.mean(resid ** 2)))
    print(f"[psi] amp-bulk linearity: slope {cf[0]:+.6e} /um, "
          f"intercept {cf[1]:.6f}, residual RMS {lin_rms:.2e}")
    gate_lin = lin_rms < 1e-3
    print(f"[psi] GATE linearity (<1e-3): {'PASS' if gate_lin else 'FAIL'}")

    psi_anode = float(np.interp(0.0, z_rel, psi))
    print(f"[psi] GATE psi(anode) = {psi_anode:.6f} (exact 1) : "
          f"{'PASS' if abs(psi_anode - 1) < 2e-3 else 'FAIL'}")

    # ── psi where the ion stops ──────────────────────────────────────────────
    # Absorbed on a wire -> psi = 0 exactly (Dirichlet). Escaped -> it freezes
    # just above the weave; take psi at the mesh topside as its residual.
    z_top_rel = S.Z_TOP - S.Z_ANODE
    psi_escape = float(np.interp(z_top_rel, z_rel, psi))
    fr = json.load(open(a.funnel))["ion_endpoints"]
    f_abs, f_esc = fr["frac_absorbed_on_mesh"], fr["frac_escaped_to_drift"]
    psi_end = f_abs * 0.0 + f_esc * psi_escape
    print(f"[psi] psi at mesh topside {psi_escape:.5f}; ion fates "
          f"{f_abs:.3f} absorbed / {f_esc:.3f} escaped -> <psi_end> "
          f"{psi_end:.5f}")

    # ── Average over the MEASURED alpha_z profile ────────────────────────────
    cal = json.load(open(os.path.expanduser(a.calib)))["point"]
    zh = cal["alpha_z_hist"]
    cnt = np.asarray(zh["counts"], float)
    ed = np.asarray(zh["edges"], float)
    zc = 0.5 * (ed[:-1] + ed[1:])                # height above the ESL [um]
    w = cnt / cnt.sum()
    z_mean = float((w * zc).sum())

    psi_at_birth = np.interp(zc, z_rel, psi)
    f_e_true = float((w * (1.0 - psi_at_birth)).sum())
    f_i_true = float((w * (psi_at_birth - psi_end)).sum())
    # Parallel-plate reference, i.e. exactly what S3 and ions.py assume.
    f_e_pp = float((w * (zc / S.AMP_GAP)).sum())
    f_i_pp = 1.0 - f_e_pp

    # Renormalised split: the fraction of the charge the readout ACTUALLY sees
    # that is slow. This is the number the digitizer's f_ion dial corresponds
    # to, because apply_ion_transit mixes prompt and rectangle to unit sum.
    seen = f_e_true + f_i_true
    f_ion_eff = f_i_true / seen if seen else float("nan")

    print("\n=== charge split on the READOUT electrode ===")
    print(f"  mean ion birth height        {z_mean:.2f} um "
          f"(S3 v2 quoted 13.84)")
    print(f"  parallel plate  f_e {f_e_pp:.4f}   f_ion {f_i_pp:.4f}")
    print(f"  true psi        f_e {f_e_true:.4f}   f_ion {f_i_true:.4f}"
          f"   (sum {seen:.4f}; {1-seen:.4f} lost to escaped ions)")
    print(f"  f_ion as the digitizer dial sees it: {f_ion_eff:.4f}")
    print(f"  SHIFT vs parallel plate: {f_ion_eff - f_i_pp:+.4f}")

    out = {
        "schema": "psi_readout/1",
        "resolution": tag, "lc_wire": lc_wire, "lc_far": lc_far,
        "n_side": a.n_side, "nz": a.nz,
        "geometry": {"pitch_um": S.PITCH, "wire_r_um": S.WIRE_R,
                     "amp_gap_um": S.AMP_GAP, "z_anode_um": S.Z_ANODE,
                     "z_under_um": S.Z_UNDER, "z_top_um": S.Z_TOP,
                     "z_cath_um": S.Z_CATH, "v_mesh": S.V_MESH},
        "gates": {"linearity_rms": lin_rms, "linearity_pass": bool(gate_lin),
                  "psi_anode": psi_anode},
        "bulk_fit": {"slope_per_um": float(cf[0]),
                     "intercept": float(cf[1]),
                     "parallel_plate_slope_per_um": -1.0 / S.AMP_GAP},
        "psi_profile": {"z_above_anode_um": z_rel.tolist(),
                        "psi": psi.tolist(),
                        "gas_fraction": gasfrac.tolist()},
        "ion_fates": {"frac_absorbed": f_abs, "frac_escaped": f_esc,
                      "psi_at_mesh_topside": psi_escape, "psi_end": psi_end},
        "split": {"z_birth_mean_um": z_mean,
                  "f_electron_parallel_plate": f_e_pp,
                  "f_ion_parallel_plate": f_i_pp,
                  "f_electron_true": f_e_true, "f_ion_true": f_i_true,
                  "fraction_seen": seen, "f_ion_effective": f_ion_eff,
                  "shift_vs_parallel_plate": f_ion_eff - f_i_pp},
        "inputs": {"calib": os.path.abspath(os.path.expanduser(a.calib)),
                   "funnel": os.path.abspath(a.funnel)},
        "seconds": time.time() - t0,
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n[psi] wrote {a.out}  ({time.time()-t0:.0f}s)")
    return 0 if (gate_lin and abs(psi_anode - 1) < 2e-3) else 1


if __name__ == "__main__":
    sys.exit(main())
