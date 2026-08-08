#!/usr/bin/env python3
"""
w2_production.py — the W2 grid re-solve, structured for lxplus condor.

Plan next-in-order item 3 (design/RESPONSE_SIM_PLAN.md): brute-force W2
re-solve of the S1 kernel grid, chosen 2026-08-08 over the (unphysical)
effective-diagonal shortcut. The solver and its acceptance battery are
`wpot_w2.py` / `w2_validate.py`; THIS file is the memory- and flop-ordered
production path on top of them, validated against the prototype assembly by
battery test 8.

JOB STRUCTURE — one condor job per Bloch family, because the expensive object
(the family eigendecomposition) is drive- and rho_s-independent:

    stage `family`   66 jobs: 64 Y-box families (ny=512) + 2 X-box families
                     (ny=16, the 1.56 mm comb box). Each job assembles its
                     family once, then propagates ALL of its box's drives
                     (2 Y or 40 X channels) x all 4 rho_s x the 61 log times,
                     and writes one mode-coefficient slab (~0.1-2 GB) to EOS.
                     rho_s costs nothing extra: M_sheet ∝ 1/rho_s, so the
                     eigenvectors are shared and only exp(-w t) rescales
                     (w2_validate test 3 proves this exact).
    stage `combine`  cheap, anywhere: gather the slabs, ifft2, write the
                     4 production npz products (kernels.run_point schema, so
                     every downstream reader works unchanged).

REAL ARITHMETIC. At the production registration (ESL phase -width/2, pad
lattice centred) sigma_s(x) and the gap mask are even, so sigma_hat and ghat
are real and every family operator is a REAL symmetric matrix in the plain
Fourier basis — no basis change needed, just a guarded cast. dsyevd beats
zheevd ~4x in flops and 2x in memory, which is what makes a 24,960-mode
family fit a 4-core / 32 GB condor slot (~5 GB per live matrix, peak ~25 GB).
The cast is guarded: asymmetric registrations raise rather than silently
truncate imaginary parts.

WHAT THIS RUN IS AND IS NOT. ny=512 (97.5 µm cells): the 100 µm gap edge is
resolved by fractional-coverage masks, not the grid — the w2_validate CG
ny-ladder measures what that costs at the prompt level and the plan records
the Richardson mitigation. The W1 grid's ny=1024 upgrade needed no extra
machinery; for W2, ny=1024 doubles the family to 49,920 modes (20 GB
matrices) and needs the y-parity/mirror-family reduction first — a recorded
follow-up, not part of this submission.

    # one family, locally or in a condor job
    python3 -m response.solver.w2_production family --box y --fam 17 --outdir X
    # after all 66 slabs exist
    python3 -m response.solver.w2_production combine --slabdir X --outdir Y
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

import numpy as np
import scipy.linalg as sla

from ..common import constants as C
from . import kernels as K
from .wpot_w2 import W2Solver

RHOS_OHM_SQ = C.RHO_S_SCAN_OHM_SQ          # (0.5, 1, 2, 5) MOhm/sq
RHO_REF = 1.0e6                            # families are solved at this rho_s
NY_Y, NY_X = 512, 16                       # NY_X = kernels.x_ny(NY_Y)
NX = 3120
D_K_M = C.KAPTON_THICK_UM * 1e-6


def _git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def make_solver(box):
    if box == "y":
        return W2Solver(RHO_REF, nx=NX, ny=NY_Y, ly_m=K.Y_BOX_M,
                        esl_phase_m=K.ESL_PHASE_CENTERED_M)
    if box == "x":
        return W2Solver(RHO_REF, nx=NX, ny=NY_X, ly_m=K.X_BOX_M,
                        esl_phase_m=K.ESL_PHASE_CENTERED_M)
    raise ValueError(box)


def drives_for(s, box):
    """
    [(name, vpad)] for every channel kernel this box owns.

    Drives are SNAPPED to the same hard cell assignment as the solver's metal
    mask: a cell is either on the pad (driven at the full potential) or it is
    not. kernels.py's fractional-coverage patterns are the right tool when the
    pad edge only sets a charge total; here the pattern is a Dirichlet clamp
    and must live exactly on the clamped set, or the drive would fight the
    floating-gap condition on the shared cells (see W2Solver's mask_snap note).
    """
    if box == "y":
        pats = [(f"Y{p}", K.y_channel_pattern(s, p)) for p in (0, 1)]
    else:
        pats = [(f"X{c}", K.x_channel_pattern(s, c))
                for c in range(C.N_PAD_PER_SUPER)]
    # STRICTLY greater: the pad x-edges fall exactly on cell centres (the
    # 680/780 µm lattice against wpot's 10 µm grid), so edge cells are covered
    # at exactly 0.5. The gap mask takes ties (frac >= 0.5 -> gap, W2Solver);
    # a >= here handed the same cells to BOTH sides and killed the first
    # canaries on the overlap guard (cluster 13353068). With > the two snaps
    # are exactly complementary: a channel's own coverage never exceeds the
    # total metal coverage, so drive > 0.5 implies metal > 0.5 implies not gap.
    return [(n, (np.round(p, 9) > 0.5).astype(float)) for n, p in pats]


def _real_guard(s):
    """The real-symmetric cast is only legal for even registrations."""
    for name, arr in (("sigma_hat", s.sigma_hat), ("ghat", s.ghat)):
        rel = float(np.abs(arr.imag).max() / max(np.abs(arr).max(), 1e-300))
        if rel > 1e-10:
            raise RuntimeError(
                f"{name} has imaginary parts ({rel:.1e}): registration is not "
                "even-symmetric, so the real-arithmetic production path is "
                "invalid here. Use the prototype (complex) solver instead.")


def _chunked_lookup(table, IA, IB, mod, out, chunk=2048):
    """out[i, j] = table[(IA[i] - IA[j]) % mod] without a full int index array."""
    for i0 in range(0, len(IA), chunk):
        d = (IA[i0:i0 + chunk, None] - IA[None, :]) % mod
        out[i0:i0 + chunk] = table[d]
    return out


def solve_family(s, drives, ifam, times, rhos=RHOS_OHM_SQ, verbose=True,
                 label=""):
    """
    One family job: assemble, eigendecompose, propagate every drive x rho_s.

    `s` is any W2Solver with an even (real-castable) registration; `drives`
    is [(name, vpad)] with every vpad vanishing on the gap support. Returns
    (slab, meta): slab complex64 (n_drive, n_rho, nt, N) of V0-hat family
    coefficients, plus everything combine needs to scatter them.
    """
    t00 = time.time()
    _real_guard(s)
    for name, vpad in drives:
        leak = float(np.abs(vpad * s.rgap).max())
        if leak > 1e-12:
            raise RuntimeError(f"drive {name} overlaps the gap support "
                               f"({leak:.1e}) — snap it to the metal mask")
    fams = s.families()
    IX, IY = fams[ifam]
    N = len(IX)
    log = (lambda m: print(f"  [{time.time()-t00:7.1f}s] {m}", flush=True)) \
        if verbose else (lambda m: None)
    log(f"{label} family {ifam}/{len(fams)}  N={N}  "
        f"(gx={s.gx}, gy={s.gy}, ny={s.ny})")

    kx, ky = s.kx[IX], s.ky[IY]
    A12, Cg, D = s._diags(np.hypot(kx, ky))
    if D is None:
        raise RuntimeError("W1-grounded substrate reached production path")

    # R: the gap-mask convolution matrix, real, chunked assembly
    R = np.empty((N, N))
    gre = np.ascontiguousarray(s.ghat.real)
    for i0 in range(0, N, 2048):
        dx = (IX[i0:i0 + 2048, None] - IX[None, :]) % s.nx
        dy = (IY[i0:i0 + 2048, None] - IY[None, :]) % s.ny
        R[i0:i0 + 2048] = gre[dy, dx]
    log(f"R assembled ({R.nbytes/1e9:.1f} GB)")

    # Gd = R diag(D) R, its spectral pseudo-inverse folded straight into C_eff
    G = R * D[None, :] @ R                    # (N,N) dgemm
    G = 0.5 * (G + G.T)
    log("Gd built")
    wg, Vg = sla.eigh(G, overwrite_a=True, check_finite=False, driver="evd")
    del G
    keep = wg > s.pinv_rtol * max(float(wg.max()), 1e-300)
    Vk = np.ascontiguousarray(Vg[:, keep])
    wk = wg[keep]
    del Vg
    log(f"Gd eigh done: rank {keep.sum()}/{N}")
    Wm = (A12[:, None] * R) @ Vk              # (N, Nk)
    del R, Vk
    Ceff = (Wm / wk[None, :]) @ Wm.T
    del Wm
    Ceff *= -1.0
    Ceff[np.diag_indices(N)] += Cg
    Ceff = 0.5 * (Ceff + Ceff.T)
    log("C_eff built")
    L = sla.cholesky(Ceff, lower=True, overwrite_a=True, check_finite=False)
    del Ceff
    log("Cholesky done")

    # M_sheet at RHO_REF, chunked; then B = L^-1 M L^-T and its eigh
    M = np.empty((N, N))
    sre = np.ascontiguousarray(s.sigma_hat.real)
    for i0 in range(0, N, 2048):
        dx = (IX[i0:i0 + 2048, None] - IX[None, :]) % s.nx
        Mc = sre[dx] * (kx[i0:i0 + 2048, None] * kx[None, :]
                        + ky[i0:i0 + 2048, None] * ky[None, :])
        Mc[IY[i0:i0 + 2048, None] != IY[None, :]] = 0.0
        M[i0:i0 + 2048] = Mc
    M = 0.5 * (M + M.T)
    log("M_sheet assembled")
    B = sla.solve_triangular(L, M, lower=True, overwrite_b=True,
                             check_finite=False)
    B = sla.solve_triangular(L, B.T, lower=True, overwrite_b=True,
                             check_finite=False)
    del M
    B = 0.5 * (B + B.T)
    w, Q = sla.eigh(B, overwrite_a=True, check_finite=False, driver="evd")
    del B
    w = np.maximum(w, 0.0)
    log(f"B eigh done (w max {w.max():.3e} 1/s at rho_ref)")

    # Propagate every drive x rho. Prompt comes from the full-grid CG (cheap,
    # matrix-free, and its identity with the family assembly is battery test 2)
    drv = drives
    times = np.asarray(times, float)
    slab = np.empty((len(drv), len(rhos), len(times), N), np.complex64)
    prompt_caps = {}
    for di, (name, vpad) in enumerate(drv):
        v0, it = s.solve_prompt_cg(vpad)
        prompt_caps[name] = float(v0.mean())
        v0f = np.fft.fft2(v0)[IY, IX]
        coef = Q.T @ (L.T @ v0f)
        for ri, rho in enumerate(rhos):
            wt = np.exp(-np.outer(w * (RHO_REF / rho), times))
            U = Q @ (wt * coef[:, None])
            slab[di, ri] = sla.solve_triangular(
                L.T, U, lower=False, check_finite=False).T.astype(np.complex64)
        log(f"drive {name} propagated ({it} CG it)")

    meta = {"box": label, "ifam": ifam, "n_fam": len(fams), "N": N,
            "nx": s.nx, "ny": s.ny, "ly_m": s.ly,
            "rho_ref": RHO_REF, "rhos": list(rhos),
            "drives": [n for n, _ in drv], "prompt_capture": prompt_caps,
            "gd_rank": int(keep.sum()), "pinv_rtol": s.pinv_rtol,
            "d_k_m": D_K_M, "d_glue_m": s.d_g,
            "boundary": "W2 nominal substrate (glue trench + FR4 to L5)",
            "git": _git_hash(), "wall_s": time.time() - t00}
    return slab, IX, IY, meta


def cmd_family(a):
    times = K.log_times(60)
    s = make_solver(a.box)
    slab, IX, IY, meta = solve_family(s, drives_for(s, a.box), a.fam, times,
                                      label=f"box={a.box}")
    os.makedirs(a.outdir, exist_ok=True)
    path = os.path.join(a.outdir, f"w2slab_{a.box}{a.fam:03d}.npz")
    np.savez_compressed(path, slab=slab, IX=IX, IY=IY, t=times,
                        meta=json.dumps(meta))
    print(f"-> {path}  ({os.path.getsize(path)/1e6:.0f} MB, "
          f"{meta['wall_s']/60:.0f} min)")
    return 0


def _gather(slabdir, box, n_fam, ri, di, nt, ny, nx):
    """Assemble one (drive, rho) kernel from every family slab of a box."""
    out_hat = np.zeros((nt, ny, nx), complex)
    for f in range(n_fam):
        with np.load(os.path.join(slabdir, f"w2slab_{box}{f:03d}.npz")) as z:
            IX, IY = z["IX"], z["IY"]
            out_hat[:, IY, IX] = z["slab"][di, ri].astype(complex)
    return np.real(np.fft.ifft2(out_hat, axes=(1, 2))).astype(np.float32)


def cmd_combine(a):
    times = K.log_times(60)
    nt = len(times)
    sy, sx = make_solver("y"), make_solver("x")
    # family counts from the slabs actually present
    n_y = len(sy.families())
    n_x = len(sx.families())
    for f in range(n_y):
        p = os.path.join(a.slabdir, f"w2slab_y{f:03d}.npz")
        if not os.path.exists(p):
            raise FileNotFoundError(f"missing slab {p} — combine needs all")
    for f in range(n_x):
        p = os.path.join(a.slabdir, f"w2slab_x{f:03d}.npz")
        if not os.path.exists(p):
            raise FileNotFoundError(f"missing slab {p}")
    os.makedirs(a.outdir, exist_ok=True)
    d_g = C.GLUE_THICK_UM * 1e-6
    for ri, rho in enumerate(RHOS_OHM_SQ):
        gy = {p: _gather(a.slabdir, "y", n_y, ri, p, nt, sy.ny, sy.nx)
              for p in (0, 1)}
        gx = {c: _gather(a.slabdir, "x", n_x, ri, c, nt, sx.ny, sx.nx)
              for c in range(C.N_PAD_PER_SUPER)}
        ytot = K.sum_over_rows(sy, gy)
        xtot = K.sum_over_columns(sx, gx)
        tbar = (ytot + xtot).mean(axis=(1, 2))
        xfrac = xtot.mean(axis=(1, 2)) / tbar
        res = {"rho_s_ohm_sq": rho, "d_kapton_m": D_K_M,
               "boundary_model": "W2",
               "sub_layers_um_eps": [[t * 1e6, e]
                                     for t, e in sy.sub_layers],
               "d_glue_m": d_g, "glue_eps_r": C.GLUE_EPS_R,
               "sum_rule_expect": K.prompt_sum_rule(D_K_M),
               "channel_capture_prompt": float(tbar[0]),
               "channel_capture_late": float(tbar[-1]),
               "x_fraction_prompt": float(xfrac[0]),
               "x_fraction_late": float(xfrac[-1]),
               "git": _git_hash(), "nx": sy.nx, "ny": sy.ny}
        res["sharing"] = K.sharing_report(sy, gy, sx, gx, times, verbose=False)
        tag = f"rho{rho/1e6:g}M_dk{D_K_M*1e6:g}um_g{d_g*1e6:.0f}um"
        path = os.path.join(a.outdir, f"greens_comb_w2_{tag}.npz")
        np.savez_compressed(
            path, t=times, x=sy.x, y_Y=sy.y, y_X=sx.y,
            G_Y_even=gy[0], G_Y_odd=gy[1],
            G_X=np.stack([gx[c] for c in range(C.N_PAD_PER_SUPER)]),
            x_cols=K.pad_x(np.arange(C.N_PAD_PER_SUPER)),
            meta=json.dumps(res))
        print(f"rho {rho/1e6:g}M: capture prompt {tbar[0]:.4f} late "
              f"{tbar[-1]:.4f}  X/Y prompt {xfrac[0]:.3f}"
              f"  -> {path} ({os.path.getsize(path)/1e6:.0f} MB)", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("family", help="solve one Bloch-family job")
    f.add_argument("--box", choices=("y", "x"), required=True)
    f.add_argument("--fam", type=int, required=True)
    f.add_argument("--outdir", required=True)
    g = sub.add_parser("combine", help="assemble products from all slabs")
    g.add_argument("--slabdir", required=True)
    g.add_argument("--outdir", required=True)
    a = ap.parse_args()
    return cmd_family(a) if a.cmd == "family" else cmd_combine(a)


if __name__ == "__main__":
    raise SystemExit(main())
