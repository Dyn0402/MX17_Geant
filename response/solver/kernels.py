#!/usr/bin/env python3
"""
kernels.py — T2b: the response kernels of the REAL readout channels.

design/RESPONSE_SIM_PLAN.md task T2b, "first priority": replace the solid-line
drive patterns in wpot.py with the comb patterns that the board actually has.
T2 (response/common/channel_map.py) read the checkerboard straight out of the
gerber connector stubs:

    a Y channel = one pad ROW,    pads at columns i with (i + j) EVEN
    an X channel = one pad COLUMN, pads at rows    j with (i + j) ODD

so a channel is 256 pads on a 1.56 mm pitch — a comb — not 512 pads on a
0.78 mm pitch. Every sharing number in plan §9 has to come from these.


HOW MANY DISTINCT KERNELS THERE ACTUALLY ARE
============================================
The stack is periodic in x with the 31.2 mm superperiod (39 ESL pitches =
40 pad pitches) and uniform in y. Combined with the checkerboard that collapses
the 1024 channels to a small, exact set:

  Y channels (rows).  A row's comb sits at every other COLUMN, so translating
  a Y channel by one row shifts its comb by one pad pitch in x and flips which
  columns it owns. The row's own y position is irrelevant (the sheet is uniform
  in y — pure translation). Only the column PARITY survives, so there are
  exactly **2 distinct Y kernels**, even-row and odd-row, each a function of
  the deposit position (x0, Δy) over the full superperiod in x.

  X channels (columns). A column's comb runs along y, where everything is
  translation invariant, so the row parity is again just a translation. What
  does NOT drop out is the column's own x position relative to the ESL strips:
  that is the 31.2 mm beat of plan §1. There are **40 distinct X kernels**, one
  per column phase.

So the beat lives in the X view. The Y view integrates over all x phases in its
own comb and only ever sees the beat through where the DEPOSIT lands.

Practical consequence for the solve: an X channel's drive is exactly periodic
in y with period 1.56 mm, so its box can BE 1.56 mm in y — the periodic images
are the rest of the comb, not an artifact. Those 40 solves are cheap. The Y
kernel needs the full ±25 mm of y because it is localised there, and it is the
expensive one.


THE SUM RULE (this file's main correctness test)
================================================
Drive EVERY pad at V_w. That is a uniform Dirichlet plane, which excites only
k = 0, where the sheet has no in-plane gradient and therefore never relaxes:

    Psi(t) = S(0)/C(0) = (eps_r/d) / (eps_r/d + 1/g)     for ALL t

= 0.875 at d = 75 µm, g = 150 µm. So for a unit charge anywhere on the ESL,

    sum over all X channels + sum over all Y channels = 0.875, at every time,

with the remaining 0.125 induced on the mesh. It exercises the comb assembly,
the Bloch blocks, the 40 X phases and the row sum all at once against a number
known in closed form, and it is where the plan §9 "X/Y charge balance
0.49/0.51" prediction comes from.

    python3 -m response.solver.kernels --quick
    python3 -m response.solver.kernels --outdir ~/x17/response_sim/s1
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

import numpy as np

from ..common import constants as C
from . import wpot as W
from .wpot import WeightingSolver

# ── Absolute geometry, in the physical active-area frame ─────────────────────
# Pad centres: 512 on a 0.78 mm pitch, centred on the active area, so the
# centre of the board falls between two pads (512 is even).
PAD_ORIGIN_M = -(C.PAD_N - 1) / 2 * C.PAD_PITCH_M          # -199.29 mm

# ESL registration (plan §12 item 3, unknown -> "nominal aligned at centre"):
# one 550 µm strip centred on x = 0. esl_sigma_profile puts the strip at
# [phase, phase + width), so the centred phase is -width/2.
ESL_PHASE_CENTERED_M = -C.ESL_WIDTH_M / 2

# The Y box: 64 pad rows = 49.92 mm, which is BOTH commensurate with the
# 0.78 mm row pitch (so the row sum is an exact stride) and ~ the ±25 mm the
# plan's export grid asks for.
Y_BOX_M = 64 * C.PAD_PITCH_M                               # 49.92 mm
# The X box: exactly one comb period. The periodic images ARE the comb.
X_BOX_M = 2 * C.PAD_PITCH_M                                # 1.56 mm


def pad_x(i):
    """x centre [m] of pad column i, folded into the solver's 31.2 mm box."""
    return np.mod(PAD_ORIGIN_M + np.asarray(i) * C.PAD_PITCH_M, C.SUPERPERIOD_M)


def log_times(n=60, t_min=1e-10, t_max=1e-5):
    """The plan §2 time axis: t = 0 (prompt) plus n log-spaced points."""
    return np.concatenate([[0.0], np.geomspace(t_min, t_max, n)])


# ── Drive patterns ───────────────────────────────────────────────────────────

def _rect_cov(coords, centre, width, box):
    """
    Fraction of each grid CELL covered by a rectangle of the given width —
    i.e. an exactly area-preserving (anti-aliased) indicator, not a hard mask.

    This is not cosmetic. A 0.68 mm pad on a 0.78 mm lattice cannot be
    represented by a hard mask at any affordable resolution: the mask's area
    jumps as samples cross the pad edge, so the total induced charge came out
    0.76 at ny=64 and 0.57 at ny=256 — a 33 % swing with nothing physical
    changing. With fractional coverage the pad area is exact at every ny and
    only the edge *shape* converges, which it does fast.

    A happy side effect: at width = the pad PITCH, adjacent cells' fractions
    sum to 1 identically, so the checkerboard still partitions the plane
    exactly and the closed-form sum rule stays exact.
    """
    h = box / len(coords)
    d = ((np.asarray(coords) - centre + box / 2) % box) - box / 2
    lo = np.clip(d - h / 2, -width / 2, width / 2)
    hi = np.clip(d + h / 2, -width / 2, width / 2)
    return (hi - lo) / h


def _pads_at(solver, xs_m, ys_m, size_m=None, v=1.0):
    """
    Union of square pads at the given centres, on the solver's (ny, nx) grid.

    size_m defaults to the true 0.68 mm pad, leaving the 100 µm inter-pad gap
    at ground. Passing the 0.78 mm PITCH instead makes the pads tile the plane,
    which is the configuration the closed-form sum rule applies to.
    """
    size_m = C.PAD_SIZE_UM * 1e-6 if size_m is None else size_m
    cx = np.zeros(solver.nx)
    cy = np.zeros(solver.ny)
    for x0 in np.atleast_1d(xs_m):
        cx += _rect_cov(solver.x, x0, size_m, solver.lx)
    for y0 in np.atleast_1d(ys_m):
        cy += _rect_cov(solver.y, y0, size_m, solver.ly)
    return v * np.outer(cy, cx)


def y_channel_pattern(solver, col_parity, v=1.0, size_m=None):
    """
    One Y channel (a pad ROW, measuring y), placed at y = 0.

    Its pads are the columns of the given parity — 20 of them inside the
    31.2 mm box, on a 1.56 mm pitch. `col_parity` is (row index) mod 2, since
    a row owns the columns i with (i + row) even.
    """
    cols = np.arange(col_parity, C.PAD_N, 2)
    xs = np.unique(np.round(pad_x(cols), 12))
    return _pads_at(solver, xs, [0.0], v=v, size_m=size_m)


def x_channel_pattern(solver, col, v=1.0, size_m=None):
    """
    One X channel (a pad COLUMN, measuring x), at column index `col`.

    In a y box of exactly 1.56 mm this is a single pad; the periodic images
    supply the other 255. In a taller box it is drawn as an explicit comb.

    Row registration is NOT free to choose. Column i owns the rows j with
    (i + j) ODD, while the Y kernels are built with their row at y = 0 owning
    the columns of the matching parity. So column i's pads sit on the rows of
    parity (i + 1) mod 2 — i.e. even columns are offset by one row pitch from
    odd ones. Get this wrong and the X and Y channels overlap on some pads and
    miss others, which the sum rule catches immediately.
    """
    x0 = pad_x(col)
    n = int(round(solver.ly / (2 * C.PAD_PITCH_M)))
    y_off = ((col + 1) % 2) * C.PAD_PITCH_M
    ys = (np.arange(n) - n // 2) * (2 * C.PAD_PITCH_M) + y_off
    return _pads_at(solver, [x0], ys, v=v, size_m=size_m)


# ── Kernel construction ──────────────────────────────────────────────────────

def _solver(rho_s, d_k, nx, ny, ly, **kw):
    return WeightingSolver(rho_s_ohm_sq=rho_s, d_kapton_m=d_k, nx=nx, ny=ny,
                           ly_m=ly, phase_m=ESL_PHASE_CENTERED_M,
                           tau_drain_s=None, **kw)


def y_kernels(rho_s, d_k, times, nx=3120, ny=512, verbose=True, size_m=None):
    """
    The 2 distinct Y-channel kernels, G_Y[parity](t, y0, x0).

    Value = charge induced on that channel by a unit point charge sitting on
    the ESL at (x0, y0) since t = 0, with the channel's row at y = 0
    (reciprocity, plan §3).
    """
    s = _solver(rho_s, d_k, nx, ny, Y_BOX_M)
    out = {}
    for par in (0, 1):
        t0 = time.time()
        out[par] = s.solve(y_channel_pattern(s, par, size_m=size_m),
                           times).astype(np.float32)
        if verbose:
            print(f"    Y kernel parity {par}: {time.time()-t0:6.1f} s", flush=True)
    return s, out


def x_ny(ny_y):
    """
    y samples for the X box that put the two boxes on the SAME y grid.

    The Y box is 32 comb periods long, so matching the sample spacing means
    ny_X = ny_Y / 32. Without this the sum rule compares arrays whose y indices
    mean different things.
    """
    if ny_y % 32:
        raise ValueError(f"ny_Y={ny_y} must be a multiple of 32 "
                         "(the Y box is 32 x 1.56 mm)")
    return ny_y // 32


def x_kernels(rho_s, d_k, times, nx=3120, ny=32, cols=None, verbose=True,
              size_m=None):
    """
    The distinct X-channel kernels, G_X[col](t, y0, x0), on a 1.56 mm y box.

    `cols` defaults to the 40 column phases; columns c and c+40 give the same
    kernel because the whole stack is 31.2 mm periodic.
    """
    s = _solver(rho_s, d_k, nx, ny, X_BOX_M)
    cols = range(C.N_PAD_PER_SUPER) if cols is None else cols
    out = {}
    t0 = time.time()
    for c in cols:
        out[c] = s.solve(x_channel_pattern(s, c, size_m=size_m),
                         times).astype(np.float32)
    if verbose:
        print(f"    {len(out)} X kernels: {time.time()-t0:6.1f} s", flush=True)
    return s, out


def _to_y0_origin(g, ny):
    """Roll a solver array so that y index 0 means y = 0 (solver puts it at ny//2)."""
    return np.roll(g, -(ny // 2), axis=1)


def sum_over_rows(sy, gy):
    """
    Sum a Y kernel over ALL rows: sum_r G_Y[parity(r)](x0, y0 - 0.78 mm * r).

    Rows alternate parity, so the result is periodic in y0 with 1.56 mm. The
    box is an exact multiple of both pitches, so this is a strided sum with no
    interpolation. Returns (nt, ny/32, nx), y index 0 at y0 = 0.
    """
    step = int(round(C.PAD_PITCH_M / (sy.ly / sy.ny)))      # samples per row
    assert abs(step * sy.ly / sy.ny - C.PAD_PITCH_M) < 1e-12, "y box not commensurate"
    nt, ny, nx = gy[0].shape
    n_sub = 2 * step                                        # one 1.56 mm period
    tot = np.zeros((nt, n_sub, nx), dtype=np.float64)
    for par in (0, 1):
        # Channel at row r contributes G_par(y0 - r*pitch); summing over the
        # rows of one parity is a stride-n_sub sum after shifting by par rows.
        # np.roll(g, +s)[q] = g[q - s], and q - par*step is what we want.
        g = np.roll(_to_y0_origin(gy[par], ny), +par * step, axis=1)
        tot += g.reshape(nt, ny // n_sub, n_sub, nx).sum(axis=1)
    return tot


def sum_over_columns(sx, gx):
    """Sum the 40 X-channel phase kernels: charge on ALL X channels."""
    tot = np.sum([g.astype(np.float64) for g in gx.values()], axis=0)
    return _to_y0_origin(tot, sx.ny)


# ── Observables ──────────────────────────────────────────────────────────────

def prompt_sum_rule(d_k, gap=C.AMP_GAP_M, eps_r=C.KAPTON_EPS_R,
                    d_g=None, glue_eps_r=C.GLUE_EPS_R):
    """
    S(0)/C(0): the fraction of the charge's image that lands on the pad plane.

    The insulator is kapton + lamination glue in SERIES, so what enters is the
    series thickness, not the kapton alone:

        S(0)/C(0) = (eps_k/d_eff) / (eps_k/d_eff + 1/gap),
        d_eff     = d_k + d_g (eps_k/eps_g)

    Getting this wrong is not a small error and it is not cosmetic — every
    solve is certified against this number. Bare 50 µm kapton gives 0.913043;
    the real kapton+glue stack gives 0.881583 (the old 75 µm single-layer
    default gave 0.875000). Pass d_g=0.0 to reproduce the pre-2026-08-08
    value, which is REQUIRED when checking a legacy product: its stored
    sum_rule_expect was computed on the bare-kapton stack, so recomputing it
    with the glue would report a 4.5 % failure that is really a version skew.
    """
    d_eff = C.insulator_d_eff_m(
        d_k, eps_r,
        (C.GLUE_THICK_UM * 1e-6 if d_g is None else d_g), glue_eps_r)
    return (eps_r / d_eff) / (eps_r / d_eff + 1.0 / gap)


def nearest_column(x0_m):
    """Column index (mod the 40-phase superperiod) whose pad centre is nearest x0."""
    xs = pad_x(np.arange(C.N_PAD_PER_SUPER))
    dx = np.abs(((xs - np.mod(x0_m, C.SUPERPERIOD_M) + C.SUPERPERIOD_M / 2)
                 % C.SUPERPERIOD_M) - C.SUPERPERIOD_M / 2)
    return int(np.argmin(dx))


# ── ESL phase: where a point sits relative to the resistive strips ───────────
# The production sheet is built with phase_m = ESL_PHASE_CENTERED_M, which puts
# a strip CENTRE at x = 0 mod the 800 µm pitch. So the natural phase coordinate
# is x mod pitch measured from a strip centre:
#
#     strip  : phase < 275 µm  or  phase >= 525 µm      (550 µm wide)
#     gap    : 275 <= phase < 525 µm, gap centre at 400 µm
#
STRIP_CENTRE_PHASE_M = 0.0
GAP_CENTRE_PHASE_M = C.ESL_PITCH_M / 2                     # 400 µm


def esl_phase(x_m):
    """Position within one ESL pitch, measured from a STRIP CENTRE [m]."""
    return np.mod(np.asarray(x_m, dtype=float), C.ESL_PITCH_M)


def on_strip(x_m):
    """True where x sits on conductive ESL, under the production registration."""
    ph = esl_phase(x_m)
    half = C.ESL_WIDTH_M / 2
    return (ph < half) | (ph >= C.ESL_PITCH_M - half)


def column_nearest_phase(phase_want_m, parity=None):
    """
    Pad column whose centre is closest to a wanted ESL PHASE (circular metric).

    `nearest_column` answers "which pad is nearest this x", which is the WRONG
    selector whenever the question is "where does this deposit sit relative to
    the resistive strips". The 800/780 µm beat advances the phase by only 20 µm
    per pad, so the two pad centres straddling an arbitrary x can easily BOTH be
    on-strip: audit 2026-08-07 A2 found exactly that, and the §3 "in-gap"
    sharing row was a second on-strip deposit for the whole life of the product
    set. Selecting in phase cannot make that mistake.

    `parity` restricts the search to even (0) or odd (1) columns. With
    `row0_parity = 0` the Y view owns the even columns, so pinning the parity
    is what makes an on-strip/in-gap pair differ ONLY in ESL phase rather than
    also swapping which view owns the pad under the deposit.
    """
    cols = np.arange(C.N_PAD_PER_SUPER)
    if parity is not None:
        cols = cols[cols % 2 == int(parity)]
    dphi = np.abs((esl_phase(pad_x(cols)) - phase_want_m
                   + C.ESL_PITCH_M / 2) % C.ESL_PITCH_M - C.ESL_PITCH_M / 2)
    return int(cols[np.argmin(dphi)])


def charge_budget_y(sy, gy, x0_m, dmax=3, row0_parity=0):
    """
    Charge on the Y channel at row offset d from the deposit's row, for a
    deposit at (x0, y0 = the centre of row `row0_parity`).

    Returns {d: q(t)} for d = -dmax..dmax. The channel at row j uses the kernel
    whose comb sits on the columns of parity (j mod 2), read at Dy = -d*pitch.
    """
    ix = int(np.argmin(np.abs(sy.x - np.mod(x0_m, sy.lx))))
    step = int(round(C.PAD_PITCH_M / (sy.ly / sy.ny)))
    iy0 = sy.ny // 2
    return {d: gy[(row0_parity + d) % 2][:, (iy0 - d * step) % sy.ny, ix]
               .astype(float)
            for d in range(-dmax, dmax + 1)}


def charge_budget_x(sx, gx, x0_m, dmax=3, y0_m=0.0):
    """
    Charge on the X channel at column offset d from the column nearest x0.

    y0 is FOLDED into the 1.56 mm comb period, mirroring `charge_budget_y`'s
    x-fold above. The X box IS one comb period and its periodic images are the
    rest of the comb, so folding is exact, not an approximation. Without it
    (pre-2026-08-07, audit C4) an out-of-box y silently clamped to the box edge
    — `argmin` over `sx.y` has no way to say "out of range" — so a caller
    asking about y = 5 mm got the answer for y = 0.78 mm and no warning.
    """
    ix = int(np.argmin(np.abs(sx.x - np.mod(x0_m, sx.lx))))
    # sx.y is centred (it runs -ly/2 .. +ly/2), so fold into that window rather
    # than into [0, ly).
    y_fold = np.mod(y0_m + sx.ly / 2, sx.ly) - sx.ly / 2
    iy = int(np.argmin(np.abs(sx.y - y_fold)))
    c0 = nearest_column(x0_m)
    return {d: gx[(c0 + d) % C.N_PAD_PER_SUPER][:, iy, ix].astype(float)
            for d in range(-dmax, dmax + 1)}


# ── Driver ───────────────────────────────────────────────────────────────────

def _tau_1e(times, frac):
    """
    DEPRECATED — do not use this as a transport timescale. Kept only so old
    products' meta blocks stay interpretable.

    It returns when a neighbour's charge covers 1 - 1/e of the way from its
    prompt value to its maximum. With tau_drain_s=None (the production default)
    the sheet conserves charge and the kernel NEVER saturates, so "maximum" is
    just the value at the window edge and the whole quantity is a function of
    where the window was cut. It moved tau_X from 42 ns to 1169 ns on a d_k
    change that can move real timescales by at most 1.44, and briefly supported
    a wrong conclusion about rho_s.

    Use response/solver/spread.py instead: sigma^2(t) = sigma_p^2 + 2Dt from
    moments of the kernel, which needs no normalisation and no saturation.
    """
    f = np.asarray(frac, dtype=float)
    a, b = f[0], f.max()
    if b <= a:
        return float("nan")
    target = a + (1 - 1 / np.e) * (b - a)
    i = int(np.argmax(f >= target))
    return float(times[i])


def sharing_report(sy, gy, sx, gx, times, verbose=True):
    """
    The plan §9 charge-level sharing observables, for a deposit ON a resistive
    strip and one in an inter-strip GAP — the two cases V2/V4 showed behave
    completely differently, because gap charge has no DC path to ground.

    Both deposits are selected by ESL PHASE, not by absolute x (audit
    2026-08-07 A2 / Fix 2): the 800/780 µm beat advances the phase by only
    20 µm per pad, so picking "the pad nearest x = 15.6 mm" landed BOTH rows on
    a strip, and the "in-gap" row was a second on-strip deposit for the whole
    life of the pre-fix product set.

    The in-gap column is then forced to the SAME parity as the on-strip one, so
    the pair differs only in ESL phase and not also in which view owns the pad
    under the deposit (row 0 owns the even columns in Y, the odd ones in X).
    """
    dmax = 3
    out = {}
    # on-strip is free to take the best phase match; in-gap follows its parity.
    c_strip = column_nearest_phase(STRIP_CENTRE_PHASE_M)
    for name, phase_want in (("on-strip", STRIP_CENTRE_PHASE_M),
                             ("in-gap", GAP_CENTRE_PHASE_M)):
        c0 = column_nearest_phase(phase_want, parity=c_strip % 2)
        x0 = float(pad_x(c0))
        # Which view owns the pad under the deposit: row j owns column i when
        # (i + j) is even, so with j = 0 the Y view owns it iff c0 is even.
        owner = "Y" if c0 % 2 == 0 else "X"
        # The sheet conductivity the SOLVER actually sampled under the deposit:
        # 0 means the deposit really is in a gap. This is the self-evidence the
        # old report lacked.
        sig = float(W.esl_sigma_profile(np.array([x0]), sy.rho_s,
                                        phase_m=ESL_PHASE_CENTERED_M)[0])

        by = charge_budget_y(sy, gy, x0, dmax=dmax, row0_parity=0)
        bx = charge_budget_x(sx, gx, x0, dmax=dmax, y0_m=0.0)
        rec = {"x0_m": x0, "col": c0, "pad_owner": owner,
               "esl_phase_um": float(esl_phase(x0) * 1e6),
               "sigma_s_S_sq": sig,
               "in_gap": bool(sig == 0.0),
               "selector": "esl_phase"}
        for tag, b in (("X", bx), ("Y", by)):
            ds = list(range(-dmax, dmax + 1))
            q = np.array([b[d] for d in ds])            # (nd, nt), absolute charge
            # Normalise per TIME SLICE, by that view's total over the window.
            # Not by the d=0 or the prompt peak: with a checkerboard the peak
            # channel alternates with d and can be near zero at t=0, and
            # dividing by it turns an ordinary kernel into ratios of 20-100.
            tot = q.sum(axis=0)
            share = q / np.where(np.abs(tot) > 1e-30, tot, np.nan)
            rec[tag] = {
                "d": ds,
                "share_prompt": [float(v) for v in share[:, 0]],
                "share_late": [float(v) for v in share[:, -1]],
                "q_prompt": [float(v) for v in q[:, 0]],
                "q_late": [float(v) for v in q[:, -1]],
                "view_total_prompt": float(tot[0]),
                "view_total_late": float(tot[-1]),
                "tau_1e_s": _tau_1e(times, q[ds.index(1)]),
            }
        out[name] = rec

        if verbose:
            where = "IN GAP" if rec["in_gap"] else "ON STRIP"
            print(f"    -- deposit {name}: pad column {c0} at x0 = "
                  f"{x0*1e3:.2f} mm, ESL phase {rec['esl_phase_um']:.0f} µm "
                  f"-> sigma_s = {sig:.3e} S/sq ({where}), "
                  f"pad owned by the {owner} view")
            for tag in ("X", "Y"):
                r = rec[tag]
                bp = " ".join(f"{v:6.3f}" for v in r["share_prompt"])
                bl = " ".join(f"{v:6.3f}" for v in r["share_late"])
                tau = r["tau_1e_s"]
                print(f"       {tag}  d=-3..3  prompt {bp}   "
                      f"(view total {r['view_total_prompt']:.4f})")
                print(f"       {' '*len(tag)}          late {bl}   "
                      f"(view total {r['view_total_late']:.4f})   tau(1/e) "
                      f"{'n/a' if np.isnan(tau) else f'{tau*1e9:.0f} ns'}")
    return out


def _git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def check_sum_rule(rho_s, d_k, times, nx, ny, verbose=True):
    """
    The closed-form test: with PITCH-sized pads the checkerboard tiles the
    plane, all-pads-driven is a uniform Dirichlet plane, and the total induced
    charge must equal S(0)/C(0) at every time.

    Run separately from production because it needs the tiling pad size, not
    the real 0.68 mm one.
    """
    sy, gy = y_kernels(rho_s, d_k, times, nx=nx, ny=ny, verbose=verbose,
                       size_m=C.PAD_PITCH_M)
    sx, gx = x_kernels(rho_s, d_k, times, nx=nx, ny=x_ny(ny), verbose=verbose,
                       size_m=C.PAD_PITCH_M)

    # First, a pure pattern check with no solver in it: do the channels really
    # partition the plane? If this fails the kernel sum cannot be right either.
    step = int(round(C.PAD_PITCH_M / (sy.ly / sy.ny)))
    cover = np.zeros((sy.ny, sy.nx))
    for par in (0, 1):
        pat = y_channel_pattern(sy, par, size_m=C.PAD_PITCH_M)
        for r in range(par, int(round(sy.ly / C.PAD_PITCH_M)), 2):
            cover += np.roll(pat, r * step, axis=0)
    for c in range(C.N_PAD_PER_SUPER):
        # draw the X combs on the TALL grid, not by tiling the 1.56 mm box —
        # the two boxes' y origins differ by a row pitch, so tiling misregisters
        cover += x_channel_pattern(sy, c, size_m=C.PAD_PITCH_M)
    tile_err = float(np.abs(cover - 1.0).max())

    ytot = sum_over_rows(sy, gy)
    xtot = sum_over_columns(sx, gx)
    tot = ytot + xtot
    expect = prompt_sum_rule(d_k)
    err = float(np.abs(tot - expect).max() / expect)
    xfrac = float(xtot.mean() / tot.mean())
    return {"expect": expect, "err": err, "x_fraction": xfrac,
            "tiling_err": tile_err, "pass": err < 2e-2 and tile_err < 1e-12}


def run_point(rho_s, d_k, times, nx=3120, ny=512, outdir=None, verbose=True,
              sum_rule=True):
    """One (rho_s, d_k) scan point: both Y kernels, all 40 X kernels, checks."""
    print(f"\n=== rho_s = {rho_s/1e6:.1f} MΩ/sq   d_kapton = {d_k*1e6:.0f} µm ===",
          flush=True)
    sy, gy = y_kernels(rho_s, d_k, times, nx=nx, ny=ny, verbose=verbose)
    sx, gx = x_kernels(rho_s, d_k, times, nx=nx, ny=x_ny(ny), verbose=verbose)

    ytot = sum_over_rows(sy, gy)           # (nt, ny/32, nx), y0 = 0 at index 0
    xtot = sum_over_columns(sx, gx)
    tot = ytot + xtot
    # With the real 0.68 mm pads the plane does NOT tile, so this is a
    # prediction rather than a closed form: the missing part is the image
    # charge that lands on the 100 µm inter-pad gaps and on the mesh.
    expect = prompt_sum_rule(d_k)
    tbar = tot.mean(axis=(1, 2))
    xfrac = xtot.mean(axis=(1, 2)) / tbar
    print(f"    charge on all 1024 channels: prompt {tbar[0]:.4f} of unit "
          f"deposit  ({tbar[0]/expect:.4f} of the tiling-pad {expect:.3f};")
    print(f"      the balance is on the mesh and on the grounded 100 µm "
          f"inter-pad copper).  Late (t={times[-1]*1e6:.0f} µs): {tbar[-1]:.4f}")
    print(f"    X/Y charge balance, prompt: {xfrac[0]:.3f} / {1-xfrac[0]:.3f} "
          f"(exactly 0.5 is REQUIRED — at t=0 the sheet is a plain dielectric "
          f"and the checkerboard is 90°-symmetric)")
    print(f"    X/Y charge balance, late:   {xfrac[-1]:.3f} / {1-xfrac[-1]:.3f}"
          f"   [data, charge-integrated: 0.49/0.51]")

    res = {"rho_s_ohm_sq": rho_s, "d_kapton_m": d_k, "sum_rule_expect": expect,
           # The dielectric stack this product was solved with. Recorded
           # explicitly because a bare-kapton solve and a kapton+glue solve are
           # otherwise indistinguishable on disk, and they differ by 4.5 % in
           # the sum rule. Absent field == solved before 2026-08-08 == no glue.
           "d_glue_m": C.GLUE_THICK_UM * 1e-6,
           "glue_eps_r": C.GLUE_EPS_R,
           "d_insulator_eff_m": C.insulator_d_eff_m(d_k),
           "channel_capture_prompt": float(tbar[0]),
           "channel_capture_late": float(tbar[-1]),
           "x_fraction_prompt": float(xfrac[0]),
           "x_fraction_late": float(xfrac[-1]),
           "git": _git_hash(), "nx": nx, "ny": ny}

    res["sharing"] = sharing_report(sy, gy, sx, gx, times)

    if sum_rule:
        sr = check_sum_rule(rho_s, d_k, times[:6], nx=min(nx, 780),
                            ny=min(ny, 128), verbose=False)
        print(f"    [closed form] tiling pads, sum = S(0)/C(0) = {sr['expect']:.6f}"
              f"   err {sr['err']:.2e}  tiling {sr['tiling_err']:.1e}   {'PASS' if sr['pass'] else 'FAIL'}")
        res["sum_rule"] = sr

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        # The glue is IN the tag. Without it a kapton+glue solve at 50 µm and
        # the legacy bare-kapton solve at 50 µm produce byte-different kernels
        # under one filename, and the second silently overwrites the first on
        # EOS. `_g0um` is therefore not the same name as the pre-2026-08-08
        # product, which is the point: a fresh file, not an overwrite.
        d_g = C.GLUE_THICK_UM * 1e-6
        tag = f"rho{rho_s/1e6:g}M_dk{d_k*1e6:g}um_g{d_g*1e6:.0f}um"
        path = os.path.join(outdir, f"greens_comb_{tag}.npz")
        np.savez_compressed(
            path,
            t=times, x=sy.x, y_Y=sy.y, y_X=sx.y,
            G_Y_even=gy[0], G_Y_odd=gy[1],
            G_X=np.stack([gx[c] for c in range(C.N_PAD_PER_SUPER)]),
            x_cols=pad_x(np.arange(C.N_PAD_PER_SUPER)),
            meta=json.dumps(res))
        print(f"    -> {path}  ({os.path.getsize(path)/1e6:.0f} MB)")
    return sy, gy, sx, gx, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="coarse grid, one scan point — a few seconds")
    ap.add_argument("--rho", type=float, default=1.0, help="MΩ/sq")
    ap.add_argument("--dk", type=float, default=C.KAPTON_THICK_UM,
                    help="kapton µm (50 = confirmed as-built; the glue is "
                         "added on top of this by the solver)")
    ap.add_argument("--nx", type=int, default=3120)
    ap.add_argument("--ny", type=int, default=512)
    ap.add_argument("--nt", type=int, default=60)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--scan", action="store_true",
                    help="the full 4 x 3 rho_s x d_k grid (plan §3)")
    a = ap.parse_args()

    print("T2b — comb channel kernels (plan §3 electrode definition, T2b)\n"
          f"      a channel = {C.PAD_N//2} pads on "
          f"{2*C.PAD_PITCH_UM/1000:.2f} mm, NOT a solid line")

    if a.quick:
        a.nx, a.ny, a.nt = 3120, 256, 12

    times = log_times(a.nt)
    # Kapton is 50 µm, confirmed — not a scan axis since 2026-08-08. --scan is
    # now rho_s alone, which is the one parameter still genuinely unknown.
    grid = ([(r, C.KAPTON_THICK_UM * 1e-6) for r in C.RHO_S_SCAN_OHM_SQ]
            if a.scan else [(a.rho * 1e6, a.dk * 1e-6)])

    for rho_s, d_k in grid:
        run_point(rho_s, d_k, times, nx=a.nx, ny=a.ny, outdir=a.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
