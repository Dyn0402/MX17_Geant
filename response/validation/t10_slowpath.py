#!/usr/bin/env python3
"""
t10_slowpath.py — T10, the half of it that was never done: the ion's LATERAL
shape.

WHAT WAS OPEN (plan §0a P3, §7 step 6). Ions carry ~91 % of the induced charge
and climb the entire 150 µm amplification gap, but the production fast path
gives them the SURFACE kernel's lateral shape — the narrowest one there is —
frozen at z = 0, for the whole transit. `response/digitizer/ions.py` and
`kernel_lut.apply_ion_transit` both say so in their own docstrings. T10's
existing certification (`test_lut_vs_solver.py`) covers caching fidelity and
explicitly not this. It is the last known physics approximation in the chain.

WHAT THE PLAN ASKED FOR, AND WHAT THIS DOES INSTEAD. §3 point 6 specifies
Garfield++ `Sensor` + `ComponentGrid::LoadWeightingField` over "the S1 Psi time
slices". Those slices do not exist and never did — `WeightingSolver.solve()`
returns the z = 0 plane and `greens_comb_*.npz` is z = 0 by construction — so
read literally the task begins with an S1 re-solve.

It does not need one. In model W1 the gas gap is charge-free and bounded above
by a GROUNDED plane, so Psi at height z follows from Psi at z = 0 exactly, mode
by mode (`response/solver/zextend.py`, verified there against an independent
real-space finite-difference Laplace solve, which converges to it at order 2).
Given that, the slow path is a quadrature, and doing it directly is both exact
and free of the two failure modes the spec's route carries: resampling a
volumetric grid through `ComponentGrid`'s interpolation, and the unmerged local
`ComponentGrid.cc` region-flag patch that T6/T7 needed, whose applicability to
the weighting-field path is unverified. The deviation is deliberate and is
recorded in the plan.

Independence is preserved where it matters: this file shares no code with the
LUT. It does not window, stride, decimate, re-index by channel offset, resample
onto a uniform grid before differentiating, or run the cumsum running mean. It
reads the S1 product and integrates the extended Ramo-Shockley expression.

THE MEASUREMENT. For a charge pair born at height z_b, the electron falls to
z = 0 and the ion climbs to z = g. With g_n(z, tau) = Psi_n at the deposit's
(x, y) and height z, the induced charge on channel n is

    Q_n(t) = q [ g_n(z_b,t) - g_n(0,t) ]              <- electron, prompt
             - q Int_0^t v_i dg_n/dz(z_i(t'), t-t') dt'   <- ion, over transit

(signs fixed so Q is positive on the LUT's convention; the two "placement"
terms of the pair cancel identically at t = 0, as they must). The fast path
replaces both g_n(z, .) by g_n(0, .) scaled by the flat weighting 1 - z/g.

THE INTERNAL CHECK THAT MAKES THIS TRUSTWORTHY: run the slow path with every
k != 0 mode of the kernel deleted, and it must reproduce the fast path to
round-off, because 1 - z/g IS the k -> 0 limit of the z factor. It runs on
every deposit, exercises the whole trajectory quadrature and the charge
bookkeeping against a path that shares none of it, and fails loudly on a sign,
a factor of dt or an off-by-one in the lag.

WHAT THAT SELF-CHECK CANNOT SEE, learned the hard way on 2026-08-08: it
compares two paths that read the SAME kernel table, so any error in building
that table is invisible to it. A seconds-vs-ns mix-up in `to_uniform` flattened
the kernel's entire time dependence and the self-check still passed at 6e-15.
`to_uniform` now refuses the ambiguous call rather than extrapolating.

WHY QUADRATURE AND NOT "O(100) EVENTS". The plan asks for O(100) events per
point. The trajectories here are parameterized, so the only stochastic variable
is the birth height z_b — the Polya gain is one multiplicative constant per
avalanche and cancels out of every quantity below. Averaging over z_b by
quadrature on the MEASURED alpha_z histogram is the n -> infinity limit of that
sampling, at lower cost and with no Monte-Carlo noise to mistake for a residual.

    python3 -m response.validation.t10_slowpath --kernel <greens_comb.npz>
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from ..common import constants as C
from ..solver import kernels as K
from ..solver.zextend import ZSlicer, certify_set
from ..digitizer import ions as ION

TOL = 0.02              # plan T10: < 2 % waveform residual, fast vs slow
DEFAULT_CALIB = "~/x17/response_sim/avalanche/aval_calib.json"
DEFAULT_MESH_V = 490.0  # det3 bench point (plan §0a P1)


# ── the avalanche, as both paths see it ──────────────────────────────────────

def birth_heights(calib_pt=None, gap=C.AMP_GAP_M, n_fallback=48):
    """
    Quadrature nodes and weights for the ion/electron birth height z_b.

    From the S3 campaign's measured `alpha_z_hist` when one is supplied — that
    histogram is the thing that corrected the analytic model's 5 µm guess to a
    measured 13.8 µm mean and moved f_ion from 0.967 to 0.908, so it is the
    right input. Falls back to the exponential exp(-z/lambda) with lambda from
    the same mean, which is what the profile is.
    """
    if calib_pt and calib_pt.get("alpha_z_hist"):
        zh = calib_pt["alpha_z_hist"]
        w = np.asarray(zh["counts"], dtype=float)
        e = np.asarray(zh["edges"], dtype=float) * 1e-6
        z = 0.5 * (e[:-1] + e[1:])
        keep = w > 0
        z, w = z[keep], w[keep]
        return z, w / w.sum(), "S3 measured alpha_z_hist"
    lam = 14.1e-6
    z = (np.arange(n_fallback) + 0.5) * gap / n_fallback
    w = np.exp(-z / lam)
    return z, w / w.sum(), f"analytic exp(-z/{lam*1e6:.1f} µm)"


def ion_speed_m_ns(gap=C.AMP_GAP_M, v_mesh=DEFAULT_MESH_V,
                   mu_cm2_vs=ION.MU_ION_CM2_VS):
    """Ion drift speed [m/ns]: v = mu E = mu V / g, constant across the gap."""
    return gap / (ION.ion_transit_ns(gap, v_mesh, mu_cm2_vs))


def flat_longitudinal(t_ns, zb, wb, gap=C.AMP_GAP_M, v_ion=None):
    """
    The FAST path's delivery profile h(t), built from the same avalanche model.

    Deliberately reconstructed here rather than read from the S3 `i_elec +
    i_ion` product: the point of the comparison is to isolate the lateral
    approximation, so the two paths must differ in NOTHING else. Production's
    measured h carries the same physics plus Garfield's own flat weighting, and
    substituting it would fold an unrelated (and separately certified)
    difference into the residual.

    Each pair contributes z_b/g promptly and 1 - z_b/g spread flat over its
    transit (g - z_b)/v_i. The ion's delivery RATE is then v_i/g for every z_b —
    independent of birth height — so the ensemble's ion term is just
    (v_i/g) P(transit still running), which is what is built below.

    NOT renormalised. sum(h)*dt comes out at 1 by construction, and forcing it
    would hide a window that failed to contain the transit; the caller asserts
    on it instead.
    """
    v_ion = ion_speed_m_ns(gap) if v_ion is None else v_ion
    dt = t_ns[1] - t_ns[0]
    h = np.zeros_like(t_ns)
    h[0] += float((wb * zb / gap).sum()) / dt       # the electron delta
    for z_i, w_i in zip(zb, wb):
        T = (gap - z_i) / v_ion                     # ns
        n = min(int(np.floor(T / dt)), len(h))
        rate = w_i * v_ion / gap                    # per ns, indep. of z_i
        h[:n] += rate
        if 0 <= n < len(h):
            h[n] += rate * (T - n * dt) / dt        # partial last bin
    return h


# ── time-axis handling ───────────────────────────────────────────────────────

def to_uniform(t_src, Y, t_dst):
    """
    Piecewise-linear resample of the S1 log time axis onto a uniform grid.

    Same interpolation the LUT uses, because the S1 product IS sampled on that
    log axis and any comparison has to put both paths on one grid. Applied to
    CHARGE, never to current: G has a step at t = 0 and differentiating after
    resampling (as here) keeps that step in one sample, where the LUT's
    `_to_current` also puts it.

    BOTH AXES ARE IN SECONDS. The S1 product's `t` is seconds (0, 1e-10 ..
    1e-5) while everything else in this module works in ns, and the first
    version of this file passed the ns grid straight in. `searchsorted` then
    put every sample past the end of the source axis, the weight clipped to 1,
    and the whole time dependence collapsed to "prompt at sample 0, fully
    relaxed forever after" — silently, and invisibly to the k = 0 self-check,
    because BOTH paths read the same corrupted table. `kernel_lut.py:73`
    converts at exactly this boundary; this now refuses to guess.
    """
    t_src = np.asarray(t_src, float)
    t_dst = np.asarray(t_dst, float)
    if t_dst.size and t_src.size and t_dst.max() > t_src.max() * (1 + 1e-9):
        raise ValueError(
            f"requested t up to {t_dst.max():g} s but the S1 product only "
            f"covers {t_src.max():g} s — a unit mix-up (ns vs s), or a window "
            "longer than the solve. Extrapolating here silently flattens the "
            "kernel's whole time dependence.")
    Y = np.asarray(Y, float)
    j = np.clip(np.searchsorted(t_src, t_dst, side="right") - 1,
                0, len(t_src) - 2)
    w = np.clip((t_dst - t_src[j]) / (t_src[j + 1] - t_src[j]), 0.0, 1.0)
    w = w.reshape((-1,) + (1,) * (Y.ndim - 1))       # broadcast over (nz, ...)
    lo, hi = Y[j], Y[j + 1]
    return lo + (hi - lo) * w


# ── the two induction paths ──────────────────────────────────────────────────

class HeightTable:
    """
    One channel's kernel and its z-derivative, on the transit's own lattice.

    THE LATTICE IS THE POINT. The ion's height advances by exactly v_i*dt per
    time step, so if the birth heights are snapped to that same lattice, every
    height either path ever needs is a column of ONE table — computed once per
    channel instead of once per (channel, birth height). Snapping moves a birth
    height by at most v_i*dt/2 = 0.25 µm against a 2.5 µm histogram bin, i.e.
    well inside the input's own resolution.
    """

    def __init__(self, radial, t_src, t_ns, gap=C.AMP_GAP_M, v_ion=None,
                 sigma_m=0.0, k0_only=False):
        self.gap, self.dt = gap, float(t_ns[1] - t_ns[0])
        self.v_ion = ion_speed_m_ns(gap) if v_ion is None else v_ion
        self.dz = self.v_ion * self.dt
        self.nz = int(np.floor(gap / self.dz)) + 1      # top node <= gap
        self.z = np.arange(self.nz) * self.dz
        kw = dict(sigma_m=sigma_m, k0_only=k0_only)
        # t_src is SECONDS (the S1 product's own axis); t_ns is ns.
        t_s = np.asarray(t_ns, float) * 1e-9
        self.val = to_uniform(t_src, radial.values(self.z, **kw), t_s)
        self.dval = to_uniform(t_src, radial.values(self.z, deriv=True, **kw),
                               t_s)

    def iz(self, z_m):
        return int(np.clip(round(float(z_m) / self.dz), 0, self.nz - 1))


def snap_heights(zb, gap=C.AMP_GAP_M, v_ion=None, dt_ns=1.0):
    """
    Move birth heights onto the transit lattice, ONCE, for both paths.

    Both paths must see the same z_b or the k=0 identity is only approximate
    and the self-check degrades to O(dt/T) ~ 3e-3 — which is the same size as
    some of the effects being measured, so it has to be exact rather than
    small. The move is at most v_i*dt/2 = 0.25 µm against a 2.5 µm histogram
    bin.
    """
    v_ion = ion_speed_m_ns(gap) if v_ion is None else v_ion
    dz = v_ion * dt_ns
    nz = int(np.floor(gap / dz)) + 1
    return np.clip(np.round(np.asarray(zb, float) / dz), 0, nz - 1) * dz


def induce_slow(tab, t_ns, zb, wb):
    """
    Q_n(t) with the ion carrying its TRUE lateral shape at every height.

    Sign convention matches the LUT's (charge positive). The pair's two
    placement terms cancel identically, so the electron reduces to the step
    g(0,t) - g(z_b,t) and everything else is the ion's transit integral.
    """
    nt = len(t_ns)
    dt = tab.dt
    out = np.zeros(nt)
    for z_i, w_i in zip(zb, wb):
        jb = tab.iz(z_i)
        # --- electron: instantaneous on this grid ------------------------
        # It crosses z_b ~ 14 µm at ~100 µm/ns, i.e. in ~0.15 ns, against a
        # 1 ns grid and a sheet whose fastest relaxation is ~8 ns.
        out += w_i * (tab.val[:, 0] - tab.val[:, jb])
        # --- ion: the SAME quadrature the fast path uses ------------------
        # Rectangle bins of dt with the last one partial — not a trapezoid.
        # Not a stylistic choice: `flat_longitudinal` builds h as exactly that
        # rectangle, so matching it makes the k=0 comparison an algebraic
        # identity instead of an agreement to O(dt/T), and makes the residual
        # this file reports second-order clean (both paths carry the same
        # O(dt) quadrature error, which cancels in the difference).
        T = (tab.gap - tab.z[jb]) / tab.v_ion
        n = min(int(np.floor(T / dt)), tab.nz - 1 - jb)
        wq = np.full(n + 1, dt)
        wq[n] = T - n * dt
        acc = np.zeros(nt)
        for j in range(n + 1):
            if wq[j] <= 0.0:
                continue
            acc[j:] += wq[j] * tab.dval[:nt - j, jb + j]
        out -= w_i * tab.v_ion * acc
    return out


def f_electron_eff(tab, zb, wb):
    """
    The electron's PROMPT share of this channel's charge, as it really is.

    The fast path gives every channel one number, f_e = <z_b>/g, because the
    flat weighting says a pair born at z_b splits z_b/g : 1 - z_b/g. The true
    split is 1 - Psi_n(z_b)/Psi_n(0) — a ratio of two DIFFERENTLY FILTERED
    kernels, not a rescaling of one — so it is per channel. Returned alongside
    the flat value; their disagreement is the whole of T10's lateral effect
    expressed as one interpretable number instead of a residual.
    """
    g0 = tab.val[0, 0]
    if g0 == 0:
        return float("nan"), float("nan")
    num = sum(w * (g0 - tab.val[0, tab.iz(z)]) for z, w in zip(zb, wb))
    return float(num / g0), float((wb * zb).sum() / tab.gap)


def induce_fast(tab, t_ns, zb, wb):
    """
    Q_n(t) the production way: the surface kernel convolved with h(t).

    This is `kernel_lut.apply_longitudinal` reduced to its physics — a linear
    convolution of the z = 0 kernel with the delivery profile — with none of
    the LUT's caching machinery, which `test_lut_vs_solver` certifies
    separately.
    """
    h = flat_longitudinal(t_ns, zb, wb, gap=tab.gap, v_ion=tab.v_ion)
    return np.convolve(tab.val[:, 0], h * tab.dt)[:len(t_ns)]


# ── loading the probes out of an S1 product ──────────────────────────────────

def _y_window(y, half_mm):
    keep = np.abs(y) <= half_mm * 1e-3
    i0 = int(np.argmax(keep))
    return slice(i0, i0 + int(keep.sum()))


def build_probes(path, x0_m, dmax=3, row0_parity=0, y0_m=0.0,
                 y_half_mm=8.0, dk_g=1e-3, verbose=True, y_decimate=1):
    """
    RadialProbes for the 2*dmax+1 channels of each view at one deposit.

    Channel indexing is taken from `response.solver.kernels` (charge_budget_y /
    charge_budget_x) and NOT re-derived: four separate 2026-08-07 bugs came from
    ad-hoc summation over the S1 arrays disagreeing with the solver's own
    helpers, which were right.
    """
    out = {"X": {}, "Y": {}}
    with np.load(path) as d:
        meta = json.loads(str(d["meta"]))
        t_src = np.asarray(d["t"], float)
        x = np.asarray(d["x"], float)
        y_Y, y_X = np.asarray(d["y_Y"], float), np.asarray(d["y_X"], float)
        # `y_decimate` throws away every other y sample, which is the ONLY
        # controlled way to ask whether the answer is limited by the product's
        # own y resolution: the ny = 512 and ny = 1024 products differ in the
        # dielectric stack as well, so comparing them confounds two things.
        if y_decimate > 1:
            y_Y, y_X = y_Y[::y_decimate], y_X[::y_decimate]
        ix = int(np.argmin(np.abs(x - np.mod(x0_m, C.SUPERPERIOD_M))))

        # --- Y view: two parity kernels, read at y = -d * pad pitch ------
        ny = len(y_Y)
        step = int(round(C.PAD_PITCH_M / (y_Y[1] - y_Y[0])))
        iy0 = ny // 2
        sl = _y_window(y_Y, y_half_mm)
        margin = y_half_mm * 1e-3 - dmax * C.PAD_PITCH_M
        for par in (0, 1):
            G = np.asarray(d["G_Y_even" if par == 0 else "G_Y_odd"]
                           )[:, ::y_decimate, :][:, sl, :]
            zs = ZSlicer(G, x, y_Y[sl], periodic_y=False, y_margin_m=margin)
            for dd in range(-dmax, dmax + 1):
                if (row0_parity + dd) % 2 != par:
                    continue
                iy_full = (iy0 - dd * step) % ny
                out["Y"][dd] = zs.probe(iy_full - sl.start, ix, dk_g=dk_g)
            zs.free()
            del G, zs
            if verbose:
                print(f"      Y parity {par} probes built", flush=True)

        # --- X view: one column kernel per offset ------------------------
        c0 = K.nearest_column(x0_m)
        ly_X = len(y_X) * (y_X[1] - y_X[0])
        y_fold = np.mod(y0_m + ly_X / 2, ly_X) - ly_X / 2
        iyx = int(np.argmin(np.abs(y_X - y_fold)))
        GX = d["G_X"]
        for dd in range(-dmax, dmax + 1):
            c = (c0 + dd) % C.N_PAD_PER_SUPER
            zs = ZSlicer(np.asarray(GX[c])[:, ::y_decimate, :], x, y_X,
                         periodic_y=True)
            out["X"][dd] = zs.probe(iyx, ix, dk_g=dk_g)
        del GX
        if verbose:
            print(f"      X columns {c0}+d probes built", flush=True)
    return out, t_src, meta, c0


# ── the campaign ─────────────────────────────────────────────────────────────

_SHAPERS = {}


def _shape(cur, dt_ns):
    """
    DREAM-shaped response to each row of `cur`.

    Row by row on purpose: `DreamShaper.apply` truncates its kernel with
    `self.h[:len(i)]`, and on a 2-D input `len(i)` is the CHANNEL count, so a
    batched call would silently convolve with a 7-sample shaper.
    """
    from ..dream.shaper import DreamShaper
    sh = _SHAPERS.setdefault(dt_ns, DreamShaper(dt_ns=dt_ns))
    return np.stack([sh.apply(row, dt_ns=dt_ns) for row in np.atleast_2d(cur)])


def compare_at(probes, t_src, t_ns, zb, wb, sigma_m=0.0, gap=C.AMP_GAP_M,
               v_ion=None):
    """
    Slow vs fast for every channel of both views at one deposit.

    ON NORMALISATION, WHICH DECIDES THE VERDICT AND SO HAS TO BE ARGUED.
    A residual needs a scale, and there are two defensible ones:

      * `_view`  — the largest channel IN THAT VIEW. Right for a SHARING
        question (c1 is defined within a view), wrong for a waveform one: for
        an in-gap deposit the Y view holds 0.5 % of the event, so a 26 % Y
        residual by this measure is 0.14 % of anything the DAQ records. The
        first version of this file used it and reported 61 %, which is a ratio,
        not a residual.
      * `_event` — the largest channel in EITHER view. This is the one the
        plan's "< 2 % waveform residual" bar means, because the waveform's own
        scale is the pulse the detector actually produces.

    Both are returned. The verdict uses `_event`; `_view` and the c1 shift are
    printed next to it so the reader can see the difference rather than take
    the favourable number on trust.
    """
    raw = {}
    for view in ("X", "Y"):
        ds = sorted(probes[view])
        qs, qf, ws, wf, fes = [], [], [], [], []
        for dd in ds:
            tab = HeightTable(probes[view][dd], t_src, t_ns, gap=gap,
                              v_ion=v_ion, sigma_m=sigma_m)
            Qs = induce_slow(tab, t_ns, zb, wb)
            Qf = induce_fast(tab, t_ns, zb, wb)
            fes.append(f_electron_eff(tab, zb, wb)[0])
            qs.append(Qs[-1])
            qf.append(Qf[-1])
            # current: backward difference with the prompt step kept in one
            # sample, exactly as kernel_lut._to_current does
            ws.append(np.diff(Qs, prepend=0.0) / (t_ns[1] - t_ns[0]))
            wf.append(np.diff(Qf, prepend=0.0) / (t_ns[1] - t_ns[0]))
        dt = t_ns[1] - t_ns[0]
        raw[view] = dict(
            d=ds, qs=np.array(qs), qf=np.array(qf),
            ws=np.array(ws), wf=np.array(wf),
            sh=_shape(np.ascontiguousarray(np.array(ws)), dt),
            hf=_shape(np.ascontiguousarray(np.array(wf)), dt),
            fes=np.array(fes))

    ev_w = max(float(np.abs(r["wf"]).max()) for r in raw.values())
    ev_h = max(float(np.abs(r["hf"]).max()) for r in raw.values())
    ev_q = max(float(np.abs(r["qf"]).max()) for r in raw.values())

    res = {}
    for view, r in raw.items():
        pf, ps = np.abs(r["hf"]).max(axis=1), np.abs(r["sh"]).max(axis=1)
        ip, im = r["d"].index(1), r["d"].index(-1)
        res[view] = {
            "d": r["d"], "q_slow": r["qs"], "q_fast": r["qf"],
            "share_slow": r["qs"] / r["qs"].sum(),
            "share_fast": r["qf"] / r["qf"].sum(),
            "f_e_eff": r["fes"],
            "f_e_flat": float((wb * zb).sum() / C.AMP_GAP_M),
            "q_prompt": r["wf"][:, 0] * (t_ns[1] - t_ns[0]),
            "peak_slow": ps, "peak_fast": pf,
            # c1: the plan §9 headline, on the shaped amplitude the DAQ sees
            "c1_fast": float((pf[ip] + pf[im]) / pf.sum()),
            "c1_slow": float((ps[ip] + ps[im]) / ps.sum()),
            "resid_charge_event": float(np.abs(r["qs"] - r["qf"]).max() / ev_q),
            "resid_wave_event": float(np.abs(r["ws"] - r["wf"]).max() / ev_w),
            "resid_shaped_event": float(np.abs(r["sh"] - r["hf"]).max() / ev_h),
            "resid_prompt_event": float(
                np.abs(r["ws"][:, 0] - r["wf"][:, 0]).max() / ev_w),
            "resid_charge_view": float(np.abs(r["qs"] - r["qf"]).max()
                                       / np.abs(r["qf"]).max()),
            "resid_shaped_view": float(np.abs(r["sh"] - r["hf"]).max()
                                       / np.abs(r["hf"]).max()),
        }
        res[view]["c1_rel"] = (res[view]["c1_slow"] / res[view]["c1_fast"] - 1.0
                               if res[view]["c1_fast"] else float("nan"))
    return res


def selfcheck(probes, t_src, t_ns, zb, wb, gap=C.AMP_GAP_M, v_ion=None):
    """
    The slow path, with lateral structure deleted, MUST equal the fast path.

    Both paths are then computing the same physics by completely different
    routes: the fast one convolves the surface kernel with an analytic delivery
    profile, the slow one integrates dPsi/dz along a trajectory and adds an
    electron step. They agree only if the trajectory quadrature, the lag
    bookkeeping, the electron term's sign and the dt factors are all right, so
    this is the test that makes a NON-zero residual elsewhere believable.
    """
    worst = 0.0
    for view in ("X", "Y"):
        for dd, pr in probes[view].items():
            tab = HeightTable(pr, t_src, t_ns, gap=gap, v_ion=v_ion,
                              k0_only=True)
            Qs = induce_slow(tab, t_ns, zb, wb)
            Qf = induce_fast(tab, t_ns, zb, wb)
            sc = np.abs(Qf).max()
            if sc > 0:
                worst = max(worst, float(np.abs(Qs - Qf).max() / sc))
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True, help="greens_comb_*.npz")
    ap.add_argument("--calib", default=DEFAULT_CALIB)
    ap.add_argument("--mesh-v", type=float, default=DEFAULT_MESH_V)
    ap.add_argument("--dt-ns", type=float, default=1.0)
    ap.add_argument("--t-max-ns", type=float, default=3000.0,
                    help="match kernel_lut.CombKernelLUT.T_MAX_NS_DEFAULT")
    ap.add_argument("--dmax", type=int, default=3)
    ap.add_argument("--out", default=None, help="write a JSON result block")
    ap.add_argument("--y-decimate", type=int, default=1,
                    help="halve the product's y sampling, to see whether the "
                         "answer is resolution-limited")
    a = ap.parse_args()

    t_ns = np.arange(0.0, a.t_max_ns + a.dt_ns, a.dt_ns)
    v_ion = ion_speed_m_ns(C.AMP_GAP_M, a.mesh_v)

    # --- the avalanche both paths share ---------------------------------
    pt = None
    cal = os.path.expanduser(a.calib)
    if os.path.exists(cal):
        pts = json.load(open(cal))["points"]
        key = next((k for k in pts if k.endswith(f"@{int(a.mesh_v)}V")), None)
        pt = pts.get(key) if key else None
    zb, wb, zsrc = birth_heights(pt)
    # ONE snap, shared by both paths (see snap_heights).
    zb = snap_heights(zb, v_ion=v_ion, dt_ns=a.dt_ns)

    print("T10 slow path — does the ion's LATERAL shape matter?\n")
    print(f"  kernel        {os.path.basename(a.kernel)}")
    print(f"  birth heights {zsrc}, <z> = "
          f"{(zb*wb).sum()*1e6:.2f} µm over {len(zb)} nodes")
    print(f"  ion speed     {v_ion*1e6:.4f} µm/ns at {a.mesh_v:.0f} V "
          f"-> transit {C.AMP_GAP_M/v_ion:.0f} ns")
    print(f"  flat f_ion    {1 - float((wb*zb).sum())/C.AMP_GAP_M:.4f}"
          f"   (what the fast path applies to every channel)")
    print(f"  time grid     {a.dt_ns:g} ns to {a.t_max_ns:g} ns\n")

    h = flat_longitudinal(t_ns, zb, wb, v_ion=v_ion)
    area = float(h.sum() * a.dt_ns)
    assert abs(area - 1.0) < 1e-6, f"delivery profile area {area}, window short?"
    print(f"  delivery profile area {area:.9f}   OK (window contains the "
          f"transit)\n")

    out = {"kernel": os.path.basename(a.kernel), "mesh_V": a.mesh_v,
           "z_source": zsrc, "z_mean_um": float((zb*wb).sum()*1e6),
           "dt_ns": a.dt_ns, "t_max_ns": a.t_max_ns, "deposits": {}}
    ok = True

    # Two deposits selected by ESL PHASE, same convention as
    # kernels.sharing_report — never by absolute x (audit 2026-08-07 A2).
    c_strip = K.column_nearest_phase(K.STRIP_CENTRE_PHASE_M)
    for name, phase in (("on-strip", K.STRIP_CENTRE_PHASE_M),
                        ("in-gap", K.GAP_CENTRE_PHASE_M)):
        c0 = K.column_nearest_phase(phase, parity=c_strip % 2)
        x0 = float(K.pad_x(c0))
        print(f"  == deposit {name}: column {c0}, x0 = {x0*1e3:.3f} mm, "
              f"ESL phase {K.esl_phase(x0)*1e6:.0f} µm ==")
        probes, t_src, meta, _ = build_probes(a.kernel, x0, dmax=a.dmax,
                                              y_decimate=a.y_decimate)

        # Certified over the whole X set against their common largest
        # amplitude — see `certify_set` for why per-probe is the wrong
        # normalisation. Only X: the Y slicers are freed after their probes are
        # built, because holding their transforms alive OOMs a 16 GB host.
        err = certify_set(list(probes["X"].values()), tol=3e-5)
        print(f"      radial-bin fidelity vs the direct reduction: "
              f"{err:.2e}   OK")
        sc = selfcheck(probes, t_src, t_ns, zb, wb, v_ion=v_ion)
        good = sc < 1e-6
        ok &= good
        print(f"      k=0-only self-check (slow must == fast): {sc:.2e}"
              f"   {'OK' if good else 'FAIL'}\n")

        rec = {"col": c0, "x0_m": x0, "sigma": {}}
        for tag, sig in (("point", 0.0), ("aval 34 µm", 34e-6),
                         ("z=2mm, 213 µm", 213e-6),
                         ("z=10mm, 477 µm", 477e-6),
                         ("z=30mm, 826 µm", 826e-6)):
            r = compare_at(probes, t_src, t_ns, zb, wb, sigma_m=sig,
                           v_ion=v_ion)
            rec["sigma"][tag] = {
                v: {k: (val.tolist() if isinstance(val, np.ndarray) else val)
                    for k, val in r[v].items()} for v in ("X", "Y")}
            print(f"      smear {tag:<16} " + "  |  ".join(
                f"{v}: shaped {r[v]['resid_shaped_event']:6.2%} "
                f"charge {r[v]['resid_charge_event']:6.2%} "
                f"c1 {r[v]['c1_fast']:.3f}->{r[v]['c1_slow']:.3f}"
                for v in ("X", "Y")))

        # the per-channel detail at the production-relevant smear
        r = rec["sigma"]["aval 34 µm"]
        for view in ("X", "Y"):
            v = r[view]
            print(f"\n      {view} view, avalanche smear only "
                  f"(d = {v['d'][0]}..{v['d'][-1]}):")
            print("        q fast   " + " ".join(f"{q:8.5f}" for q in v["q_fast"]))
            print("        q slow   " + " ".join(f"{q:8.5f}" for q in v["q_slow"]))
            print("        share f  " + " ".join(f"{q:8.4f}" for q in v["share_fast"]))
            print("        share s  " + " ".join(f"{q:8.4f}" for q in v["share_slow"]))
            print("        q prompt " + " ".join(f"{q:8.5f}" for q in v["q_prompt"]))
            print(f"        f_e true " + " ".join(f"{q:8.4f}" for q in v["f_e_eff"])
                  + f"   <- the fast path uses {v['f_e_flat']:.4f} for all")
            print("        (f_e is a ratio to the PROMPT amplitude on the row "
                  "above; where that is\n         near zero the ratio is "
                  "large and the absolute error is not)")
        print()
        out["deposits"][name] = rec

    def worst(metric):
        return max(rec["sigma"][s][v][metric]
                   for rec in out["deposits"].values()
                   for s in rec["sigma"] for v in ("X", "Y"))

    for m in ("resid_wave_event", "resid_prompt_event", "resid_charge_event",
              "resid_shaped_event", "resid_shaped_view", "resid_charge_view"):
        out["worst_" + m] = worst(m)
    out["worst_c1_rel"] = max(
        abs(rec["sigma"][s][v]["c1_rel"])
        for rec in out["deposits"].values()
        for s in rec["sigma"] for v in ("X", "Y"))
    out["tol"] = TOL
    # The bar is applied to the SHAPED waveform normalised by the EVENT's
    # largest shaped amplitude. See `compare_at` for why: the within-view
    # normalisation is a sharing ratio, not a waveform residual, and on a view
    # holding 0.5 % of the event it reports 61 % for a 0.14 % difference.
    out["pass"] = bool(out["worst_resid_shaped_event"] <= TOL)

    print(f"  worst residual vs the EVENT's largest shaped amplitude "
          f"(the bar):")
    print(f"    DREAM-shaped waveform  {out['worst_resid_shaped_event']:.2%}"
          f"      raw current {out['worst_resid_wave_event']:.2%}"
          f"   channel charge {out['worst_resid_charge_event']:.2%}")
    print(f"  same, normalised WITHIN each view (a sharing ratio, not a "
          f"residual):")
    print(f"    DREAM-shaped {out['worst_resid_shaped_view']:.2%}"
          f"   channel charge {out['worst_resid_charge_view']:.2%}")
    print(f"  worst shift in c1, the §9 observable: "
          f"{out['worst_c1_rel']:+.1%}")
    print(f"  T10 bar {TOL:.0%}  ->  "
          f"{'PASS' if out['pass'] else 'FAIL — the fast path is not certified'}")
    if not ok:
        print("\n  !! the k=0 self-check did not pass; the residuals above are "
              "not trustworthy")
    if a.out:
        json.dump(out, open(a.out, "w"), indent=1)
        print(f"\n  -> {a.out}")
    return 0 if (ok and out["pass"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
