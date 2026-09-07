#!/usr/bin/env python3
"""
figures.py — figures for the S1 solver, built from the solver's own output.

Every panel here is computed on the spot; nothing is a sketch or a cartoon.

    python3 -m response.solver.figures [--outdir design/figures/response]
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..common import constants as C
from . import wpot as W
from .wpot import WeightingSolver, pad_pattern, strip_pattern_y, stack_coeffs

# Validated categorical palette (dataviz reference instance, light mode:
# blue / orange / aqua / violet — passes the adjacent-pair CVD and
# normal-vision floors).
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
SERIES = [BLUE, ORANGE, AQUA, VIOLET]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": "#e6e6e2", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "font.size": 10, "axes.titlesize": 11,
    "legend.frameon": False, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0,
})


def _save(fig, outdir, name):
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")
    return p


# ── F1: the stack transfer functions ─────────────────────────────────────────

def fig_transfer(outdir):
    """C(k) and S(k): what the layer stack does to each lateral mode."""
    k = np.logspace(1, 6, 500)
    d_k = C.KAPTON_THICK_UM * 1e-6
    # Kapton is confirmed 50 µm, so the interesting axis is no longer its
    # thickness but the lamination glue: how much adhesive is left over the pad.
    # The bare-kapton curve is the pre-2026-08-08 model, kept as the reference
    # that shows how much the glue actually matters.
    cases = [(0.0, "bare kapton 50 µm (superseded)")]
    cases += [(C.glue_over_pad_um(t) * 1e-6,
               f"+ glue {C.glue_over_pad_um(t):.1f} µm"
               + (" (production)" if t == C.GLUE_THICK_SUPPLIED_UM else ""))
              for t in C.GLUE_THICK_SUPPLIED_UM_BRACKET]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
    for i, (d_g, lab) in enumerate(cases):
        Cv, Sv = stack_coeffs(k, C.AMP_GAP_M, d_k, d_glue_m=d_g)
        ls = "-" if d_g else "--"
        ax[0].loglog(k, Cv * 1e6, color=SERIES[i], ls=ls, label=lab)
        ax[1].semilogx(k, Sv / Cv, color=SERIES[i], ls=ls, label=lab)
        ax[0].axhline(C.c_prime(kapton_m=d_k, glue_m=d_g) * 1e6,
                      color=SERIES[i], ls=":", lw=1.2)
    for a in ax:
        for feat, lab in ((C.ESL_PITCH_UM, "800 µm ESL pitch"),
                          (C.PAD_PITCH_UM, None),
                          (C.ESL_GAP_UM, "250 µm gap")):
            a.axvline(2 * np.pi / (feat * 1e-6), color=MUTED, ls="--", lw=0.9)
            if lab:
                a.text(2 * np.pi / (feat * 1e-6), a.get_ylim()[1], " " + lab,
                       rotation=90, va="top", ha="left", fontsize=7.5,
                       color=MUTED)
    ax[0].set_xlabel("lateral wavenumber k  [rad/m]")
    ax[0].set_ylabel("areal capacitance C(k)  [µF/m²]")
    ax[0].set_title("C(k): dotted lines are the k→0 limit c′")
    ax[1].set_xlabel("lateral wavenumber k  [rad/m]")
    ax[1].set_ylabel("prompt coupling  S(k)/C(k)")
    ax[1].set_title("Prompt (t=0⁺) weighting-potential transfer")
    ax[0].legend(); ax[1].legend()
    fig.suptitle("S1 — what the MX17 layer stack does to each lateral mode",
                 fontsize=12, y=1.02)
    return _save(fig, outdir, "s1_transfer_functions.png")


# ── F2: prompt map + time evolution ──────────────────────────────────────────

def fig_evolution(outdir, rho_s=1e6, d_k=75e-6, tau=6.3e-6):
    """Ψ(x,y,t) for one driven pad: the map, and cuts along x and y."""
    s = WeightingSolver(rho_s_ohm_sq=rho_s, d_kapton_m=d_k, nx=1560, ny=512,
                        ly_m=0.0512, tau_drain_s=tau)
    vpad = pad_pattern(s, x0_m=s.lx / 2, y0_m=0.0)
    times = np.array([0.0, 3e-9, 10e-9, 30e-9, 100e-9, 300e-9, 1e-6])
    psi = s.solve(vpad, times)

    ix = np.argmin(np.abs(s.x - s.lx / 2))
    iy = s.ny // 2
    xr = (s.x - s.lx / 2) * 1e3
    yr = s.y * 1e3
    win = np.abs(xr) < 3.0
    winy = np.abs(yr) < 10.0
    winm = np.abs(xr) < 1.6
    winmy = np.abs(yr) < 2.4

    fig = plt.figure(figsize=(12.5, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=0.32)

    axm = fig.add_subplot(gs[0])
    m = np.clip(psi[0][np.ix_(winmy, winm)], 1e-4, None)
    im = axm.pcolormesh(xr[winm], yr[winmy], m, cmap="magma", shading="auto",
                        norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=m.max()))
    axm.set_xlabel("x − x_pad  [mm]"); axm.set_ylabel("y − y_pad  [mm]")
    axm.set_title("Prompt Ψ(x, y, t=0⁺)\n(one driven 0.68 mm pad)")
    axm.grid(False)
    fig.colorbar(im, ax=axm, label="Ψ / V_w")

    axx = fig.add_subplot(gs[1])
    axy = fig.add_subplot(gs[2])
    cmap = plt.get_cmap("viridis")
    for it, t in enumerate(times):
        col = cmap(it / (len(times) - 1))
        lab = "prompt" if t == 0 else (f"{t*1e9:.0f} ns" if t < 1e-6
                                       else f"{t*1e6:.0f} µs")
        # Clip at the noise floor: below ~1e-4 the curves are Gibbs ringing
        # from the step-function sigma_s(x), not signal, and plotting it on a
        # log axis reads as a filled band.
        axx.plot(xr[win], np.clip(psi[it][iy][win], 1e-4, None),
                 color=col, label=lab, lw=1.6)
        axy.plot(yr[winy], np.clip(psi[it][:, ix][winy], 1e-4, None),
                 color=col, lw=1.6)

    # Mark where the insulating ESL gaps are.
    for a, arr in ((axx, xr[win]),):
        sig = W.esl_sigma_profile(s.x, rho_s)[win]
        a.fill_between(arr, 0, 1, where=(sig == 0), transform=a.get_xaxis_transform(),
                       color=MUTED, alpha=0.12, lw=0)
    axx.set_xlabel("x − x_pad  [mm]  (ACROSS the strips)")
    axy.set_xlabel("y − y_pad  [mm]  (ALONG the strips)")
    for a in (axx, axy):
        a.set_ylabel("Ψ / V_w")
        a.set_yscale("log"); a.set_ylim(9e-5, 1.4)
    axx.set_title("Across x: gaps (shaded) block DC transport")
    axy.set_title("Along y: unimpeded resistive spreading")
    axx.legend(fontsize=7.5, ncol=2)
    fig.suptitle(f"S1 — dynamic weighting potential, ρ_s = {rho_s/1e6:.0f} MΩ/sq, "
                 f"kapton {d_k*1e6:.0f} µm, τ_drain = {tau*1e6:.1f} µs",
                 fontsize=12, y=1.03)
    return _save(fig, outdir, "s1_psi_evolution.png"), (s, psi, times)


# ── F3: the X/Y anisotropy, the headline validation target ───────────────────

def fig_anisotropy(outdir, s, psi, times):
    """
    RMS spread of Ψ along x and along y versus time.

    This is the headline: the measured detector has tau_Y ~ 410 ns against
    tau_X ~ 230 ns and kY 1.8-2.9 (plan §9). The plan requires that asymmetry
    to come out of GEOMETRY ALONE. Here it does — the strips run along y, so
    resistive transport spreads charge along y while the 250 µm insulating
    gaps block it across x.
    """
    ix = np.argmin(np.abs(s.x - s.lx / 2))
    iy = s.ny // 2
    xr = (s.x - s.lx / 2)
    yr = s.y

    def rms(prof, coord):
        w = np.abs(prof)
        w = np.maximum(w - w.min(), 0)
        if w.sum() <= 0:
            return np.nan
        return np.sqrt((w * coord ** 2).sum() / w.sum())

    sx = np.array([rms(psi[i][iy], xr) for i in range(len(times))]) * 1e3
    sy = np.array([rms(psi[i][:, ix], yr) for i in range(len(times))]) * 1e3

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
    tn = np.maximum(times, 1e-10) * 1e9
    ax[0].plot(tn, sx, "o-", color=BLUE, ms=6, label="across strips (x → X view)")
    ax[0].plot(tn, sy, "s-", color=ORANGE, ms=6, label="along strips (y → Y view)")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("time since the charge landed  [ns]")
    ax[0].set_ylabel("RMS spread of Ψ  [mm]")
    ax[0].set_title("Charge spreading is anisotropic — from geometry alone")
    ax[0].legend()

    ratio = sy / sx
    ax[1].plot(tn, ratio, "o-", color=VIOLET, ms=6)
    ax[1].axhspan(1.8, 2.9, color=AQUA, alpha=0.18, lw=0)
    ax[1].text(tn[1], 2.35, " measured k_Y = 1.8–2.9\n (blind target, plan §9)",
               fontsize=8.5, color=INK2, va="center")
    ax[1].axhline(1.0, color=MUTED, ls="--", lw=1)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("time since the charge landed  [ns]")
    ax[1].set_ylabel("spread ratio  σ_y / σ_x")
    ax[1].set_title("Y/X anisotropy vs the measured band")
    fig.suptitle("S1 — the X/Y sharing asymmetry is a geometric prediction, "
                 "not a fitted parameter", fontsize=12, y=1.03)
    _save(fig, outdir, "s1_xy_anisotropy.png")
    return sx, sy, ratio


# ── F4: the 31.2 mm beat ─────────────────────────────────────────────────────

def fig_beat(outdir, rho_s=1e6, d_k=75e-6, tau=6.3e-6, n_phase=40):
    """
    The ESL/pad pitch beat. Each of the 40 pads in one 31.2 mm superperiod sits
    at a different phase relative to the 800 µm strip pattern, so each sees a
    different kernel. Plan §1 calls this out as a prediction to look for in
    data; here is its size.
    """
    # 10 µm sampling: the phase walks by 20 µm per pad, so a 20 µm grid aliases
    # against the very quantity being measured.
    s = WeightingSolver(rho_s_ohm_sq=rho_s, d_kapton_m=d_k, nx=3120, ny=64,
                        ly_m=0.0512, tau_drain_s=tau)

    def coverage(x0):
        """Analytic overlap of a pad with the conducting strips — no grid."""
        half = C.PAD_SIZE_UM * 1e-6 / 2
        u = np.linspace(x0 - half, x0 + half, 4001)
        return float((np.mod(u, C.ESL_PITCH_M) < C.ESL_WIDTH_M).mean())

    # The PROMPT kernel carries no beat at all, by construction: at t = 0+ the
    # sheet is a plain dielectric and sigma_s(x) does not enter, so the prompt
    # response is translation invariant. The beat can only appear once
    # conduction has had time to act, so it must be measured at finite t.
    t_probe = np.array([0.0, 30e-9, 100e-9, 300e-9])
    xs, cover, keep = [], [], []
    for j in range(n_phase):
        x0 = (j + 0.5) * C.PAD_PITCH_M
        v = pad_pattern(s, x0_m=x0, y0_m=0.0)
        p = s.solve(v, t_probe)
        i0 = np.argmin(np.abs(s.x - x0))
        # Ψ at the pad centre, relative to its own prompt value: how much of
        # the signal this pad still holds after the sheet has relaxed.
        keep.append([p[it, s.ny // 2, i0] / p[0, s.ny // 2, i0]
                     for it in range(1, len(t_probe))])
        cover.append(coverage(x0))
        xs.append(x0 * 1e3)
    keep = np.array(keep)

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0))
    ax[0].plot(xs, cover, "o-", color=AQUA, ms=5)
    ax[0].set_xlabel("pad centre x within one superperiod  [mm]")
    ax[0].set_ylabel("fraction of the pad over conducting paste")
    ax[0].set_title("Pad-to-strip registration walks 20 µm per pitch\n"
                    "and repeats after exactly 31.2 mm")

    for i, t in enumerate(t_probe[1:]):
        r = keep[:, i]
        ax[1].plot(xs, r / r.mean(), "o-", color=SERIES[i], ms=5,
                   label=f"t = {t*1e9:.0f} ns  "
                         f"({(r.max()-r.min())/r.mean()*100:.1f} % p-p)")
    ax[1].set_xlabel("pad centre x within one superperiod  [mm]")
    ax[1].set_ylabel("Ψ(t) / Ψ(prompt) at the pad centre, / its mean")
    ax[1].set_title("Kernel modulation appears only once the sheet conducts\n"
                    "(the prompt kernel has exactly zero beat)")
    ax[1].legend(fontsize=8.5)
    fig.suptitle("S1 — the 31.2 mm beat between the 800 µm ESL pitch and the "
                 "780 µm pad pitch", fontsize=12, y=1.03)
    _save(fig, outdir, "s1_superperiod_beat.png")
    return xs, cover, keep


# ── F5: the drain, and the A1 problem ────────────────────────────────────────

def fig_drain(outdir):
    """Assumption A1's predicted drain time against the measured band."""
    rho = np.logspace(np.log10(2e5), np.log10(2e7), 200)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for i, L_mm in enumerate((412.0, 100.0, 31.2, 11.5)):
        tau = np.array([W.drain_time_from_strip_ends(r, L_mm * 1e-3) for r in rho])
        lab = (f"grounded every {L_mm:.0f} mm"
               + ("  ← A1 (active width)" if L_mm == 412 else ""))
        ax.loglog(rho / 1e6, tau * 1e6, color=SERIES[i], label=lab)
    ax.axhspan(5.3, 7.3, color=MUTED, alpha=0.25, lw=0)
    ax.text(0.22, 6.2, " measured τ_g = 5.3–7.3 µs", fontsize=9, color=INK,
            va="center")
    ax.axvspan(0.5, 5.0, color=AQUA, alpha=0.10, lw=0)
    ax.text(1.4, 2e5, "ρ_s scan range", fontsize=8.5, color=INK2, ha="center")
    ax.set_xlabel("ESL surface resistivity ρ_s  [MΩ/sq]")
    ax.set_ylabel("predicted global drain time τ  [µs]")
    ax.set_title("A1 (strips grounded only at the y-ends) misses the measured\n"
                 "drain time by ~3 orders of magnitude", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper left")
    return _save(fig, outdir, "s1_drain_vs_A1.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="design/figures/response")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    print("S1 figures:")
    fig_transfer(a.outdir)
    _, (s, psi, times) = fig_evolution(a.outdir)
    fig_anisotropy(a.outdir, s, psi, times)
    fig_beat(a.outdir)
    fig_drain(a.outdir)


if __name__ == "__main__":
    main()
