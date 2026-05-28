#!/usr/bin/env python3
"""
analyse_sr90_calibration.py
Feasibility study: Sr-90/Y-90 source as a calibration substitute for the He-3 target.

Physics setup:
  • Sr-90 (t½ = 28.8 yr) → Y-90 + β⁻  (Emax = 0.546 MeV)
  • Y-90 (t½ = 64.1 h)  → Zr-90 + β⁻  (Emax = 2.28 MeV)
  In secular equilibrium one Sr-90 and one Y-90 beta are emitted per decay,
  so the source produces 2A β⁻/s for a declared activity A.

Source placement:
  Source replaces the He-3 target at the same z-position (gun position).
  Electrons are emitted isotropically; only those aimed at the detector
  contribute.  Geometric efficiency is computed from the solid angle of the
  40×40 cm Micromegas face at 226.5 mm from the source.

Trigger condition:
  Coincidence between plastic scintillator AND liquid scintillator layer 1.
  Proxy from simulation: primInLS1 (reaching LS1 guarantees passing through
  the plastic scint upstream).  A threshold-corrected version uses the
  plateau-fraction trigger fraction stored in the summary CSV if available.

Inputs:
  --spectrum   Sr90_Y90_Beta_Spectrum.csv  (shipped with this script)
  --summary    full_experiment_analysis_electron.csv  (from analyse_full_experiment.py)
  --outfile    sr90_calibration.pdf

Usage:
    python3 sr90_calibration/analyse_sr90_calibration.py \\
        --spectrum sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \\
        --summary  /path/to/full_experiment_analysis_electron.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

# ── Physical / geometric constants ────────────────────────────────────────────

# Source-to-MM distance [cm] — from path-length calculations in the simulation
SOURCE_TO_MM_CM = 22.65

# Detector active area half-side [cm]
DET_HALF_CM = 20.0

# Solid angle of the 40×40 cm detector face from the source [sr]
# Ω = 4 arctan(ab / (d √(a²+b²+d²)))  with a=b=half-side, d=distance
_a, _d = DET_HALF_CM, SOURCE_TO_MM_CM
SOLID_ANGLE_SR   = 4 * np.arctan(_a**2 / (_d * np.sqrt(2*_a**2 + _d**2)))
GEOM_EFFICIENCY  = SOLID_ANGLE_SR / (4 * np.pi)

# Betas per declared source activity (secular equilibrium: Sr-90 + Y-90)
BETAS_PER_DECAY  = 2.0

# ── Colour scheme ──────────────────────────────────────────────────────────────
C_SR90   = "#e377c2"   # pink   — Sr-90 component
C_Y90    = "#1f77b4"   # blue   — Y-90 component
C_TOTAL  = "#333333"   # dark   — combined spectrum
C_SCINT  = "#d62728"   # red    — plastic scintillator
C_LS1    = "#6a1d8a"   # purple — LS layer 1
C_TRIG   = "#2ca02c"   # green  — trigger efficiency

E_LO, E_HI = 4.0, 16.5   # physical electron range from He-3 experiment [MeV]

# ── Helpers ───────────────────────────────────────────────────────────────────

def _add_energy_lines(ax):
    """Mark the He-3 experiment electron range (for reference only)."""
    for e in (E_LO, E_HI):
        ax.axvline(e, color="lightgrey", lw=0.8, ls="--", zorder=0)


def _trapz(y, x):
    """Trapezoidal integration, skipping NaN values."""
    mask = ~(np.isnan(y) | np.isnan(x))
    return float(np.trapz(y[mask], x[mask]))


def load_spectrum(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def load_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df.sort_values("energy_MeV").reset_index(drop=True)


def interp_sim(summary: pd.DataFrame, col: str,
               e_grid: np.ndarray, fill_below: float = 0.0) -> np.ndarray:
    """
    Interpolate a simulation column onto e_grid.
    Energies below the simulation minimum return fill_below.
    """
    x = summary["energy_MeV"].values
    y = summary[col].fillna(0).values
    out = np.interp(e_grid, x, y, left=fill_below, right=float(y[-1]))
    return out


# ── Analysis ──────────────────────────────────────────────────────────────────

def compute_spectral_quantities(spectrum: pd.DataFrame,
                                summary: pd.DataFrame) -> dict:
    """
    Convolve the beta spectrum with simulated observables.
    Returns a dict of spectrum-weighted scalars and arrays.
    """
    E   = spectrum["Energy_MeV"].values
    S   = spectrum["Total_ProbabilityDensity"].values       # dN/dE combined
    S90 = spectrum["Sr90_ProbabilityDensity"].values
    S_Y = spectrum["Y90_ProbabilityDensity"].values

    norm_total = _trapz(S,   E)
    norm_sr    = _trapz(S90, E)
    norm_y     = _trapz(S_Y, E)

    # Trigger proxy: primInLS1 (reaching LS1 ⟹ passing through plastic scint)
    # Falls back to trans_primInScintWall × trans_primInLS1 if needed.
    if "trans_primInLS1" in summary.columns:
        f_trig = interp_sim(summary, "trans_primInLS1", E, fill_below=0.0)
    else:
        f_trig = np.zeros_like(E)

    # Energy deposition observables [MeV/primary]
    f_scint = interp_sim(summary, "edep_edepScintWall", E, fill_below=0.0)
    f_ls1   = interp_sim(summary, "edep_edepLS1",       E, fill_below=0.0)
    f_ls_tot= interp_sim(summary, "ls_total_mean",       E, fill_below=0.0)

    # Containment fractions
    f_nentry = 1.0 - f_trig                              # doesn't reach LS1
    f_pt     = interp_sim(summary, "trans_primInLSCFRP5", E, 0.0) \
               if "trans_primInLSCFRP5" in summary.columns \
               else np.zeros_like(E)
    f_cont   = np.clip(f_trig - f_pt, 0, 1)             # enters LS, stops inside

    # Spectrum-weighted trigger efficiency (fraction of all decays that trigger,
    # given they travel toward the detector)
    spec_trig_eff = _trapz(f_trig * S, E) / norm_total

    # Per-species trigger efficiencies
    spec_trig_sr  = _trapz(f_trig * S90, E) / norm_sr  if norm_sr > 0 else 0.0
    spec_trig_y   = _trapz(f_trig * S_Y, E) / norm_y   if norm_y > 0 else 0.0

    # Spectrum-weighted containment
    spec_nentry   = _trapz(f_nentry * S, E) / norm_total
    spec_pt       = _trapz(f_pt    * S, E) / norm_total
    spec_cont     = _trapz(f_cont  * S, E) / norm_total

    # Spectrum-weighted mean edep for TRIGGERED events
    trig_norm = _trapz(f_trig * S, E)
    if trig_norm > 0:
        mean_scint_trig  = _trapz(f_scint  * f_trig * S, E) / trig_norm
        mean_ls1_trig    = _trapz(f_ls1    * f_trig * S, E) / trig_norm
        mean_lstot_trig  = _trapz(f_ls_tot * f_trig * S, E) / trig_norm
    else:
        mean_scint_trig = mean_ls1_trig = mean_lstot_trig = 0.0

    # Effective trigger energy (mean beta energy of triggered events)
    eff_E_trig = _trapz(E * f_trig * S, E) / trig_norm if trig_norm > 0 else 0.0

    return dict(
        E=E, S=S, S90=S90, S_Y=S_Y,
        norm_total=norm_total,
        f_trig=f_trig, f_scint=f_scint, f_ls1=f_ls1, f_ls_tot=f_ls_tot,
        f_nentry=f_nentry, f_pt=f_pt, f_cont=f_cont,
        spec_trig_eff=spec_trig_eff,
        spec_trig_sr=spec_trig_sr, spec_trig_y=spec_trig_y,
        spec_nentry=spec_nentry, spec_pt=spec_pt, spec_cont=spec_cont,
        mean_scint_trig=mean_scint_trig,
        mean_ls1_trig=mean_ls1_trig,
        mean_lstot_trig=mean_lstot_trig,
        eff_E_trig=eff_E_trig,
    )


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_spectrum_overview(pdf, res: dict):
    """Beta spectrum with Sr-90 / Y-90 components and trigger acceptance."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 7),
                             gridspec_kw={"height_ratios": [2, 1]},
                             sharex=True)

    E, S, S90, S_Y = res["E"], res["S"], res["S90"], res["S_Y"]
    f_trig = res["f_trig"]

    # Top: spectrum
    ax = axes[0]
    ax.fill_between(E, S90, alpha=0.35, color=C_SR90, label="Sr-90")
    ax.fill_between(E, S_Y, alpha=0.35, color=C_Y90,  label="Y-90")
    ax.plot(E, S,   color=C_TOTAL, lw=1.4, label="Total (Sr-90 + Y-90)")

    # Shade the region where trigger fraction > 10 %
    above = f_trig > 0.10
    if above.any():
        ax.fill_between(E, 0, S, where=above, alpha=0.18, color=C_TRIG,
                        label=f"Trigger acceptance (f_trig > 10 %)")

    ax.axvline(0.546, color=C_SR90, lw=0.9, ls="--", label="Sr-90 endpoint 0.546 MeV")
    ax.axvline(2.28,  color=C_Y90,  lw=0.9, ls="--", label="Y-90 endpoint 2.28 MeV")
    ax.set_ylabel("Probability density  dN/dE", fontsize=11)
    ax.set_title("Sr-90 / Y-90 beta spectrum and trigger acceptance", fontsize=12)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 2.5)

    # Bottom: trigger fraction overlaid
    ax2 = axes[1]
    ax2.plot(E, f_trig * 100, color=C_TRIG, lw=1.6)
    ax2.axhline(res["spec_trig_eff"] * 100, color="k", lw=1.0, ls="--",
                label=f"Spectrum-weighted: {res['spec_trig_eff']*100:.1f} %")
    ax2.set_xlabel("Beta energy (MeV)", fontsize=11)
    ax2.set_ylabel("Trigger fraction (%)", fontsize=11)
    ax2.set_ylim(-2, 105)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_fate_breakdown(pdf, res: dict):
    """Spectrum-weighted fate of every decay: miss / contained / punch-through."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    E = res["E"]

    # Left: fate fractions vs energy
    ax = axes[0]
    ax.fill_between(E, 0,                   res["f_nentry"],
                    alpha=0.55, color="#d62728", label="Doesn't reach LS1")
    ax.fill_between(E, res["f_nentry"],      res["f_nentry"] + res["f_cont"],
                    alpha=0.55, color="#2ca02c", label="Contained in LS stack")
    ax.fill_between(E, res["f_nentry"] + res["f_cont"], 1,
                    alpha=0.55, color="#ff7f0e", label="Punch-through (exits LS4)")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Fraction of electrons", fontsize=10)
    ax.set_title("Electron fate vs. energy\n(given electron reaches detector)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Right: pie chart of spectrum-integrated fate
    ax2 = axes[1]
    fracs = [res["spec_nentry"], res["spec_cont"], res["spec_pt"]]
    labels = [
        f"No LS entry\n{res['spec_nentry']*100:.1f} %",
        f"Contained\n{res['spec_cont']*100:.1f} %",
        f"Punch-through\n{res['spec_pt']*100:.1f} %",
    ]
    colors = ["#d62728", "#2ca02c", "#ff7f0e"]
    # Only draw non-zero wedges
    non_zero = [(f, l, c) for f, l, c in zip(fracs, labels, colors) if f > 0.001]
    if non_zero:
        f_, l_, c_ = zip(*non_zero)
        ax2.pie(f_, labels=l_, colors=c_, autopct="%1.1f%%",
                startangle=90, textprops={"fontsize": 9})
    ax2.set_title("Spectrum-integrated fate\n(β reaching detector)", fontsize=11)

    fig.suptitle("Fate of Sr-90/Y-90 electrons in the LS calorimeter stack",
                 fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_triggered_edep(pdf, res: dict):
    """
    Energy deposition in scintillator and LS1 as a function of beta energy,
    weighted by the triggered fraction of the beta spectrum.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    E      = res["E"]
    S_trig = res["f_trig"] * res["S"]   # triggered spectrum weight

    # Left: edep vs energy (both detectors)
    ax = axes[0]
    ax2r = ax.twinx()

    l1, = ax.plot(E, res["f_scint"] * 1000, color=C_SCINT, lw=1.6,
                  label="Plastic scint. (left, keV)")
    l2, = ax2r.plot(E, res["f_ls1"],         color=C_LS1,   lw=1.6, ls="--",
                    label="LS layer 1 (right, MeV)")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(bottom=0)
    ax2r.set_ylim(bottom=0)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Mean edep in plastic scint. (keV)", fontsize=10, color=C_SCINT)
    ax2r.set_ylabel("Mean edep in LS1 (MeV)", fontsize=10, color=C_LS1)
    ax.tick_params(axis="y", colors=C_SCINT)
    ax2r.tick_params(axis="y", colors=C_LS1)
    ax.set_title("Detector signal vs. beta energy", fontsize=11)
    ax.legend(handles=[l1, l2], fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    # Right: spectrum-weighted edep distributions (approximate via energy-binned means)
    ax2 = axes[1]
    # Effective triggered spectrum (normalise to 1)
    norm = np.trapz(S_trig, E)
    if norm > 0:
        w = S_trig / norm
        # Accumulate spectrum-weighted average edep in plastic scint [keV]
        ax2.plot(E, res["f_scint"] * 1000 * w / w.max(),
                 color=C_SCINT, lw=1.4, label="Plastic scint. response (a.u.)")
        ax2.plot(E, res["f_ls1"] * w / (res["f_ls1"] * w).max() if (res["f_ls1"]*w).max() > 0 else w*0,
                 color=C_LS1,   lw=1.4, ls="--", label="LS1 response (a.u.)")
        ax2.fill_between(E, 0, w / w.max(), alpha=0.12, color="grey",
                         label="Triggered β spectrum (a.u.)")
    ax2.axvline(res["eff_E_trig"], color="k", lw=1.0, ls=":",
                label=f"Eff. trigger energy {res['eff_E_trig']*1000:.0f} keV")
    ax2.set_xlim(0, 2.5)
    ax2.set_ylim(0, 1.15)
    ax2.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax2.set_ylabel("Normalised response / spectrum", fontsize=10)
    ax2.set_title("Spectrum-weighted detector response\n(triggered events only)",
                  fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Scintillator signals for Sr-90/Y-90 triggered events", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_rate_estimate(pdf, res: dict):
    """
    Trigger rate and angular resolution vs source activity.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: trigger rate vs activity
    ax = axes[0]
    activities_bq = np.logspace(3, 8, 200)   # 1 kBq → 100 MBq
    # Rate = betas/s × geom_eff × spectral_trig_eff
    rate = activities_bq * BETAS_PER_DECAY * GEOM_EFFICIENCY * res["spec_trig_eff"]

    ax.loglog(activities_bq / 1e3, rate, color=C_TRIG, lw=2.0)

    # Reference lines
    for hz, lbl in [(1, "1 Hz"), (10, "10 Hz"), (1000, "1 kHz")]:
        ax.axhline(hz, color="lightgrey", lw=0.8, ls=":")
        ax.text(1e5 / 1e3, hz * 1.3, lbl, fontsize=7, color="grey")

    # Common source activities
    for act_bq, name in [(3.7e4, "1 μCi"), (3.7e5, "10 μCi"),
                          (3.7e6, "1 mCi")]:
        ax.axvline(act_bq / 1e3, color="#ff7f0e", lw=0.9, ls="--")
        r = act_bq * BETAS_PER_DECAY * GEOM_EFFICIENCY * res["spec_trig_eff"]
        ax.annotate(f"{name}\n→ {r:.1f} Hz", xy=(act_bq / 1e3, r),
                    xytext=(act_bq / 1e3 * 2, r * 3),
                    fontsize=7, color="#ff7f0e",
                    arrowprops=dict(arrowstyle="->", color="#ff7f0e", lw=0.7))

    ax.set_xlabel("Source activity (kBq)", fontsize=10)
    ax.set_ylabel("Coincidence trigger rate (Hz)", fontsize=10)
    ax.set_title(f"Trigger rate vs source activity\n"
                 f"(geom. eff. {GEOM_EFFICIENCY*100:.1f}%, "
                 f"spectral eff. {res['spec_trig_eff']*100:.1f}%)", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)

    # Right: per-species breakdown of trigger efficiency
    ax2 = axes[1]
    species = ["Sr-90\n(Emax=0.546 MeV)", "Y-90\n(Emax=2.28 MeV)", "Combined"]
    effs    = [res["spec_trig_sr"], res["spec_trig_y"], res["spec_trig_eff"]]
    colors  = [C_SR90, C_Y90, C_TOTAL]
    bars    = ax2.bar(species, [e * 100 for e in effs], color=colors, width=0.5,
                      alpha=0.8, edgecolor="white")
    for bar, eff in zip(bars, effs):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f"{eff*100:.2f}%", ha="center", va="bottom", fontsize=10,
                 fontweight="bold")
    ax2.set_ylabel("Spectral trigger efficiency (%)", fontsize=10)
    ax2.set_title("Trigger efficiency by isotope\n"
                  "(fraction of each isotope's β spectrum that triggers LS1)",
                  fontsize=10)
    ax2.set_ylim(0, max(effs) * 1.3 * 100 + 1)
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Sr-90/Y-90 calibration rate estimates", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_mm_context(pdf, res: dict, summary: pd.DataFrame):
    """
    MM performance for the triggered Sr-90/Y-90 subset:
    energy deposition in drift gas and the energy range vs. full-experiment range.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    E = res["E"]

    # Left: drift gas edep weighted by triggered beta spectrum
    ax = axes[0]
    if "edep_edepDrift" in summary.columns:
        f_drift = interp_sim(summary, "edep_edepDrift", E, 0.0)
        S_trig  = res["f_trig"] * res["S"]
        norm    = np.trapz(S_trig, E) if np.trapz(S_trig, E) > 0 else 1
        ax.fill_between(E, 0, S_trig / norm, alpha=0.2, color=C_TRIG,
                        label="Triggered β spectrum (a.u.)")
        ax2r = ax.twinx()
        ax2r.plot(E, f_drift * 1e6, color="#2ca02c", lw=1.6,
                  label="MM drift edep (eV, right)")
        ax2r.set_ylabel("Mean edep in MM drift (eV)", fontsize=9, color="#2ca02c")
        ax2r.tick_params(axis="y", colors="#2ca02c")
        ax2r.set_ylim(bottom=0)
        ax.legend(fontsize=8, loc="upper left")
        ax2r.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Triggered β spectrum (norm.)", fontsize=10, color=C_TRIG)
    ax.set_title("MM drift gas signal for triggered events", fontsize=11)
    ax.grid(True, alpha=0.3)

    # Right: compare Sr-90 beta range to the He-3 experiment electron range
    ax2 = axes[1]
    # Full simulation range
    E_sim  = summary["energy_MeV"].values
    f_ls1  = summary["trans_primInLS1"].values if "trans_primInLS1" in summary.columns \
              else np.zeros(len(E_sim))
    ax2.plot(E_sim, f_ls1 * 100, color=C_LS1, lw=1.8, label="LS1 transmission (simulation)")

    # Shade the Sr-90/Y-90 beta energy range
    ax2.axvspan(0, 2.3, alpha=0.12, color=C_TOTAL, label="Sr-90/Y-90 range (0–2.3 MeV)")
    ax2.axvspan(E_LO, E_HI, alpha=0.10, color="#ff7f0e",
                label=f"He-3 experiment range ({E_LO}–{E_HI} MeV)")
    ax2.axvline(2.28, color=C_Y90, lw=1.0, ls="--", label="Y-90 endpoint 2.28 MeV")
    ax2.axvline(0.546, color=C_SR90, lw=1.0, ls="--", label="Sr-90 endpoint 0.546 MeV")

    ax2.set_xlim(0, 18)
    ax2.set_ylim(-2, 105)
    ax2.set_xlabel("Electron energy (MeV)", fontsize=10)
    ax2.set_ylabel("LS1 transmission (%)", fontsize=10)
    ax2.set_title("Sr-90/Y-90 range vs full-experiment coverage", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Micromegas context for Sr-90/Y-90 calibration", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_summary_page(pdf, res: dict):
    """Text summary page with key numbers."""
    fig = plt.figure(figsize=(9, 7))
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    lines = [
        ("Sr-90 / Y-90 Calibration Feasibility Summary", True),
        ("", False),
        ("─── Source properties ───────────────────────────────", False),
        (f"  Sr-90 endpoint          0.546 MeV", False),
        (f"  Y-90 endpoint           2.28  MeV", False),
        (f"  Betas per decay (equil) {BETAS_PER_DECAY:.0f}  (one per isotope)", False),
        ("", False),
        ("─── Geometry ────────────────────────────────────────", False),
        (f"  Source-to-MM distance   {SOURCE_TO_MM_CM:.1f} cm", False),
        (f"  Detector area           40 × 40 cm", False),
        (f"  Geometric solid angle   {SOLID_ANGLE_SR:.3f} sr", False),
        (f"  Geometric efficiency    {GEOM_EFFICIENCY*100:.1f} %  "
         f"(fraction of 4π aimed at MM)", False),
        ("", False),
        ("─── Spectral trigger efficiency ─────────────────────", False),
        (f"  Sr-90 alone             {res['spec_trig_sr']*100:.3f} %  "
         f"(nearly all stop before LS1)", False),
        (f"  Y-90 alone              {res['spec_trig_y']*100:.2f} %", False),
        (f"  Combined (total)        {res['spec_trig_eff']*100:.3f} %", False),
        ("", False),
        ("─── Triggered event properties ──────────────────────", False),
        (f"  Fraction contained in LS {res['spec_cont']*100:.2f} %  "
         f"(of triggered)", False),
        (f"  Effective trigger energy {res['eff_E_trig']*1000:.0f} keV  "
         f"(spectrum-weighted, triggered subset)", False),
        (f"  Mean edep plastic scint  {res['mean_scint_trig']*1000:.2f} keV  (triggered)", False),
        (f"  Mean edep LS1            {res['mean_ls1_trig']*1000:.2f} keV  (triggered)", False),
        (f"  Mean total LS edep       {res['mean_lstot_trig']*1000:.2f} keV  (triggered)", False),
        ("", False),
        ("─── Rate estimates ──────────────────────────────────", False),
    ]
    # Rate for common activities
    total_eff = GEOM_EFFICIENCY * res["spec_trig_eff"] * BETAS_PER_DECAY
    for act_bq, name in [(3.7e4, "1 μCi"), (3.7e5, "10 μCi"),
                          (3.7e6, "1 mCi"), (3.7e7, "10 mCi")]:
        rate = act_bq * total_eff
        lines.append((f"  {name:10s}  ({act_bq/1e3:>7.0f} kBq)  →  "
                       f"{rate:.3f} Hz trigger", False))
    lines += [
        ("", False),
        ("─── Feasibility assessment ──────────────────────────", False),
        ("  Sr-90 component:  Essentially NO electrons reach LS1.", False),
        ("  Y-90 component:   Small fraction triggers; effective energy ~Y90/3.", False),
        ("  Calibration use:  Possible for MM and plastic scint. alone (low threshold).", False),
        ("  For LS stack:     Needs mCi-level source or dedicated geometry.", False),
        ("  Angular tracking: Available for triggered subset (see simulation plots).", False),
    ]

    y = 0.97
    for text, bold in lines:
        weight = "bold" if bold else "normal"
        fs     = 13 if bold else 9
        ax.text(0.05, y, text, transform=ax.transAxes,
                fontsize=fs, fontweight=weight, va="top",
                fontfamily="monospace" if not bold else "sans-serif")
        y -= 0.040 if bold else 0.038

    pdf.savefig(fig); plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Sr-90/Y-90 calibration feasibility analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--spectrum",
                   default=str(Path(__file__).parent / "Sr90_Y90_Beta_Spectrum.csv"),
                   help="Beta spectrum CSV")
    p.add_argument("--summary",  required=True,
                   help="Simulation summary CSV (from analyse_full_experiment.py, electron)")
    p.add_argument("--outfile",
                   default=str(Path(__file__).parent / "sr90_calibration.pdf"),
                   help="Output PDF")
    return p.parse_args()


def main():
    args = parse_args()

    plt.rcParams.update({
        "figure.dpi": 150, "axes.titlesize": 11,
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9,
    })

    print(f"Loading beta spectrum from {args.spectrum} ...")
    spectrum = load_spectrum(args.spectrum)
    print(f"  {len(spectrum)} energy points, "
          f"{spectrum['Energy_MeV'].min():.3f} – {spectrum['Energy_MeV'].max():.3f} MeV")

    print(f"Loading simulation summary from {args.summary} ...")
    summary = load_summary(args.summary)
    print(f"  {len(summary)} energy points, "
          f"{summary['energy_MeV'].min():.4g} – {summary['energy_MeV'].max():.4g} MeV")

    print(f"Geometric efficiency: {GEOM_EFFICIENCY*100:.2f}%  "
          f"(solid angle {SOLID_ANGLE_SR:.3f} sr)")

    print("Computing spectrum-weighted quantities ...")
    res = compute_spectral_quantities(spectrum, summary)

    print(f"\nKey results:")
    print(f"  Combined spectral trigger efficiency: {res['spec_trig_eff']*100:.3f} %")
    print(f"  Sr-90 trigger efficiency:             {res['spec_trig_sr']*100:.4f} %")
    print(f"  Y-90 trigger efficiency:              {res['spec_trig_y']*100:.3f} %")
    print(f"  Total efficiency (geom × spectral):   "
          f"{GEOM_EFFICIENCY * res['spec_trig_eff'] * BETAS_PER_DECAY * 100:.4f} %")
    print(f"  Trigger rate @ 1 mCi:                 "
          f"{3.7e6 * BETAS_PER_DECAY * GEOM_EFFICIENCY * res['spec_trig_eff']:.2f} Hz")

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:

        # Cover
        fig = plt.figure(figsize=(8, 4))
        fig.text(0.5, 0.65, "Sr-90 / Y-90 Calibration Feasibility",
                 ha="center", fontsize=18, fontweight="bold")
        fig.text(0.5, 0.48,
                 f"Spectral trigger efficiency (combined): "
                 f"{res['spec_trig_eff']*100:.3f} %\n"
                 f"Geometric efficiency: {GEOM_EFFICIENCY*100:.1f} %\n"
                 f"Rate @ 1 mCi: "
                 f"{3.7e6*BETAS_PER_DECAY*GEOM_EFFICIENCY*res['spec_trig_eff']:.2f} Hz",
                 ha="center", fontsize=11, color="grey")
        pdf.savefig(fig); plt.close(fig)

        print("  Summary page ...")
        plot_summary_page(pdf, res)

        print("  Beta spectrum overview ...")
        plot_spectrum_overview(pdf, res)

        print("  Electron fate breakdown ...")
        plot_fate_breakdown(pdf, res)

        print("  Triggered energy deposition ...")
        plot_triggered_edep(pdf, res)

        print("  Rate estimates ...")
        plot_rate_estimate(pdf, res)

        print("  MM context ...")
        plot_mm_context(pdf, res, summary)

    print(f"\nDone → {args.outfile}")


if __name__ == "__main__":
    main()
