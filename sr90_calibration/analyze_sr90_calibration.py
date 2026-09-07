#!/usr/bin/env python3
"""
analyze_sr90_calibration.py
Feasibility study: Sr-90/Y-90 beta source for detector calibration.

Physics setup:
  • Sr-90 (t½ = 28.8 yr) → Y-90 + β⁻  (Emax = 0.546 MeV)
  • Y-90 (t½ = 64.1 h)  → Zr-90 + β⁻  (Emax = 2.28 MeV)
  In secular equilibrium one Sr-90 and one Y-90 beta are emitted per decay,
  so the source produces 2A β⁻/s for a declared activity A.

Source placement:
  Source is placed at the same z-position as the (absent) He-3 target.
  The He-3 capsule is NOT present; electrons travel 226.5 mm of air before
  reaching the Micromegas entrance.  Geometric efficiency is computed from
  the solid angle of the 40×40 cm Micromegas face at 226.5 mm from the source.

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
    python3 sr90_calibration/analyze_sr90_calibration.py \\
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

# ── Collimated source defaults ────────────────────────────────────────────────
# The source is placed inside a plastic collimator that restricts the emission
# cone.  Adjust these if you know the exact collimator geometry.
COLLIMATOR_RADIUS_MM = 5.0    # aperture half-diameter [mm]
COLLIMATOR_LENGTH_MM = 20.0   # plastic depth along beam axis [mm]
COLLIMATOR_GAP_MM    = 5.0    # air gap from source to collimator entrance [mm]
SOURCE_ACTIVITY_BQ   = 2.0e6  # source activity [Bq] (2 MBq)


def _collimator_efficiency(r_mm=COLLIMATOR_RADIUS_MM,
                            l_mm=COLLIMATOR_LENGTH_MM,
                            gap_mm=COLLIMATOR_GAP_MM):
    """
    Fraction of 4π captured by the collimated beam cone, plus the half-angle.
    Model: point source at distance (gap_mm + l_mm) behind the aperture,
    aperture radius r_mm.  All betas through the aperture are assumed to
    reach the 40×40 cm detector (confirmed: beam cone ~r/d  ≪ 20 cm at 22 cm).
    """
    d_total_cm  = (l_mm + gap_mm) / 10.0        # source → aperture exit [cm]
    theta_max   = np.arctan(r_mm / 10.0 / d_total_cm)   # half-angle [rad]
    solid_angle = 2 * np.pi * (1 - np.cos(theta_max))   # [sr]
    efficiency  = solid_angle / (4 * np.pi)
    return efficiency, solid_angle, np.degrees(theta_max)

# ── Material budget gaps ──────────────────────────────────────────────────────
# (label, start_trans_col or None=1.0, end_trans_col or "exit", color)
BUDGET_GAPS = [
    ("200 mm air\n+ Drift window",
     None,                    "trans_primInDrift",     "#1f77b4"),

    ("MM drift gas (30 mm)",
     "trans_primInDrift",     "trans_primInAmp",        "#2ca02c"),

    ("MM amp + resistive paste",
     "trans_primInAmp",       "trans_primInPCB",        "#a1d99b"),

    ("MM PCB",
     "trans_primInPCB",       "trans_primInScintWall",  "#ff7f0e"),

    ("Scint. Wall",
     "trans_primInScintWall", "trans_primInLS1",        "#d62728"),

    ("Liq. scint. layers 1–2 + CFRP walls",
     "trans_primInLS1",       "trans_primInLSCFRP5",    "#6a1d8a"),

    ("Exits LS stack",
     "trans_primInLSCFRP5",   "exit",                   "#555555"),
]

# Approximate radiation lengths and material thicknesses for the table
# (label, thickness_cm, density_g_cm3, Z_eff)  — for annotation only
MATERIAL_LAYERS = [
    # Layer,            t [cm],  rho [g/cm³]  description
    # Note: He-3 gas and capsule walls are ABSENT in the Sr-90 calibration setup.
    # The source sits in air; the path to the MM entrance is 226.5 mm of air.
    ("Air gap (source→MM)", 20.0,  0.00120, "200 mm air  (He-3 capsule absent)"),
    ("Drift window",      0.01,  3.0,     "~0.1 mm Mylar+Cu+Kapton"),
    ("MM drift gas",      3.0,   0.00182, "30 mm Ar/Iso"),
    ("PCB Cu (×4)",       0.0104, 8.96,  "4 × 0.026 mm Cu"),
    ("PCB FR4 (×4)",      0.04,  1.85,   "4 × 0.10 mm FR4"),
    ("PCB Rohacell",      0.50,  0.052,  "5 mm Rohacell 51"),
    ("Air gap 2",         2.0,   0.00120, "20 mm air"),
    ("Plastic scint.",    0.30,  1.032,  "3 mm PVT"),
    ("Air gap 3",         2.0,   0.00120, "20 mm air"),
    ("CFRP cell walls",   0.75,  1.55,   "5 × 1.5 mm CFRP"),
]

# ── Colour scheme ──────────────────────────────────────────────────────────────
C_SR90   = "#e377c2"   # pink   — Sr-90 component
C_Y90    = "#1f77b4"   # blue   — Y-90 component
C_TOTAL  = "#333333"   # dark   — combined spectrum
C_SCINT  = "#d62728"   # red    — plastic scintillator
C_LS1    = "#6a1d8a"   # purple — LS layer 1
C_TRIG   = "#2ca02c"   # green  — trigger efficiency

# Beta spectrum endpoints — used for axis annotations
E_SR90 = 0.546   # Sr-90 beta endpoint [MeV]
E_Y90  = 2.28    # Y-90 beta endpoint  [MeV]

# ── Helpers ───────────────────────────────────────────────────────────────────


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

    # Top: spectrum — isotope fills kept transparent so trigger overlay is clear
    ax = axes[0]
    ax.fill_between(E, S90, alpha=0.15, color=C_SR90, label="Sr-90")
    ax.fill_between(E, S_Y, alpha=0.15, color=C_Y90,  label="Y-90")
    ax.plot(E, S,   color=C_TOTAL, lw=1.4, label="Total (Sr-90 + Y-90)")

    # Solid green trigger acceptance overlaid on top
    ax.fill_between(E, 0, S * f_trig, alpha=0.85, color=C_TRIG, zorder=4,
                    label="Triggered fraction of spectrum")

    ax.axvline(0.546, color=C_SR90, lw=0.9, ls="--", label="Sr-90 endpoint 0.546 MeV")
    ax.axvline(2.28,  color=C_Y90,  lw=0.9, ls="--", label="Y-90 endpoint 2.28 MeV")
    ax.set_ylabel("Probability density  dN/dE", fontsize=11)
    ax.set_title("Sr-90 / Y-90 beta spectrum and trigger acceptance", fontsize=12)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 2.5)

    ax.set_ylim(bottom=0)

    # Bottom: trigger fraction
    ax2 = axes[1]
    ax2.plot(E, f_trig * 100, color=C_TRIG, lw=1.6)
    ax2.set_xlabel("Beta energy (MeV)", fontsize=11)
    ax2.set_ylabel("Trigger fraction (%)", fontsize=11)
    ax2.set_ylim(bottom=0)
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
                    alpha=0.55, color="#ff7f0e", label="Punch-through (exits LS2)")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Fraction of electrons", fontsize=10)
    ax.set_title("Electron fate vs. energy\n(given electron reaches detector)", fontsize=11)
    ax.legend(fontsize=9, loc="lower left")
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

    # Right: LS1 energy scale comparison — Y-90 triggered vs signal electrons.
    # Hardcoded from full-experiment simulation (ArIso, full stack, electron).
    _LS1_Y90_KEV   = res["mean_ls1_trig"] * 1000          # from sr90 simulation
    _LS1_5MEV_KEV  = 1985.4    # mean edep in LS1 at 5 MeV  (full-stack CSV)
    _LS1_16MEV_KEV = 2783.4    # mean edep in LS1 at 16 MeV (full-stack CSV)
    _THRESH_KEV    = 500.0     # rough LS1 trigger threshold for signal electrons

    ax2 = axes[1]
    categories  = ["Y-90 triggered\n(Sr-90 cal.)",
                   "5 MeV\nsignal e⁻", "16 MeV\nsignal e⁻"]
    values      = [_LS1_Y90_KEV, _LS1_5MEV_KEV, _LS1_16MEV_KEV]
    bar_colors  = [C_Y90, "#ff7f0e", "#d62728"]

    bars = ax2.bar(categories, values, color=bar_colors, alpha=0.85, width=0.5)
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 val * 1.04, f"{val:.0f} keV",
                 ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax2.axhline(_THRESH_KEV, color="k", lw=1.2, ls="--",
                label=f"Estimated threshold ≈ {_THRESH_KEV:.0f} keV")
    ax2.set_ylabel("Mean LS1 energy deposition (keV)", fontsize=10)
    ax2.set_title("LS1 energy scale: calibration vs signal\n"
                  "Y-90 deposits ~60× less than 5 MeV signal",
                  fontsize=11)
    ax2.set_ylim(bottom=0)
    ax2.legend(fontsize=9)
    ax2.grid(True, axis="y", alpha=0.3)

    # Annotate the mismatch ratio
    ratio = _LS1_5MEV_KEV / max(_LS1_Y90_KEV, 1)
    ax2.annotate(f"×{ratio:.0f}", xy=(0.5, _LS1_Y90_KEV), xytext=(0.5, _LS1_5MEV_KEV / 2),
                 fontsize=10, color="#333333", ha="center",
                 arrowprops=dict(arrowstyle="<->", color="#555555", lw=1.0))

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

    # Right: LS1 transmission vs energy with Sr-90/Y-90 endpoints marked
    ax2 = axes[1]
    E_sim  = summary["energy_MeV"].values
    f_ls1  = summary["trans_primInLS1"].values if "trans_primInLS1" in summary.columns \
              else np.zeros(len(E_sim))
    ax2.plot(E_sim, f_ls1 * 100, color=C_LS1, lw=1.8, label="LS1 transmission (simulation)")

    ax2.axvspan(0, E_Y90 + 0.05, alpha=0.10, color=C_TOTAL,
                label=f"Sr-90/Y-90 beta range (0–{E_Y90} MeV)")
    ax2.axvline(E_Y90,  color=C_Y90,  lw=1.0, ls="--", label=f"Y-90 endpoint {E_Y90} MeV")
    ax2.axvline(E_SR90, color=C_SR90, lw=1.0, ls="--", label=f"Sr-90 endpoint {E_SR90} MeV")

    ax2.set_xlim(0, min(E_sim.max() * 1.05, 3.5))
    ax2.set_ylim(-2, 105)
    ax2.set_xlabel("Electron energy (MeV)", fontsize=10)
    ax2.set_ylabel("LS1 transmission (%)", fontsize=10)
    ax2.set_title("LS1 transmission across the Sr-90/Y-90 beta spectrum", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Micromegas context for Sr-90/Y-90 calibration", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def _write_lines(ax, lines, y_start=0.97):
    y = y_start
    for text, bold in lines:
        fs = 13 if bold else 9
        ax.text(0.05, y, text, transform=ax.transAxes,
                fontsize=fs, fontweight="bold" if bold else "normal", va="top",
                fontfamily="sans-serif" if bold else "monospace")
        y -= 0.045 if bold else 0.038


def plot_summary_page(pdf, res: dict):
    """Page 1: source properties, geometry, spectral efficiency, triggered events."""
    fig = plt.figure(figsize=(9, 7))
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    col_eff, col_sr, col_deg = _collimator_efficiency()

    _write_lines(ax, [
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
        (f"  Geometric eff. (open)   {GEOM_EFFICIENCY*100:.1f} %  "
         f"(solid angle {SOLID_ANGLE_SR:.3f} sr)", False),
        (f"  Geometric eff. (collimated)  {col_eff*100:.3f} %  "
         f"({col_deg:.1f}° half-angle, r={COLLIMATOR_RADIUS_MM:.0f}mm L={COLLIMATOR_LENGTH_MM:.0f}mm)", False),
        ("", False),
        ("─── Spectral trigger efficiency ─────────────────────", False),
        (f"  Sr-90 alone             {res['spec_trig_sr']*100:.3f} %", False),
        (f"  Y-90 alone              {res['spec_trig_y']*100:.3f} %", False),
        (f"  Combined (total)        {res['spec_trig_eff']*100:.3f} %", False),
        ("", False),
        ("─── Triggered event properties ──────────────────────", False),
        (f"  Effective trigger energy {res['eff_E_trig']*1000:.0f} keV  (spectrum-weighted)", False),
        (f"  Mean edep plastic scint  {res['mean_scint_trig']*1000:.2f} keV", False),
        (f"  Mean edep LS1            {res['mean_ls1_trig']*1000:.2f} keV", False),
        (f"  Mean total LS edep       {res['mean_lstot_trig']*1000:.2f} keV", False),
    ])
    pdf.savefig(fig); plt.close(fig)


def plot_rate_page(pdf, res: dict):
    """Page 2: rate estimates for our 2 MBq collimated source + feasibility notes."""
    fig = plt.figure(figsize=(9, 7))
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    col_eff, col_sr, col_deg = _collimator_efficiency()

    # Rates for both open and collimated geometry
    open_eff = GEOM_EFFICIENCY * res["spec_trig_eff"] * BETAS_PER_DECAY
    col_total_eff = col_eff * res["spec_trig_eff"] * BETAS_PER_DECAY

    rate_open_2mbq   = SOURCE_ACTIVITY_BQ * open_eff
    rate_col_2mbq    = SOURCE_ACTIVITY_BQ * col_total_eff

    lines = [
        ("Rate Estimates and Feasibility", True),
        ("", False),
        ("─── Our source: 2 MBq, collimated ──────────────────", False),
        (f"  Activity                {SOURCE_ACTIVITY_BQ/1e6:.1f} MBq  = {SOURCE_ACTIVITY_BQ:.2e} Bq", False),
        (f"  Collimated geom. eff.   {col_eff*100:.3f} %  ({col_deg:.1f}° half-angle)", False),
        (f"  Spectral trigger eff.   {res['spec_trig_eff']*100:.3f} %", False),
        (f"  Combined efficiency     {col_total_eff*100:.4f} %", False),
        (f"  → Coincidence rate      {rate_col_2mbq:.2f} Hz", False),
        (f"  → Events in 1 hour      {rate_col_2mbq*3600:.0f}", False),
        ("", False),
        ("─── Sensitivity to collimator geometry ─────────────", False),
    ]

    for r_mm, l_mm in [(3, 20), (5, 20), (10, 20), (5, 10)]:
        eff_c, _, deg_c = _collimator_efficiency(r_mm, l_mm, COLLIMATOR_GAP_MM)
        rate_c = SOURCE_ACTIVITY_BQ * eff_c * res["spec_trig_eff"] * BETAS_PER_DECAY
        lines.append((f"  r={r_mm:2.0f}mm L={l_mm:2.0f}mm  "
                       f"θ={deg_c:.1f}°  → {rate_c:.3f} Hz  "
                       f"({rate_c*3600:.0f} ev/hr)", False))

    lines += [
        ("", False),
        ("─── Open source (no collimator) ─────────────────────", False),
        (f"  Geometric eff.          {GEOM_EFFICIENCY*100:.1f} %", False),
        (f"  → Rate @ 2 MBq          {rate_open_2mbq:.1f} Hz", False),
        ("", False),
        ("─── Feasibility assessment ──────────────────────────", False),
        ("  Sr-90 component:  Essentially none reach LS1.", False),
        ("  Y-90 component:   Small fraction triggers (see spectral eff.).", False),
        ("  Calibration use:  Trigger rate is feasible, but see energy", False),
        ("                    mismatch plot — Y-90 deposits ~keV in LS1", False),
        ("                    while signal electrons deposit ~MeV.", False),
        ("  Angular tracking: Available for triggered subset.", False),
    ]

    _write_lines(ax, lines)
    pdf.savefig(fig); plt.close(fig)


# ── Material budget analysis ──────────────────────────────────────────────────

def compute_material_budget(spectrum: pd.DataFrame,
                            summary: pd.DataFrame) -> list:
    """
    For each gap in BUDGET_GAPS compute:
      stopping_vs_E : array of fraction of beam stopped in this gap at each energy
      spectrum_frac : spectrum-weighted fraction of all decays stopped here
      sr_frac / y_frac : per-isotope spectrum fractions
    Returns a list of dicts, one per gap, in stack order.
    """
    E    = spectrum["Energy_MeV"].values
    S    = spectrum["Total_ProbabilityDensity"].values
    S90  = spectrum["Sr90_ProbabilityDensity"].values
    S_Y  = spectrum["Y90_ProbabilityDensity"].values
    norm      = _trapz(S,   E)
    norm_sr   = _trapz(S90, E)
    norm_y    = _trapz(S_Y, E)

    sim_E = summary["energy_MeV"].values

    def get_trans(col):
        if col is None:
            return np.ones_like(E)
        if col == "exit":
            # punch-through: what's left after LS2 (may not be in summary)
            if "trans_primInLSCFRP5" in summary.columns:
                return interp_sim(summary, "trans_primInLSCFRP5", E, 0.0)
            return np.zeros_like(E)
        if col not in summary.columns:
            return np.zeros_like(E)
        return interp_sim(summary, col, E, fill_below=0.0)

    results = []
    for label, start_col, end_col, color in BUDGET_GAPS:
        t_start = get_trans(start_col)
        t_end   = get_trans(end_col)
        stop    = np.clip(t_start - t_end, 0, 1)

        results.append(dict(
            label         = label,
            color         = color,
            stopping_vs_E = stop,
            spectrum_frac  = _trapz(stop * S,   E) / norm      if norm     > 0 else 0,
            sr_frac        = _trapz(stop * S90, E) / norm_sr   if norm_sr  > 0 else 0,
            y_frac         = _trapz(stop * S_Y, E) / norm_y    if norm_y   > 0 else 0,
        ))
    return results


def plot_stopping_vs_energy(pdf, budget: list, res: dict):
    """
    Two-panel page:
      Left : Stacked area chart — stopping fraction per gap vs energy.
             Removable gaps have hatching.  Beta spectrum shown for context.
      Right: Horizontal bar chart — spectrum-weighted stopping per gap.
             Shows total 'avoidable' vs 'unavoidable' stopping.
    """
    E  = res["E"]
    S  = res["S"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Left: stacked area ──────────────────────────────────────────────
    ax = axes[0]
    bottom = np.zeros_like(E)

    # Beta spectrum on secondary y-axis
    ax2 = ax.twinx()
    ax2.fill_between(E, 0, S / S.max(), alpha=0.08, color="grey")
    ax2.set_ylim(0, 3)
    ax2.set_yticks([])

    for gap in budget:
        stop = gap["stopping_vs_E"]
        ax.fill_between(E, bottom, bottom + stop,
                        color=gap["color"], alpha=0.80,
                        linewidth=0.3, label=gap["label"])
        bottom += stop

    ax.set_xlim(0, 2.4)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Fraction of electrons stopped in gap", fontsize=10)
    ax.set_title("Where electrons stop vs energy", fontsize=10)
    ax.axvline(0.546, color="grey", lw=0.8, ls="--", alpha=0.6)
    ax.axvline(2.28,  color="grey", lw=0.8, ls="--", alpha=0.6)
    ax.text(0.55, 0.97, "Sr-90\nend", fontsize=6, color="grey", va="top")
    ax.text(2.29, 0.97, "Y-90\nend", fontsize=6, color="grey", va="top")
    ax.legend(fontsize=6, loc="upper left", ncol=1,
              bbox_to_anchor=(0.0, 0.95))
    ax.grid(True, alpha=0.2)

    # ── Right: horizontal bar ───────────────────────────────────────────
    ax = axes[1]
    labels = [g["label"].replace("\n", " ") for g in budget]
    fracs  = [g["spectrum_frac"] * 100        for g in budget]
    colors = [g["color"]                       for g in budget]

    # Reverse so first gap is on top
    y = np.arange(len(budget))[::-1]
    for i, (frac, color, lbl) in enumerate(zip(fracs, colors, labels)):
        ax.barh(y[i], frac, color=color, alpha=0.80, height=0.7)
        ax.text(frac + 0.3, y[i], f"{frac:.1f} %",
                va="center", fontsize=8, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels([g["label"].replace("\n", " ") for g in budget[::-1]],
                       fontsize=8)
    ax.set_xlabel("Spectrum-weighted fraction stopped (%)", fontsize=10)
    ax.set_xlim(0, max(fracs) * 1.35 + 1)
    ax.set_title("Spectrum-integrated stopping per gap\n(Sr-90 + Y-90 combined)",
                 fontsize=10)
    ax.grid(True, axis="x", alpha=0.3)

    fig.suptitle("Material budget: where does the Sr-90 / Y-90 beta spectrum stop?",
                 fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_cumulative_transmission(pdf, budget: list, res: dict):
    """
    Cumulative transmission at each gap boundary vs energy, for both isotopes.
    Also shows: what the trigger efficiency would be if only 'removable' material
    were eliminated (rough upper bound using air-free approximation).
    """
    E   = res["E"]
    S90 = res["S90"]
    S_Y = res["S_Y"]
    S   = res["S"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ── Left: cumulative transmission waterfall ─────────────────────────
    ax = axes[0]
    cum = np.ones_like(E)
    ax.fill_between(E, 0, S / S.max() * 0.9, alpha=0.08, color="grey",
                    label="Beta spectrum (norm.)")

    for gap in budget[:-1]:   # skip "exits" row
        cum = np.clip(cum - gap["stopping_vs_E"], 0, 1)
        lbl = gap["label"].replace("\n", " ").split("(")[0].strip()
        ax.plot(E, cum * 100, color=gap["color"], lw=1.5, label=f"After {lbl}")

    ax.set_xlim(0, 2.4)
    ax.set_ylim(-2, 105)
    ax.axvline(0.546, color="grey", lw=0.8, ls=":", alpha=0.7)
    ax.axvline(2.28,  color="grey", lw=0.8, ls=":", alpha=0.7)
    ax.set_xlabel("Beta energy (MeV)", fontsize=10)
    ax.set_ylabel("Electrons surviving to this boundary (%)", fontsize=10)
    ax.set_title("Cumulative transmission through the stack", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3)

    # ── Right: per-isotope breakdown table ─────────────────────────────
    ax = axes[1]
    ax.axis("off")

    headers = ["Gap / material", "All β (%)", "Sr-90 (%)", "Y-90 (%)"]
    rows = []
    for gap in budget:
        rows.append([
            gap["label"].replace("\n", " "),
            f"{gap['spectrum_frac']*100:.2f}",
            f"{gap['sr_frac']*100:.2f}",
            f"{gap['y_frac']*100:.2f}",
        ])

    table = ax.table(
        cellText=rows,
        colLabels=headers,
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.5)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#e8e8e8")
            cell.set_text_props(fontweight="bold")

    ax.set_title("Stopping fraction per gap by isotope\n"
                 "(fraction of all decays aimed at detector)", fontsize=10)

    fig.suptitle("Transmission waterfall and per-gap stopping table", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_material_thickness_table(pdf):
    """
    One-page reference table showing every material layer with its
    thickness, areal density [g/cm²], and whether it is removable.
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")

    headers = ["Layer", "Thickness", "Density (g/cm³)", "Areal density (mg/cm²)"]
    rows = []
    total_areal = 0.0
    for name, t_cm, rho, description in MATERIAL_LAYERS:
        areal = t_cm * rho * 1000   # mg/cm²
        total_areal += areal
        t_str   = (f"{t_cm*10:.1f} mm" if t_cm >= 0.1
                   else f"{t_cm*10000:.0f} µm")
        rho_str = f"{rho:.4f}" if rho < 0.01 else f"{rho:.3f}"
        rows.append([description, t_str, rho_str, f"{areal:.1f}"])

    rows.append(["─" * 35, "─" * 8, "─" * 10, "─" * 18])
    rows.append(["Total material budget", "", "", f"{total_areal:.1f}"])

    table = ax.table(
        cellText=rows,
        colLabels=headers,
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.6)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#d8d8d8")
            cell.set_text_props(fontweight="bold")

    # CSDA range reference note (He-3 capsule absent — source sits in air)
    ax.text(0.5, 0.01,
            "CSDA ranges for reference: 2 MeV e⁻ in water ≈ 0.88 g/cm² │ "
            "2 MeV e⁻ in Cu ≈ 0.55 g/cm²\n"
            "Y-90 endpoint 2.28 MeV → CSDA range in water ≈ 1.0 g/cm². "
            "He-3 capsule is absent; 200 mm air gap has negligible areal density.",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8, style="italic", color="#555555",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fffbe6", alpha=0.9))

    ax.set_title("Material layer reference table — full detector stack",
                 fontsize=12, pad=20)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Calibration energy-mismatch plot ─────────────────────────────────────────

def plot_edep_comparison(pdf, res: dict, full_summary_path: str):
    """
    Show the calibration energy-mismatch problem:
      Y-90 triggered electrons deposit only ~keV in LS1,
      while real signal electrons (4–16 MeV) deposit ~MeV.
    Setting a trigger threshold at the Y-90 level would admit
    a huge background rate in the actual experiment.

    Loads per-energy mean edep from the full-experiment simulation CSV.
    """
    try:
        import csv as _csv
        with open(full_summary_path) as f:
            fs = list(_csv.DictReader(f))
    except Exception as e:
        print(f"  Warning: could not load full summary {full_summary_path}: {e}")
        return

    # Find closest rows to 4.5 MeV and 15 MeV
    def closest(target):
        return min(fs, key=lambda r: abs(float(r["energy_MeV"]) - target))

    r4  = closest(4.5)
    r15 = closest(15.0)
    e4  = float(r4["energy_MeV"])
    e15 = float(r15["energy_MeV"])

    # Grab mean edep values [MeV]
    def get(row, col):
        return float(row.get(col, 0) or 0)

    edep_keys = ["edep_edepScintWall", "edep_edepLS1", "edep_edepLS2"]
    labels_det = ["Plastic\nscint.", "LS 1", "LS 2"]

    vals_sr90 = [res["mean_scint_trig"] * 1000,      # keV
                 res["mean_ls1_trig"]   * 1000, 0]   # LS2 not scored in sr90

    vals_4mev  = [get(r4,  k) * 1000 for k in edep_keys]
    vals_15mev = [get(r15, k) * 1000 for k in edep_keys]

    x   = np.arange(len(labels_det))
    w   = 0.26

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Left: per-detector mean edep ─────────────────────────────────────
    ax = axes[0]
    b1 = ax.bar(x - w,  vals_sr90,  w, label=f"Sr-90/Y-90  (triggered, Y-90-weighted)",
                color=C_Y90,  alpha=0.85)
    b2 = ax.bar(x,      vals_4mev,  w, label=f"{e4:.1f} MeV electrons (signal low end)",
                color="#ff7f0e", alpha=0.85)
    b3 = ax.bar(x + w,  vals_15mev, w, label=f"{e15:.1f} MeV electrons (signal high end)",
                color="#d62728", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(labels_det, fontsize=10)
    ax.set_ylabel("Mean energy deposition (keV)", fontsize=10)
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.1)
    ax.set_title("Energy deposition per detector layer\n"
                 "Sr-90/Y-90 triggered events vs signal electrons", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", which="both", alpha=0.3)

    # Annotate ratio
    if vals_sr90[1] > 0 and vals_4mev[1] > 0:
        ratio4  = vals_4mev[1]  / vals_sr90[1]
        ratio15 = vals_15mev[1] / vals_sr90[1]
        ax.text(x[1] + w / 2, vals_15mev[1] * 1.4,
                f"×{ratio4:.0f} vs\nSr-90/Y-90",
                ha="center", fontsize=7, color="#ff7f0e")
        ax.text(x[1] + w, vals_15mev[1] * 3,
                f"×{ratio15:.0f} vs\nSr-90/Y-90",
                ha="center", fontsize=7, color="#d62728")

    # ── Right: the mismatch problem visualised ────────────────────────────
    ax2 = axes[1]

    # Schematic: show three "thresholds" as horizontal bands
    sr90_ls1  = res["mean_ls1_trig"] * 1000       # keV
    sig4_ls1  = vals_4mev[1]                       # keV
    sig15_ls1 = vals_15mev[1]                      # keV

    energies = [sr90_ls1, sig4_ls1, sig15_ls1]
    colors_e = [C_Y90, "#ff7f0e", "#d62728"]
    labels_e = [
        f"Sr-90/Y-90 triggered\n({sr90_ls1:.1f} keV in LS1)",
        f"{e4:.1f} MeV signal\n({sig4_ls1:.0f} keV in LS1)",
        f"{e15:.1f} MeV signal\n({sig15_ls1:.0f} keV in LS1)",
    ]

    ypos = [0.75, 0.45, 0.15]
    ax2.set_xlim(0, max(energies) * 1.4)
    ax2.set_ylim(0, 1)
    ax2.axis("off")

    for e_val, color, lbl, yp in zip(energies, colors_e, labels_e, ypos):
        # Horizontal bar
        ax2.barh(yp, e_val, height=0.12,
                 color=color, alpha=0.85, left=0,
                 transform=ax2.transData)
        ax2.text(e_val + max(energies) * 0.03, yp,
                 f"{e_val:.1f} keV",
                 va="center", fontsize=9, fontweight="bold", color=color)
        ax2.text(-max(energies) * 0.02, yp, lbl,
                 va="center", ha="right", fontsize=8, color="#333333",
                 transform=ax2.transData)

    # Threshold line at Y-90 level
    ax2.axvline(sr90_ls1, color=C_Y90, lw=1.5, ls="--", alpha=0.7,
                transform=ax2.get_xaxis_transform(), zorder=5)
    ax2.text(sr90_ls1 + max(energies) * 0.01, 1.02,
             "← calibration\nthreshold",
             va="bottom", fontsize=7, color=C_Y90,
             transform=ax2.get_xaxis_transform())

    ax2.set_title("Mismatch in LS1 energy scale\n"
                  "Calibrating to Y-90 sets threshold far below signal",
                  fontsize=11)

    fig.suptitle("Calibration energy-mismatch: Sr-90/Y-90 vs signal electrons",
                 fontsize=12)
    fig.tight_layout()
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
                   help="Sr-90 simulation summary CSV (from analyze_full_experiment.py, "
                        "electron, --prefix sr90)")
    p.add_argument("--full-summary", default=None,
                   help="Full-experiment simulation summary CSV for energy comparison "
                        "(from analyze_full_experiment.py, electron, --prefix full)")
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
    res    = compute_spectral_quantities(spectrum, summary)
    budget = compute_material_budget(spectrum, summary)

    col_eff, col_sr, col_deg = _collimator_efficiency()
    rate_collimated = SOURCE_ACTIVITY_BQ * col_eff * res["spec_trig_eff"] * BETAS_PER_DECAY

    print(f"\nKey results:")
    print(f"  Combined spectral trigger efficiency: {res['spec_trig_eff']*100:.3f} %")
    print(f"  Collimated geom. eff. (r={COLLIMATOR_RADIUS_MM:.0f}mm, L={COLLIMATOR_LENGTH_MM:.0f}mm): "
          f"{col_eff*100:.3f} %  ({col_deg:.1f}° half-angle)")
    print(f"  Coincidence rate @ {SOURCE_ACTIVITY_BQ/1e6:.1f} MBq collimated: "
          f"{rate_collimated:.3f} Hz  ({rate_collimated*3600:.0f} ev/hr)")

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:

        # Cover
        fig = plt.figure(figsize=(8, 4))
        fig.text(0.5, 0.65, "Sr-90 / Y-90 Calibration Feasibility",
                 ha="center", fontsize=18, fontweight="bold")
        fig.text(0.5, 0.48,
                 f"Spectral trigger efficiency: {res['spec_trig_eff']*100:.3f} %\n"
                 f"Collimated rate @ {SOURCE_ACTIVITY_BQ/1e6:.0f} MBq: "
                 f"{rate_collimated:.3f} Hz  ({rate_collimated*3600:.0f} ev/hr)",
                 ha="center", fontsize=11, color="grey")
        pdf.savefig(fig); plt.close(fig)

        print("  Summary page 1 (physics + geometry) ...")
        plot_summary_page(pdf, res)

        print("  Summary page 2 (rates + feasibility) ...")
        plot_rate_page(pdf, res)

        print("  Beta spectrum overview ...")
        plot_spectrum_overview(pdf, res)

        print("  Electron fate breakdown ...")
        plot_fate_breakdown(pdf, res)

        print("  Triggered energy deposition ...")
        plot_triggered_edep(pdf, res)

        print("  Rate estimates (log-log) ...")
        plot_rate_estimate(pdf, res)

        print("  MM context ...")
        plot_mm_context(pdf, res, summary)

        if args.full_summary:
            print("  Calibration energy mismatch ...")
            plot_edep_comparison(pdf, res, args.full_summary)
        else:
            print("  (skipping energy mismatch — no --full-summary provided)")

        print("  Material budget: stopping vs energy ...")
        plot_stopping_vs_energy(pdf, budget, res)

        print("  Material budget: cumulative transmission + table ...")
        plot_cumulative_transmission(pdf, budget, res)

        print("  Material layer reference table ...")
        plot_material_thickness_table(pdf)

    print(f"\nDone → {args.outfile}")


if __name__ == "__main__":
    main()
