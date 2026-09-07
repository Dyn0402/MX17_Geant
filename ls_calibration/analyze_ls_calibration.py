#!/usr/bin/env python3
"""
analyze_ls_calibration.py
Analysis for kLSCalib and kBackScintCalib simulation modes.

Each mode runs with the Sr-90/Y-90 beta spectrum sampled directly by the gun.
This script reads the resulting event CSV(s), plots the measured edep spectra,
and assesses calibration quality.

Usage:
    # LS calibration
    python analyze_ls_calibration.py \\
        --indir /path/to/lscalib/csvs --mode ls --outpdf lscalib_ls.pdf

    # Back scint calibration
    python analyze_ls_calibration.py \\
        --indir /path/to/backscint/csvs --mode backscint --outpdf lscalib_back.pdf

CSV column layout (from kLSCalib / kBackScintCalib):
    eventID, edepDrift_eV, ..., primaryKE_MeV,
    edepLS1_eV, edepLSCFRP_eV, primInLS1,
    edepBackScint_eV, primInBackScint
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpdf

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Constants ─────────────────────────────────────────────────────────────────
# Sr-90 → Sr-90 endpoint: 0.546 MeV
# Y-90  → Y-90  endpoint: 2.280 MeV
SR90_ENDPOINT_MEV = 0.546
Y90_ENDPOINT_MEV  = 2.280


# ── Data loading ──────────────────────────────────────────────────────────────

def load_csvs(indir, pattern="*_events.csv"):
    """Load and concatenate all matching CSVs from indir."""
    files = sorted(glob.glob(os.path.join(indir, pattern)))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {indir}")
    frames = []
    for f in files:
        try:
            df = pd.read_csv(f)
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not load {f}: {e}", file=sys.stderr)
    if not frames:
        raise RuntimeError("No CSV files loaded successfully.")
    df = pd.concat(frames, ignore_index=True)
    # Normalise column names
    df.columns = [c.strip() for c in df.columns]
    print(f"  Loaded {len(frames)} file(s), {len(df):,} events total.")
    return df


def get_col(df, *candidates, default=0.0):
    """Return first matching column as numpy array, or array of default."""
    for c in candidates:
        # exact match
        if c in df.columns:
            return df[c].values.astype(float)
        # case-insensitive
        lc = {col.lower(): col for col in df.columns}
        if c.lower() in lc:
            return df[lc[c.lower()]].values.astype(float)
    return np.full(len(df), default)


# ── Plot helpers ──────────────────────────────────────────────────────────────

def _vline(ax, x, label, color="grey", lw=1.2):
    ax.axvline(x, color=color, lw=lw, ls="--", alpha=0.7, label=label)


def plot_edep_spectrum(pdf, df, mode):
    """Main edep histogram for the active detector."""
    if mode == "ls":
        edep_eV = get_col(df, "edepLS1_eV", "edepls1_ev", "edepls1")
        det_label = "Liquid scintillator (LS1, 2 cm LAB)"
        color = "#1f77b4"
    else:
        edep_eV = get_col(df, "edepBackScint_eV", "edepbackscint_ev", "edepbackscint")
        det_label = "Back plastic scintillator (PVT, 2 cm)"
        color = "#d62728"

    edep_keV = edep_eV * 1e-3

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"Energy deposit spectrum — {det_label}", fontsize=12)

    for ax, xmax, title in [
        (axes[0], max(edep_keV.max()*1.05, 100), "Full range"),
        (axes[1], 500,                            "Zoom: 0–500 keV"),
    ]:
        bins = np.linspace(0, xmax, 200)
        ax.hist(edep_keV, bins=bins, color=color, alpha=0.75, density=True)
        ax.set_xlabel("Energy deposit [keV]", fontsize=10)
        ax.set_ylabel("Normalised counts / keV", fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, xmax)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_edep_vs_primary_ke(pdf, df, mode):
    """2D: primary KE (sampled beta energy) vs edep in detector."""
    ke   = get_col(df, "primaryKE_MeV", "primaryke_mev")
    if mode == "ls":
        edep = get_col(df, "edepLS1_eV", "edepls1_ev") * 1e-3
        ylabel = "LS edep [keV]"
        color  = "#1f77b4"
    else:
        edep = get_col(df, "edepBackScint_eV", "edepbackscint_ev") * 1e-3
        ylabel = "Back scint edep [keV]"
        color  = "#d62728"

    # Only events where edep > 0
    mask = edep > 0
    ke_m, edep_m = ke[mask], edep[mask]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hexbin(ke_m, edep_m, gridsize=80, mincnt=1, cmap="Blues", bins="log")
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label("counts (log scale)")
    _vline(ax, SR90_ENDPOINT_MEV, f"Sr-90 Q = {SR90_ENDPOINT_MEV} MeV")
    _vline(ax, Y90_ENDPOINT_MEV,  f"Y-90 Q = {Y90_ENDPOINT_MEV} MeV", color="orange")
    ax.set_xlabel("Primary beta KE [MeV]", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title("Sampled beta energy vs. detector edep", fontsize=11)
    ax.set_xlim(0, Y90_ENDPOINT_MEV * 1.05)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_mean_edep_vs_ke(pdf, df, mode):
    """Mean edep vs primary KE (binned)."""
    ke   = get_col(df, "primaryKE_MeV", "primaryke_mev")
    if mode == "ls":
        edep = get_col(df, "edepLS1_eV", "edepls1_ev") * 1e-3
        ylabel = "⟨LS edep⟩ [keV]"
        color  = "#1f77b4"
        wall   = get_col(df, "edepLSCFRP_eV", "edeplscfrp_ev") * 1e-3
        wall_label = "⟨LS wall edep⟩ [keV]"
    else:
        edep = get_col(df, "edepBackScint_eV", "edepbackscint_ev") * 1e-3
        ylabel = "⟨Back scint edep⟩ [keV]"
        color  = "#d62728"
        wall   = None
        wall_label = None

    Ebins = np.linspace(0, Y90_ENDPOINT_MEV, 50)
    idx   = np.digitize(ke, Ebins) - 1

    Ec, mean_e, std_e, eff = [], [], [], []
    mean_w = []
    for i in range(len(Ebins)-1):
        sel = edep[idx == i]
        if len(sel) < 20:
            continue
        Ec.append(0.5*(Ebins[i]+Ebins[i+1]))
        mean_e.append(np.mean(sel))
        std_e.append(np.std(sel))
        eff.append(np.mean(sel > 0))
        if wall is not None:
            mean_w.append(np.mean(wall[idx == i]))

    Ec = np.array(Ec)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Mean edep and detection efficiency vs. primary beta energy", fontsize=11)

    axes[0].errorbar(Ec, mean_e, yerr=std_e, fmt="o-", color=color, ms=3, lw=1.3,
                     elinewidth=0.8, label=ylabel.replace("⟩","").replace("⟨",""))
    if wall is not None and mean_w:
        axes[0].plot(Ec, mean_w, "s--", color="#888888", ms=3, lw=1, label=wall_label)
    _vline(axes[0], SR90_ENDPOINT_MEV, f"Sr-90 Q")
    _vline(axes[0], Y90_ENDPOINT_MEV,  f"Y-90 Q", color="orange")
    axes[0].set_xlabel("Primary beta KE [MeV]")
    axes[0].set_ylabel("Energy deposit [keV]")
    axes[0].set_xlim(0, Y90_ENDPOINT_MEV*1.05)
    axes[0].set_ylim(bottom=0)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(Ec, np.array(eff)*100, "o-", color=color, ms=3, lw=1.3)
    _vline(axes[1], SR90_ENDPOINT_MEV, f"Sr-90 Q")
    _vline(axes[1], Y90_ENDPOINT_MEV,  f"Y-90 Q", color="orange")
    axes[1].set_xlabel("Primary beta KE [MeV]")
    axes[1].set_ylabel("Detection efficiency [%]")
    axes[1].set_xlim(0, Y90_ENDPOINT_MEV*1.05)
    axes[1].set_ylim(0, 105)
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_summary_page(pdf, df, mode, indir):
    """Text summary."""
    if mode == "ls":
        edep = get_col(df, "edepLS1_eV", "edepls1_ev") * 1e-3
        det = "Liquid scintillator (LS1)"
    else:
        edep = get_col(df, "edepBackScint_eV", "edepbackscint_ev") * 1e-3
        det = "Back plastic scintillator (PVT)"

    ke  = get_col(df, "primaryKE_MeV", "primaryke_mev")
    eff = float(np.mean(edep > 0)) * 100
    mean_all = float(np.mean(edep))
    mean_trig = float(np.mean(edep[edep > 0])) if eff > 0 else 0.0

    # Y-90-only betas (above Sr-90 Q)
    y90_mask = ke > SR90_ENDPOINT_MEV
    if y90_mask.sum() > 0:
        mean_y90 = float(np.mean(edep[y90_mask]))
    else:
        mean_y90 = 0.0

    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    fig.suptitle("Calibration Simulation — Summary", fontsize=14, y=0.97)

    if mode == "ls":
        geom = ("Source gun → 100 mm air → LS module\n"
                "  LS: 2 mm CFRP | 600 µm inner CFRP | 40 µm Al | 2 cm LAB | (mirror)")
    else:
        geom = ("Source gun → 100 mm air → back scint bar\n"
                "  Bar: 30×20×2 cm PVT, wrapped in 20 µm Al foil + 200 µm mylar tape")

    lines = (
        f"Mode              : {mode} calibration\n"
        f"Detector          : {det}\n"
        f"Input             : {indir}\n"
        f"Events            : {len(df):,}\n"
        f"Spectrum          : Sr-90/Y-90 beta spectrum sampled at gun\n"
        f"\n"
        f"─── Geometry ──────────────────────────────────────────\n"
        f"  {geom}\n"
        f"\n"
        f"─── Results (spectrum-integrated) ─────────────────────\n"
        f"  Detection efficiency (edep > 0) : {eff:.1f} %\n"
        f"  ⟨edep⟩  all events              : {mean_all:.1f} keV\n"
        f"  ⟨edep⟩  triggered events        : {mean_trig:.1f} keV\n"
        f"  ⟨edep⟩  Y-90 betas only         : {mean_y90:.1f} keV  (E > {SR90_ENDPOINT_MEV} MeV)\n"
        f"\n"
        f"─── Endpoints ─────────────────────────────────────────\n"
        f"  Sr-90 endpoint : {SR90_ENDPOINT_MEV} MeV\n"
        f"  Y-90  endpoint : {Y90_ENDPOINT_MEV} MeV\n"
    )

    ax.text(0.03, 0.93, lines, transform=ax.transAxes,
            fontsize=10, va="top", ha="left", fontfamily="monospace",
            bbox=dict(facecolor="#f5f5f5", edgecolor="#cccccc", boxstyle="round,pad=0.5"))

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir",   required=True,
                    help="Directory containing *_events.csv files")
    ap.add_argument("--mode",    choices=["ls", "backscint"], default="ls",
                    help="Which detector: 'ls' (liquid scint) or 'backscint' (plastic)")
    ap.add_argument("--pattern", default="*_events.csv",
                    help="Glob pattern for input files")
    ap.add_argument("--outpdf",  default=None,
                    help="Output PDF path (default: <indir>/<mode>_calib_analysis.pdf)")
    args = ap.parse_args()

    outpdf = args.outpdf or os.path.join(args.indir, f"{args.mode}_calib_analysis.pdf")
    det_name = "Liquid scintillator" if args.mode == "ls" else "Back plastic scintillator"

    print(f"Loading {args.mode.upper()} calibration data from {args.indir} ...")
    df = load_csvs(args.indir, pattern=args.pattern)

    print(f"Writing PDF → {outpdf}")
    with mpdf.PdfPages(outpdf) as pdf:
        plot_summary_page(pdf, df, args.mode, args.indir)
        plot_edep_spectrum(pdf, df, args.mode)
        plot_mean_edep_vs_ke(pdf, df, args.mode)
        plot_edep_vs_primary_ke(pdf, df, args.mode)

    print("Done.")


if __name__ == "__main__":
    main()
