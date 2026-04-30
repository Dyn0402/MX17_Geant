#!/usr/bin/env python3
"""
plot_results.py
Generates publication-quality sensitivity comparison plots from summary.csv.

Produces:
  1. sensitivity_vs_energy.pdf  -- <Nprim_drift> vs energy, one panel per particle,
                                    curves per gas, log-log axes.
  2. gas_comparison_bar.pdf     -- bar chart of mean Nprim at selected energies.
  3. efficiency_heatmap.pdf     -- detection efficiency (fraction events with Nprim>0)
                                   as heatmap: gas x energy, one per particle type.
  4. nprim_distributions.pdf    -- histograms of Nprim for selected (gas,particle,E).
     (requires --rootdir to read EventTrees directly)

Usage:
    python3 scripts/plot_results.py --summary /path/to/summary.csv [options]
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")  # headless on lxplus
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("ERROR: numpy/pandas/matplotlib not found. Source setup_lxplus.sh first.")
    sys.exit(1)

# ---- Styling ----
GAS_COLORS = {
    "ArCF4":    "#1f77b4",
    "HeEth":    "#2ca02c",
    "ArCO2":    "#d62728",
    "ArCF4Iso": "#ff7f0e",
    "NeIso":    "#9467bd",
    "NeCF4":    "#8c564b",
    "ArCF4CO2": "#17becf",
    "PureCF4":  "#e377c2",
}
GAS_LABELS = {
    "ArCF4":    "Ar/CF₄ 90/10",
    "HeEth":    "He/C₂H₆ 96.5/3.5",
    "ArCO2":    "Ar/CO₂ 70/30",
    "ArCF4Iso": "Ar/CF₄/iC₄H₁₀ 88/10/2",
    "NeIso":    "Ne/iC₄H₁₀ 95/5",
    "NeCF4":    "Ne/CF₄ 90/10",
    "ArCF4CO2": "Ar/CF₄/CO₂ 45/40/15",
    "PureCF4":  "Pure CF₄",
}
PARTICLE_TITLES = {
    "gamma":    "Photons (γ)",
    "electron": "Electrons (e⁻)",
    "neutron":  "Neutrons (n)",
}

# Shielding comparison styling (by Al thickness in mm)
AL_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 4: "#d62728"}
AL_STYLES = {0: "-",       1: "--",      4: ":"}
AL_LABELS = {0: "No shielding (0 mm Al)", 1: "1 mm Al", 4: "4 mm Al"}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="Path to summary.csv")
    p.add_argument("--outdir",  default=None,
                   help="Output directory for plots (default: summary dir)")
    p.add_argument("--rootdir", default=None,
                   help="Directory with merged ROOT files (for distribution plots)")
    return p.parse_args()


def load_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "al_mm" not in df.columns:
        df["al_mm"] = 0.0  # backward compat: old summaries have no Al column
    return df


def _draw_sensitivity_panel(ax, sub: pd.DataFrame, particle: str,
                            yscale: str, xlim=None):
    """Draw one sensitivity panel onto ax. Helper for fig_sensitivity_vs_energy."""
    mean_vals = []

    for gas in sorted(sub["gas"].unique()):
        g = sub[sub["gas"] == gas].sort_values("energy_MeV")
        if xlim is not None:
            g = g[(g["energy_MeV"] >= xlim[0]) & (g["energy_MeV"] <= xlim[1])]
        if g.empty:
            continue

        clip_low = 0.01 if yscale == "log" else 0.0
        means_raw = g["nPrimDrift_mean"]
        means = means_raw.clip(lower=clip_low)
        mean_vals.extend(means.values)

        # Error band: standard error on the mean (SEM = std / sqrt(n))
        # Zeros are included in both mean and std, so this is the SEM on the
        # true per-event mean (which captures efficiency correctly).
        sem = g["nPrimDrift_std"] / np.sqrt(g["n_events"])
        band_lo = (means_raw - sem).clip(lower=clip_low)
        band_hi = (means_raw + sem).clip(lower=clip_low)
        ax.fill_between(g["energy_MeV"].values, band_lo.values, band_hi.values,
                        alpha=0.25, color=GAS_COLORS.get(gas, "grey"))

        ax.plot(g["energy_MeV"], means, "o-", color=GAS_COLORS.get(gas, "grey"),
                label=GAS_LABELS.get(gas, gas), linewidth=1.8, markersize=4)

    ax.set_xscale("log")
    ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_xlabel("Particle energy [MeV]", fontsize=11)
    ax.set_ylabel("Mean primary ion pairs in drift gap", fontsize=11)
    ax.set_title(PARTICLE_TITLES.get(particle, particle), fontsize=12)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.3)

    # Set y-limits from mean line data only, ignoring the shaded bands
    if mean_vals:
        lo, hi = min(mean_vals), max(mean_vals)
        if yscale == "log":
            ax.set_ylim(bottom=max(0.01, lo * 0.3), top=hi * 3)
        else:
            ax.set_ylim(bottom=0, top=hi * 1.15)
    elif yscale == "log":
        ax.set_ylim(bottom=0.01)


def _make_sensitivity_fig(df: pd.DataFrame, title: str,
                          xlims: dict, yscale: str):
    """Build and return a sensitivity figure. xlims maps particle -> (xmin, xmax) or None."""
    particles = [p for p in ["gamma", "electron", "neutron"] if p in df["particle"].unique()]
    fig, axes = plt.subplots(1, len(particles), figsize=(6 * len(particles), 5), sharey=False)
    if len(particles) == 1:
        axes = [axes]
    for ax, particle in zip(axes, particles):
        sub = df[df["particle"] == particle].copy()
        _draw_sensitivity_panel(ax, sub, particle, yscale, xlim=xlims.get(particle))
    fig.suptitle(title, fontsize=13)
    plt.tight_layout()
    return fig


def fig_sensitivity_vs_energy(df: pd.DataFrame, outdir: Path):
    """
    Write sensitivity_vs_energy.pdf with three pages:
      1. Full range, log-log
      2. Zoomed range, log-log  (γ: 1 keV–20 MeV; n: 50 keV–20 MeV; e⁻: unchanged)
      3. Zoomed range, log-linear
    """
    base_title = "Micromegas Drift Gap Primary Ionization\n(3 cm gas, W-value corrected)"
    zoom_xlims = {
        "gamma":    (1e-3, 20.0e3),   # 1 keV – 20 GeV
        "electron": None,           # unchanged
        "neutron":  (0.05, 20.0e3),   # 50 keV – 20 GeV
    }
    full_xlims = {"gamma": None, "electron": None, "neutron": None}

    pages = [
        (full_xlims, "log",    base_title),
        (zoom_xlims, "log",    base_title + " — zoomed"),
        (zoom_xlims, "linear", base_title + " — zoomed, linear y"),
    ]

    out = outdir / "sensitivity_vs_energy.pdf"
    with PdfPages(out) as pdf:
        for xlims, yscale, title in pages:
            fig = _make_sensitivity_fig(df, title, xlims, yscale)
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {out} (3 pages)")


def fig_efficiency(df: pd.DataFrame, outdir: Path):
    """Heatmap of detection efficiency (Nprim > 0) per particle."""
    particles = [p for p in ["gamma", "electron", "neutron"] if p in df["particle"].unique()]
    gases_sorted = sorted(df["gas"].unique(),
                          key=lambda g: list(GAS_LABELS.keys()).index(g)
                          if g in GAS_LABELS else 99)

    fig, axes = plt.subplots(1, len(particles),
                              figsize=(5 * len(particles), max(3, len(gases_sorted) * 0.6 + 1)))
    if len(particles) == 1:
        axes = [axes]

    for ax, particle in zip(axes, particles):
        sub = df[df["particle"] == particle]
        energies = sorted(sub["energy_MeV"].unique())
        data = np.zeros((len(gases_sorted), len(energies)))
        for i, gas in enumerate(gases_sorted):
            for j, e in enumerate(energies):
                row = sub[(sub["gas"] == gas) & (np.isclose(sub["energy_MeV"], e))]
                if not row.empty:
                    data[i, j] = row["efficiency"].values[0]

        im = ax.imshow(data, aspect="auto", vmin=0, vmax=1,
                       cmap="RdYlGn", origin="upper")
        ax.set_xticks(range(len(energies)))
        ax.set_xticklabels([f"{e:.3g}" for e in energies],
                           rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(len(gases_sorted)))
        ax.set_yticklabels([GAS_LABELS.get(g, g) for g in gases_sorted], fontsize=8)
        ax.set_xlabel("Energy [MeV]", fontsize=10)
        ax.set_title(PARTICLE_TITLES.get(particle, particle), fontsize=11)
        plt.colorbar(im, ax=ax, label="Detection efficiency")

        # Annotate cells
        for i in range(len(gases_sorted)):
            for j in range(len(energies)):
                val = data[i, j]
                color = "white" if val < 0.4 or val > 0.8 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6, color=color)

    fig.suptitle("Detection Efficiency (fraction events with Nprim > 0)", fontsize=12)
    plt.tight_layout()
    out = outdir / "efficiency_heatmap.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig_relative_sensitivity(df: pd.DataFrame, outdir: Path):
    """
    Normalised sensitivity: for each (particle, energy), normalise Nprim across gases.
    Shows which gas is most/least sensitive at each energy point.
    """
    particles = [p for p in ["gamma", "electron", "neutron"] if p in df["particle"].unique()]
    fig, axes = plt.subplots(len(particles), 1,
                              figsize=(10, 4 * len(particles)), sharex=False)
    if len(particles) == 1:
        axes = [axes]

    for ax, particle in zip(axes, particles):
        sub = df[df["particle"] == particle].copy()
        energies = sorted(sub["energy_MeV"].unique())
        gases_here = sorted(sub["gas"].unique())

        # Build pivot
        pivot = sub.pivot_table(index="energy_MeV", columns="gas",
                                values="nPrimDrift_mean")
        # Normalise each row to max
        row_max = pivot.max(axis=1)
        norm = pivot.div(row_max, axis=0)

        for gas in gases_here:
            if gas not in norm.columns:
                continue
            ax.plot(norm.index, norm[gas], "o-",
                    color=GAS_COLORS.get(gas, "grey"),
                    label=GAS_LABELS.get(gas, gas), linewidth=1.6, markersize=4)

        ax.set_xscale("log")
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("Relative Nprim (1 = best gas)", fontsize=10)
        ax.set_title(PARTICLE_TITLES.get(particle, particle), fontsize=11)
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, which="both", alpha=0.3)
        ax.axhline(1.0, color="k", lw=0.8, ls="--")

    axes[-1].set_xlabel("Energy [MeV]", fontsize=11)
    fig.suptitle("Relative Gas Sensitivity (normalised to best gas per energy point)",
                 fontsize=12)
    plt.tight_layout()
    out = outdir / "relative_sensitivity.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig_nprim_distributions_from_root(rootdir: Path, outdir: Path, df_summary: pd.DataFrame):
    """
    Read EventTree directly to plot Nprim histograms at selected points.
    Only runs if --rootdir is given and uproot is available.
    """
    try:
        import uproot
    except ImportError:
        print("  Skipping distribution plots (uproot not available)")
        return

    # Select interesting points: one energy per particle per gas
    select = [
        ("gamma",    1.0),
        ("gamma",    0.06),
        ("electron", 5.0),
        ("neutron",  1.0),
        ("neutron",  2.5e-8),
    ]

    gases_sorted = [g for g in GAS_LABELS if g in df_summary["gas"].unique()]

    pdf_path = outdir / "nprim_distributions.pdf"
    with PdfPages(pdf_path) as pdf:
        for (particle, energy) in select:
            fig, axes = plt.subplots(1, len(gases_sorted),
                                     figsize=(4 * len(gases_sorted), 4),
                                     sharey=False)
            if len(gases_sorted) == 1:
                axes = [axes]

            for ax, gas in zip(axes, gases_sorted):
                e_str = f"{energy:.6g}".replace(".", "p")
                fnames = list(rootdir.glob(f"{gas}_{particle}_{e_str}MeV_merged.root"))
                if not fnames:
                    fnames = list(rootdir.glob(f"{gas}_{particle}_{e_str}MeV_t*.root"))
                if not fnames:
                    ax.text(0.5, 0.5, "No file", transform=ax.transAxes, ha="center")
                    ax.set_title(GAS_LABELS.get(gas, gas), fontsize=9)
                    continue

                try:
                    with uproot.open(fnames[0]) as f:
                        arr = f["EventTree"]["nPrimDrift"].array(library="np")
                except Exception as e:
                    ax.text(0.5, 0.5, f"Error:\n{e}", transform=ax.transAxes,
                            ha="center", fontsize=7)
                    continue

                arr = arr[arr > 0]  # non-zero events only
                if len(arr) == 0:
                    ax.text(0.5, 0.5, "No ionization", transform=ax.transAxes, ha="center")
                else:
                    ax.hist(arr, bins=50, color=GAS_COLORS.get(gas, "grey"),
                            edgecolor="none", alpha=0.8, density=True)
                    ax.axvline(np.mean(arr), color="k", lw=1.5, ls="--",
                               label=f"μ={np.mean(arr):.1f}")
                    ax.legend(fontsize=8)

                ax.set_title(GAS_LABELS.get(gas, gas), fontsize=9)
                ax.set_xlabel("N primary ion pairs (drift)", fontsize=9)
                ax.set_ylabel("Probability density", fontsize=9)

            title = f"{PARTICLE_TITLES.get(particle, particle)} @ {energy:.4g} MeV"
            fig.suptitle(title, fontsize=11)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {pdf_path}")


def fig_shielding_comparison(df: pd.DataFrame, outdir: Path):
    """
    For each gas that has runs with Al shielding, produce a comparison plot:
    one panel per particle type, one curve per Al thickness (0, 1, 4 mm).
    Saved to shielding_comparison.pdf.
    """
    shielded_gases = df[df["al_mm"] > 0]["gas"].unique()
    if len(shielded_gases) == 0:
        return

    out = outdir / "shielding_comparison.pdf"
    with PdfPages(out) as pdf:
        for gas in sorted(shielded_gases):
            sub = df[df["gas"] == gas].copy()
            al_thicknesses = sorted(sub["al_mm"].unique())

            particles = [p for p in ["gamma", "electron", "neutron"]
                         if p in sub["particle"].unique()]
            fig, axes = plt.subplots(1, len(particles),
                                     figsize=(6 * len(particles), 5), sharey=False)
            if len(particles) == 1:
                axes = [axes]

            for ax, particle in zip(axes, particles):
                psub = sub[sub["particle"] == particle].copy()
                clip_low = 0.01

                for al_mm in al_thicknesses:
                    al_key = int(round(al_mm))
                    g = psub[np.isclose(psub["al_mm"], al_mm)].sort_values("energy_MeV")
                    if g.empty:
                        continue

                    means = g["nPrimDrift_mean"].clip(lower=clip_low)
                    sem   = g["nPrimDrift_std"] / np.sqrt(g["n_events"])
                    band_lo = (g["nPrimDrift_mean"] - sem).clip(lower=clip_low)
                    band_hi = (g["nPrimDrift_mean"] + sem).clip(lower=clip_low)

                    color = AL_COLORS.get(al_key, "grey")
                    style = AL_STYLES.get(al_key, "-")
                    label = AL_LABELS.get(al_key, f"{al_mm:g} mm Al")

                    ax.fill_between(g["energy_MeV"].values, band_lo.values, band_hi.values,
                                    alpha=0.2, color=color)
                    ax.plot(g["energy_MeV"], means, style, color=color,
                            label=label, linewidth=1.8, markersize=4)

                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.set_xlabel("Particle energy [MeV]", fontsize=11)
                ax.set_ylabel("Mean primary ion pairs in drift gap", fontsize=11)
                ax.set_title(PARTICLE_TITLES.get(particle, particle), fontsize=12)
                ax.legend(fontsize=9, loc="best")
                ax.grid(True, which="both", alpha=0.3)

                all_means = psub["nPrimDrift_mean"].clip(lower=clip_low)
                if not all_means.empty:
                    ax.set_ylim(bottom=max(0.01, all_means.min() * 0.3),
                                top=all_means.max() * 3)

            fig.suptitle(f"Al Shielding Effect — {GAS_LABELS.get(gas, gas)}\n"
                         f"(3 cm drift gap, Mylar window + 2 cm air gap + Al)",
                         fontsize=12)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {out} ({len(shielded_gases)} gas(es))")


def main():
    args = parse_args()
    summary_path = Path(args.summary)
    outdir = Path(args.outdir) if args.outdir else summary_path.parent

    print(f"Reading summary: {summary_path}")
    df = load_summary(str(summary_path))
    print(f"  {len(df)} rows, gases: {sorted(df['gas'].unique())}, "
          f"particles: {sorted(df['particle'].unique())}")

    # No-shielding subset for the standard sensitivity/efficiency plots
    df_noshield = df[df["al_mm"] == 0].copy()

    print("\nGenerating plots...")
    fig_sensitivity_vs_energy(df_noshield, outdir)
    fig_efficiency(df_noshield, outdir)
    fig_relative_sensitivity(df_noshield, outdir)

    fig_shielding_comparison(df, outdir)

    if args.rootdir:
        rootdir = Path(args.rootdir)
        print("Generating Nprim distribution plots from ROOT files...")
        fig_nprim_distributions_from_root(rootdir, outdir, df_noshield)

    print(f"\nAll plots saved to: {outdir}")


if __name__ == "__main__":
    main()
