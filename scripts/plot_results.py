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
    # Mixtures
    "ArCF4":    "#1f77b4",
    "HeEth":    "#2ca02c",
    "ArCO2":    "#d62728",
    "ArCF4Iso": "#ff7f0e",
    "NeIso":    "#9467bd",
    "NeCF4":    "#8c564b",
    "ArCF4CO2": "#17becf",
    "PureCF4":  "#e377c2",
    # Pure gases
    "PureAr":     "#5ba3d9",
    "PureHe":     "#55b96a",
    "PureNe":     "#e8943a",
    "PureEthane": "#b07dd4",
    "PureIso":    "#a0633a",
    "PureCO2":    "#cc4444",
}
GAS_LABELS = {
    # Mixtures
    "ArCF4":    "Ar/CF₄ 90/10",
    "HeEth":    "He/C₂H₆ 96.5/3.5",
    "ArCO2":    "Ar/CO₂ 70/30",
    "ArCF4Iso": "Ar/CF₄/iC₄H₁₀ 88/10/2",
    "NeIso":    "Ne/iC₄H₁₀ 95/5",
    "NeCF4":    "Ne/CF₄ 90/10",
    "ArCF4CO2": "Ar/CF₄/CO₂ 45/40/15",
    "PureCF4":  "Pure CF₄",
    # Pure gases
    "PureAr":     "Pure Ar",
    "PureHe":     "Pure He",
    "PureNe":     "Pure Ne",
    "PureEthane": "Pure C₂H₆",
    "PureIso":    "Pure iC₄H₁₀",
    "PureCO2":    "Pure CO₂",
}

# Gas group membership (PureCF4 appears in both for comparison)
MIXTURE_GASES = frozenset({
    "ArCF4", "HeEth", "ArCO2", "ArCF4Iso", "NeIso", "NeCF4", "ArCF4CO2", "PureCF4",
})
PURE_GASES = frozenset({
    "PureAr", "PureHe", "PureNe", "PureEthane", "PureIso", "PureCO2", "PureCF4",
})
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
        df["al_mm"] = 0.0
    if "primInDrift_fraction" not in df.columns:
        # Older summaries: use efficiency as best available proxy
        df["primInDrift_fraction"] = df["efficiency"]
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
    For each gas with Al shielding runs, produce a 4-row comparison per page:
      Row 0: mean primary ion pairs, log y-axis
      Row 1: primary-in-drift fraction, y fixed 0–1
      Row 2: mean primary ion pairs, linear y-axis
      Row 3: primary-in-drift fraction, y-max scaled to data
    X-axis: linear for electrons, log for all others.
    Neutron columns in rows 1–3 have x lower-limit set to 1e-2 MeV.
    """
    shielded_gases = df[df["al_mm"] > 0]["gas"].unique()
    if len(shielded_gases) == 0:
        return

    def _xaxis(ax, particle, first_row):
        if particle == "electron":
            ax.set_xscale("linear")
        else:
            ax.set_xscale("log")
            if particle == "neutron" and not first_row:
                ax.set_xlim(left=1e-2)

    out = outdir / "shielding_comparison.pdf"
    with PdfPages(out) as pdf:
        for gas in sorted(shielded_gases):
            sub = df[df["gas"] == gas].copy()
            al_thicknesses = sorted(sub["al_mm"].unique())

            particles = [p for p in ["gamma", "electron", "neutron", "proton", "muon"]
                         if p in sub["particle"].unique()]
            ncols = len(particles)
            fig, axes_grid = plt.subplots(4, ncols,
                                          figsize=(6 * ncols, 20),
                                          sharey=False, squeeze=False)

            for col, particle in enumerate(particles):
                psub = sub[sub["particle"] == particle].copy()
                ax0, ax1, ax2, ax3 = (axes_grid[r][col] for r in range(4))
                clip_low = 0.01

                all_means_raw = []
                all_fracs     = []

                for al_mm in al_thicknesses:
                    al_key = int(round(al_mm))
                    g = psub[np.isclose(psub["al_mm"], al_mm)].sort_values("energy_MeV")
                    if g.empty:
                        continue
                    color = AL_COLORS.get(al_key, "grey")
                    style = AL_STYLES.get(al_key, "-")
                    label = AL_LABELS.get(al_key, f"{al_mm:g} mm Al")

                    raw   = g["nPrimDrift_mean"]
                    sem   = g["nPrimDrift_std"] / np.sqrt(g["n_events"])
                    frac  = g["primInDrift_fraction"]
                    sem_f = np.sqrt(frac * (1 - frac) / g["n_events"])

                    all_means_raw.extend(raw.values)
                    all_fracs.extend(frac.values)

                    # rows 0 & 2: mean primaries (log-clipped vs raw)
                    for ax, mn, blo, bhi in [
                        (ax0, raw.clip(lower=clip_low),
                              (raw - sem).clip(lower=clip_low),
                              (raw + sem).clip(lower=clip_low)),
                        (ax2, raw,
                              (raw - sem).clip(lower=0),
                              (raw + sem)),
                    ]:
                        ax.fill_between(g["energy_MeV"].values, blo.values,
                                        bhi.values, alpha=0.2, color=color)
                        ax.plot(g["energy_MeV"], mn, style, color=color,
                                label=label, linewidth=1.8, markersize=4)

                    # rows 1 & 3: fraction
                    for ax in (ax1, ax3):
                        ax.fill_between(g["energy_MeV"].values,
                                        (frac - sem_f).clip(0).values,
                                        (frac + sem_f).clip(upper=1).values,
                                        alpha=0.2, color=color)
                        ax.plot(g["energy_MeV"], frac, style, color=color,
                                label=label, linewidth=1.8, markersize=4)

                # ---- axis formatting ----
                title = PARTICLE_TITLES.get(particle, particle)
                for row, ax in enumerate([ax0, ax1, ax2, ax3]):
                    _xaxis(ax, particle, first_row=(row == 0))
                    ax.set_xlabel("Particle energy [MeV]", fontsize=10)
                    ax.set_title(title, fontsize=11)
                    ax.legend(fontsize=8, loc="best")
                    ax.grid(True, which="both", alpha=0.3)

                # row 0: log-y mean primaries
                ax0.set_yscale("log")
                ax0.set_ylabel("Mean primary ion pairs in drift gap", fontsize=10)
                valid = [v for v in all_means_raw if v > 0]
                if valid:
                    ax0.set_ylim(bottom=max(clip_low, min(valid) * 0.3),
                                 top=max(valid) * 3)

                # row 1: fraction, fixed 0–1
                ax1.set_ylim(0, 1.05)
                ax1.set_ylabel("Fraction of primaries reaching drift", fontsize=10)
                ax1.axhline(1.0, color="k", lw=0.6, ls="--", alpha=0.4)

                # row 2: linear-y mean primaries
                ax2.set_yscale("linear")
                ax2.set_ylabel("Mean primary ion pairs in drift gap", fontsize=10)
                if all_means_raw:
                    ax2.set_ylim(bottom=0, top=max(all_means_raw) * 1.15)

                # row 3: fraction, y-max from data
                frac_max = max(all_fracs) if all_fracs else 1.0
                ax3.set_ylim(0, max(frac_max * 1.15, 0.02))
                ax3.set_ylabel("Fraction of primaries reaching drift", fontsize=10)

            fig.suptitle(f"Al Shielding Effect — {GAS_LABELS.get(gas, gas)}\n"
                         f"(3 cm drift gap, Mylar window + 2 cm air gap + Al)",
                         fontsize=12)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {out} ({len(shielded_gases)} gas(es))")


def fig_efficiency_vs_energy(df: pd.DataFrame, outdir: Path):
    """
    Line plot of primary-in-drift fraction vs energy, one panel per particle,
    one curve per gas.  Mirrors the sensitivity_vs_energy layout but shows
    what fraction of events produce any ionisation in the drift.
    """
    particles = [p for p in ["gamma", "electron", "neutron", "proton", "muon"]
                 if p in df["particle"].unique()]
    fig, axes = plt.subplots(1, len(particles),
                             figsize=(6 * len(particles), 5), sharey=False,
                             squeeze=False)

    for ax, particle in zip(axes[0], particles):
        sub = df[df["particle"] == particle].copy()
        for gas in sorted(sub["gas"].unique()):
            g = sub[sub["gas"] == gas].sort_values("energy_MeV")
            if g.empty:
                continue
            frac = g["primInDrift_fraction"]
            sem  = np.sqrt(frac * (1 - frac) / g["n_events"])
            color = GAS_COLORS.get(gas, "grey")
            ax.fill_between(g["energy_MeV"].values,
                            (frac - sem).clip(0).values,
                            (frac + sem).clip(upper=1).values,
                            alpha=0.2, color=color)
            ax.plot(g["energy_MeV"], frac, "o-", color=color,
                    label=GAS_LABELS.get(gas, gas), linewidth=1.8, markersize=4)

        ax.set_xscale("log")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Particle energy [MeV]", fontsize=11)
        ax.set_ylabel("Fraction of primaries reaching drift", fontsize=11)
        ax.set_title(PARTICLE_TITLES.get(particle, particle), fontsize=12)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, which="both", alpha=0.3)
        ax.axhline(1.0, color="k", lw=0.6, ls="--", alpha=0.4)

    fig.suptitle("Primary Transmission to Drift Gap\n"
                 "(fraction of events with primary reaching drift volume)",
                 fontsize=13)
    plt.tight_layout()
    out = outdir / "efficiency_vs_energy.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


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

    # Split into mixture and pure-gas subsets
    df_mix  = df_noshield[df_noshield["gas"].isin(MIXTURE_GASES)].copy()
    df_pure = df_noshield[df_noshield["gas"].isin(PURE_GASES)].copy()

    mix_outdir  = outdir / "mixtures"
    pure_outdir = outdir / "pure_gases"
    mix_outdir.mkdir(exist_ok=True)
    pure_outdir.mkdir(exist_ok=True)

    if not df_mix.empty:
        print("\nGenerating mixture gas plots...")
        fig_sensitivity_vs_energy(df_mix, mix_outdir)
        fig_efficiency(df_mix, mix_outdir)
        fig_efficiency_vs_energy(df_mix, mix_outdir)
        fig_relative_sensitivity(df_mix, mix_outdir)
    else:
        print("\n(no mixture gas data — skipping mixture plots)")

    if not df_pure.empty:
        print("\nGenerating pure gas plots...")
        fig_sensitivity_vs_energy(df_pure, pure_outdir)
        fig_efficiency(df_pure, pure_outdir)
        fig_efficiency_vs_energy(df_pure, pure_outdir)
        fig_relative_sensitivity(df_pure, pure_outdir)
    else:
        print("\n(no pure gas data — skipping pure gas plots)")

    fig_shielding_comparison(df, outdir)

    if args.rootdir:
        rootdir = Path(args.rootdir)
        print("Generating Nprim distribution plots from ROOT files...")
        fig_nprim_distributions_from_root(rootdir, outdir, df_noshield)

    print(f"\nAll plots saved to: {outdir}")


if __name__ == "__main__":
    main()
