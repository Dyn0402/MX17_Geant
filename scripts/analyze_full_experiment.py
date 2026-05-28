#!/usr/bin/env python3
"""
analyze_full_experiment.py
Analysis of full-experiment electron transmission and calorimetry.

Reads ROOT files produced by mm_sim -m full and generates:
  1. Transmission fraction per layer vs energy
  2. Mean energy deposition per layer vs energy
  3. LS calorimetry: total LS edep vs initial energy (linearity + resolution)
  4. Fully-contained fraction vs energy
  5. Layer-by-layer edep profile at selected energies
  6. Angular distributions in the Micromegas drift gas (multiple scattering)
  7. Angular resolution (RMS) vs energy

Usage
    python3 scripts/analyze_full_experiment.py \\
        --indir /eos/user/d/dneff/mx17_geant_sim_results/full \\
        --gas ArIso --particle electron \\
        --outfile full_experiment_analysis.pdf
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import pandas as pd

# ── Layer definitions ────────────────────────────────────────────────────────
# (edep_branch, primIn_branch, display_label, color)
LAYERS = [
    ("edepHe3Gas",    "primInHe3Gas",    "He-3 gas",        "#4db8ff"),
    ("edepDrift",     "primInDrift",     "MM drift gas",    "#2ca02c"),
    ("edepAmp",       "primInAmp",       "MM amp gas",      "#98df8a"),
    ("edepResistPaste", None,            "Resistive paste", "#bcbd22"),
    ("edepPCB",       "primInPCB",       "PCB stack",       "#ff7f0e"),
    ("edepScintWall", "primInScintWall", "Plastic scint.",  "#d62728"),
    ("edepLS1",       "primInLS1",       "Liq. scint. 1",   "#9467bd"),
    ("edepLS2",       "primInLS2",       "Liq. scint. 2",   "#c5b0d5"),
    ("edepLS3",       "primInLS3",       "Liq. scint. 3",   "#8c564b"),
    ("edepLS4",       "primInLS4",       "Liq. scint. 4",   "#e377c2"),
]

LS_EDEP_BRANCHES = ["edepLS1", "edepLS2", "edepLS3", "edepLS4"]

# Energies (MeV) used as snapshots for bar-chart and angular plots
SNAPSHOT_ENERGIES_MEV = [0.5, 1.0, 3.0, 5.0, 10.0]

# For angular analysis: max events to process per energy (speed limit)
MAX_EVENTS_ANGULAR = 3000

# ── File I/O ─────────────────────────────────────────────────────────────────

def parse_energy_from_tag(tag: str) -> float:
    """Extract energy in MeV from a job tag like full_ArIso_electron_0p5MeV."""
    m = re.search(r"_(\d+(?:p\d+)?)MeV", tag)
    if not m:
        return None
    return float(m.group(1).replace("p", "."))


def group_files(indir: str, gas: str, particle: str):
    """
    Return dict  energy_MeV -> [list of ROOT file paths].
    Handles both merged files and per-thread _t<N>.root files.
    """
    pattern = Path(indir) / f"full_{gas}_{particle}_*.root"
    all_files = sorted(Path(indir).glob(f"full_{gas}_{particle}_*.root"))
    if not all_files:
        sys.exit(f"No files found matching {pattern}")

    groups = defaultdict(list)
    for f in all_files:
        # Strip _tN suffix to get the base tag
        base = re.sub(r"_t\d+$", "", f.stem)
        energy = parse_energy_from_tag(base)
        if energy is not None:
            groups[energy].append(str(f))

    return dict(sorted(groups.items()))


def read_event_tree(files: list) -> pd.DataFrame:
    """Concatenate EventTree from one or more ROOT files into a DataFrame."""
    dfs = []
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                if "EventTree" not in f:
                    continue
                dfs.append(f["EventTree"].arrays(library="pd"))
        except Exception as e:
            print(f"  Warning: could not read {fpath}: {e}")
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def read_cluster_tree(files: list, max_events: int = MAX_EVENTS_ANGULAR) -> pd.DataFrame:
    """
    Read ClusterTree (DriftGas, primary track only) for angular analysis.
    Limits to the first max_events events for speed.
    """
    rows = []
    n_events_seen = set()

    for fpath in files:
        if len(n_events_seen) >= max_events:
            break
        try:
            with uproot.open(fpath) as f:
                if "ClusterTree" not in f:
                    continue
                tree = f["ClusterTree"]
                branches = ["eventID", "trackID", "x", "y", "z", "volume"]
                arrays = tree.arrays(branches, library="np")

                # Decode volume byte-strings (Char_t[32] → str)
                vols = arrays["volume"]
                if vols.dtype.kind == "S":  # byte strings
                    vol_str = np.array([v.decode("ascii").rstrip("\x00")
                                        for v in vols])
                else:
                    vol_str = vols.astype(str)

                # Filter: DriftGas and primary track (trackID==1)
                mask = (vol_str == "DriftGas") & (arrays["trackID"] == 1)
                evids = arrays["eventID"][mask]
                xs    = arrays["x"][mask]
                ys    = arrays["y"][mask]
                zs    = arrays["z"][mask]

                for evid, x, y, z in zip(evids, xs, ys, zs):
                    if evid not in n_events_seen:
                        if len(n_events_seen) >= max_events:
                            break
                        n_events_seen.add(evid)
                    rows.append({"eventID": evid, "x": x, "y": y, "z": z})
        except Exception as e:
            print(f"  Warning: could not read clusters from {fpath}: {e}")

    if not rows:
        return pd.DataFrame(columns=["eventID", "x", "y", "z"])
    return pd.DataFrame(rows)


# ── Analysis functions ────────────────────────────────────────────────────────

def compute_summary(groups: dict) -> pd.DataFrame:
    """
    For each energy, load EventTree and compute per-layer statistics.
    Returns a DataFrame with one row per energy.
    """
    records = []
    n_total = len(groups)
    for idx, (energy, files) in enumerate(groups.items(), 1):
        print(f"  [{idx}/{n_total}] E = {energy:.4g} MeV  ({len(files)} file(s))")
        df = read_event_tree(files)
        if df.empty:
            continue

        n = len(df)
        row = {"energy_MeV": energy, "n_events": n}

        for (edep_br, primin_br, label, _color) in LAYERS:
            # Mean edep over all events [eV → MeV]
            if edep_br in df.columns:
                row[f"edep_{edep_br}"] = df[edep_br].mean() / 1e6
                row[f"edep_std_{edep_br}"] = df[edep_br].std() / 1e6
            # Transmission fraction
            if primin_br and primin_br in df.columns:
                row[f"trans_{primin_br}"] = df[primin_br].mean()

        # Total LS edep per event [eV → MeV]
        ls_cols = [b for b in LS_EDEP_BRANCHES if b in df.columns]
        if ls_cols:
            ls_total = df[ls_cols].sum(axis=1) / 1e6  # MeV per event
            row["ls_total_mean"]  = ls_total.mean()
            row["ls_total_std"]   = ls_total.std()
            row["ls_total_median"]= ls_total.median()

            # Containment: reached LS1 AND deposited >50% of (energy - typical upstream loss)
            # Proxy: events where ls_total > 0.5 * energy (loose; upstream losses ~few %)
            row["frac_in_ls"]  = (ls_total > 0).mean()
            row["frac_any_ls4"] = df["primInLS4"].mean() if "primInLS4" in df.columns else np.nan

            # Energy resolution: for events that reach LS1, σ/μ of LS total
            if "primInLS1" in df.columns:
                mask_in = df["primInLS1"].astype(bool)
                ls_in = ls_total[mask_in]
                if len(ls_in) > 10:
                    row["ls_res_mean"] = ls_in.mean()
                    row["ls_res_std"]  = ls_in.std()
                    row["ls_res_fwhm"] = 2.355 * ls_in.std() / ls_in.mean() if ls_in.mean() > 0 else np.nan

        records.append(row)

    return pd.DataFrame(records).sort_values("energy_MeV").reset_index(drop=True)


def compute_angles(cluster_df: pd.DataFrame) -> np.ndarray:
    """
    For each event in cluster_df, fit a straight line to the primary track
    in the DriftGas (x vs z, y vs z) and return the polar angle in mrad.
    Returns array of angles [mrad], one per event with >= 5 clusters.
    """
    angles = []
    for evid, grp in cluster_df.groupby("eventID"):
        if len(grp) < 5:
            continue
        z = grp["z"].values
        x = grp["x"].values
        y = grp["y"].values
        z_span = z.max() - z.min()
        if z_span < 0.5:  # less than 0.5 mm of track — skip
            continue
        # Linear fit: x = ax*z + bx, y = ay*z + by
        ax = np.polyfit(z, x, 1)[0]  # dx/dz
        ay = np.polyfit(z, y, 1)[0]  # dy/dz
        theta = np.degrees(np.arctan(np.sqrt(ax**2 + ay**2))) * 1000 / (180/np.pi)
        # Simpler: theta in mrad = 1000 * arctan(sqrt(ax^2+ay^2)) radians
        theta_rad = np.arctan(np.sqrt(ax**2 + ay**2))
        angles.append(np.degrees(theta_rad))  # store in degrees for readability
    return np.array(angles)


def find_snapshot_files(groups: dict, target_energies: list) -> dict:
    """
    For each target energy, find the closest available energy in groups.
    Returns dict: target_energy -> (actual_energy, files).
    """
    avail = np.array(sorted(groups.keys()))
    result = {}
    for e_target in target_energies:
        idx = np.argmin(np.abs(avail - e_target))
        e_actual = avail[idx]
        result[e_target] = (e_actual, groups[e_actual])
    return result


# ── Plotting ──────────────────────────────────────────────────────────────────

STYLE = dict(linewidth=1.6, marker="none")

def _set_xaxis(ax, xlabel="Electron energy (MeV)"):
    ax.set_xlabel(xlabel, fontsize=11)


def plot_transmission(pdf: PdfPages, summary: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9, 5.5))

    en = summary["energy_MeV"].values
    for (edep_br, primin_br, label, color) in LAYERS:
        col = f"trans_{primin_br}"
        if primin_br and col in summary.columns:
            ax.plot(en, summary[col].values, label=label, color=color, **STYLE)

    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel("Transmission fraction", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Fraction of electrons reaching each layer", fontsize=12)
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    ax.axhline(1.0, color="grey", lw=0.6, ls="--")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_edep(pdf: PdfPages, summary: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9, 5.5))

    en = summary["energy_MeV"].values
    for (edep_br, primin_br, label, color) in LAYERS:
        col = f"edep_{edep_br}"
        if col in summary.columns:
            vals = summary[col].values
            # Skip layers with negligible edep
            if np.nanmax(vals) < 1e-8:
                continue
            ax.plot(en, vals, label=label, color=color, **STYLE)

    ax.set_yscale("log")
    ax.set_ylabel("Mean energy deposition (MeV/electron)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Mean energy deposited per layer (averaged over all events)", fontsize=12)
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_calorimetry(pdf: PdfPages, summary: pd.DataFrame):
    """Three-panel calorimetry figure."""
    if "ls_total_mean" not in summary.columns:
        return

    en   = summary["energy_MeV"].values
    mean = summary["ls_total_mean"].values
    std  = summary["ls_total_std"].values

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: total LS edep vs initial energy
    ax = axes[0]
    ax.errorbar(en, mean, yerr=std, fmt="o", ms=3, lw=1,
                color="#9467bd", ecolor="#c5b0d5", capsize=2, label="Mean ± 1σ")
    # Ideal line (slope=1 if no upstream losses)
    ax.plot(en, en, "k--", lw=0.8, label="E_dep = E_initial")
    ax.set_yscale("log")
    ax.set_xlabel("Initial electron energy (MeV)", fontsize=10)
    ax.set_ylabel("Total LS edep (MeV)", fontsize=10)
    ax.set_title("LS calorimeter response", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Panel B: linearity — edep/E_initial
    ax = axes[1]
    linearity = mean / en
    ax.plot(en, linearity, "o-", ms=3, lw=1.4, color="#9467bd")
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Initial electron energy (MeV)", fontsize=10)
    ax.set_ylabel("E_dep(LS) / E_initial", fontsize=10)
    ax.set_title("Calorimeter linearity", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)

    # Panel C: energy resolution σ/μ
    ax = axes[2]
    if "ls_res_std" in summary.columns:
        res_mean = summary["ls_res_mean"].values
        res_std  = summary["ls_res_std"].values
        # Only for events that reach LS
        resolution = res_std / res_mean
        ax.plot(en, resolution * 100, "o-", ms=3, lw=1.4, color="#d62728")
        ax.set_xlabel("Initial electron energy (MeV)", fontsize=10)
        ax.set_ylabel("σ/μ  (%)", fontsize=10)
        ax.set_title("LS energy resolution\n(events reaching LS1)", fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_ylim(bottom=0)

    fig.suptitle("Liquid scintillator calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_containment(pdf: PdfPages, summary: pd.DataFrame):
    """
    Containment analysis: fraction reaching LS1 and fraction reaching LS4.
    Overlay of both requirements highlights the 'fully contained' window.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    en = summary["energy_MeV"].values

    # Left: transmission to first and last LS layer
    ax = axes[0]
    for col, label, color in [
        ("trans_primInLS1", "Reaches LS1 (entry)",  "#9467bd"),
        ("trans_primInLS4", "Reaches LS4 (last layer)", "#e377c2"),
    ]:
        if col in summary.columns:
            ax.plot(en, summary[col].values, label=label, color=color, lw=1.8)
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel("Fraction of electrons", fontsize=10)
    _set_xaxis(ax)
    ax.set_title("LS layer entry fractions", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    # Right: "fully contained" proxy
    # Good calorimetry requires: enters LS1 AND is stopped before exiting LS4.
    # Proxy: enters LS1 but does NOT reach LS4 (stopped in LS1-LS3),
    # OR enters LS4 but deposits significant energy there.
    # We approximate: fraction_enters_LS1 - fraction_exits_past_LS4 (i.e. enters LS4 with enough energy to punch through)
    # A simpler observable: enters LS1 AND the LS total edep / E_initial is high.
    ax = axes[1]
    if "ls_total_mean" in summary.columns and "trans_primInLS1" in summary.columns:
        frac_enters = summary["trans_primInLS1"].values
        linearity   = summary["ls_total_mean"].values / en

        # "Fully contained" proxy: enters LS1 AND linearity > 0.7 (rough)
        # We can't get this event-by-event from summary alone, so show the two
        # components and note what they imply.
        ax.plot(en, frac_enters,   label="Enters LS1",    color="#9467bd", lw=1.8)
        ax.plot(en, linearity, label="LS edep / E₀", color="#ff7f0e",
                lw=1.8, linestyle="--")
        ax.axhline(1.0, color="grey", lw=0.7, ls=":")
        ax.set_ylim(-0.05, 1.2)
        ax.set_ylabel("Fraction / ratio", fontsize=10)
        _set_xaxis(ax)
        ax.set_title("Containment quality\n(high linearity = well contained)", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, which="both", alpha=0.3)

    fig.suptitle("Containment analysis for LS calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_layer_profiles(pdf: PdfPages, groups: dict, snapshot_energies: list):
    """Bar chart of edep in each layer at selected energies."""
    snapshots = find_snapshot_files(groups, snapshot_energies)

    layer_labels = [label for (_, _, label, _) in LAYERS]
    layer_colors = [color for (_, _, _, color) in LAYERS]
    layer_edep_cols = [f"edep_{edep_br}" for (edep_br, _, _, _) in LAYERS]

    # Build summary for snapshot energies only
    snap_records = {}
    for e_target, (e_actual, files) in snapshots.items():
        df = read_event_tree(files)
        if df.empty:
            continue
        row = {}
        for (edep_br, _, label, _) in LAYERS:
            if edep_br in df.columns:
                row[label] = df[edep_br].mean() / 1e6
            else:
                row[label] = 0.0
        snap_records[e_actual] = row

    if not snap_records:
        return

    fig, ax = plt.subplots(figsize=(11, 5))
    n_layers = len(LAYERS)
    n_energies = len(snap_records)
    bar_width = 0.7 / n_energies
    x = np.arange(n_layers)

    energies_sorted = sorted(snap_records.keys())
    colors_e = cm.viridis(np.linspace(0.2, 0.85, n_energies))

    for i, e in enumerate(energies_sorted):
        vals = [snap_records[e].get(label, 0) for (_, _, label, _) in LAYERS]
        offset = (i - n_energies / 2 + 0.5) * bar_width
        ax.bar(x + offset, vals, bar_width * 0.9,
               label=f"{e:.3g} MeV", color=colors_e[i], alpha=0.85)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(layer_labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Mean energy deposition (MeV/electron)", fontsize=10)
    ax.set_title("Energy deposition profile at selected electron energies", fontsize=12)
    ax.legend(title="E₀", fontsize=9, title_fontsize=9)
    ax.grid(True, axis="y", which="both", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_angular_distributions(pdf: PdfPages, groups: dict,
                                snapshot_energies: list):
    """
    Reconstruct straight-line tracks in the MM drift gas.
    Plot angle distribution histograms for selected energies.
    """
    snapshots = find_snapshot_files(groups, snapshot_energies)
    n_snap = len(snapshots)
    colors = cm.plasma(np.linspace(0.15, 0.85, n_snap))

    fig, axes = plt.subplots(1, n_snap, figsize=(3.2 * n_snap, 4),
                              sharey=False, sharex=False)
    if n_snap == 1:
        axes = [axes]

    for ax, ((e_target, (e_actual, files)), color) in zip(
            axes, zip(snapshots.items(), colors)):
        print(f"  Angular analysis at {e_actual:.3g} MeV ...")
        cdf = read_cluster_tree(files, max_events=MAX_EVENTS_ANGULAR)
        if cdf.empty:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        angles_deg = compute_angles(cdf)
        if len(angles_deg) == 0:
            ax.text(0.5, 0.5, "No tracks", transform=ax.transAxes, ha="center")
            continue

        angles_mrad = np.deg2rad(angles_deg) * 1000  # mrad
        rms_mrad = np.sqrt(np.mean(angles_mrad**2))

        ax.hist(angles_mrad, bins=60, density=True, color=color,
                alpha=0.8, edgecolor="none")
        ax.axvline(rms_mrad, color="k", lw=1.2, ls="--",
                   label=f"RMS = {rms_mrad:.1f} mrad")
        ax.set_xlabel("Polar angle (mrad)", fontsize=9)
        ax.set_ylabel("Density", fontsize=9) if ax == axes[0] else None
        ax.set_title(f"{e_actual:.3g} MeV", fontsize=10)
        ax.legend(fontsize=8)
        ax.set_xlim(left=0)

    fig.suptitle("Track angle distribution in MM drift gas\n"
                 "(spread from multiple scattering in upstream material)",
                 fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_angular_resolution(pdf: PdfPages, groups: dict,
                             step: int = 10):
    """
    RMS polar angle reconstructed in the drift gas vs electron energy.
    Processes every `step`-th energy to keep runtime manageable.
    """
    energies_all = sorted(groups.keys())
    sel_energies  = energies_all[::step]

    rms_mrad_vals = []
    en_vals       = []

    n_total = len(sel_energies)
    for idx, energy in enumerate(sel_energies, 1):
        print(f"  Angular RMS [{idx}/{n_total}]  E = {energy:.4g} MeV")
        files = groups[energy]
        cdf   = read_cluster_tree(files, max_events=MAX_EVENTS_ANGULAR)
        if cdf.empty:
            continue
        angles_deg  = compute_angles(cdf)
        if len(angles_deg) < 20:
            continue
        angles_mrad = np.deg2rad(angles_deg) * 1000
        rms_mrad_vals.append(np.sqrt(np.mean(angles_mrad**2)))
        en_vals.append(energy)

    if not en_vals:
        return

    en_vals       = np.array(en_vals)
    rms_mrad_vals = np.array(rms_mrad_vals)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(en_vals, rms_mrad_vals, "o-", ms=4, lw=1.6, color="#1f77b4",
            label="RMS polar angle")

    # Highland formula reference (crude estimate for electrons through ~few g/cm2)
    # theta_rms ~ (13.6 MeV / (beta*p)) * sqrt(x/X0) * [1 + 0.038*ln(x/X0)]
    # rough upstream x/X0 ≈ 0.05 (He-3 target + walls + air + MM dead layers)
    x_over_X0 = 0.05
    p_vals = en_vals  # ultra-relativistic approximation p ≈ E for electrons
    theta_highland_rad = (13.6e-3 / p_vals) * np.sqrt(x_over_X0) * (
        1 + 0.038 * np.log(x_over_X0))
    ax.plot(en_vals, theta_highland_rad * 1000, "k--", lw=1,
            label=f"Highland (x/X₀≈{x_over_X0})")

    ax.set_ylabel("RMS polar angle (mrad)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Angular resolution from MS in upstream material\n"
                 "(reconstructed from straight-line fit in MM drift gas)",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(bottom=0)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Full-experiment analysis",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--indir",    required=True,
                   help="Directory containing full-experiment ROOT files")
    p.add_argument("--gas",      default="ArIso",
                   help="Gas tag used in filenames")
    p.add_argument("--particle", default="electron",
                   help="Particle type used in filenames")
    p.add_argument("--outfile",  default="full_experiment_analysis.pdf",
                   help="Output PDF path")
    p.add_argument("--angular-step", type=int, default=10,
                   help="Process every Nth energy point for angular RMS plot")
    p.add_argument("--no-angular", action="store_true",
                   help="Skip angular analysis (faster)")
    return p.parse_args()


def main():
    args = parse_args()

    plt.rcParams.update({
        "figure.dpi":      150,
        "axes.titlesize":  11,
        "axes.labelsize":  10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })

    print(f"Scanning {args.indir} for {args.gas}/{args.particle} files ...")
    groups = group_files(args.indir, args.gas, args.particle)
    n_energies = len(groups)
    print(f"Found {n_energies} energy points: "
          f"{min(groups):.4g} – {max(groups):.4g} MeV")

    # ── Event-level summary ───────────────────────────────────────────────
    print("\nLoading EventTree summaries ...")
    summary = compute_summary(groups)
    if summary.empty:
        sys.exit("No data could be loaded.")

    summary_csv = Path(args.outfile).with_suffix(".csv")
    summary.to_csv(summary_csv, index=False)
    print(f"Summary table saved → {summary_csv}")

    # ── Plots ─────────────────────────────────────────────────────────────
    print(f"\nWriting plots → {args.outfile}")
    with PdfPages(args.outfile) as pdf:

        # Cover page
        fig = plt.figure(figsize=(8, 4))
        fig.text(0.5, 0.65, "Full-experiment stack analysis",
                 ha="center", va="center", fontsize=18, fontweight="bold")
        fig.text(0.5, 0.5,
                 f"Gas: {args.gas}   Particle: {args.particle}\n"
                 f"Energy range: {min(groups):.4g} – {max(groups):.4g} MeV "
                 f"({n_energies} points)\n"
                 f"Input: {args.indir}",
                 ha="center", va="center", fontsize=10, color="grey")
        pdf.savefig(fig)
        plt.close(fig)

        print("  Plot 1: Transmission fractions ...")
        plot_transmission(pdf, summary)

        print("  Plot 2: Energy deposition per layer ...")
        plot_edep(pdf, summary)

        print("  Plot 3: LS calorimetry ...")
        plot_calorimetry(pdf, summary)

        print("  Plot 4: Containment analysis ...")
        plot_containment(pdf, summary)

        print("  Plot 5: Layer edep profiles at selected energies ...")
        plot_layer_profiles(pdf, groups, SNAPSHOT_ENERGIES_MEV)

        if not args.no_angular:
            print("  Plot 6: Angular distributions ...")
            plot_angular_distributions(pdf, groups, SNAPSHOT_ENERGIES_MEV)

            print(f"  Plot 7: Angular resolution vs energy "
                  f"(every {args.angular_step}th point) ...")
            plot_angular_resolution(pdf, groups, step=args.angular_step)
        else:
            print("  Skipping angular analysis (--no-angular)")

    print(f"\nDone.  Output: {args.outfile}")


if __name__ == "__main__":
    main()
