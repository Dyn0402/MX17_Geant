#!/usr/bin/env python3
"""
analyze_full_experiment.py
Analysis of full-experiment electron/positron transmission and calorimetry.

Color scheme: same detector layer → same color (grouped by subsystem).
Line style  : solid = electron, dashed = positron.

Each transmission/edep plot is produced twice:
  - Full version  : all 10 layers
  - Simple version: MM drift gas, Plastic scint., Liq. scint. 1, Liq. scint. 2

Usage
    python3 scripts/analyze_full_experiment.py \\
        --indir /eos/user/d/dneff/mx17_geant_sim_results/full \\
        --gas ArIso \\
        --particles electron positron \\
        --outfile full_experiment_analysis.pdf
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import pandas as pd

# ── Color scheme (grouped by detector subsystem) ─────────────────────────────
# He-3 target:   blue
# MM elements:   green family  (dark → light)
# PCB:           orange
# Plastic scint: red
# LS elements:   purple family (dark → light)

LAYERS = [
    # (edep_branch, primIn_branch, label, color)
    ("edepHe3Gas",      "primInHe3Gas",    "He-3 gas",        "#1f77b4"),
    ("edepDrift",       "primInDrift",     "MM drift gas",    "#2ca02c"),
    ("edepAmp",         "primInAmp",       "MM amp gas",      "#74c476"),
    ("edepResistPaste", None,              "Resistive paste", "#a1d99b"),
    ("edepPCB",         "primInPCB",       "PCB stack",       "#ff7f0e"),
    ("edepScintWall",   "primInScintWall", "Plastic scint.",  "#d62728"),
    ("edepLS1",         "primInLS1",       "Liq. scint. 1",   "#6a1d8a"),
    ("edepLS2",         "primInLS2",       "Liq. scint. 2",   "#9b59b6"),
    ("edepLS3",         "primInLS3",       "Liq. scint. 3",   "#bb8fce"),
    ("edepLS4",         "primInLS4",       "Liq. scint. 4",   "#d7bde2"),
]

SIMPLIFIED_LAYERS = [
    ("edepDrift",     "primInDrift",     "MM drift gas",   "#2ca02c"),
    ("edepScintWall", "primInScintWall", "Plastic scint.", "#d62728"),
    ("edepLS1",       "primInLS1",       "Liq. scint. 1",  "#6a1d8a"),
    ("edepLS4",       "primInLS4",       "Liq. scint. 4",  "#d7bde2"),
]

# Coarse grouping: aggregate detector subsystems
# (label, color, [edep branches to sum], trans_col for entry transmission)
COARSE_GROUPS = [
    ("He-3 target",        "#1f77b4",
     ["edepHe3Gas"],
     "primInHe3Gas"),
    ("Micromegas (all)",   "#2ca02c",
     ["edepDrift", "edepAmp", "edepResistPaste"],
     "primInDrift"),
    ("PCB stack",          "#ff7f0e",
     ["edepPCB"],
     "primInPCB"),
    ("Plastic scint.",     "#d62728",
     ["edepScintWall"],
     "primInScintWall"),
    ("Liquid scint. (all)","#6a1d8a",
     ["edepLS1", "edepLS2", "edepLS3", "edepLS4"],
     "primInLS1"),
]

# Cumulative path from gun (He-3 centre, z = 0) to entry face of each
# scored volume [mm].  Used for survival-vs-path plots.
# He-3 gas radius 25 mm + exit wall (Al 0.5 + CFRP 0.9) + 200 mm air
# + MM dead layers (~0.1 mm) → drift gas entry ≈ 226.5 mm.
SURVIVAL_COLS = [
    # (trans_col,        z_entry_mm, label)
    ("primInHe3Gas",       25.0,   "He-3 exit"),
    ("primInDrift",       226.5,   "MM drift"),
    ("primInAmp",         256.5,   "MM amp"),
    ("primInPCB",         256.8,   "PCB"),
    ("primInScintWall",   282.6,   "Plastic scint."),
    ("primInLS1",         307.3,   "LS1"),
    ("primInLS2",         323.8,   "LS2"),
    ("primInLS3",         340.3,   "LS3"),
    ("primInLS4",         356.8,   "LS4"),
]
LS4_EXIT_MM = 371.8   # end of LS4 (LS_CFRP_5 skipped — not scored)

SURVIVAL_ENERGIES_MEV = [1.0, 3.0, 5.0, 10.0, 16.5]

# Schematic stack for the diagram strip.
# (z_start_mm, z_end_mm, face_color, text_color, label)
# Air gaps are omitted — they take the figure background color.
DETECTOR_STACK = [
    ( 0.0,   25.0,  "#1f77b4", "white",  "He-3 gas"),   # He-3 gas (half-diameter)
    (25.0,   25.5,  "#aec7e8", "black",  ""),            # Al capsule wall
    (25.5,   26.4,  "#555555", "white",  ""),            # CFRP capsule wall
    # 26.4 – 226.4  Air gap 1 (background)
    (226.4,  226.5, "#c7e9c0", "black",  ""),            # MM dead layers (Mylar+cathode)
    (226.5,  256.5, "#2ca02c", "white",  "Drift"),       # Drift gas  30 mm
    (256.5,  256.55,"#808080", "white",  ""),            # Micromesh
    (256.55, 257.0, "#74c476", "black",  ""),            # Amp gas + resistive paste
    (257.0,  262.6, "#ff7f0e", "black",  "PCB"),         # PCB stack  5.6 mm
    # 262.6 – 282.6  Air gap 2 (background)
    (282.6,  285.8, "#d62728", "white",  "Scint."),      # Plastic scint. wall  3.2 mm
    # 285.8 – 305.8  Air gap 3 (background)
    (305.8,  307.3, "#555555", "white",  ""),            # CFRP  1.5 mm
    (307.3,  322.3, "#6a1d8a", "white",  "LS 1"),        # Liq. scint. 1  15 mm
    (322.3,  323.8, "#555555", "white",  ""),            # CFRP
    (323.8,  338.8, "#9b59b6", "white",  "LS 2"),        # Liq. scint. 2
    (338.8,  340.3, "#555555", "white",  ""),            # CFRP
    (340.3,  355.3, "#bb8fce", "black",  "LS 3"),        # Liq. scint. 3
    (355.3,  356.8, "#555555", "white",  ""),            # CFRP
    (356.8,  371.8, "#d7bde2", "black",  "LS 4"),        # Liq. scint. 4
]

# Line style per particle
PARTICLE_STYLES = {
    "electron": {"ls": "-",  "lw": 1.8, "label": "e⁻"},
    "positron": {"ls": "--", "lw": 1.8, "label": "e⁺"},
}

LS_EDEP_BRANCHES  = ["edepLS1", "edepLS2", "edepLS3", "edepLS4"]
SNAPSHOT_ENERGIES = [0.5, 1.0, 3.0, 5.0, 10.0]
MAX_EVENTS_ANG    = 3000


# ── File I/O ──────────────────────────────────────────────────────────────────

def parse_energy(tag: str):
    m = re.search(r"_(\d+(?:p\d+)?)MeV", tag)
    return float(m.group(1).replace("p", ".")) if m else None


def group_files(indir: str, gas: str, particle: str) -> dict:
    """energy_MeV → [root file paths]"""
    groups = defaultdict(list)
    for f in sorted(Path(indir).glob(f"full_{gas}_{particle}_*.root")):
        base   = re.sub(r"_t\d+$", "", f.stem)
        energy = parse_energy(base)
        if energy is not None:
            groups[energy].append(str(f))
    if not groups:
        print(f"  No files found for {gas}/{particle} in {indir}")
    return dict(sorted(groups.items()))


def read_event_tree(files: list) -> pd.DataFrame:
    dfs = []
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                if "EventTree" in f:
                    dfs.append(f["EventTree"].arrays(library="pd"))
        except Exception as e:
            print(f"  Warning: {fpath}: {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def read_cluster_tree(files: list, max_events: int = MAX_EVENTS_ANG) -> pd.DataFrame:
    rows, seen = [], set()
    for fpath in files:
        if len(seen) >= max_events:
            break
        try:
            with uproot.open(fpath) as f:
                if "ClusterTree" not in f:
                    continue
                t = f["ClusterTree"]
                a = t.arrays(["eventID", "trackID", "x", "y", "z", "volume"],
                              library="np")
                vols = a["volume"]
                vol_str = (np.array([v.decode("ascii").rstrip("\x00") for v in vols])
                           if vols.dtype.kind == "S" else vols.astype(str))
                mask = (vol_str == "DriftGas") & (a["trackID"] == 1)
                for ev, x, y, z in zip(a["eventID"][mask], a["x"][mask],
                                        a["y"][mask], a["z"][mask]):
                    if ev not in seen:
                        if len(seen) >= max_events:
                            break
                        seen.add(ev)
                    rows.append({"eventID": ev, "x": x, "y": y, "z": z})
        except Exception as e:
            print(f"  Warning (clusters): {fpath}: {e}")
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=["eventID","x","y","z"])


# ── Summary computation ───────────────────────────────────────────────────────

def compute_summary(groups: dict) -> pd.DataFrame:
    records = []
    n = len(groups)
    for idx, (energy, files) in enumerate(groups.items(), 1):
        print(f"    [{idx}/{n}] {energy:.4g} MeV")
        df = read_event_tree(files)
        if df.empty:
            continue
        n_ev = len(df)
        row  = {"energy_MeV": energy, "n_events": n_ev}

        for (edep_br, primin_br, label, _) in LAYERS:
            if edep_br in df.columns:
                row[f"edep_{edep_br}"]     = df[edep_br].mean() / 1e6
                row[f"edep_std_{edep_br}"] = df[edep_br].std()  / 1e6
            if primin_br and primin_br in df.columns:
                row[f"trans_{primin_br}"] = df[primin_br].mean()

        ls_cols = [b for b in LS_EDEP_BRANCHES if b in df.columns]
        if ls_cols:
            ls_tot = df[ls_cols].sum(axis=1) / 1e6
            row["ls_total_mean"]   = ls_tot.mean()
            row["ls_total_std"]    = ls_tot.std()
            row["ls_total_median"] = ls_tot.median()
            row["frac_any_ls4"]    = (df["primInLS4"].mean()
                                      if "primInLS4" in df.columns else np.nan)
            if "primInLS1" in df.columns:
                ls_in = ls_tot[df["primInLS1"].astype(bool)]
                if len(ls_in) > 10:
                    row["ls_res_mean"] = ls_in.mean()
                    row["ls_res_std"]  = ls_in.std()

        records.append(row)
    return pd.DataFrame(records).sort_values("energy_MeV").reset_index(drop=True)


def compute_angles(cdf: pd.DataFrame) -> np.ndarray:
    angles = []
    for _, grp in cdf.groupby("eventID"):
        if len(grp) < 5:
            continue
        z, x, y = grp["z"].values, grp["x"].values, grp["y"].values
        if z.max() - z.min() < 0.5:
            continue
        ax_ = np.polyfit(z, x, 1)[0]
        ay_ = np.polyfit(z, y, 1)[0]
        angles.append(np.degrees(np.arctan(np.sqrt(ax_**2 + ay_**2))))
    return np.array(angles)


def find_snapshots(groups: dict, targets: list) -> dict:
    avail = np.array(sorted(groups.keys()))
    return {t: (avail[np.argmin(np.abs(avail - t))],
                groups[avail[np.argmin(np.abs(avail - t))]]) for t in targets}


# ── Legend helpers ────────────────────────────────────────────────────────────

def _layer_legend(ax, layer_list, loc="lower right", ncol=2):
    """Color-patch legend for layers."""
    handles = [mpatches.Patch(color=color, label=label)
               for (_, _, label, color) in layer_list]
    return ax.legend(handles=handles, loc=loc, fontsize=7,
                     ncol=ncol, title="Layer", title_fontsize=7)


def _particle_legend(ax, summaries, loc="upper right"):
    """Line-style legend for particles."""
    handles = [mlines.Line2D([], [], color="k",
                             ls=PARTICLE_STYLES[p]["ls"],
                             lw=PARTICLE_STYLES[p]["lw"],
                             label=PARTICLE_STYLES[p]["label"])
               for p in summaries if p in PARTICLE_STYLES]
    if handles:
        return ax.legend(handles=handles, loc=loc, fontsize=9,
                         title="Particle", title_fontsize=8)
    return None


def _add_both_legends(ax, layer_list, summaries,
                      layer_loc="lower right", part_loc="upper right"):
    leg1 = _layer_legend(ax, layer_list, loc=layer_loc)
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc=part_loc)


E_LO, E_HI = 4.0, 16.5  # physical e0 range [MeV]

def _add_energy_lines(ax):
    """Mark the physical electron energy range with vertical lines."""
    ax.axvline(E_LO, color="grey", lw=1.0, ls="--", zorder=0)
    ax.axvline(E_HI, color="grey", lw=1.0, ls="--", zorder=0)
    # Shade the region between the two lines
    ylim = ax.get_ylim()
    ax.axvspan(E_LO, E_HI, alpha=0.06, color="grey", zorder=0)
    ax.set_ylim(ylim)  # restore in case axvspan shifted it

def _set_xaxis(ax, xlabel="Electron energy (MeV)"):
    ax.set_xlabel(xlabel, fontsize=11)
    _add_energy_lines(ax)


# ── Plot functions ────────────────────────────────────────────────────────────

def _plot_transmission_on(ax, summaries: dict, layer_list: list):
    for (_, primin_br, label, color) in layer_list:
        if not primin_br:
            continue
        col = f"trans_{primin_br}"
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6})
            ax.plot(summary["energy_MeV"].values, summary[col].values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_ylim(-0.02, 1.05)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("Transmission fraction", fontsize=11)
    ax.grid(True, alpha=0.3)


def plot_transmission(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _plot_transmission_on(ax, summaries, LAYERS)
    _set_xaxis(ax)
    ax.set_title("Fraction reaching each layer — all layers", fontsize=12)
    _add_both_legends(ax, LAYERS, summaries)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_transmission_simple(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    _plot_transmission_on(ax, summaries, SIMPLIFIED_LAYERS)
    _set_xaxis(ax)
    ax.set_title("Fraction reaching each layer — key layers", fontsize=12)
    _add_both_legends(ax, SIMPLIFIED_LAYERS, summaries,
                      layer_loc="lower right", part_loc="center right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def _plot_edep_on(ax, summaries: dict, layer_list: list):
    for (edep_br, _, label, color) in layer_list:
        col = f"edep_{edep_br}"
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            vals = summary[col].values
            if np.nanmax(vals) < 1e-8:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6})
            ax.plot(summary["energy_MeV"].values, vals,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_yscale("log")
    ax.set_ylabel("Mean energy deposition (MeV / primary)", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)


def plot_edep(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _plot_edep_on(ax, summaries, LAYERS)
    _set_xaxis(ax)
    ax.set_title("Mean energy deposited per layer (all events) — all layers",
                 fontsize=12)
    _add_both_legends(ax, LAYERS, summaries, layer_loc="upper left")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_edep_simple(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    _plot_edep_on(ax, summaries, SIMPLIFIED_LAYERS)
    _set_xaxis(ax)
    ax.set_title("Mean energy deposited per layer (all events) — key layers",
                 fontsize=12)
    _add_both_legends(ax, SIMPLIFIED_LAYERS, summaries, layer_loc="upper left")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_calorimetry(pdf: PdfPages, summaries: dict):
    if not any("ls_total_mean" in s.columns for s in summaries.values()):
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: LS edep vs initial energy
    ax = axes[0]
    e_max = 0.0
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en   = summary["energy_MeV"].values
        mean = summary["ls_total_mean"].values
        std  = summary["ls_total_std"].values
        ax.errorbar(en, mean, yerr=std, fmt="o", ms=2.5, lw=1,
                    ls=pst["ls"], color="#9467bd" if particle == "electron" else "#e377c2",
                    ecolor="lightgrey", capsize=1.5, label=pst["label"])
        e_max = max(e_max, en.max())
    ax.plot([0, e_max], [0, e_max], "k--", lw=0.8, label="E_dep = E₀", zorder=0)
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("Total LS edep (MeV)", fontsize=10)
    ax.set_title("LS calorimeter response", fontsize=11)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel B: linearity
    ax = axes[1]
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst  = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en   = summary["energy_MeV"].values
        lin  = summary["ls_total_mean"].values / en
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, lin, ls=pst["ls"], lw=pst["lw"], color=color,
                label=pst["label"])
    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_ylim(0, 1.15)
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("E_dep(LS) / E₀", fontsize=10)
    ax.set_title("Calorimeter linearity", fontsize=11)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel C: energy resolution
    ax = axes[2]
    for particle, summary in summaries.items():
        if "ls_res_std" not in summary.columns:
            continue
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        en    = summary["energy_MeV"].values
        res   = summary["ls_res_std"].values / summary["ls_res_mean"].values
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, res * 100, ls=pst["ls"], lw=pst["lw"], color=color,
                label=pst["label"])
    _add_energy_lines(ax)
    ax.set_xlabel("Initial energy (MeV)", fontsize=10)
    ax.set_ylabel("σ/μ  (%)", fontsize=10)
    ax.set_title("LS energy resolution\n(events reaching LS1)", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    fig.suptitle("Liquid scintillator calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def plot_containment(pdf: PdfPages, summaries: dict):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    for particle, summary in summaries.items():
        pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        for col, color in [("trans_primInLS1", "#6a1d8a"),
                            ("trans_primInLS4", "#d7bde2")]:
            if col in summary.columns:
                layer = "LS1" if "LS1" in col else "LS4"
                ax.plot(summary["energy_MeV"].values, summary[col].values,
                        color=color, ls=pst["ls"], lw=pst["lw"],
                        label=f"{layer} ({pst['label']})")
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel("Fraction of primaries", fontsize=10)
    _set_xaxis(ax)
    ax.set_title("LS layer entry fractions", fontsize=11)
    ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for particle, summary in summaries.items():
        if "ls_total_mean" not in summary.columns:
            continue
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8, "label": particle})
        en    = summary["energy_MeV"].values
        color = "#9467bd" if particle == "electron" else "#e377c2"
        ax.plot(en, summary["trans_primInLS1"].values if "trans_primInLS1" in summary.columns else np.nan,
                color=color, ls=pst["ls"], lw=pst["lw"], label=f"Enters LS1 ({pst['label']})")
        ax.plot(en, summary["ls_total_mean"].values / en,
                color=color, ls=":" , lw=1.2, label=f"LS edep/E₀ ({pst['label']})")
    ax.axhline(1.0, color="grey", lw=0.7, ls=":")
    ax.set_ylim(-0.05, 1.2)
    ax.set_ylabel("Fraction / ratio", fontsize=10)
    _set_xaxis(ax)
    ax.set_title("Containment quality", fontsize=11)
    ax.legend(fontsize=8, ncol=2); ax.grid(True, alpha=0.3)

    fig.suptitle("Containment analysis for LS calorimetry", fontsize=13, y=1.01)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def plot_layer_profiles(pdf: PdfPages, all_groups: dict,
                        snapshot_energies: list):
    """Bar chart of edep per layer at snapshot energies, one panel per particle."""
    n_particles = len(all_groups)
    if n_particles == 0:
        return
    fig, axes = plt.subplots(1, n_particles,
                              figsize=(11 * n_particles / 2 + 3, 5),
                              sharey=True)
    if n_particles == 1:
        axes = [axes]

    layer_labels = [label for (_, _, label, _) in LAYERS]
    x = np.arange(len(LAYERS))
    colors_e = cm.viridis(np.linspace(0.2, 0.85, len(snapshot_energies)))

    for ax, (particle, groups) in zip(axes, all_groups.items()):
        snaps = find_snapshots(groups, snapshot_energies)
        bar_width = 0.7 / len(snaps)
        for i, (e_target, (e_actual, files)) in enumerate(sorted(snaps.items())):
            df = read_event_tree(files)
            if df.empty:
                continue
            vals = []
            for (edep_br, _, _, _) in LAYERS:
                vals.append(df[edep_br].mean() / 1e6 if edep_br in df.columns else 0.0)
            offset = (i - len(snaps) / 2 + 0.5) * bar_width
            ax.bar(x + offset, vals, bar_width * 0.9,
                   label=f"{e_actual:.3g} MeV", color=colors_e[i], alpha=0.85)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(layer_labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("Mean edep (MeV/primary)", fontsize=10)
        ax.set_title(f"Layer profile — {PARTICLE_STYLES.get(particle,{}).get('label', particle)}",
                     fontsize=11)
        ax.legend(title="E₀", fontsize=8, title_fontsize=8)
        ax.grid(True, axis="y", which="both", alpha=0.3)

    fig.suptitle("Energy deposition profile at selected energies", fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_angular_distributions(pdf: PdfPages, all_groups: dict,
                                snapshot_energies: list):
    """One row of histograms per particle."""
    n_particles = len(all_groups)
    n_snaps     = len(snapshot_energies)
    colors      = cm.plasma(np.linspace(0.15, 0.85, n_snaps))

    fig, axes = plt.subplots(n_particles, n_snaps,
                              figsize=(3.0 * n_snaps, 3.5 * n_particles),
                              sharey="row", sharex=False)
    if n_particles == 1:
        axes = axes[np.newaxis, :]
    if n_snaps == 1:
        axes = axes[:, np.newaxis]

    for row, (particle, groups) in enumerate(all_groups.items()):
        snaps = find_snapshots(groups, snapshot_energies)
        for col, ((e_target, (e_actual, files)), color) in enumerate(
                zip(sorted(snaps.items()), colors)):
            ax = axes[row, col]
            print(f"  Angular: {particle} {e_actual:.3g} MeV ...")
            cdf    = read_cluster_tree(files)
            angles = compute_angles(cdf)
            if len(angles) == 0:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                continue
            rms    = np.sqrt(np.mean(angles**2))
            ax.hist(angles, bins=60, density=True, color=color, alpha=0.8)
            ax.axvline(rms, color="k", lw=1.1, ls="--",
                       label=f"RMS={rms:.3f}°")
            ax.set_title(f"{e_actual:.3g} MeV", fontsize=9)
            ax.legend(fontsize=7)
            ax.set_xlim(left=0)
            if col == 0:
                plab = PARTICLE_STYLES.get(particle, {}).get("label", particle)
                ax.set_ylabel(f"{plab}\nDensity", fontsize=9)
            if row == n_particles - 1:
                ax.set_xlabel("Polar angle (degrees)", fontsize=8)

    fig.suptitle("Track angle in MM drift gas (multiple scattering spread)",
                 fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_angular_resolution(pdf: PdfPages, all_groups: dict, step: int = 10):
    """RMS polar angle vs energy for all particles on one axes."""
    fig, ax = plt.subplots(figsize=(8, 4.5))

    for particle, groups in all_groups.items():
        pst   = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.6, "label": particle})
        color = "#1f77b4" if particle == "electron" else "#ff7f0e"
        energies_sel = sorted(groups.keys())[::step]
        en_vals, rms_vals = [], []
        n = len(energies_sel)
        for idx, energy in enumerate(energies_sel, 1):
            print(f"  Angular RMS {particle} [{idx}/{n}] {energy:.4g} MeV")
            cdf    = read_cluster_tree(groups[energy])
            angles = compute_angles(cdf)
            if len(angles) < 20:
                continue
            rms_vals.append(np.sqrt(np.mean(angles**2)))
            en_vals.append(energy)
        if en_vals:
            ax.plot(en_vals, rms_vals, ls=pst["ls"], lw=pst["lw"],
                    color=color, marker="o", ms=3, label=pst["label"])

    ax.set_ylabel("RMS polar angle (degrees)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Angular resolution from MS in upstream material\n"
                 "(straight-line fit in MM drift gas)", fontsize=11)
    ax.legend(fontsize=9); ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Coarse grouped plots ──────────────────────────────────────────────────────

def _coarse_legend(ax, loc="lower right"):
    handles = [mpatches.Patch(color=color, label=label)
               for (label, color, _, _) in COARSE_GROUPS]
    return ax.legend(handles=handles, loc=loc, fontsize=8,
                     title="Group", title_fontsize=8)


def plot_transmission_coarse(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    for (label, color, _edep_brs, trans_col) in COARSE_GROUPS:
        col = f"trans_{trans_col}"
        for particle, summary in summaries.items():
            if col not in summary.columns:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
            ax.plot(summary["energy_MeV"].values, summary[col].values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_ylim(-0.02, 1.05)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("Transmission fraction", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Transmission by detector group", fontsize=12)
    ax.grid(True, alpha=0.3)
    leg1 = _coarse_legend(ax, loc="lower right")
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc="center right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def plot_edep_coarse(pdf: PdfPages, summaries: dict):
    fig, ax = plt.subplots(figsize=(8, 5))
    for (label, color, edep_brs, _trans_col) in COARSE_GROUPS:
        for particle, summary in summaries.items():
            cols = [f"edep_{b}" for b in edep_brs if f"edep_{b}" in summary.columns]
            if not cols:
                continue
            total = sum(summary[c] for c in cols)
            if total.max() < 1e-8:
                continue
            pst = PARTICLE_STYLES.get(particle, {"ls": "-", "lw": 1.8})
            ax.plot(summary["energy_MeV"].values, total.values,
                    color=color, ls=pst["ls"], lw=pst["lw"])
    ax.set_yscale("log")
    ax.set_ylabel("Mean energy deposition (MeV / primary)", fontsize=11)
    _set_xaxis(ax)
    ax.set_title("Energy deposition by detector group", fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    leg1 = _coarse_legend(ax, loc="upper left")
    ax.add_artist(leg1)
    _particle_legend(ax, summaries, loc="lower right")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ── Survival vs path-length plots ─────────────────────────────────────────────

def _draw_stack_diagram(ax, xlim):
    """
    Thin schematic of the detector stack sharing the path-length x-axis.
    Air gaps take the figure background colour (drawn as nothing).
    Very thin layers get a minimum display width so they remain visible.
    """
    MIN_VIS_MM = 1.8   # minimum drawn width [mm] for thin layers

    ax.set_facecolor(ax.get_figure().get_facecolor())

    for (z0, z1, fc, tc, label) in DETECTOR_STACK:
        w_real = z1 - z0
        w_draw = max(w_real, MIN_VIS_MM)
        zc     = (z0 + z1) / 2                     # keep block centred on true position
        rect   = mpatches.Rectangle(
            (zc - w_draw / 2, 0.04), w_draw, 0.92,
            facecolor=fc, edgecolor="white", linewidth=0.5, zorder=2)
        ax.add_patch(rect)

        if label and w_draw >= 3.0:
            rot = 0 if w_draw >= 14 else 90
            fs  = 7.0 if w_draw >= 14 else 5.5
            ax.text(zc, 0.5, label,
                    ha="center", va="center",
                    fontsize=fs, color=tc, rotation=rot,
                    fontweight="bold", zorder=3, clip_on=True)

    # Beam-direction arrow at the left
    ax.annotate("", xy=(xlim[0] + 12, 0.5), xytext=(xlim[0] + 1, 0.5),
                arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0),
                zorder=4)
    ax.text(xlim[0] + 0.5, 0.5, "beam", fontsize=5.5, color="#666666",
            ha="left", va="center", style="italic")

    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1)
    ax.set_autoscale_on(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(left=False, bottom=False,
                   labelleft=False, labelbottom=False)


def plot_survival(pdf: PdfPages, all_summaries: dict):
    """
    Top row: thin schematic of the detector stack (shared x-axis).
    Then one row per particle: survival S(z) and loss rate −dS/dz.
    """
    n_par  = len(all_summaries)
    XLIM   = (0.0, LS4_EXIT_MM + 6)
    DIAG_H = 0.22   # height ratio of diagram row relative to data rows

    fig = plt.figure(figsize=(12, 1.1 + 4.5 * n_par))
    gs  = fig.add_gridspec(
        n_par + 1, 2,
        height_ratios=[DIAG_H] + [1.0] * n_par,
        hspace=0.06, wspace=0.30,
    )

    ax_diag = fig.add_subplot(gs[0, :])
    _draw_stack_diagram(ax_diag, xlim=XLIM)

    snap_colors = cm.plasma(np.linspace(0.1, 0.85, len(SURVIVAL_ENERGIES_MEV)))

    for row_idx, (particle, summary) in enumerate(all_summaries.items()):
        ax_s = fig.add_subplot(gs[row_idx + 1, 0], sharex=ax_diag)
        ax_l = fig.add_subplot(gs[row_idx + 1, 1], sharex=ax_diag)

        plab     = PARTICLE_STYLES.get(particle, {}).get("label", particle)
        avail_en = summary["energy_MeV"].values

        for e_target, color in zip(SURVIVAL_ENERGIES_MEV, snap_colors):
            idx      = np.argmin(np.abs(avail_en - e_target))
            e_actual = avail_en[idx]
            row_data = summary.iloc[idx]

            z_pts, s_pts = [0.0], [1.0]
            for (trans_col, z_entry, _lbl) in SURVIVAL_COLS:
                col = f"trans_{trans_col}"
                if col in summary.columns:
                    val = row_data[col]
                    if not np.isnan(val):
                        z_pts.append(z_entry)
                        s_pts.append(float(val))

            z_pts = np.array(z_pts)
            s_pts = np.array(s_pts)
            lbl   = f"{e_actual:.3g} MeV"

            ax_s.plot(z_pts, s_pts, color=color, lw=1.6,
                      marker="o", ms=4, label=lbl)

            dz   = np.diff(z_pts)
            ds   = np.diff(s_pts)
            loss = -ds / dz
            z_mid = 0.5 * (z_pts[:-1] + z_pts[1:])
            ax_l.plot(z_mid, loss, color=color, lw=1.5,
                      marker="o", ms=4, label=lbl)

        # Shared boundary guides on both panels
        for (_, z_entry, _) in SURVIVAL_COLS:
            ax_s.axvline(z_entry, color="lightgrey", lw=0.7, ls=":", zorder=0)
            ax_l.axvline(z_entry, color="lightgrey", lw=0.7, ls=":", zorder=0)
        ax_s.axvline(LS4_EXIT_MM, color="lightgrey", lw=0.7, ls=":", zorder=0)
        ax_l.axvline(LS4_EXIT_MM, color="lightgrey", lw=0.7, ls=":", zorder=0)

        ax_s.set_ylabel("Survival fraction", fontsize=10)
        ax_s.set_title(f"{plab} — Survival", fontsize=11)
        ax_s.set_ylim(-0.02, 1.08)
        ax_s.axhline(1.0, color="grey", lw=0.6, ls=":")
        ax_s.legend(title="E₀", fontsize=8, title_fontsize=8, loc="upper right")
        ax_s.grid(True, alpha=0.3)

        ax_l.set_ylabel("Loss rate  −dS/dz  (mm⁻¹)", fontsize=10)
        ax_l.set_title(f"{plab} — Loss rate", fontsize=11)
        ax_l.set_ylim(bottom=0)
        ax_l.legend(title="E₀", fontsize=8, title_fontsize=8)
        ax_l.grid(True, alpha=0.3)

        # x labels only on the bottom data row
        is_bottom = (row_idx == n_par - 1)
        if is_bottom:
            ax_s.set_xlabel("Path from gun (mm)", fontsize=10)
            ax_l.set_xlabel("Path from gun (mm)", fontsize=10)
        else:
            plt.setp(ax_s.get_xticklabels(), visible=False)
            plt.setp(ax_l.get_xticklabels(), visible=False)

    ax_diag.set_xlim(*XLIM)

    fig.suptitle("Primary particle survival vs cumulative path through the stack",
                 fontsize=12)
    fig.savefig(pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Full-experiment stack analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--indir",    required=True)
    p.add_argument("--gas",      default="ArIso")
    p.add_argument("--particles", nargs="+", default=["electron", "positron"],
                   help="Particle types to load (space-separated)")
    p.add_argument("--outfile",  default="full_experiment_analysis.pdf")
    p.add_argument("--angular-step", type=int, default=10,
                   help="Process every Nth energy for angular RMS plot")
    p.add_argument("--no-angular", action="store_true",
                   help="Skip angular analysis (faster)")
    return p.parse_args()


def main():
    args = parse_args()

    plt.rcParams.update({
        "figure.dpi": 150, "axes.titlesize": 11,
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9,
    })

    all_groups    = {}
    all_summaries = {}

    for particle in args.particles:
        print(f"\nLoading {particle} ...")
        groups = group_files(args.indir, args.gas, particle)
        if not groups:
            continue
        print(f"  {len(groups)} energy points: "
              f"{min(groups):.4g} – {max(groups):.4g} MeV")
        all_groups[particle] = groups
        print(f"  Computing summaries ...")
        all_summaries[particle] = compute_summary(groups)

    if not all_summaries:
        sys.exit("No data loaded.")

    # Save summary CSVs
    for particle, summary in all_summaries.items():
        csv = Path(args.outfile).with_name(
            Path(args.outfile).stem + f"_{particle}.csv")
        summary.to_csv(csv, index=False)
        print(f"Summary saved → {csv}")

    e_range = {p: f"{min(g):.4g}–{max(g):.4g} MeV"
               for p, g in all_groups.items()}

    print(f"\nWriting {args.outfile} ...")
    with PdfPages(args.outfile) as pdf:

        # Cover
        fig = plt.figure(figsize=(8, 4))
        fig.text(0.5, 0.68, "Full-experiment stack analysis",
                 ha="center", fontsize=18, fontweight="bold")
        fig.text(0.5, 0.52,
                 f"Gas: {args.gas}\n" +
                 "\n".join(f"{p}: {r}" for p, r in e_range.items()) +
                 f"\nInput: {args.indir}",
                 ha="center", fontsize=10, color="grey")
        pdf.savefig(fig); plt.close(fig)

        print("  Transmission (all layers) ...")
        plot_transmission(pdf, all_summaries)
        print("  Transmission (simplified: drift / plastic scint. / LS1 / LS4) ...")
        plot_transmission_simple(pdf, all_summaries)
        print("  Transmission (coarse: detector groups) ...")
        plot_transmission_coarse(pdf, all_summaries)

        print("  Edep (all layers) ...")
        plot_edep(pdf, all_summaries)
        print("  Edep (simplified) ...")
        plot_edep_simple(pdf, all_summaries)
        print("  Edep (coarse: detector groups) ...")
        plot_edep_coarse(pdf, all_summaries)

        print("  Calorimetry ...")
        plot_calorimetry(pdf, all_summaries)

        print("  Containment ...")
        plot_containment(pdf, all_summaries)

        print("  Survival vs path length ...")
        plot_survival(pdf, all_summaries)

        print("  Layer profiles ...")
        plot_layer_profiles(pdf, all_groups, SNAPSHOT_ENERGIES)

        if not args.no_angular:
            print("  Angular distributions ...")
            plot_angular_distributions(pdf, all_groups, SNAPSHOT_ENERGIES)
            print(f"  Angular RMS (every {args.angular_step}th point) ...")
            plot_angular_resolution(pdf, all_groups, step=args.angular_step)
        else:
            print("  Skipping angular analysis (--no-angular)")

    print(f"\nDone → {args.outfile}")


if __name__ == "__main__":
    main()
