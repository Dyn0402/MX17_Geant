#!/usr/bin/env python3
# -- coding: utf-8 --
"""
Created on April 23 4:36 PM 2026
Created in PyCharm
Created as MX17_Geant/plot_track_stats.py

@author: Dylan Neff, dylan
"""

import os
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist

# --- Configuration ---
BASE_PATH = "/media/ucla/x17/mx17_geant_sim_results/"
# FILE_PATH = "ArCF4Iso_neutron_10000MeV_t0.root"
# FILE_PATH = "ArCF4Iso_neutron_1MeV_t0.root"
# FILE_PATH = "ArCO2_neutron_1p07MeV_t0.root"
# FILE_PATH = "HeEth_neutron_1p07MeV_t0.root"
SUMMARY_CSV = "summary_stats.csv"
TREE_NAME = "ClusterTree"
VOL_LIMITS = 100.0  # mm
SPAN_THRESHOLD = 3  # mm: tracks smaller than this are "Points"
ENTRY_STOP = 1000000

# --- Updated Geometry Configuration (from G4 Construction) ---
DET_XY_HALF = 200.0   # 40cm / 2
Z_ACTIVE_MIN = -15.07 # Cathode interface
Z_ACTIVE_MAX = 15.11  # Anode interface


def calculate_max_span(df):
    """Diameter of the cluster cloud."""
    if len(df) < 2: return 0.0
    return np.max(pdist(df[['x', 'y', 'z']].values))


def process_file(file_name, show_diagnostic=True):
    """Analyzes a single ROOT file and appends results to CSV."""
    full_path = os.path.join(BASE_PATH, file_name)
    if not os.path.exists(full_path):
        print(f"File not found: {full_path}")
        return

    with uproot.open(full_path) as file:
        df = file[TREE_NAME].arrays(library="pd", entry_stop=ENTRY_STOP)

    df['particle'] = df['particle'].astype(str).astype('category')
    grouped = df.groupby(['eventID', 'trackID'], observed=True)

    # 1. Aggregate Track Data
    tracks = grouped.agg(
        particle=('particle', 'first'),
        total_edep=('edep', 'sum'),
        final_ke=('ke', 'last'),
        last_x=('x', 'last'), last_y=('y', 'last'), last_z=('z', 'last')
    )
    tracks['max_span'] = grouped.apply(calculate_max_span, include_groups=False)

    # 2. Logic: Point vs Track & Escape
    tracks['is_point'] = tracks['max_span'] < SPAN_THRESHOLD
    # 2. Precise Escape Logic
    # We define escape as leaving the gas volume while still having energy
    margin = 0.01  # Small tolerance for float precision

    # Check if the last cluster position hit any of the 6 faces of the gas volume
    escaped_x = tracks['last_x'].abs() >= (DET_XY_HALF - margin)
    escaped_y = tracks['last_y'].abs() >= (DET_XY_HALF - margin)
    escaped_z_front = tracks['last_z'] <= (Z_ACTIVE_MIN + margin)
    escaped_z_back = tracks['last_z'] >= (Z_ACTIVE_MAX - margin)

    tracks['escaped'] = (tracks['final_ke'] > 1e-6) & (escaped_x | escaped_y | escaped_z_front | escaped_z_back)
    tracks['status'] = tracks['escaped'].map({True: 'Escaped', False: 'Contained'})

    # 3. Energy Budget Metrics
    total_e = tracks['total_edep'].sum()
    point_e = tracks[tracks['is_point']]['total_edep'].sum()
    track_e = total_e - point_e
    p_pct = (point_e / total_e * 100) if total_e > 0 else 0
    t_pct = 100 - p_pct

    # 4. Save to CSV
    summary_row = pd.DataFrame([{
        'file_name': file_name,
        'pct_point_energy': p_pct,
        'avg_span': tracks['max_span'].mean(),
        'weighted_avg_span': np.average(tracks['max_span'], weights=tracks['total_edep']) if total_e > 0 else 0,
        'escape_fraction': tracks['escaped'].mean()
    }])
    header = not os.path.exists(SUMMARY_CSV)
    summary_row.to_csv(SUMMARY_CSV, mode='a', index=False, header=header)

    # 5. Diagnostic Plots
    if show_diagnostic:
        sns.set_theme(style="ticks")
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"Analysis: {file_name}", fontsize=16, fontweight='bold')

        # Plot A: Energy Deposition per Particle
        sns.boxplot(data=tracks, x='particle', y='total_edep', hue='particle', legend=False, ax=axes[0, 0],
                    palette='viridis')
        axes[0, 0].set_title("Energy Deposition")

        # Plot B: Escape Statistics (Stacked Bar)
        escape_counts = tracks.groupby(['particle', 'status'], observed=True).size().unstack(fill_value=0)
        escape_counts.plot(kind='bar', stacked=True, ax=axes[0, 1], color=['#2ecc71', '#e74c3c'])
        axes[0, 1].set_title("Containment by Particle Type")
        axes[0, 1].set_ylabel("Track Count")

        # Plot C: Max Span Distribution (Log Scale)
        sns.histplot(data=tracks, x='max_span', hue='particle', element="step", log_scale=True, ax=axes[1, 0])
        axes[1, 0].axvline(SPAN_THRESHOLD, color='red', linestyle='--', label='Point Threshold')
        axes[1, 0].set_title("Spatial Extent (Max Span)")
        axes[1, 0].legend()

        # Plot D: Span vs Edep Scatter
        sns.scatterplot(data=tracks, x='max_span', y='total_edep', hue='particle', alpha=0.4, ax=axes[1, 1])
        axes[1, 1].set_xscale('log')
        axes[1, 1].set_title("Span vs. Energy Deposition")

        # Overlay text for Energy Budget
        stats_text = (f"Energy Budget:\n"
                      f"  Point-like: {p_pct:.1f}%\n"
                      f"  Track-like: {t_pct:.1f}%")
        fig.text(0.5, 0.02, stats_text, ha='center', fontsize=12,
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.5'))

        plt.tight_layout(rect=[0, 0.05, 1, 0.95])
        plt.show()


def plot_comparison_results():
    """Reads CSV and compares all processed files."""
    if not os.path.exists(SUMMARY_CSV):
        print("No CSV found.")
        return

    df = pd.read_csv(SUMMARY_CSV).drop_duplicates(subset=['file_name'], keep='last')

    fig, ax1 = plt.subplots(figsize=(12, 6))
    sns.barplot(data=df, x='file_name', y='pct_point_energy', palette="viridis", ax=ax1)
    ax1.set_title("Multi-File Comparison: Point Energy Fraction", fontsize=14)
    ax1.set_ylabel("% Total Energy deposited as Points (<0.1mm)")
    plt.xticks(rotation=30, ha='right')

    ax2 = ax1.twinx()
    sns.lineplot(data=df, x='file_name', y='escape_fraction', marker='o', color='red', label='Escape Fraction', ax=ax2)
    ax2.set_ylabel("Fraction of Tracks Escaping Volume")
    ax2.grid(False)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example Usage:
    # process_file("ArCO2_neutron_1p07MeV_t0.root", show_diagnostic=True)
    process_file("ArCO2_electron_1MeV_t0.root", show_diagnostic=True)
    # plot_comparison_results()