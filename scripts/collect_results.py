#!/usr/bin/env python3
"""
collect_results.py
After all Condor jobs finish:
  1. Runs `hadd` to merge per-thread ROOT files into one per job tag.
  2. Reads the merged EventTree from each ROOT file.
  3. Writes a flat summary CSV: one row per (gas, particle, energy) combination
     with mean/median/std of Nprimary_drift and Edep_drift.
  4. Optionally writes a cluster-positions CSV (can be large).

Usage:
    python3 scripts/collect_results.py --indir /path/to/mm_results [options]

Requirements (on lxplus after sourcing setup_lxplus.sh):
    uproot, numpy, pandas  (all in LCG_106)
    hadd (comes with ROOT)
"""

import argparse
import gc
import os
import sys
import re
import subprocess
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    import uproot
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: uproot/numpy/pandas not found.")
    print("Source setup_lxplus.sh first, or: pip install uproot numpy pandas")
    sys.exit(1)

try:
    from tqdm import tqdm
    _TQDM = True
except ImportError:
    _TQDM = False
    print("Tip: install tqdm for progress bars:  pip install tqdm")


def _tqdm_write(msg):
    """Print without disrupting an active tqdm bar."""
    if _TQDM:
        tqdm.write(msg)
    else:
        print(msg)


# Regex to parse job tag: gas_particle_energyMeV[_Al{mm}mm]
TAG_RE = re.compile(
    r"^(?P<gas>ArCF4CO2|ArCF4Iso|ArCF4|HeEth|ArCO2|NeIso|NeCF4|PureCF4)"
    r"_(?P<particle>gamma|electron|neutron|proton|muon)"
    r"_(?P<energy>[0-9ep.+-]+)MeV"
    r"(_Al(?P<al_mm>[0-9p]+)mm)?"
    r"(_t\d+)?\.root$"
)


def parse_tag(fname: str):
    m = TAG_RE.match(fname)
    if not m:
        return None
    energy_str = m.group("energy").replace("p", ".")
    try:
        energy = float(energy_str)
    except ValueError:
        return None
    al_mm_raw = m.group("al_mm")
    al_mm = float(al_mm_raw.replace("p", ".")) if al_mm_raw else 0.0
    return m.group("gas"), m.group("particle"), energy, al_mm


def find_root_files(indir: Path):
    """Group ROOT files by (gas, particle, energy) tag."""
    groups = defaultdict(list)
    for f in sorted(indir.glob("*.root")):
        parsed = parse_tag(f.name)
        if parsed:
            groups[parsed].append(f)
    return groups


def hadd_group(key, files, merged_dir: Path, dry_run=False, quiet=False) -> Path:
    """Merge thread files into one file using hadd."""
    gas, particle, energy, al_mm = key
    e_str  = f"{energy:.6g}".replace(".", "p")
    al_str = f"_Al{al_mm:g}mm".replace(".", "p") if al_mm > 0 else ""
    outname = merged_dir / f"{gas}_{particle}_{e_str}MeV{al_str}_merged.root"

    if outname.exists():
        return outname  # already merged

    if len(files) == 1:
        # Only one thread file -- symlink or copy
        if not dry_run:
            outname.symlink_to(files[0].resolve())
        return outname

    cmd = ["hadd", "-f", str(outname)] + [str(f) for f in sorted(files)]
    if not quiet:
        print(f"  hadd -> {outname.name}  ({len(files)} files)")
    if not dry_run:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            _tqdm_write(f"  WARNING: hadd failed for {outname.name}\n{result.stderr[:500]}")
            return None
    return outname


def read_event_tree(root_file: Path) -> pd.DataFrame:
    """Read EventTree into a pandas DataFrame."""
    try:
        with uproot.open(root_file) as f:
            if "EventTree" not in f:
                print(f"  WARNING: No EventTree in {root_file.name}")
                return pd.DataFrame()
            tree = f["EventTree"]
            df = tree.arrays(
                ["eventID", "edepDrift", "edepAmp",
                 "nPrimDrift", "nPrimAmp",
                 "nClusDrift", "nClusAmp",
                 "primInDrift", "primInAmp"],
                library="pd"
            )
        return df
    except Exception as e:
        print(f"  ERROR reading {root_file.name}: {e}")
        return pd.DataFrame()


def summarize(df: pd.DataFrame, gas: str, particle: str, energy: float,
              al_mm: float = 0.0) -> dict:
    """Compute summary statistics for one (gas, particle, energy, al_mm) point."""
    n = len(df)
    if n == 0:
        return {}

    # Efficiency: fraction of events with any drift ionization
    eff = (df["nPrimDrift"] > 0).mean()

    def stats(col):
        nonzero = df[col][df[col] > 0]
        return {
            "mean":   df[col].mean(),
            "median": df[col].median(),
            "std":    df[col].std(),
            "p10":    df[col].quantile(0.10),
            "p90":    df[col].quantile(0.90),
            "mean_nonzero":  nonzero.mean() if len(nonzero) else 0.0,
        }

    s_np  = stats("nPrimDrift")
    s_ed  = stats("edepDrift")
    s_npa = stats("nPrimAmp")

    row = {
        "gas":       gas,
        "particle":  particle,
        "energy_MeV": energy,
        "al_mm":     al_mm,
        "n_events":  n,
        "efficiency": eff,
        "primInDrift_fraction": float(df["primInDrift"].mean()),
        # Drift
        "nPrimDrift_mean":         s_np["mean"],
        "nPrimDrift_median":       s_np["median"],
        "nPrimDrift_std":          s_np["std"],
        "nPrimDrift_p10":          s_np["p10"],
        "nPrimDrift_p90":          s_np["p90"],
        "nPrimDrift_mean_nonzero": s_np["mean_nonzero"],
        "edepDrift_eV_mean":       s_ed["mean"],
        "edepDrift_eV_median":     s_ed["median"],
        "edepDrift_eV_std":        s_ed["std"],
        # Amp
        "nPrimAmp_mean":   s_npa["mean"],
        "nPrimAmp_median": s_npa["median"],
        "nClusDrift_mean": df["nClusDrift"].mean(),
    }
    return row


def collect_clusters(root_file: Path, gas: str, particle: str,
                     energy: float, max_events: int = 500) -> pd.DataFrame:
    """
    Read ClusterTree (positions) for up to max_events events.
    Used for spatial distribution plots.
    """
    try:
        with uproot.open(root_file) as f:
            if "ClusterTree" not in f:
                return pd.DataFrame()
            tree = f["ClusterTree"]
            total = tree.num_entries
            stop = min(total, max_events * 500)  # approximate cluster count
            df = tree.arrays(
                ["eventID", "x", "y", "z", "edep", "nPrimary",
                 "volume", "particle", "ke"],
                entry_stop=stop,
                library="pd"
            )
        df["gas"]      = gas
        df["beam"]     = particle
        df["energy"]   = energy
        return df
    except Exception as e:
        print(f"  ERROR reading ClusterTree from {root_file.name}: {e}")
        return pd.DataFrame()


def parse_args():
    p = argparse.ArgumentParser(description="Collect and summarise Micromegas sim results")
    p.add_argument("--indir",  required=True,
                   help="Directory containing job output ROOT files")
    p.add_argument("--outdir", default=None,
                   help="Output directory for summary files (default: indir/summary)")
    p.add_argument("--no-hadd",  action="store_true",
                   help="Skip hadd (if files already merged)")
    p.add_argument("--clusters", action="store_true",
                   help="Also write cluster position CSV (can be large)")
    p.add_argument("--dry-run",  action="store_true",
                   help="Print what would be done without doing it")
    p.add_argument("--workers", type=int, default=4,
                   help="Number of parallel workers for hadd and file reading (default: 4)")
    return p.parse_args()


def main():
    args = parse_args()
    indir  = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "summary"
    outdir.mkdir(parents=True, exist_ok=True)
    merged_dir = outdir / "merged"
    merged_dir.mkdir(exist_ok=True)

    print(f"Input dir : {indir}")
    print(f"Output dir: {outdir}")

    # Find files
    groups = find_root_files(indir)
    if not groups:
        print("ERROR: No matching ROOT files found in", indir)
        sys.exit(1)
    print(f"Found {len(groups)} unique (gas, particle, energy) combinations.")

    # hadd
    merged = {}
    if not args.no_hadd:
        print(f"\n--- Merging thread files (workers={args.workers}) ---")
        keys_files = sorted(groups.items())
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            future_to_key = {
                pool.submit(hadd_group, key, files, merged_dir, args.dry_run, _TQDM): key
                for key, files in keys_files
            }
            it = as_completed(future_to_key)
            if _TQDM:
                it = tqdm(it, total=len(future_to_key), desc="hadd", unit="job")
            for fut in it:
                key = future_to_key[fut]
                mf = fut.result()
                if mf:
                    merged[key] = mf
    else:
        for key, files in groups.items():
            # Assume single file or already merged
            merged[key] = files[0]

    if args.dry_run:
        print("\nDry run complete.")
        return

    # Read and summarise in parallel
    print(f"\n--- Computing summaries (workers={args.workers}) ---")
    summary_rows = []
    out_clus = outdir / "clusters_sample.csv"
    if args.clusters and out_clus.exists():
        out_clus.unlink()
    cluster_rows_written = 0

    valid_keys = [k for k in sorted(merged.keys())
                  if merged[k] and merged[k].exists()]

    def _summarize_one(key):
        gas, particle, energy, al_mm = key
        root_file = merged[key]
        df = read_event_tree(root_file)
        if df.empty:
            return key, None, pd.DataFrame()
        row = summarize(df, gas, particle, energy, al_mm)
        del df
        gc.collect()
        cdf = collect_clusters(root_file, gas, particle, energy) if args.clusters else pd.DataFrame()
        return key, row, cdf

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        future_to_key = {pool.submit(_summarize_one, key): key for key in valid_keys}
        it = as_completed(future_to_key)
        if _TQDM:
            it = tqdm(it, total=len(future_to_key), desc="summarise", unit="file")
        for fut in it:
            key = future_to_key[fut]
            gas, particle, energy, al_mm = key
            al_tag = f"  Al={al_mm:g}mm" if al_mm > 0 else ""
            try:
                rkey, row, cdf = fut.result()
                results[rkey] = (row, cdf)
                if not _TQDM:
                    if row:
                        print(f"  {gas:12s} {particle:10s} {energy:10.4g} MeV{al_tag}"
                              f"  n={row['n_events']}  eff={row['efficiency']:.3f}"
                              f"  <Nprim_drift>={row['nPrimDrift_mean']:.1f}")
                    else:
                        print(f"  {gas:12s} {particle:10s} {energy:10.4g} MeV{al_tag}  EMPTY")
            except Exception as exc:
                _tqdm_write(f"  ERROR {gas} {particle} {energy:.4g} MeV{al_tag}: {exc}")

    # Collect results in sorted order
    for key in valid_keys:
        if key not in results:
            continue
        row, cdf = results[key]
        if row:
            summary_rows.append(row)
        if args.clusters and not cdf.empty:
            write_header = cluster_rows_written == 0
            cdf.to_csv(out_clus, mode="a", index=False,
                       header=write_header, float_format="%.5g")
            cluster_rows_written += len(cdf)

    # Write summary CSV
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        # Sort sensibly
        part_order = {"gamma": 0, "electron": 1, "neutron": 2, "proton": 3}
        summary_df["_pord"] = summary_df["particle"].map(part_order).fillna(99)
        summary_df = summary_df.sort_values(["gas", "al_mm", "_pord", "energy_MeV"])
        summary_df = summary_df.drop(columns=["_pord"])

        out_csv = outdir / "summary.csv"
        summary_df.to_csv(out_csv, index=False, float_format="%.5g")
        print(f"\nSummary written: {out_csv}  ({len(summary_df)} rows)")
    else:
        print("WARNING: No summary rows produced.")

    if args.clusters and cluster_rows_written:
        print(f"Cluster sample written: {out_clus}  ({cluster_rows_written} rows)")

    print("\nDone! Next step:")
    print(f"  python3 scripts/plot_results.py --summary {outdir}/summary.csv")


if __name__ == "__main__":
    main()
