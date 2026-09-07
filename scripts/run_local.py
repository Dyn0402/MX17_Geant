#!/usr/bin/env python3
"""
run_local.py
Run mm_sim locally in parallel using multiple threads.
Each thread launches a subprocess for one energy point; no GIL issue.

Usage:
    python scripts/run_local.py --mode sr90 --outdir /tmp/sr90_out --workers 8
    python scripts/run_local.py --mode lscalib --outdir /tmp/lscalib_out --nevents 50000
"""

import argparse
import multiprocessing
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# Energy grids
ENERGIES_SR90     = sorted(set(
    [round(e * 0.05, 3) for e in range(2, 21)] +   # 0.1–1.0 MeV step 0.05
    [round(e * 0.1,  2) for e in range(11, 31)]     # 1.1–3.0 MeV step 0.1
))
ENERGIES_LSCALIB  = sorted(set(
    [round(e * 0.1, 2) for e in range(1, 11)]  +
    [round(1.0 + e * 0.1, 2) for e in range(1, 14)] +
    [2.28, 2.5, 3.0]
))
ENERGIES_FULL     = sorted(set(
    [round(e * 0.5, 1) for e in range(1, 9)]   +   # 0.5–4.0 MeV
    [round(4.0 + e * 1.0, 1) for e in range(1, 14)] +  # 5–17 MeV
    [16.5, 18.0]
))


def energy_grid(mode):
    if mode in ("sr90", "sr90nomm"):
        return ENERGIES_SR90
    if mode == "lscalib":
        return ENERGIES_LSCALIB
    return ENERGIES_FULL


def run_one(exe, mode, particle, energy_mev, nevents, gas, outdir, seed, timeout):
    tag     = f"{mode}_{particle}_e{energy_mev:.3f}"
    outfile = os.path.join(outdir, tag)
    cmd = [
        exe,
        "-m", mode,
        "-p", particle,
        "-e", str(energy_mev),
        "-n", str(nevents),
        "-g", gas,
        "-o", outfile,
        "-s", str(seed),
        "-t", "1",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            return (energy_mev, False, result.stderr[-500:])
        return (energy_mev, True, "")
    except subprocess.TimeoutExpired:
        return (energy_mev, False, "TIMEOUT")
    except Exception as e:
        return (energy_mev, False, str(e))


def main():
    default_exe = os.path.join(os.path.dirname(__file__), "..", "build", "mm_sim")

    ap = argparse.ArgumentParser()
    ap.add_argument("--exe",       default=default_exe, help="Path to mm_sim binary")
    ap.add_argument("--mode",      default="lscalib",
                    choices=["vacuum","full","sr90","sr90nomm","lscalib"])
    ap.add_argument("--outdir",    required=True,       help="Output directory")
    ap.add_argument("--particle",  default="electron")
    ap.add_argument("--gas",       default="ArCF4")
    ap.add_argument("--nevents",   type=int, default=10000)
    ap.add_argument("--workers",   type=int,
                    default=max(1, multiprocessing.cpu_count()-1))
    ap.add_argument("--max-jobs",  type=int, default=0,  help="Limit number of jobs (0=all)")
    ap.add_argument("--seed-base", type=int, default=12345)
    ap.add_argument("--timeout",   type=int, default=3600, help="Per-job timeout [s]")
    ap.add_argument("--dry-run",   action="store_true")
    args = ap.parse_args()

    exe = os.path.abspath(args.exe)
    if not os.path.isfile(exe):
        print(f"ERROR: executable not found: {exe}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    energies = energy_grid(args.mode)
    if args.max_jobs > 0:
        energies = energies[:args.max_jobs]

    print(f"Mode: {args.mode}  |  Particle: {args.particle}  |  "
          f"Gas: {args.gas}  |  Events: {args.nevents}")
    print(f"Energy grid: {len(energies)} points  "
          f"({energies[0]:.3f}–{energies[-1]:.3f} MeV)")
    print(f"Workers: {args.workers}  |  Output: {args.outdir}")

    if args.dry_run:
        for i, E in enumerate(energies):
            print(f"  [dry-run] E={E:.3f} MeV  seed={args.seed_base+i}")
        return

    done = failed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(run_one, exe, args.mode, args.particle, E,
                      args.nevents, args.gas, args.outdir,
                      args.seed_base + i, args.timeout): E
            for i, E in enumerate(energies)
        }
        for fut in as_completed(futs):
            E = futs[fut]
            energy, ok, msg = fut.result()
            if ok:
                done += 1
                print(f"  OK  E={energy:.3f} MeV  ({done}/{len(energies)})")
            else:
                failed += 1
                print(f"  FAIL  E={energy:.3f} MeV: {msg}")

    print(f"\nDone: {done} OK, {failed} failed.")


if __name__ == "__main__":
    main()
