#!/usr/bin/env python3
"""
submit_condor_lscalib.py
Submit kLSCalib or kBackScintCalib jobs to HTCondor on lxplus.

HTCondor rules observed here:
  - executable, log, output, error must be on AFS (not /eos directly).
  - The wrapper SCRIPT may reference /eos paths (exe, spectrum, outfile).
  - A single 'queue ... from (...)' avoids the deprecated multi-queue warning.

Usage:
    python scripts/submit_condor_lscalib.py \\
        --exe      /eos/home-d/dneff/sim/mm_sim \\
        --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \\
        --outdir   /eos/home-d/dneff/sim/lscalib_out \\
        --mode ls --njobs 20 --nevents 100000

    python scripts/submit_condor_lscalib.py \\
        --exe      /eos/home-d/dneff/sim/mm_sim \\
        --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \\
        --outdir   /eos/home-d/dneff/sim/backscint_out \\
        --mode backscint --njobs 20 --nevents 100000
"""

import argparse
import os
import stat
import subprocess
import sys
from pathlib import Path


# ── Wrapper script ────────────────────────────────────────────────────────────

def write_wrapper(job_dir: Path, exe: str, mode: str, spectrum: str,
                  outdir: str, nevents: int, src_dist: float) -> Path:
    """Write a single shared wrapper that takes (outfile, seed) as $1 $2."""
    g4_mode = "lscalib" if mode == "ls" else "backscintcalib"
    wrapper = job_dir / f"run_{mode}_calib.sh"

    content = f"""\
#!/bin/bash
set -e
OUTFILE=$1
SEED=$2

source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh

{exe} \\
  -m {g4_mode} \\
  -p electron \\
  -n {nevents} \\
  -g ArCF4 \\
  -o "$OUTFILE" \\
  -s "$SEED" \\
  --spectrum {spectrum} \\
  --src-dist {src_dist:.1f}

echo "Done: $OUTFILE  (seed $SEED)"
"""
    wrapper.write_text(content)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


# ── Submit file ───────────────────────────────────────────────────────────────

def write_submit(job_dir: Path, wrapper: Path, jobs: list,
                 flavour: str = "longlunch") -> Path:
    """
    jobs: list of (outfile, seed, tag)  — all outfile paths are on EOS.
    The submit file itself (and logs) live in job_dir on AFS.
    """
    log_dir = job_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    sub = job_dir / "lscalib.sub"
    lines = [
        f"executable            = {wrapper}",
        f"output                = {log_dir}/$(tag).out",
        f"error                 = {log_dir}/$(tag).err",
        f"log                   = {log_dir}/condor.log",
        f'+JobFlavour           = "{flavour}"',
        "request_cpus          = 1",
        "request_memory        = 512",
        'requirements          = (OpSysAndVer =?= "AlmaLinux9")',
        "",
        "arguments             = $(outfile) $(seed)",
        "",
        "queue outfile,seed,tag from (",
    ]
    for outfile, seed, tag in jobs:
        lines.append(f"  {outfile}, {seed}, {tag}")
    lines.append(")")

    sub.write_text("\n".join(lines) + "\n")
    return sub


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Submit LS or back-scint calibration jobs to HTCondor.")
    ap.add_argument("--exe",       required=True,
                    help="Path to mm_sim executable (EOS or AFS)")
    ap.add_argument("--spectrum",  required=True,
                    help="Path to Sr90_Y90_Beta_Spectrum.csv (EOS or AFS)")
    ap.add_argument("--outdir",    required=True,
                    help="EOS directory for simulation output ROOT/CSV files")
    ap.add_argument("--mode",      choices=["ls", "backscint"], default="ls",
                    help="'ls'=liquid scint only  |  'backscint'=back plastic scint only")
    ap.add_argument("--jobdir",
                    default=os.path.join(os.environ.get("HOME", "."),
                                        "condor", "mx17_lscalib"),
                    help="AFS directory for wrapper, submit file, and logs "
                         "(default: ~/condor/mx17_lscalib)")
    ap.add_argument("--njobs",     type=int, default=10,
                    help="Number of parallel jobs (default: 10)")
    ap.add_argument("--nevents",   type=int, default=100000,
                    help="Events per job (default: 100 000)")
    ap.add_argument("--seed-base", type=int, default=42000)
    ap.add_argument("--src-dist",  type=float, default=100.0,
                    help="Source-to-detector air gap [mm] (default: 100)")
    ap.add_argument("--flavour",   default="longlunch",
                    help="HTCondor job flavour (default: longlunch)")
    ap.add_argument("--dry-run",   action="store_true",
                    help="Print what would be submitted without actually submitting")
    args = ap.parse_args()

    job_dir = Path(args.jobdir)
    outdir  = args.outdir
    det_tag = "ls" if args.mode == "ls" else "backscint"

    # Validate paths
    if not Path(args.exe).is_file():
        print(f"WARNING: executable not found: {args.exe}", file=sys.stderr)
    if not Path(args.spectrum).is_file():
        print(f"WARNING: spectrum file not found: {args.spectrum}", file=sys.stderr)

    if not args.dry_run:
        job_dir.mkdir(parents=True, exist_ok=True)
        Path(outdir).mkdir(parents=True, exist_ok=True)
    else:
        job_dir.mkdir(parents=True, exist_ok=True)  # AFS only — safe locally

    # Single wrapper shared by all jobs; (outfile, seed) passed as arguments
    wrapper = write_wrapper(job_dir, args.exe, args.mode, args.spectrum,
                            outdir, args.nevents, args.src_dist)

    # Build job list
    jobs = []
    for i in range(args.njobs):
        seed    = args.seed_base + i
        tag     = f"{det_tag}_calib_job{i:03d}"
        outfile = os.path.join(outdir, tag)
        jobs.append((outfile, seed, tag))

    sub = write_submit(job_dir, wrapper, jobs, flavour=args.flavour)

    total = args.njobs * args.nevents
    print(f"Mode        : {args.mode} calibration ({args.mode})")
    print(f"Jobs        : {args.njobs}  ×  {args.nevents:,} events  =  {total:,} total")
    print(f"Output (EOS): {outdir}")
    print(f"Jobs (AFS)  : {job_dir}")
    print(f"Submit file : {sub}")
    print()
    print("After jobs complete, merge with:")
    if outdir.startswith("/eos"):
        print(f"  hadd {outdir}/{det_tag}_calib_all.root "
              f"{outdir}/{det_tag}_calib_job*.root")
        print("or for CSV:")
        print(f"  head -1 {outdir}/{det_tag}_calib_job000_events.csv > merged.csv")
        print(f"  tail -n +2 -q {outdir}/{det_tag}_calib_job*_events.csv >> merged.csv")
    print()

    if args.dry_run:
        print("Dry run — not submitting.")
        print(f"Submit file contents:\n{sub.read_text()}")
        return

    r = subprocess.run(["condor_submit", str(sub)], capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(r.returncode)


if __name__ == "__main__":
    main()
