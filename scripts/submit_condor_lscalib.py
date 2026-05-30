#!/usr/bin/env python3
"""
submit_condor_lscalib.py
Submit kLSCalib or kBackScintCalib jobs to HTCondor on lxplus.

Each job runs N events with the Sr-90/Y-90 beta spectrum sampled at the gun
(--spectrum path to the CSV file).  Multiple parallel jobs are submitted to
accumulate statistics; merge with hadd afterwards.

Usage:
    python scripts/submit_condor_lscalib.py \\
        --exe /eos/home-d/dneff/sim/mm_sim \\
        --spectrum /eos/home-d/dneff/sim/Sr90_Y90_Beta_Spectrum.csv \\
        --outdir /eos/home-d/dneff/sim/lscalib_out \\
        --mode ls \\
        --njobs 20 --nevents 100000
"""

import argparse
import os
import subprocess
import sys


def write_wrapper(jobdir, exe, mode, spectrum, outfile, nevents, seed, src_dist):
    name    = os.path.basename(outfile)
    wrapper = os.path.join(jobdir, "wrappers", f"{name}.sh")
    os.makedirs(os.path.dirname(wrapper), exist_ok=True)

    g4_mode = "lscalib" if mode == "ls" else "backscintcalib"

    with open(wrapper, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh\n")
        f.write(f"{exe} \\\n")
        f.write(f"  -m {g4_mode} \\\n")
        f.write(f"  -p electron \\\n")
        f.write(f"  -n {nevents} \\\n")
        f.write(f"  -g ArCF4 \\\n")   # gas not used in calib modes but required
        f.write(f"  -o {outfile} \\\n")
        f.write(f"  -s {seed} \\\n")
        f.write(f"  --spectrum {spectrum} \\\n")
        f.write(f"  --src-dist {src_dist:.1f}\n")
    os.chmod(wrapper, 0o755)
    return wrapper


def write_jdl(jobdir, wrappers, log_dir):
    jdl = os.path.join(jobdir, "lscalib_jobs.jdl")
    with open(jdl, "w") as f:
        f.write("universe = vanilla\n")
        f.write("+JobFlavour = \"longlunch\"\n")
        f.write("request_memory = 512MB\n\n")
        for w in wrappers:
            tag = os.path.splitext(os.path.basename(w))[0]
            f.write(f"executable = {w}\n")
            f.write(f"log    = {log_dir}/{tag}.log\n")
            f.write(f"output = {log_dir}/{tag}.out\n")
            f.write(f"error  = {log_dir}/{tag}.err\n")
            f.write("queue 1\n\n")
    return jdl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe",       required=True,  help="Path to mm_sim on EOS")
    ap.add_argument("--spectrum",  required=True,  help="Path to Sr90_Y90 CSV on EOS")
    ap.add_argument("--outdir",    required=True,  help="Output directory on EOS")
    ap.add_argument("--mode",      choices=["ls","backscint"], default="ls",
                    help="'ls'=liquid scint only,  'backscint'=back plastic scint only")
    ap.add_argument("--jobdir",    default=None,   help="Job scripts dir (default: outdir/jobs)")
    ap.add_argument("--njobs",     type=int, default=10, help="Number of parallel jobs")
    ap.add_argument("--nevents",   type=int, default=100000, help="Events per job")
    ap.add_argument("--seed-base", type=int, default=42000)
    ap.add_argument("--src-dist",  type=float, default=100.0,
                    help="Source-to-detector air gap [mm]")
    ap.add_argument("--dry-run",   action="store_true")
    args = ap.parse_args()

    jobdir  = args.jobdir or os.path.join(args.outdir, "jobs")
    log_dir = os.path.join(jobdir, "logs")
    for d in [args.outdir, jobdir, log_dir]:
        os.makedirs(d, exist_ok=True)

    det_tag = "ls" if args.mode == "ls" else "backscint"
    wrappers = []
    for i in range(args.njobs):
        seed    = args.seed_base + i
        outfile = os.path.join(args.outdir, f"{det_tag}_calib_job{i:03d}")
        w = write_wrapper(jobdir, args.exe, args.mode, args.spectrum,
                          outfile, args.nevents, seed, args.src_dist)
        wrappers.append(w)

    jdl = write_jdl(jobdir, wrappers, log_dir)

    total = args.njobs * args.nevents
    print(f"Mode     : {args.mode} calibration")
    print(f"Jobs     : {args.njobs}  ×  {args.nevents:,} events  =  {total:,} total")
    print(f"Output   : {args.outdir}")
    print(f"Spectrum : {args.spectrum}")
    print(f"JDL      : {jdl}")
    print()
    print("After jobs complete, merge with:")
    print(f"  hadd {args.outdir}/{det_tag}_calib_all.root {args.outdir}/{det_tag}_calib_job*.root")
    print("or for CSV:")
    print(f"  head -1 {args.outdir}/{det_tag}_calib_job000_events.csv > merged.csv")
    print(f"  tail -n +2 -q {args.outdir}/{det_tag}_calib_job*_events.csv >> merged.csv")

    if not args.dry_run:
        r = subprocess.run(["condor_submit", jdl], capture_output=True, text=True)
        print(r.stdout)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            sys.exit(r.returncode)
    else:
        print("\nDry run — not submitting.")


if __name__ == "__main__":
    main()
