#!/usr/bin/env python3
"""
collect.py — merge the S3 seed slices into aval_calib.json and plot them.

The campaign (mx17_aval.sub) splits every voltage point across 8 independent
seed slices so it parallelises. This merges them back, fits the Polya, and
produces the figures.

    # pull the results back from lxplus first
    rsync -av lxplus:/afs/cern.ch/user/d/dneff/work/mx17_response/avalanche/results/ \\
              ~/x17/response_sim/avalanche/

    python3 -m response.avalanche.collect ~/x17/response_sim/avalanche
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np


def polya_fit(gains):
    """
    Fit P(g) proportional to (g/gbar)^theta exp(-(1+theta) g/gbar).

    theta is recovered from the relative variance, which for a Polya is
    exactly  var/mean^2 = 1/(1+theta).  That is a moment estimator, not a
    likelihood fit — it is unbiased enough for a calibration constant and it
    cannot fail to converge, which matters in an unattended pipeline.
    """
    g = np.asarray(gains, dtype=float)
    g = g[g > 0]
    if len(g) < 10:
        return None
    mean = g.mean()
    rel_var = g.var(ddof=1) / mean ** 2
    theta = 1.0 / rel_var - 1.0
    return {"gain_mean": float(mean), "theta": float(theta),
            "rel_var": float(rel_var), "n": int(len(g)),
            "gain_median": float(np.median(g))}


Z_BINS, Z_RANGE = 60, (0.0, 150.0)


def _moments(seq):
    """(N, sum, sum of squares) of a list, via numpy so it is one pass."""
    a = np.asarray(seq, dtype=float)
    return int(a.size), float(a.sum()), float(np.square(a).sum())


def reduce_file(path):
    """
    Parse ONE result file and immediately throw away everything that is not a
    sufficient statistic.

    This matters more than it looks. `r_end_um`, `t_end_ns` and `z_ion_um` hold
    one entry per electron/ion — 4.5 MILLION each in a 200-avalanche slice, and
    that is what makes each file 264 MB and the campaign 19 GB. The original
    version json.load()ed every file and kept them all, which reached 6.5 GB
    resident after 15 of 56 files and could not have finished on lxplus.

    Everything downstream needs from those arrays is moments and one histogram,
    all of which are exactly accumulable, so the merged numbers are unchanged —
    this is a reorganisation, not an approximation. `gains` is kept whole: it
    is one entry per avalanche (200), not per electron.
    """
    d = json.load(open(path))
    r, c = d["results"], d["config"]
    n_r, s_r, s_r2 = _moments(r["r_end_um"])
    n_t, s_t, s_t2 = _moments(r["t_end_ns"])
    hz, edges = np.histogram(np.asarray(r["z_ion_um"], dtype=float),
                             bins=Z_BINS, range=Z_RANGE)
    return {
        "gas": c["gas_file"], "volt": c["voltage_V"], "nev": c["nev"],
        "gains": list(r["gains"]),
        "r": (n_r, s_r, s_r2), "t": (n_t, s_t, s_t2),
        "zhist": hz, "zedges": edges,
        "i_elec": np.asarray(r["i_elec"], dtype=float),
        "i_ion": np.asarray(r["i_ion"], dtype=float),
        "signal_dt_ns": r["signal_dt_ns"],
        "field_model": d["provenance"]["field_model"],
    }


def load(indir):
    """Group every result file by (gas, voltage), reduced on the way in."""
    points = defaultdict(list)
    files = sorted(glob.glob(os.path.join(indir, "aval_*.json")))
    for i, f in enumerate(files, 1):
        try:
            s = reduce_file(f)
        except Exception as e:                       # partial transfer
            print(f"  skipping unreadable {os.path.basename(f)}: {e}")
            continue
        points[(s["gas"], s["volt"])].append(s)
        print(f"  [{i}/{len(files)}] {os.path.basename(f)}", flush=True)
    return points


def _pooled(triples):
    """Pooled (mean, std) from a list of (N, sum, sumsq)."""
    n = sum(t[0] for t in triples)
    if not n:
        return None, None
    s = sum(t[1] for t in triples)
    s2 = sum(t[2] for t in triples)
    mean = s / n
    var = max(s2 / n - mean ** 2, 0.0)
    return mean, np.sqrt(var)


def merge(slices):
    """Merge seed slices of one (gas, voltage) point."""
    gains = [g for s in slices for g in s["gains"]]
    nev = sum(s["nev"] for s in slices)
    # Signals are already per-avalanche averages; combine weighting by nev.
    i_el = sum(s["i_elec"] * s["nev"] for s in slices)
    i_ion = sum(s["i_ion"] * s["nev"] for s in slices)

    r_mean, r_std = _pooled([s["r"] for s in slices])
    t_mean, t_std = _pooled([s["t"] for s in slices])
    n_r = sum(s["r"][0] for s in slices)
    s_r2 = sum(s["r"][2] for s in slices)
    zh = sum(s["zhist"] for s in slices)

    return {
        "n_slices": len(slices), "nev_total": nev,
        "polya": polya_fit(gains),
        # sigma0 is the RMS radius sqrt(<r^2>), not the spread about the mean
        "sigma0_um": float(np.sqrt(s_r2 / n_r)) if n_r else None,
        "sigma0_rms_um": float(r_std) if r_std is not None else None,
        "t_arrival_mean_ns": float(t_mean) if t_mean is not None else None,
        "t_arrival_rms_ns": float(t_std) if t_std is not None else None,
        "signal_dt_ns": slices[0]["signal_dt_ns"],
        "i_elec": (i_el / nev).tolist(),
        "i_ion": (i_ion / nev).tolist(),
        "alpha_z_hist": {"counts": zh.tolist(),
                         "edges": slices[0]["zedges"].tolist()},
        "field_model": slices[0]["field_model"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--figdir", default="design/figures/response")
    a = ap.parse_args()

    points = load(a.indir)
    if not points:
        print(f"no aval_*.json under {a.indir} — has the campaign landed yet?")
        return 1
    print(f"{len(points)} (gas, voltage) points from "
          f"{sum(len(v) for v in points.values())} slice files")

    calib, rows = {}, []
    for (gas, volt), sl in sorted(points.items()):
        m = merge(sl)
        key = f"{gas}@{volt:.0f}V"
        rows.append((volt, m))
        p = m["polya"]
        print(f"  {key:<44s} nev={m['nev_total']:5d}  "
              f"gain={p['gain_mean']:9.1f}  theta={p['theta']:5.2f}  "
              f"sigma0={m['sigma0_um']:.1f} µm")
        calib[key] = m

    out = a.out or os.path.join(a.indir, "aval_calib.json")
    json.dump({"schema": "aval_calib/1", "points": calib}, open(out, "w"))
    print(f"wrote {out}")

    # ── figures ──────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 0
    os.makedirs(a.figdir, exist_ok=True)
    rows.sort()
    V = [r[0] for r in rows]
    G = [r[1]["polya"]["gain_mean"] for r in rows]
    TH = [r[1]["polya"]["theta"] for r in rows]
    S0 = [r[1]["sigma0_um"] for r in rows]

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    ax[0].semilogy(V, G, "o-", color="#2a78d6", ms=6)
    ax[0].set_xlabel("mesh voltage [V]"); ax[0].set_ylabel("mean gain")
    ax[0].set_title("Gain vs voltage (150 µm gap)")
    ax[1].plot(V, TH, "s-", color="#eb6834", ms=6)
    ax[1].set_xlabel("mesh voltage [V]"); ax[1].set_ylabel("Polya θ")
    ax[1].set_title("Polya shape parameter")
    ax[2].plot(V, S0, "^-", color="#1baf7a", ms=6)
    ax[2].set_xlabel("mesh voltage [V]")
    ax[2].set_ylabel("transverse avalanche σ₀ at the ESL [µm]")
    ax[2].set_title("Avalanche footprint")
    for x in ax:
        x.grid(True, color="#e6e6e2")
        x.spines["top"].set_visible(False); x.spines["right"].set_visible(False)
    fig.suptitle("S3 — avalanche calibration, Ar/iC₄H₁₀ 95/5, uniform-field "
                 "first pass", y=1.03)
    p = os.path.join(a.figdir, "s3_avalanche_calib.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
