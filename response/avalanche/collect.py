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


def load(indir):
    """Group every result file by (gas, voltage)."""
    points = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(indir, "aval_*.json"))):
        try:
            d = json.load(open(f))
        except Exception as e:                       # partial transfer
            print(f"  skipping unreadable {os.path.basename(f)}: {e}")
            continue
        c = d["config"]
        points[(c["gas_file"], c["voltage_V"])].append(d)
    return points


def merge(slices):
    """Merge seed slices of one (gas, voltage) point."""
    gains, r_end, t_end, z_ion = [], [], [], []
    i_el, i_ion, dt, nev = None, None, None, 0
    for d in slices:
        r = d["results"]
        gains += r["gains"]
        r_end += r["r_end_um"]
        t_end += r["t_end_ns"]
        z_ion += r["z_ion_um"]
        n = d["config"]["nev"]
        nev += n
        # Signals are already per-avalanche averages; combine weighting by nev.
        ie = np.asarray(r["i_elec"]) * n
        ii = np.asarray(r["i_ion"]) * n
        i_el = ie if i_el is None else i_el + ie
        i_ion = ii if i_ion is None else i_ion + ii
        dt = r["signal_dt_ns"]

    out = {
        "n_slices": len(slices), "nev_total": nev,
        "polya": polya_fit(gains),
        "sigma0_um": float(np.sqrt(np.mean(np.square(r_end)))) if r_end else None,
        "sigma0_rms_um": float(np.std(r_end)) if r_end else None,
        "t_arrival_mean_ns": float(np.mean(t_end)) if t_end else None,
        "t_arrival_rms_ns": float(np.std(t_end)) if t_end else None,
        "signal_dt_ns": dt,
        "i_elec": (i_el / nev).tolist() if i_el is not None else None,
        "i_ion": (i_ion / nev).tolist() if i_ion is not None else None,
        "alpha_z_hist": None,
        "field_model": slices[0]["provenance"]["field_model"],
    }
    if z_ion:
        h, e = np.histogram(z_ion, bins=60, range=(0, 150))
        out["alpha_z_hist"] = {"counts": h.tolist(), "edges": e.tolist()}
    out["_gains"] = gains
    out["_r_end"] = r_end
    return out


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
        m.pop("_gains"), m.pop("_r_end")
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
