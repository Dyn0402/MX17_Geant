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

    THE FIT IS CONDITIONAL ON SURVIVAL (`g > 0`), and that conditioning is now
    carried explicitly rather than silently dropped — audit A7. The digitizer
    applies this Polya to every mesh-surviving electron, so if some seeds
    produced no avalanche at all the per-electron charge would be biased high
    by 1/P(g>0). `survival` records P(g>0) so the decomposition is closed:

        E[g]  =  survival * E[g | g > 0]

    MEASURED 2026-08-07 over all 56 raw slices (all 7 voltages, 6400 seed
    electrons): survival is EXACTLY 1.0 everywhere — not one seed failed to
    multiply — so the conditional Polya IS the unconditional one and the
    shipped calib was never biased. The field is kept because "it happens to be
    1 here" and "it is 1 by construction" are different statements, and only
    the first is true.
    """
    g = np.asarray(gains, dtype=float)
    n_all = len(g)
    g = g[g > 0]
    if len(g) < 10:
        return None
    mean = g.mean()
    rel_var = g.var(ddof=1) / mean ** 2
    theta = 1.0 / rel_var - 1.0
    surv = len(g) / n_all if n_all else float("nan")
    return {"gain_mean": float(mean), "theta": float(theta),
            "rel_var": float(rel_var), "n": int(len(g)),
            "gain_median": float(np.median(g)),
            # P(g>0) and its binomial error, and the UNCONDITIONAL mean the
            # digitizer actually needs. The two routes must agree by identity;
            # asserted in merge().
            "survival": float(surv),
            "survival_err": float(np.sqrt(max(surv * (1 - surv), 0.0) / n_all))
            if n_all else float("nan"),
            "n_seeds": int(n_all),
            "gain_mean_uncond": float(surv * mean)}


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
        # Penning is a CAMPAIGN AXIS, not a constant, and it must survive the
        # reduction or the grouping below cannot separate arms that differ
        # only by it. Carried 2026-08-10 after the slope-hunt near-miss; see
        # `penning_tag`.
        "penning_mode": c.get("penning_mode"),
        "penning_rp": c.get("penning_rp"),
        # Amplification gap is a campaign axis as of 2026-08-11 (the 135 vs
        # 150 µm effective-gap scan) and has to survive for the same reason
        # Penning did: without it the two gap arms share (gas, voltage,
        # Penning) exactly and merge into one silently well-formed point.
        # Every pre-existing product ran at 150.0, so an absent key defaults
        # there rather than to None — which keeps old raw directories
        # re-collectable to byte-identical keys.
        "gap_um": float(c.get("gap_um", 150.0)),
        "gains": list(r["gains"]),
        "r": (n_r, s_r, s_r2), "t": (n_t, s_t, s_t2),
        "zhist": hz, "zedges": edges,
        "i_elec": np.asarray(r["i_elec"], dtype=float),
        "i_ion": np.asarray(r["i_ion"], dtype=float),
        "signal_dt_ns": r["signal_dt_ns"],
        "field_model": d["provenance"]["field_model"],
        "seed_z0_um": c.get("seed_z0_um"),
        "campaign_label": d["provenance"].get("campaign_label"),
    }


def penning_tag(mode, rp):
    """Identity of the Penning setting, for grouping and for the point key.

    `auto` and a manual rP are physically different configurations of the SAME
    gas at the SAME voltage, so they are different points.
    """
    if mode == "manual":
        return f"rP{rp:.2f}" if rp is not None else "rPmanual"
    return mode or "penningUNKNOWN"


def arm_tags(s):
    """Every campaign axis beyond (gas, voltage), as ordered key fragments.

    One place to add the next axis. A fragment is emitted for EVERY slice, so
    two arms differing on any axis can never share a group; `main` then decides
    which fragments actually reach the string key (only the ambiguous ones, so
    old products keep their exact key format).

    ⚠️ WHENEVER A CAMPAIGN VARIES SOMETHING NEW, IT GOES HERE. Twice now a
    campaign has been within hours of merging away the very axis it was built
    to resolve: Penning (2026-08-10, caught mid-run) and the amplification gap
    (2026-08-11, caught before submission). Both look perfectly well-formed in
    the merged file, which is what makes the failure mode dangerous.
    """
    tags = [penning_tag(s.get("penning_mode"), s.get("penning_rp"))]
    # Gap: suffix only off the 150 µm nominal, so every product written before
    # 2026-08-11 keeps its key unchanged when re-collected.
    gap = float(s.get("gap_um", 150.0))
    tags.append(None if abs(gap - 150.0) < 1e-9 else f"gap{gap:.0f}um")
    # Field model: `uniform_field` vs a meshfield map is a physics arm too, and
    # the slope hunt shipped both in one product. They stayed apart only
    # because its uniform arm happened to be the sole `auto` one — an accident,
    # not a guarantee.
    fm = str(s.get("field_model") or "")
    tags.append("uniform" if fm == "uniform_field" else None)
    return tuple(tags)


def load(indir):
    """
    Group every result file by (gas, voltage, Penning), reduced on the way in.

    ⚠️ Penning belongs in this key and its absence was a live bug, caught
    2026-08-10 with the 144-slice slope hunt still running. That campaign
    exists precisely to compare rP arms — {0.30, 0.50, 0.65, 0.80} plus an
    `auto` arm — at three shared voltages on one gas. Grouped by (gas, voltage)
    alone, every arm collapses into ONE point per voltage, silently averaged
    over the axis the campaign was built to resolve, and the merged file looks
    perfectly well-formed.

    This is the same failure as the T7 voltage-label incident (§0a): a grouping
    key that ignores what the campaign actually varied. There the 56 slices all
    shared one field map regardless of their voltage label; here the arms would
    all share one key regardless of their Penning label. The lesson generalises
    — whenever a new campaign axis is added, it has to be added HERE too.
    """
    points = defaultdict(list)
    files = sorted(glob.glob(os.path.join(indir, "aval_*.json")))
    for i, f in enumerate(files, 1):
        try:
            s = reduce_file(f)
        except Exception as e:                       # partial transfer
            print(f"  skipping unreadable {os.path.basename(f)}: {e}")
            continue
        points[(s["gas"], s["volt"], arm_tags(s))].append(s)
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

    # THE ACCEPTANCE CHECK (audit A7), asserted here where both routes to the
    # mean charge per drifted electron are in scope:
    #     survival * mean(conditional Polya)  ==  mean over ALL seeds incl zeros
    # It is an identity, so any disagreement is an accounting bug, not physics.
    pol = polya_fit(gains)
    if pol is not None:
        direct = float(np.mean(np.asarray(gains, dtype=float)))
        via = pol["gain_mean_uncond"]
        if abs(via - direct) > 1e-6 * max(abs(direct), 1.0):
            raise AssertionError(
                f"survival decomposition broken: survival*E[g|g>0] = {via:.6g} "
                f"but the direct mean over all seeds is {direct:.6g}")

    r_mean, r_std = _pooled([s["r"] for s in slices])
    t_mean, t_std = _pooled([s["t"] for s in slices])
    n_r = sum(s["r"][0] for s in slices)
    s_r2 = sum(s["r"][2] for s in slices)
    zh = sum(s["zhist"] for s in slices)

    return {
        "n_slices": len(slices), "nev_total": nev,
        "polya": pol,
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
        "campaign_label": slices[0]["campaign_label"],
        # Height above the anode where seeding happened -- sigma0/t_arrival
        # above are measured FROM here, so a consumer must not separately add
        # diffusion for whatever drift leg sits between this point and the
        # mesh (hit for real: Stage B mirrored this as a hardcoded 180 um
        # constant and double-counted it; response/digitizer commit b343c82).
        "seed_z0_um": slices[0]["seed_z0_um"],
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
    print(f"{len(points)} (gas, voltage, arm) points from "
          f"{sum(len(v) for v in points.values())} slice files")

    # The key stays `gas@voltage` while that is unambiguous, so every product
    # written before 2026-08-10 keeps its exact key format and no consumer
    # breaks. An arm suffix appears ONLY where it has to — i.e. where two arms
    # would otherwise collide, which is exactly the case the grouping fix above
    # exists for. Note the arms are SEPARATE regardless (that is `load`'s job);
    # this only decides how much of the distinction the string key shows.
    n_tags = defaultdict(set)
    for gas, volt, tags in points:
        n_tags[(gas, volt)].add(tags)
    # Per (gas, voltage), take the FEWEST leading fragments that already make
    # the keys unique. Minimal rather than "every fragment that differs", so
    # that adding an axis to `arm_tags` cannot rename a key in a product where
    # the older axes were already sufficient — re-collecting the 2026-08-09
    # slope-hunt raw still writes `...@460V@auto`, not `...@460V@auto@uniform`.
    live = {}
    for gv, tagsets in n_tags.items():
        n = 0
        while n < len(next(iter(tagsets))) and len(
                {t[:n] for t in tagsets}) < len(tagsets):
            n += 1
        live[gv] = set(range(n))

    calib, rows = {}, []
    for (gas, volt, tags), sl in sorted(points.items(),
                                        key=lambda kv: str(kv[0])):
        ptag = tags[0]
        suffix = "".join(f"@{t}" for i, t in enumerate(tags)
                         if t is not None and i in live[(gas, volt)])
        m = merge(sl)
        # Machine-readable voltage/gas, not just baked into the string key --
        # a consumer needing "what voltage was this point run at" had to
        # parse "Ar_iC4H10_95_5_Saclay_160m.gas@530V" to get it (hit for
        # real: response/digitizer/digitize.py had to hardcode 490V rather
        # than read it, for aval_calib_meshfield_pooled.json, which predates
        # this fix and still only carries the voltage in free text).
        m["voltage_V"] = float(volt)
        m["gas_file"] = gas
        # Machine-readable Penning too, for the same reason the voltage is:
        # a consumer asking "which rP was this" must not have to parse a key.
        m["penning_mode"] = sl[0].get("penning_mode")
        m["penning_rp"] = sl[0].get("penning_rp")
        m["penning_tag"] = ptag
        # Machine-readable gap for the same reason again — the 135/150 scan's
        # consumer must not have to parse "@gap135um" out of a string key.
        m["gap_um"] = float(sl[0].get("gap_um", 150.0))
        key = f"{gas}@{volt:.0f}V{suffix}"
        rows.append((gas, ptag, volt, m))
        p = m["polya"]
        print(f"  {key:<44s} nev={m['nev_total']:5d}  "
              f"gain={p['gain_mean']:9.1f}  theta={p['theta']:5.2f}  "
              f"surv={p['survival']:.4f}+-{p['survival_err']:.4f}  "
              f"sigma0={m['sigma0_um']:.1f} µm")
        calib[key] = m

    out = a.out or os.path.join(a.indir, "aval_calib_v3.json")
    # v3: same content as v2 plus the survival block (audit A7). Bumped so a
    # consumer can require it rather than discover its absence at runtime.
    json.dump({"schema": "aval_calib/3", "points": calib}, open(out, "w"))
    print(f"wrote {out}")

    # ── figures ──────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 0
    os.makedirs(a.figdir, exist_ok=True)
    # Sort by (gas, volt) explicitly -- (volt, m) used to fall back to
    # comparing the dict m when two different gases share a voltage (hit for
    # real: a 95/5 + 90/10 campaign both including 530V raised
    # "'<' not supported between instances of 'dict' and 'dict'";  the JSON
    # above was already written by that point, so no data was lost, but the
    # figure never got made).
    # Sort and series-split by (gas, Penning) -- one curve per campaign ARM.
    # Grouping by gas alone would zigzag a single line through every rP arm at
    # each shared voltage, which is unreadable and, worse, looks like scatter.
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    series = sorted(set((r[0], r[1]) for r in rows))
    gases = sorted(set(r[0] for r in rows))
    colors = plt.cm.tab10.colors

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    for i, (gas, ptag) in enumerate(series):
        grows = [r for r in rows if r[0] == gas and r[1] == ptag]
        V = [r[2] for r in grows]
        G = [r[3]["polya"]["gain_mean"] for r in grows]
        TH = [r[3]["polya"]["theta"] for r in grows]
        S0 = [r[3]["sigma0_um"] for r in grows]
        c = colors[i % len(colors)]
        label = os.path.splitext(gas)[0]
        if len({t for g, t in series if g == gas}) > 1:
            label = f"{label} {ptag}"
        ax[0].semilogy(V, G, "o-", color=c, ms=6, label=label)
        ax[1].plot(V, TH, "s-", color=c, ms=6, label=label)
        ax[2].plot(V, S0, "^-", color=c, ms=6, label=label)
    ax[0].set_xlabel("mesh voltage [V]"); ax[0].set_ylabel("mean gain")
    ax[0].set_title("Gain vs voltage (150 µm gap)")
    ax[1].set_xlabel("mesh voltage [V]"); ax[1].set_ylabel("Polya θ")
    ax[1].set_title("Polya shape parameter")
    ax[2].set_xlabel("mesh voltage [V]")
    ax[2].set_ylabel("transverse avalanche σ₀ at the ESL [µm]")
    ax[2].set_title("Avalanche footprint")
    for x in ax:
        x.grid(True, color="#e6e6e2")
        x.spines["top"].set_visible(False); x.spines["right"].set_visible(False)
        if len(series) > 1:
            x.legend(fontsize=7)
    fig.suptitle("S3 — avalanche calibration", y=1.03)
    p = os.path.join(a.figdir, "s3_avalanche_calib.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
