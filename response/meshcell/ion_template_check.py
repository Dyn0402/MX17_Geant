#!/usr/bin/env python3
"""
ion_template_check.py — validate the S3 v2 measured i_ion template against an
INDEPENDENT first-principles reconstruction.

Work-list item 3 of design/report/HANDOFF_S3_ION_2026-08-09.md:

    "Validate the S3 v2 i_ion template: which ion species/mobility (Ar+ vs
     iC4H10+ vs cluster ions differ ~x2), kinematic consistency of
     172-ns-to-half against the calib's own alpha_z_hist + amp-gap field."

and the reason it is asked at all: the schema-1 calib shipped i_elec/i_ion as
2000 zeros and silently NaN'd the LUT. The v2 arrays are populated, but
populated is not validated.

WHAT THIS DOES. It rebuilds the ion current from four inputs that share NO code
with the Garfield run that produced the template, and compares the delivery
quantiles:

  1. E_z(z) — plane-averaged from the T6 production field map
     (meshfield_production.txt), the same map S3's meshfield branch drifts in.
  2. K0(E/N) — Garfield's own Ar+/Ar mobility table, evaluated at the FIELD
     THE IONS ACTUALLY SEE rather than at zero field. This is the whole game:
     the amp gap runs at 31.0 kV/cm = 123.8 Td, where K0 = 1.212, not the
     zero-field 1.53 that `ions.py`'s analytic model uses. The analytic
     rectangle is therefore ~21 % too FAST, and the measured template being
     SLOWER than it is a feature, not a discrepancy.
  3. psi(z) — the true readout weighting potential through the woven mesh,
     from psi_readout.py. Charge induced per step is -dpsi, not -dz/gap.
  4. The birth-height distribution (calib alpha_z_hist) and the absorption-
     height distribution (funnel_ion_endpoints.json), so both ends of the
     ion's path are measured rather than assumed to be 0 and exactly `gap`.

SPECIES. The emitter hardcodes IonMobility_Ar+_Ar.txt for every gas. In
Ar/iC4H10 95/5 charge transfer moves the charge to the lowest-IP species, so
the real drifter is an isobutane / cluster ion. Garfield ships no iC4H10+-in-Ar
table; using the two it does ship, Blanc's law over the mixture gives
K0 ~ 1.46 against Ar+'s 1.53, i.e. the real ion is a few per cent SLOWER. The
handoff's "x2 mobility" worry is not there -- but the sign matters and is
reported below, because a slower ion makes the modelled rise WORSE, not better.

Run (seconds, laptop; needs psi_readout.json to exist first):

    ~/PycharmProjects/nTof_x17/.venv/bin/python ion_template_check.py \
        --field-map /media/dylan/data/x17/response_sim/meshfield/meshfield_production.txt
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Loschmidt number [cm^-3] and the bench gas state. The gas tables are Saclay
# 160 m, i.e. essentially sea level; 1 atm / 293 K is the reference the
# mobility table's K0 is reduced to.
N_LOSCHMIDT = 2.6868e19
T_GAS_K = 293.0
P_GAS_PA = 101325.0
K_BOLTZ = 1.380649e-23


def load_mobility(path):
    rows = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        a, b = ln.split()[:2]
        rows.append((float(a), float(b)))
    return np.array(rows)


def ez_profile(path):
    """Plane-averaged E_z(z) over the gas nodes of the exported map."""
    import pandas as pd
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                     names=["x", "y", "z", "ex", "ey", "ez", "v", "flag"])
    g = df[df.flag == 1].groupby("z")["ez"].mean()
    return np.column_stack([g.index.values * 1e4, g.values])   # um, V/cm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field-map", default="/media/dylan/data/x17/"
                    "response_sim/meshfield/meshfield_production.txt")
    ap.add_argument("--mobility", default="/home/dylan/garfield/Data/"
                    "IonMobility_Ar+_Ar.txt")
    ap.add_argument("--mobility-quencher", default="/home/dylan/garfield/Data/"
                    "IonMobility_C8Hn+_iC4H10.txt")
    ap.add_argument("--quencher-frac", type=float, default=0.05)
    ap.add_argument("--psi", default=os.path.join(HERE, "psi_readout.json"))
    ap.add_argument("--calib", default=os.path.join(
        HERE, "..", "avalanche", "aval_calib_meshfield_pooled.json"))
    ap.add_argument("--funnel", default=os.path.join(
        HERE, "funnel_ion_endpoints.json"))
    ap.add_argument("--out", default=os.path.join(HERE,
                                                  "ion_template_check.json"))
    a = ap.parse_args()

    n_gas = P_GAS_PA / (K_BOLTZ * T_GAS_K) / 1e6              # cm^-3
    scale = N_LOSCHMIDT / n_gas                               # K = K0 * this
    tab = load_mobility(a.mobility)
    ez = ez_profile(a.field_map)

    psi_doc = json.load(open(a.psi))
    zr = np.array(psi_doc["psi_profile"]["z_above_anode_um"])
    pv = np.array(psi_doc["psi_profile"]["psi"])
    z_anode = psi_doc["geometry"]["z_anode_um"]

    cal = json.load(open(os.path.expanduser(a.calib)))["point"]
    zh = cal["alpha_z_hist"]
    cnt = np.array(zh["counts"], float)
    ed = np.array(zh["edges"], float)
    zc, w = 0.5 * (ed[:-1] + ed[1:]), cnt / cnt.sum()

    fun = json.load(open(a.funnel))["ion_endpoints"]
    ze = np.array([r["ze_above_anode_um"] for r in fun["rows"]
                   if r["outcome"] == "absorbed_on_mesh"])

    def v_um_ns(h):
        """Ion speed at height h above the ESL [um/ns]."""
        E = np.maximum(np.interp(z_anode + h, ez[:, 0], ez[:, 1]), 1.0)
        td = E / n_gas / 1e-17
        return np.interp(td, tab[:, 0], tab[:, 1]) * scale * E * 1e4 / 1e9

    # ── the reconstruction ───────────────────────────────────────────────────
    grid = np.linspace(0.0, float(ze.max()) + 1.0, 3401)
    dz = grid[1] - grid[0]
    vg = v_um_ns(grid)
    t_of_h = np.concatenate([[0.0], np.cumsum(dz / (0.5 * (vg[1:] + vg[:-1])))])
    psi_g = np.interp(grid, zr, pv)

    nb, tmax = 3000, 600.0
    tb = np.linspace(0.0, tmax, nb + 1)
    q = np.zeros(nb)
    zes = np.percentile(ze, np.linspace(1, 99, 99))
    for wi, h0 in zip(w, zc):
        if wi <= 0:
            continue
        t0 = np.interp(h0, grid, t_of_h)
        for hend in zes:
            sel = (grid >= h0) & (grid <= hend)
            if sel.sum() < 2:
                continue
            tt = np.interp(grid[sel], grid, t_of_h) - t0
            dq = -np.diff(psi_g[sel])
            idx = np.clip(np.searchsorted(tb, tt[1:]) - 1, 0, nb - 1)
            np.add.at(q, idx, dq * wi / len(zes))
    c = np.cumsum(q)
    c /= c[-1]
    tc = 0.5 * (tb[1:] + tb[:-1])
    # AVALANCHE DEVELOPMENT TIME. This reconstruction starts every ion's clock
    # at its own birth; the Garfield template's clock starts when the SEED is
    # launched, 180 um above the mesh. So the template is later than this by
    # the time the avalanche took to reach each birth height — at most the
    # calib's own measured mean electron arrival time. Adding it is not a fit:
    # t_arrival_mean_ns is read from the same calib being checked, and it is
    # the only free constant in the comparison. Reported both ways, because a
    # reader should see that the SHAPE agrees before any offset is applied.
    t_off = float(cal["t_arrival_mean_ns"])
    qs = (0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
    recon_raw = {f"q{int(f*100)}": float(np.interp(f, c, tc)) for f in qs}
    recon = {k: v + t_off for k, v in recon_raw.items()}

    # ── the template it is being checked against ─────────────────────────────
    dt = cal["signal_dt_ns"]
    ii = np.array(cal["i_ion"], float)
    tt = (np.arange(len(ii)) + 1) * dt
    ci = np.cumsum(ii) / ii.sum()
    meas = {f"q{int(f*100)}": float(np.interp(f, ci, tt)) for f in qs}

    # ── the analytic rectangle, for contrast ─────────────────────────────────
    t_rect = (150e-4 ** 2) / (1.5 * 490.0) * 1e9
    anal = {k: float(t_rect * int(k[1:]) / 100.0) for k in meas}

    print(f"gas N = {n_gas:.4e} /cm3   K = K0 x {scale:.4f}")
    e_amp = float(np.interp(z_anode + 75.0, ez[:, 0], ez[:, 1]))
    td_amp = e_amp / n_gas / 1e-17
    k0_amp = float(np.interp(td_amp, tab[:, 0], tab[:, 1]))
    print(f"amp gap: E = {e_amp:.0f} V/cm = {td_amp:.1f} Td  ->  "
          f"K0 = {k0_amp:.3f}  (zero-field {tab[0,1]:.3f}, "
          f"{100*(1-k0_amp/tab[0,1]):.0f} % lower)")
    print(f"ion path: born <{float((w*zc).sum()):.1f}> um, absorbed "
          f"<{ze.mean():.1f}> um above the ESL")

    tabq = load_mobility(a.mobility_quencher)
    k0_q = float(np.interp(td_amp, tabq[:, 0], tabq[:, 1]))
    k_blanc = 1.0 / ((1 - a.quencher_frac) / k0_amp + a.quencher_frac / k0_q)
    print(f"species: Ar+/Ar K0 = {k0_amp:.3f}; Blanc over {a.quencher_frac:.0%}"
          f" quencher (C8Hn+ K0 = {k0_q:.3f}) -> {k_blanc:.3f}, i.e. the real "
          f"ion is {100*(k0_amp/k_blanc-1):.0f} % SLOWER")

    print("\nion charge delivered, cumulative [ns]")
    print(f"{'':>26}" + "".join(f"{k:>9}" for k in meas))
    for name, d in (("recon, ion clock", recon_raw),
                    ("recon + t_aval offset", recon),
                    ("S3 v2 measured template", meas),
                    ("analytic 306 ns rect", anal)):
        print(f"{name:>26}" + "".join(f"{d[k]:9.1f}" for k in meas))
    dev = [100 * (recon[k] - meas[k]) / meas[k] for k in meas]
    print(f"{'recon vs measured [%]':>26}" + "".join(f"{x:+9.1f}"
                                                     for x in dev))
    worst = max(abs(x) for x in dev)
    ok = worst < 5.0
    print(f"\n(t_aval offset = {t_off:.2f} ns, the calib's own measured "
          f"t_arrival_mean_ns — the one constant in this comparison)")
    print(f"GATE (all quantiles within 5 %): {'PASS' if ok else 'FAIL'} "
          f"(worst {worst:.1f} %)")

    json.dump({"schema": "ion_template_check/1",
               "gas": {"n_cm3": n_gas, "T_K": T_GAS_K, "K_scale": scale},
               "amp": {"E_Vcm": e_amp, "E_over_N_Td": td_amp,
                       "K0_at_field": k0_amp, "K0_zero_field": float(tab[0, 1]),
                       "K0_blanc_mixture": k_blanc},
               "path": {"z_birth_mean_um": float((w * zc).sum()),
                        "z_absorb_mean_um": float(ze.mean())},
               "t_avalanche_offset_ns": t_off,
               "quantiles_ns": {"recon": recon, "recon_raw": recon_raw,
                                "measured": meas,
                                "analytic_rect": anal},
               "deviation_pct": dict(zip(meas, dev)),
               "gate_pass": bool(ok), "worst_dev_pct": worst,
               "inputs": {"field_map": a.field_map, "mobility": a.mobility,
                          "psi": os.path.abspath(a.psi),
                          "calib": os.path.abspath(
                              os.path.expanduser(a.calib))}},
              open(a.out, "w"), indent=1)
    print(f"wrote {a.out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
