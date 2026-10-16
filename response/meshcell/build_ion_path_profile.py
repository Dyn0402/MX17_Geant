#!/usr/bin/env python3
"""
build_ion_path_profile.py — pack what an ion sees on its way to the mesh into a
small JSON, so the digitizer can rebuild the ion current for ANY mobility table
without the 100 MB field map or a Garfield rerun.

Consumed by response/digitizer/ions.py:reconstructed_longitudinal (the
`--ion-model reconstructed` path). The inputs are the same ones
ion_template_check.py validated against the S3 v2 template (3.8 %):

  * E_z(h)   plane-averaged over the GAS nodes of the T6 production map, at
             V_mesh = 490 V. Laplace is linear in the boundary voltage, so the
             consumer rescales by V/490. Note the bulk value is 31.0 kV/cm, not
             490/150 = 32.7: the woven mesh is transparent enough that its
             equivalent electrode plane sits ~8 um above its underside
             (mesh_plane_check.py), so the effective gap is ~158 um.
  * psi(h)   readout weighting potential through the woven mesh (psi_readout).
  * h_end    where the ions stop: percentiles of the absorbed-on-mesh heights
             (funnel_ion_endpoints). The ~5 % that escape through the mesh stop
             at psi ~ 0.003, i.e. they have delivered everything, so treating
             them as absorbed changes nothing.

    ~/PycharmProjects/nTof_x17/.venv/bin/python build_ion_path_profile.py \
        --field-map meshfield_production.txt
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field-map", default=os.path.join(HERE, "meshfield_production.txt"))
    ap.add_argument("--psi", default=os.path.join(HERE, "psi_readout.json"))
    ap.add_argument("--funnel", default=os.path.join(HERE, "funnel_ion_endpoints.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "ion_path_profile.json"))
    ap.add_argument("--dz-um", type=float, default=0.5)
    a = ap.parse_args()

    m = np.loadtxt(a.field_map, comments="#")
    gas = m[:, 7] == 1
    z_um, ez = m[gas, 2] * 1e4, m[gas, 5]
    zu, inv = np.unique(z_um, return_inverse=True)
    ez_mean = np.bincount(inv, weights=ez) / np.bincount(inv)

    psi_doc = json.load(open(a.psi))
    geo = psi_doc["geometry"]
    z_anode = geo["z_anode_um"]
    zr = np.asarray(psi_doc["psi_profile"]["z_above_anode_um"])
    pv = np.asarray(psi_doc["psi_profile"]["psi"])

    fun = json.load(open(a.funnel))["ion_endpoints"]
    ze = np.array([r["ze_above_anode_um"] for r in fun["rows"]
                   if r["outcome"] == "absorbed_on_mesh"])
    h_end = np.percentile(ze, np.linspace(1, 99, 99))

    h = np.arange(0.0, float(ze.max()) + 1.0 + 1e-9, a.dz_um)
    out = {
        "schema": "ion_path_profile/1",
        "note": ("plane-averaged gas E_z and readout psi vs height above the "
                 "ESL, plus ion absorption-height percentiles; for "
                 "ions.reconstructed_longitudinal"),
        "v_mesh_ref": float(geo["v_mesh"]),
        "amp_gap_um": float(geo["amp_gap_um"]),
        "h_um": h.round(4).tolist(),
        "ez_vcm": np.interp(z_anode + h, zu, ez_mean).round(3).tolist(),
        "psi": np.interp(h, zr, pv).round(7).tolist(),
        "h_end_pctl_um": h_end.round(3).tolist(),
        "frac_escaped": fun["frac_escaped_to_drift"],
        "sources": {"field_map": os.path.basename(a.field_map),
                    "psi": os.path.basename(a.psi),
                    "funnel": os.path.basename(a.funnel)},
    }
    json.dump(out, open(a.out, "w"), indent=0)
    print(f"wrote {a.out}: {len(h)} heights, bulk E_z "
          f"{np.interp(75, h, out['ez_vcm']):.0f} V/cm, "
          f"median absorption {np.median(ze):.1f} um")


if __name__ == "__main__":
    main()
