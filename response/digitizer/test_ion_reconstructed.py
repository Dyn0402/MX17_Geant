#!/usr/bin/env python3
"""
test_ion_reconstructed.py — the reconstructed ion current must reproduce the
measured one when it is handed the same physics.

With the Ar+ table (what the S3 emitter hardcodes) ions.reconstructed_longitudinal
has to land on the meshfield calib's own i_ion template: every charge-delivery
quantile within 5 % (ion_template_check.py got 3.8 % from the same inputs), and
the charge split exactly the calib's. Then two consistency checks of the
mobility handling:

  * a flat K0 equal to the Ar+ table's value at the bulk field must land on
    the table result (the field is uniform below the mesh, so to <1 %);
  * ion_speed_scale = s on the Ar+ table is exactly a 1/s time compression of
    the ion clock (no change to the electron spike or the avalanche offset).

Needs the Garfield Ar+ table (looked up in the Garfield Data directories).

    python3 -m response.digitizer.test_ion_reconstructed
"""
from __future__ import annotations

import json
import os

import numpy as np

from . import ions as ION

CALIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                     "avalanche", "aval_calib_meshfield_pooled.json")
QS = (0.10, 0.25, 0.50, 0.75, 0.90, 0.99)


def ion_quantiles(h, dt, f_e_area, ie):
    """Quantiles of the ION part of a unit-area profile h [ns]."""
    dt_ns = dt * 1e9
    ii = h.copy()
    ii[:len(ie)] -= ie / (ie.sum() * dt) * f_e_area
    c = np.cumsum(ii) / ii.sum()
    t = (np.arange(len(ii)) + 1) * dt_ns
    return np.array([np.interp(f, c, t) for f in QS])


def main():
    pt = json.load(open(CALIB))["point"]
    ie = np.asarray(pt["i_elec"], float)

    hm, dt, im = ION.measured_longitudinal(pt)
    hr, dtr, ir = ION.reconstructed_longitudinal(pt, "IonMobility_Ar+_Ar.txt")
    assert dt == dtr
    assert abs(ir["f_ion"] - im["f_ion"]) < 1e-12, (ir["f_ion"], im["f_ion"])
    assert abs(hr.sum() * dt - 1.0) < 1e-9

    qm = ion_quantiles(hm, dt, im["f_electron"], ie)
    qr = ion_quantiles(hr, dt, ir["f_electron"], ie)
    dev = (qr - qm) / qm * 100
    print(f"mobility {ir['mobility']}  bulk E {ir['E_bulk_Vcm']:.0f} V/cm, "
          f"K0 {ir['K0_bulk']:.3f}")
    print("        " + "".join(f"{f'q{int(f*100)}':>8}" for f in QS))
    print("measured" + "".join(f"{x:8.1f}" for x in qm))
    print("recon   " + "".join(f"{x:8.1f}" for x in qr))
    print("dev %   " + "".join(f"{x:+8.1f}" for x in dev))
    assert np.all(np.abs(dev) < 5.0), "reconstruction misses the measured template"

    # a flat K0 at the Ar+ table's bulk value: the field is uniform across
    # almost the whole gap, so this must land on the table result
    hf, _, rf = ION.reconstructed_longitudinal(pt, f"K0={ir['K0_bulk']:.5f}")
    qf = ion_quantiles(hf, dt, rf["f_electron"], ie)
    d2 = (qf - qr) / qr * 100
    print(f"\nflat K0={ir['K0_bulk']:.3f} vs the Ar+ table:  "
          + " ".join(f"{x:+.2f}%" for x in d2))
    assert np.all(np.abs(d2) < 1.0), "flat-K0 path disagrees with the table"

    # ion_speed_scale is a pure compression of the ion clock
    h2, _, r2 = ION.reconstructed_longitudinal(pt, "IonMobility_Ar+_Ar.txt",
                                               ion_speed_scale=2.0)
    q2 = ion_quantiles(h2, dt, r2["f_electron"], ie)
    t_off = pt["t_arrival_mean_ns"]
    ratio = (qr[1:] - t_off) / (q2[1:] - t_off)
    print(f"ion_speed_scale 2 compresses the ion clock by "
          + " ".join(f"{x:.3f}" for x in ratio))
    assert np.all(np.abs(ratio - 2.0) < 0.05)
    print("\nPASS")


if __name__ == "__main__":
    main()
