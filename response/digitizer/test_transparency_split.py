"""
The mesh transparency must be applied EXACTLY ONCE (fix 2026-08-11).

This guards the double-count described in `split_calib_survival`: the S3
calibration's `survival` field means different things under the two field
models, and multiplying the external `DEFAULT_TRANSPARENCY` by a meshfield
calib's `survival` thins the electrons at the mesh twice — which is what Stage
B did from 2026-08-08 (meshfield calib became production) to 2026-08-11,
costing a factor 0.873 * 0.9559 / 0.955 = 0.874 on every simulated charge.
See design/report/TRANSPARENCY_DOUBLE_COUNT_2026-08-11.md.

Runs anywhere: the LUT is stubbed, because building a real one peaks around
3.7 GB and is forbidden on the laptop (RESPONSE_SIM_PLAN machine roles).

    python3 -m response.digitizer.test_transparency_split
"""

from __future__ import annotations

import glob
import json
import os
import sys

from . import digitize as D


HERE = os.path.dirname(os.path.abspath(__file__))
AVAL = os.path.join(HERE, os.pardir, "avalanche")
PRODUCTION_CALIB = os.path.join(AVAL, "aval_calib_meshfield_pooled.json")


class _StubLUT:
    """Enough of CombKernelLUT for __init__; nothing here induces."""

    def __init__(self, *a, **kw):
        import numpy as np
        self.y_Y = np.array([-0.006, 0.0, 0.006])

    def describe(self):
        return {"rho_s_MOhm_sq": None, "d_kapton_um": None}


def _dig(calib_path, **kw):
    real, D.CombKernelLUT = D.CombKernelLUT, _StubLUT
    try:
        return D.Digitizer("stub-kernel", calib_path, with_ions=False, **kw)
    finally:
        D.CombKernelLUT = real


def _mesh_factor(dig):
    """The product that multiplies every drifting electron's survival."""
    return dig.transparency * dig.aval_survival


def check_production():
    dig = _dig(PRODUCTION_CALIB)
    surv = json.load(open(PRODUCTION_CALIB))["point"]["polya"]["survival"]

    assert dig.transparency == surv, (
        f"meshfield calib must SUPPLY the transparency, got {dig.transparency}")
    assert dig.aval_survival == 1.0, (
        f"P(g>0) must be 1.0 once `survival` is read as the mesh term, "
        f"got {dig.aval_survival}")
    # The whole point: one factor, not two.
    assert abs(_mesh_factor(dig) - surv) < 1e-12, (
        f"mesh thinning is {_mesh_factor(dig)}, must be exactly {surv}")
    # And it must agree with T6's independent 3-D G7 measurement, which is the
    # cross-check the double-count was hiding.
    assert abs(_mesh_factor(dig) - 0.955) < 0.02, (
        f"{_mesh_factor(dig)} disagrees with T6 G7 0.955 -0.045/+0.005")

    old = 0.873 * surv
    print(f"  production (meshfield)   eps={dig.transparency:.6f} "
          f"P(g>0)={dig.aval_survival}  net={_mesh_factor(dig):.6f}   "
          f"was {old:.6f}  ->  charge x{_mesh_factor(dig)/old:.4f}")


def check_uniform_field():
    """A uniform_field calib has no mesh in it, so it NEEDS the external one."""
    doc = {"schema": "aval_calib/2",
           "points": {"Ar_iC4H10_95_5@490V": {
               "field_model": "uniform_field",
               "polya": {"gain_mean": 44600.0, "theta": 1.1},
               "sigma0_um": 30.0, "voltage_V": 490.0}}}
    path = os.path.join(HERE, "_tmp_uniform_calib.json")
    with open(path, "w") as f:
        json.dump(doc, f)
    try:
        dig = _dig(path)
        assert dig.calib_transparency is None
        assert dig.transparency == D.DEFAULT_TRANSPARENCY
        assert dig.aval_survival == 1.0      # absent in schema <= 2
        # T6's 3-D value, not the superseded 1-D 0.873 (runbook, 2026-08-08).
        assert D.DEFAULT_TRANSPARENCY == 0.955
        print(f"  uniform_field            eps={dig.transparency:.6f} "
              f"[external, T6 3-D]  net={_mesh_factor(dig):.6f}")
    finally:
        os.remove(path)


def check_uniform_field_with_survival():
    """There, `survival` really is P(g>0) and must stay on that factor."""
    doc = {"schema": "aval_calib/3",
           "points": {"Ar_iC4H10_95_5@490V": {
               "field_model": "uniform_field",
               "polya": {"gain_mean": 44600.0, "theta": 1.1, "survival": 0.97},
               "sigma0_um": 30.0, "voltage_V": 490.0}}}
    path = os.path.join(HERE, "_tmp_uniform_surv_calib.json")
    with open(path, "w") as f:
        json.dump(doc, f)
    try:
        dig = _dig(path)
        assert dig.transparency == D.DEFAULT_TRANSPARENCY
        assert dig.aval_survival == 0.97
        print(f"  uniform_field + survival eps={dig.transparency:.6f} "
              f"P(g>0)={dig.aval_survival}  net={_mesh_factor(dig):.6f}")
    finally:
        os.remove(path)


def check_override():
    """An explicit transparency still wins, and says so in the provenance."""
    dig = _dig(PRODUCTION_CALIB, transparency=0.5)
    assert dig.transparency == 0.5
    assert dig.describe()["transparency_source"] == "caller override"
    print("  caller override          honoured, and labelled in describe()")


def check_every_shipped_calib():
    """No shipped calibration may end up thinned twice, or raise."""
    seen = 0
    for path in sorted(glob.glob(os.path.join(AVAL, "*.json"))):
        try:
            pt, _ = D.load_calib(path, 490.0)
        except Exception:
            continue                       # not a calib container; not our test
        eps, pmul = D.split_calib_survival(pt)
        net = (eps if eps is not None else D.DEFAULT_TRANSPARENCY) * pmul
        assert 0.80 < net <= 1.0, f"{os.path.basename(path)}: net {net}"
        # The double-count signature: a meshfield calib whose measured mesh
        # term has been multiplied by a second one lands well under 0.90.
        assert not (str(pt.get("field_model", "")).startswith("meshfield")
                    and net < 0.90), (
            f"{os.path.basename(path)} looks double-counted: {net}")
        seen += 1
    print(f"  {seen} shipped calibrations, none double-counted")


def main():
    print("Mesh transparency applied exactly once\n")
    check_production()
    check_uniform_field()
    check_uniform_field_with_survival()
    check_override()
    check_every_shipped_calib()
    print("\nPASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
