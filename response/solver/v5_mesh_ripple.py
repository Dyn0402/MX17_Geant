#!/usr/bin/env python3
"""
v5_mesh_ripple.py — V5: is treating the woven mesh as a solid plane safe?

THE ASSUMPTION UNDER TEST (plan §3, geometry model W1): the mesh at z = +g is
modelled as a flat grounded plane. A real micromesh is a weave of 19 um wires
on a 67 um pitch, so the potential on that boundary is not uniform — it is
modulated at the weave pitch, with the wires held at the mesh potential and the
holes floating toward the drift field.

WHY IT IS SAFE, AND THE POINT OF THIS CHECK. A boundary modulation of spatial
period p is a Fourier component with transverse wavenumber k = 2*pi/p, and
Laplace's equation makes every such component decay as exp(-k z) away from the
boundary. The ESL sits a full amplification gap away, z = g = 150 um, and the
weave pitch is p = 67 um, so the fundamental is suppressed by

    exp(-2*pi*g/p)

which is a very large exponent. The plan asks for < 1 %. The check is worth
doing rather than asserting because the answer depends on g/p, and g/p is only
2.24 — "many pitches away" is not obviously true at a glance, and a mesh with a
coarser weave or a narrower gap would fail it.

WHAT THIS DOES NOT COVER. This is the ripple in the WEIGHTING geometry, i.e.
whether the S1 solve may treat the upper boundary as flat. It says nothing
about mesh transparency or about electron focusing through the holes, which are
real, first-order, and handled by T6 (response/meshcell). Those are transport,
not boundary shape.

    python3 -m response.solver.v5_mesh_ripple
"""

from __future__ import annotations

import numpy as np

from ..common import constants as C

# From shared/MX17ModuleGeometry.hh via response/meshcell/mesh_transparency.C:
# 19 um wires on a 67 um pitch (400 lpi / 18 um weave scaled by 5.5 %).
WEAVE_PITCH_M = 67e-6
WIRE_DIAM_M = 19e-6

TOL = 0.01          # plan §V5: ripple at the ESL < 1 %


def ripple_at(z_m, pitch_m=WEAVE_PITCH_M, n_harm=4):
    """
    Relative amplitude of the weave ripple a distance z below the mesh.

    Each harmonic n of the boundary modulation decays as exp(-2*pi*n*z/p). The
    fundamental dominates so completely here that the sum is reported mainly to
    show the higher harmonics are not sneaking anything in.
    """
    n = np.arange(1, n_harm + 1)
    return np.exp(-2.0 * np.pi * n * z_m / pitch_m)


def main():
    g = C.AMP_GAP_M
    p = WEAVE_PITCH_M

    print("V5 — mesh-as-plane, weave ripple at the ESL\n")
    print(f"  weave pitch p      {p*1e6:.0f} µm")
    print(f"  wire diameter      {WIRE_DIAM_M*1e6:.0f} µm")
    print(f"  amplification gap  {g*1e6:.0f} µm     g/p = {g/p:.2f}")
    print()

    h = ripple_at(g, p)
    for i, v in enumerate(h, start=1):
        print(f"  harmonic n={i}   exp(-2*pi*{i}*g/p) = {v:.3e}")
    total = float(h.sum())
    print(f"\n  total ripple at the ESL (z = 0, i.e. g below the mesh): "
          f"{total:.3e}")

    ok = total < TOL
    print(f"  V5 requires < {TOL:.0%}  ->  {'PASS' if ok else 'FAIL'}")

    # The scale on which the assumption WOULD break, so the result is a
    # statement about this detector rather than a universal one.
    g_fail = -p * np.log(TOL) / (2.0 * np.pi)
    print(f"\n  the flat-mesh assumption reaches the 1 % limit at "
          f"g = {g_fail*1e6:.0f} µm for this weave;")
    print(f"  MX17 runs at {g*1e6:.0f} µm, a factor {g/g_fail:.1f} inside it.")

    # Sanity: the ripple must fall monotonically and be ~1 at the mesh itself.
    at_mesh = float(ripple_at(1e-9, p).sum())
    assert at_mesh > 3.0, "ripple should be O(n_harm) right at the boundary"
    assert ripple_at(g / 2, p)[0] > h[0], "ripple must decay with distance"

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
