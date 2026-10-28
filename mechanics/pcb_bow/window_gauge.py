#!/usr/bin/env python3
"""window_gauge.py -- the 60 um mylar window (400.2 mm square flange aperture) as an overpressure gauge.

Same solver (plate.py) in its membrane limit: 60 um PET, nu = 0.38, clamped, immovable edges, no pre-tension
(a taut-glued foil sags less at low pressure).  Cross-checked against the Vlassak & Nix square-membrane law
p = E t w^3 / ((1 - nu)(0.8 + 0.062 nu)^3 a_half^4).  The sag goes as p^(1/3): a centre-sag measurement good to 10 %
gives the overpressure to 30 %.  Writes results/window_gauge.json.
"""
import json
import os
import numpy as np
from plate import Plate, layer_props

HERE = os.path.dirname(os.path.abspath(__file__))
NU = 0.38
out = dict(p_mbar=[1, 2, 3, 5, 7, 10], rows={})
for E in (3.8e9, 4.5e9, 5.0e9):
    L = layer_props([(E, 60e-6)], nu=NU)
    p = Plate(0.4002, 31, L, edge='clamped', inplane='immovable', nu=NU)
    ws = p.sweep([-100.0 * x for x in out['p_mbar']], substeps=3)
    sag = [float(-w[15, 15] * 1e3) for w in ws]
    vn = [float((100 * x * 0.2001 ** 4 * (1 - NU) * (0.8 + 0.062 * NU) ** 3 / (E * 60e-6)) ** (1 / 3) * 1e3)
          for x in out['p_mbar']]
    out['rows'][f'{E / 1e9:.1f}GPa'] = dict(solver_mm=sag, vlassak_nix_mm=vn)
    print(f'E={E / 1e9:.1f} GPa  ' + '  '.join(f'{x:g} mbar: {s:5.1f} ({v:5.1f})' for x, s, v in
                                                zip(out['p_mbar'], sag, vn)))
json.dump(out, open(f'{HERE}/results/window_gauge.json', 'w'), indent=1)
