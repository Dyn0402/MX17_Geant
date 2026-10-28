#!/usr/bin/env python3
"""validate.py -- checks of plate.py against closed-form plate results before any MX17 number is quoted.

1. Linear Kirchhoff, square, nu = 0.3 (Timoshenko & Woinowsky-Krieger tables 8, 35): w0 = alpha q a^4 / D with
   alpha = 0.00406 (simply supported), 0.00126 (clamped); grid convergence.
2. Rotational-spring edge: K -> 0 gives SS, K -> inf gives clamped.
3. Two-layer interface: k = 0 gives D1 + D2, k -> inf gives the composite section (layer_props of the bonded stack).
4. Mesh-tension moment m on all edges: SS square obeys D lap(w) = -m with w = 0 (Marcus), centre |w| = 0.07367 m a^2/D
   (the Prandtl-torsion constant); clamped edges suppress it entirely.
5. Large deflection: a thin plate with immovable clamped edges approaches the square-membrane bulge law
   p = C E t w^3 / a_half^4 (Vlassak & Nix 1992 energy solution, C = (0.8 + 0.062 nu)^-3 / (1 - nu) = 2.60 at
   nu = 0.3; approximate to a few %), and at moderate w/t the plate result is bracketed by linear and membrane.
"""
import json
import os
import time
import numpy as np
from plate import Plate, layer_props

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
E, t, a = 20e9, 1.7e-3, 0.402
L1 = layer_props([(E, t)])
D = L1['D']
q = 100.0

print('1. linear coefficients alpha = w0 D / (q a^4)')
for N in (31, 41, 61, 81):
    row = {}
    for edge, ref in (('ss', 0.00406), ('clamped', 0.00126)):
        p = Plate(a, N, L1, edge=edge)
        w = p.solve(-q, linear=True)
        al = -w[N // 2, N // 2] * D / (q * a ** 4)
        row[edge] = al
        print(f'   N={N:3d} {edge:8s} alpha={al:.5f} (ref {ref}) err {100 * (al / ref - 1):+.2f} %')
    out[f'linear_N{N}'] = row

print('2. rotational spring edge')
N = 61
for K in (0.0, 1.0, 10.0, 100.0, 1e3, 1e4, 1e6):
    p = Plate(a, N, L1, edge='spring', K_rot=K)
    w = p.solve(-q, linear=True)
    al = -w[N // 2, N // 2] * D / (q * a ** 4)
    print(f'   K={K:8.1f} N m/m  K a/D={K * a / D:8.2f}  alpha={al:.5f}')
    out[f'spring_K{K:g}'] = al

print('3. two-layer interface (PCB 1.7 mm + Rohacell 51 5 mm + 25 um mylar)')
L2 = layer_props([(70e6, 5e-3), (4.5e9, 25e-6)])
comp = layer_props([(E, t), (70e6, 5e-3), (4.5e9, 25e-6)])
for kk in (0.0, 1e7, 1e8, 1e9, 1e10, 1e12):
    p = Plate(a, 41, L1, L2, k=kk, edge='clamped')
    w = p.solve(-q, linear=True)
    Deff = 0.00126 * q * a ** 4 / (-w[20, 20])
    print(f'   k={kk:8.1e}  D_eff={Deff:7.3f} N m   (D1+D2={L1["D"] + L2["D"]:.3f}, composite={comp["D"]:.3f})')
    out[f'interface_k{kk:g}'] = Deff
out['D_pcb'] = L1['D']; out['D_sum'] = L1['D'] + L2['D']; out['D_composite'] = comp['D']

print('4. mesh-tension moment m = 1 N m/m')
for edge in ('ss', 'clamped'):
    p = Plate(a, 61, L1, edge=edge)
    w = p.solve(0.0, m=1.0, linear=True)
    # SS square under uniform edge moment: Marcus analogy, D lap(w) = -m, w = 0 on the edge -> Prandtl-torsion
    # constant: |w0| = 0.07367 m a^2 / D
    print(f'   {edge:8s} centre w = {w[30, 30] * 1e3:+.4f} mm   (SS exact: {-0.07367 * a ** 2 / D * 1e3:+.4f} mm, '
          f'clamped exact 0)')
    out[f'mesh_moment_{edge}_mm'] = w[30, 30] * 1e3

print('5. large deflection, clamped immovable')
for tt in (0.2e-3, 1.7e-3):
    Lt = layer_props([(E, tt)])
    for qq in (50.0, 300.0, 1000.0):
        p = Plate(a, 41, Lt, edge='clamped', inplane='immovable')
        t0 = time.time()
        w = p.solve(-qq, steps=8)
        w0 = -w[20, 20]
        wl = 0.00126 * qq * a ** 4 / Lt['D']
        wm = (qq * (a / 2) ** 4 / ((0.8 + 0.062 * 0.3) ** -3 / 0.7 * E * tt)) ** (1 / 3)
        print(f'   t={tt * 1e3:.1f} mm q={qq:6.0f} Pa: w0={w0 * 1e3:7.3f} mm  linear {wl * 1e3:8.3f}  membrane '
              f'{wm * 1e3:7.3f}  w0/t={w0 / tt:5.2f}  ({time.time() - t0:.1f} s)')
        out[f'nl_t{tt * 1e3:.1f}_q{qq:g}'] = dict(w0=w0, linear=wl, membrane=wm)
json.dump(out, open(f'{HERE}/results/validate.json', 'w'), indent=1)
