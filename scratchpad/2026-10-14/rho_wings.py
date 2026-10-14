#!/usr/bin/env python3
"""Shaped point-source kernel ratios (Y wings, Y0/X0) vs rho_s from the W2 kernels. Run from repo root: PYTHONPATH=. python scratchpad/2026-10-14/rho_wings.py"""
import sys, numpy as np
sys.argv = ['x', 'unused', '']
from response.comsol import compare_kernels3d as CK
from response.dream.shaper import DreamShaper
cols = ['X0','X1','Y0','Y1','Y2','Y3','Y4']
dt = 1.0; tg = np.arange(0, 3000, dt); sh = DreamShaper(dt_ns=dt)
for rho in ('0.5M','1M','2M','5M'):
    tk, G = CK.chain_G(f'/media/ucla/x17_resp/s1_w2_ny512/greens_comb_w2_rho{rho}_dk50um_g19um.npz', cols)
    pk = {c: CK.shaped_from_charge(tk, G[c], tg[tg <= tk[-1]*1e9], sh, dt) for c in cols}
    p = {c: np.max(np.abs(v)) for c, v in pk.items()}
    t = {c: tg[np.argmax(np.abs(v))] for c, v in pk.items()}
    print(f'rho {rho:5s} X0 {p["X0"]:.3f} Y0 {p["Y0"]:.3f} | Y2/Y0 {p["Y2"]/p["Y0"]:.3f} (t {t["Y2"]:.0f}) Y4/Y0 {p["Y4"]/p["Y0"]:.3f} | X1/X0 {p["X1"]/p["X0"]:.3f}  Y share of view? Y0/(X0+Y0) {p["Y0"]/(p["X0"]+p["Y0"]):.3f}')
