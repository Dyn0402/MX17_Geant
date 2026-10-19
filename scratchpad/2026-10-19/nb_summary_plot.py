"""Summary figure (2026-10-19), numbers from nb_side.log / nb_prompt.log / nb_copy_chip_edge (interior) / noise_corr_inject.log.
Left: the one-sided X copy a per detector/run (signal) with the noise test's distance-1 excess. Right: the early-window
in-phase content alpha_|k| vs strip distance (weekend det3, clean side for X), data vs v3 / v4c, log scale."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
runs = ['det3 6-27\nFEU7 490V', 'det3 6-28\nFEU7 490V', 'det3 6-22\nFEU3 495V', 'det3 6-23\nFEU3 495V\ndrift 600', 'det2 6-22\nFEU6', 'det6 6-26\nFEU3']
a = [0.092, 0.117, 0.166, 0.172, 0.013, 0.015]; ea = [0.019, 0.011, 0.013, 0.012, 0.011, 0.021]
nx = [-0.033, 0.031, -0.021, np.nan, 0.004, 0.090]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
x = np.arange(len(runs)); ax[0].bar(x - .18, a, .36, yerr=ea, label='signal copy a (X)', color='C3')
ax[0].bar(x + .18, nx, .36, color='0.6', label='noise: corr excess at channel distance 1')
ax[0].axhline(0, color='k', lw=.6); ax[0].set_xticks(x, runs, fontsize=8); ax[0].legend(fontsize=8)
ax[0].set_title('One-sided X copy: present in signal, absent in noise\n(an electronic copy a = 0.1 / 0.2 would add +0.06..0.12 / +0.13..0.22 to the noise excess)', fontsize=9)
k = np.arange(1, 6)
dX = [0.334, 0.0495, 0.0124, 0.0055, 0.0034]; dY = [(0.4298 + 0.4486) / 2, (0.1155 + 0.1285) / 2, (0.0227 + 0.0290) / 2, (0.0061 + 0.0095) / 2, (0.0024 + 0.0028) / 2]
v4X = [0.402, 0.098, 0.0035, 0.0002, 0.0001]; v4Y = [0.481, 0.103, 0.0082, 0.0001, 0.0001]
v3X = [0.364, 0.082, 0.0038, 0.0005, 0.0003]; v3Y = [0.418, 0.089, 0.0072, 0.0002, 0.0001]
for y, st, lab in ((dX, 'C0o-', 'data X (+ side)'), (dY, 'C1o-', 'data Y'), (v4X, 'C0s--', 'v4c X'), (v4Y, 'C1s--', 'v4c Y'), (v3X, 'C0^:', 'v3 X'), (v3Y, 'C1^:', 'v3 Y')):
    ax[1].semilogy(k, np.clip(y, 1e-4, None), st, label=lab)
ax[1].set_xlabel('strip distance |k|'); ax[1].set_ylabel('in-phase content α (−200..+300 ns) / central'); ax[1].grid(alpha=.3, which='both')
ax[1].set_title('Prompt content vs distance (weekend det3, n ≈ 900): data has a tail to ±5 the sim lacks,\nand a narrower core (X ±1, ±2) than v4c', fontsize=9); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig('nb_summary.png', dpi=90)
