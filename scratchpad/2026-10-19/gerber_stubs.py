"""Where do the X and Y bus lines sit relative to their pads? (2026-10-19)
For each layer: the 390 um stubs from a pad centre, their direction (dx, dy) counted, and the bus-line positions
(long 0.1 mm segments) relative to the pad grid. Uses response.common.channel_map.parse_draws.
Usage (repo root): python3 -m scratchpad... or python3 scratchpad/2026-10-19/gerber_stubs.py"""
import sys, os, numpy as np
from collections import Counter
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from response.common.channel_map import parse_draws, STUB_LEN_MM
from response.common import constants as C
G = 'design/gerbers/readout_pcb'; pitch = C.PAD_PITCH_UM * 1e-3; origin = -(C.PAD_N - 1) / 2 * pitch
for tag, fn in (('Y', 'DFS3498A_L3-TrackY.gbr'), ('X', 'DFS3498A_L4-TrackX.gbr')):
    seg = parse_draws(f'{G}/{fn}'); L = np.hypot(seg[:, 2] - seg[:, 0], seg[:, 3] - seg[:, 1])
    st = seg[np.abs(L - STUB_LEN_MM) < 1e-4]
    # orient each stub from its on-grid end (the pad) to the other end (the bus)
    def ongrid(x, y):
        return (np.abs((x - origin) / pitch - np.rint((x - origin) / pitch)) < 1e-3) & (np.abs((y - origin) / pitch - np.rint((y - origin) / pitch)) < 1e-3)
    a = ongrid(st[:, 0], st[:, 1]); b = ongrid(st[:, 2], st[:, 3])
    p = np.where(a[:, None], st[:, :2], st[:, 2:]); q = np.where(a[:, None], st[:, 2:], st[:, :2])
    d = np.round((q - p) * 1000).astype(int)
    print(f'{tag} ({fn}): {len(st)} stubs; start on grid {a.sum()}, end on grid {b.sum()}')
    print('   stub direction (dx, dy) um:', Counter(map(tuple, d)).most_common(6))
    i = np.rint((p[:, 0] - origin) / pitch).astype(int); j = np.rint((p[:, 1] - origin) / pitch).astype(int)
    for nm, k in (('col parity', i % 2), ('row parity', j % 2)):
        for par in (0, 1):
            m = k == par; print(f'   {nm} {par}: ', Counter(map(tuple, d[m])).most_common(3))
    # long bus segments: their fixed coordinate relative to the pad grid
    long_ = seg[L > 2]
    if tag == 'X': c = long_[np.abs(long_[:, 0] - long_[:, 2]) < 1e-4][:, 0]
    else:          c = long_[np.abs(long_[:, 1] - long_[:, 3]) < 1e-4][:, 1]
    frac = np.round(((c - origin) / pitch) % 1, 3)
    print(f'   long bus segments: {len(c)}; position in pad pitches mod 1:', Counter(frac).most_common(4))
