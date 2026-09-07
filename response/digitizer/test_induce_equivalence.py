"""
test_induce_equivalence.py — the vectorised induction against the original
per-avalanche loop it replaced.

The paths must agree to float noise; the loop is kept ONLY as the oracle.

The speed-up printed here is NOT the production one. This probes 400
avalanches, where the vectorised path's fixed setup cost dominates and it can
come out slower; the 5.8x that motivated it appears at real per-event
occupancies. Correctness is what this test is for.

The kernel and calib used to be hardcoded to one machine's absolute paths, so
this never ran anywhere else and reported nothing when it failed to load
(audit, 2026-08-08).
"""
import argparse, os, time
import numpy as np
from response.common import constants as C
from response.solver import kernels as K
from response.digitizer.digitize import Digitizer

ap = argparse.ArgumentParser()
ap.add_argument("--kernel", required=True)
ap.add_argument("--calib", required=True)
ap.add_argument("--ion-model", choices=("measured", "analytic"),
                default="measured")
a = ap.parse_args()
dig = Digitizer(os.path.expanduser(a.kernel), os.path.expanduser(a.calib),
                seed=3, ion_model=a.ion_model)

def reference(dig, x, y, t, q, n_samp):
    """The original per-avalanche loop, kept only as the equivalence oracle."""
    out = {}
    pitch = C.PAD_PITCH_M
    ds = np.arange(-dig.n_side, dig.n_side + 1)
    ix = dig.lut.ix(x)
    # The column comes from the LUT, matching production (audit C2): the X band
    # defines d = 0 per x SAMPLE, so rounding the true x onto the pad lattice
    # here would make the oracle disagree with the code under test on ~1.3 % of
    # avalanches for a reason that is not a bug in either.
    col = dig.lut.col_at(x, ix=ix)
    row = np.rint((y - K.PAD_ORIGIN_M) / pitch).astype(int)
    k0 = np.rint(t / (dig.lut.dt * 1e9)).astype(int)
    inside = ((col >= -dig.n_side) & (col < C.PAD_N + dig.n_side)
              & (row >= -dig.n_side) & (row < C.PAD_N + dig.n_side))
    ok = (np.isfinite(q) & (q > 0) & (k0 >= 0) & (k0 < n_samp) & inside)
    dy_step = dig.lut.y_Y[1] - dig.lut.y_Y[0]
    dyx_step = dig.lut.y_X[1] - dig.lut.y_X[0]
    for i in np.flatnonzero(ok):
        rows = row[i] + ds
        dy = y[i] - (K.PAD_ORIGIN_M + rows * pitch)
        iy = np.clip(np.rint((dy - dig.lut.y_Y[0]) / dy_step).astype(int),
                     0, len(dig.lut.y_Y) - 1)
        curY = dig.lut.I_Y[rows % 2, iy, ix[i], :]
        y_rel = np.mod(y[i] - K.PAD_ORIGIN_M, 2 * pitch)
        jy = int(np.clip(round((y_rel - dig.lut.y_X[0]) / dyx_step),
                         0, len(dig.lut.y_X) - 1))
        curX = dig.lut.I_X[:, jy, ix[i], :]
        cols = col[i] + ds
        for j in range(len(ds)):
            for key, cur in ((("Y", int(rows[j])), curY[j]),
                             (("X", int(cols[j])), curX[j])):
                n = min(len(cur), n_samp - k0[i])
                if n <= 0: continue
                buf = out.setdefault(key, np.zeros(n_samp))
                buf[k0[i]:k0[i]+n] += q[i] * cur[:n]
    return out

rng = np.random.default_rng(11)
N = 400
x0 = float(K.pad_x(K.nearest_column(C.SUPERPERIOD_M/2)))
y0 = K.PAD_ORIGIN_M
x = x0 + rng.normal(0, 800e-6, N)
y = y0 + rng.normal(0, 800e-6, N)
t = rng.uniform(0, 700, N)
q = rng.gamma(2.6, 1e4, N)
n_samp = 1200

t0=time.time(); a = dig.induce(x, y, t, q, n_samp); ta=time.time()-t0
t0=time.time(); b = reference(dig, x, y, t, q, n_samp); tb=time.time()-t0

keys = set(a) | set(b)
worst = 0.0; scale = max(np.abs(v).max() for v in b.values())
for k in keys:
    va = a.get(k, np.zeros(n_samp)); vb = b.get(k, np.zeros(n_samp))
    worst = max(worst, np.abs(va - vb).max())
print(f"channels: vectorised {len(a)}  reference {len(b)}  union {len(keys)}")
print(f"max abs difference = {worst:.3e}   (signal scale {scale:.3e})")
rel = worst / scale if scale > 0 else float("nan")
print(f"relative           = {rel:.3e}")
print(f"time: vectorised {ta:.3f}s   loop {tb:.3f}s   speedup x{tb/ta:.1f}")
TOL = 1e-6
good = np.isfinite(rel) and rel <= TOL
print("\n" + ("PASS" if good else f"FAIL (relative {rel:.2e} > {TOL:.0e})"))
raise SystemExit(0 if good else 1)
