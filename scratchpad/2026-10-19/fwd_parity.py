"""Forward model of the parity-split copy (2026-10-19 evening): uniform Gaussian clouds over the checkerboard, every X and Y
line's induced charge from the bus_wpot.py weighting maps at the ESL, then the data's own classification and observables.
X line k (a pad row) = the solved row-0 map shifted by (780 k, 780 k). Y lines (pad columns) use the bus-free map rotated
by 90 degrees (their bus is on the deeper layer: ~0, as the solve found). Peak X line k, peak Y line c; parity
(k + c) % 2 picks the pad type under the cloud (checked against the A_x/A_y split, as for the data).
Observables per parity: copy = <(S_{k+1} - S_{k-1}) / S_k> (sign: + = toward the bus side), excess/Y, X and Y +-1 shares.
Usage: fwd_parity.py MAP_ESL.npy [MAP_ESL_NOBUS.npy] [--sigma 150,250,350,450] [--n 20000] [--shift-y 0,60] [--json OUT]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('mapx'); ap.add_argument('mapy', nargs='?', default='bus_par/d20_nobus_esl.npy')
ap.add_argument('--sigma', default='150,250,350,450'); ap.add_argument('--n', type=int, default=20000); ap.add_argument('--d', type=float, default=20.)
ap.add_argument('--json', default=None)
ap.add_argument('--shift-y', default='0', help='um, comma list: charge landing over a Y pad is moved by this along -y (bus side)')
a = ap.parse_args(); P = 780.; d = a.d
EX = np.load(a.mapx).astype(float); EY0 = np.load(a.mapy).astype(float)
nx, ny = EX.shape; nrow = int(round(ny * d / P / 2))
xc = (np.arange(nx) + .5) * d - P / 2; yc = (np.arange(ny) + .5) * d - nrow * P
from scipy.interpolate import RegularGridInterpolator
def mk(E):
    # periodic in x (1560): pad by wrapping; zero outside the solved y range
    Ex = np.concatenate([E[-10:], E, E[:10]]); xx = np.concatenate([xc[-10:] - 2 * P, xc, xc[:10] + 2 * P])
    f = RegularGridInterpolator((xx, yc), Ex, bounds_error=False, fill_value=0.)
    return lambda x, y: f(np.stack([((x + P / 2) % (2 * P)) - P / 2, y], -1))
fX, fY0 = mk(EX), mk(EY0)
# sub-cloud quadrature: deposit = 2D Gaussian sampled on a 7x7 grid of +-3 sigma
g = np.linspace(-3, 3, 13); GX, GY = np.meshgrid(g, g, indexing='ij'); W = np.exp(-.5 * (GX ** 2 + GY ** 2)); W /= W.sum()
K = np.arange(-4, 5)
def ypad(px, py):
    i = np.rint((px - P / 2) / P).astype(int); j = np.rint(py / P).astype(int)
    inside = (np.abs(px - (P / 2 + P * i)) <= 340) & (np.abs(py - P * j) <= 340)
    return inside & ((i + j) % 2 == 1)
def signals(cx, cy, s, dsh=0.):
    px = cx[:, None, None] + s * GX[None]; py = cy[:, None, None] + s * GY[None]
    if dsh: py = py - dsh * ypad(px, py)
    SX = np.stack([(fX(px - P * k, py - P * k) * W).sum((1, 2)) for k in K], 1)
    # Y column c: the row map rotated (along <- y, across <- x), pads on the Y sites (-390 + 780 c, 780 c + 1560 m)
    SY = np.stack([(fY0(py - P * c + P / 2, px + P / 2 - P * c) * W).sum((1, 2)) for c in K], 1)
    return SX, SY
rng = np.random.default_rng(1); RES = []
for dsh in [float(v) for v in a.shift_y.split(',')]:
  for s in [float(v) for v in a.sigma.split(',')]:
      cx = rng.uniform(-P, P, a.n); cy = rng.uniform(-P / 2, P / 2, a.n)
      SX, SY = signals(cx, cy, s, dsh)
      kx = SX[:, 1:-1].argmax(1) + 1; ky = SY[:, 1:-1].argmax(1) + 1; ii = np.arange(a.n)
      c0 = SX[ii, kx]; up = SX[ii, kx + 1]; dn = SX[ii, kx - 1]; y0 = SY[ii, ky]
      par = (K[kx] + K[ky]) % 2
      copy = (dn - up) / c0; exc = (dn - up) / y0; xs = (dn + up) / c0; ys = (SY[ii, ky - 1] + SY[ii, ky + 1]) / y0
      line = [f'shift {dsh:3.0f} sigma {s:3.0f}:']
      for p in (0, 1):
          m = par == p
          line.append(f'par{p} (Ax/Ay {np.median(c0[m] / y0[m]):.2f}) copy {copy[m].mean():+.3f} exc/Y {exc[m].mean():+.3f} Xsh {xs[m].mean():.2f} Ysh {ys[m].mean():.2f}')
      line.append(f'all copy {copy.mean():+.3f}')
      yp = 0 if ys[par == 0].mean() < ys[par == 1].mean() else 1      # the Y-pad class, read as in the data
      se = lambda v: float(v.std() / np.sqrt(len(v)))
      RES.append(dict(shift=dsh, sigma=s, ypad=[float(copy[par == yp].mean()), se(copy[par == yp])],
                      xpad=[float(copy[par != yp].mean()), se(copy[par != yp])], ysh_ypad=float(ys[par == yp].mean())))
      print(' | '.join(line))

if a.json:
    import json; json.dump(dict(map=a.mapx, res=RES), open(a.json, 'w'), indent=0)
