"""a1 vs stripe offset straight from an X kernel (2026-10-21).
Inputs: the production slice (gx_slice_rho5M.npz from lxplus: G (40, 9 times, 16 y, 3120 x)) or an X-only diagnosis
product (w2_production.py combine-x: G_X (40, 61 times, ny, nx)). Gaussian clouds of width sigma (x and y, y periodic
over the 1.56 mm comb box) centred over a column, kept where that column is the peak; A = (G[c+1] - G[c-1]) / G[c] in
+x order, averaged over centres, per column c (40 phases) -> vs delta(c) = ((x_c + 400) % 800) - 400 um
(sim_delta.py's definition, stripes centred at k 800 um in the kernel frame). Also the across-stripe late content:
(G[c+-2](t) - G[c+-2](0)) / G[c](t) for a cloud at the column centre, mean over columns (late minus prompt, x 1e3).
Usage: kernel_delta.py <npz> [<npz> ...]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from response.common import constants as C
from response.solver import kernels as K
P = C.PAD_PITCH_M; L = C.SUPERPERIOD_M
TW = np.array([0., 100e-9, 200e-9, 300e-9, 500e-9, 1e-6])
xc = K.pad_x(np.arange(40)); delta = ((xc * 1e6 + 400) % 800) - 400; o = np.argsort(delta)
def conv2(g, sig, dx, dy):
    nx, ny = g.shape[-1], g.shape[-2]
    kx = np.arange(nx); kx = np.minimum(kx, nx - kx) * dx; wx = np.exp(-.5 * (kx / sig) ** 2); wx /= wx.sum()
    ky = np.arange(ny); ky = np.minimum(ky, ny - ky) * dy; ly = ny * dy
    wy = sum(np.exp(-.5 * ((ky + m * ly) / sig) ** 2) for m in range(-6, 7)); wy /= wy.sum()
    F = np.fft.fft(np.fft.fft(g, axis=-1) * np.fft.fft(wx), axis=-2) * np.fft.fft(wy)[:, None]
    return np.real(np.fft.ifft(np.fft.ifft(F, axis=-2), axis=-1))
def run(path):
    d = np.load(path); G = (d['G'] if 'G' in d else d['G_X']); t = d['t']
    it = [int(np.argmin(np.abs(t - w))) for w in TW]; G = G[:, it].astype(float); t = t[it]
    x = d['x']; dx = x[1] - x[0]; dy = d['y_X'][1] - d['y_X'][0]
    res = {}
    for sig in (250e-6, 450e-6):
        Gs = conv2(G, sig, dx, dy); A = np.zeros((len(t), 40)); late = np.zeros(len(t))
        for c in range(40):
            u = ((x - xc[c] + L / 2) % L) - L / 2; m = np.abs(u) < P / 2; i0 = int(np.argmin(np.abs(u)))
            for j in range(len(t)):
                g0, gm, gp = Gs[c, j][:, m], Gs[(c - 1) % 40, j][:, m], Gs[(c + 1) % 40, j][:, m]
                pk = (g0 > gm) & (g0 > gp); A[j, c] = np.mean(((gp - gm) / g0)[pk])
                g2 = .5 * (Gs[(c - 2) % 40, j, :, i0] + Gs[(c + 2) % 40, j, :, i0] - Gs[(c - 2) % 40, 0, :, i0] - Gs[(c + 2) % 40, 0, :, i0])
                late[j] += np.mean(g2 / Gs[c, j, :, i0]) / 40
        res[f'{sig*1e6:.0f}'] = dict(t_ns=(t * 1e9).tolist(), delta=delta[o].tolist(), A=A[:, o].tolist(), late2=late.tolist(),
                                     odd_amp=[float(np.mean(np.sign(delta) * A[j])) for j in range(len(t))])
    return res
if __name__ == '__main__':
    OUT = {}
    sel = o[::3]
    for p in sys.argv[1:]:
        r = run(p); OUT[os.path.basename(p)] = r
        for sg, v in r.items():
            print(f'{os.path.basename(p)}  sigma {sg} um   delta: ' + ' '.join(f'{delta[i]:+5.0f}' for i in sel) + '   <sgn(d) A>  late+-2(x1e3)')
            for j, tn in enumerate(v['t_ns']):
                Aj = np.array(v['A'][j]); Ad = dict(zip(o, Aj))
                print(f'   t {tn:6.0f} ns        ' + ' '.join(f'{Ad[i]:+.3f}' for i in sel) + f'   {v["odd_amp"][j]:+.4f}   {1e3*v["late2"][j]:+6.1f}')
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kernel_delta.json'), 'w'), indent=1)
