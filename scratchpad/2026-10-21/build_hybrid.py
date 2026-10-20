"""Hybrid kernel product for the digitizer (2026-10-21, det3 X asymmetry seventh round): the production W2 product
(Y kernels, meta, schema) with G_X replaced by an X-only diagnosis product (w2_production.py combine-x, e.g. conducting
ESL gaps at nx = 1560), linearly upsampled (periodic) onto the production 10 um x grid. The Y view keeps the
production (insulating-gap) kernels: an X-view test only. Usage: build_hybrid.py <production npz> <gx npz> <out npz>"""
import sys, json, numpy as np
prod, gxp, out = sys.argv[1:4]
g = np.load(gxp); mg = json.loads(str(g['meta']))
with np.load(prod) as d:
    D = {k: d[k] for k in d.files}
x0, x1 = D['x'], g['x']; L = x0[1] * len(x0)
assert abs(x1[1] * len(x1) - L) < 1e-12 and np.allclose(D['t'], g['t']) and np.allclose(D['y_X'], g['y_X'])
assert np.allclose(D['x_cols'], g['x_cols'])
GX = g['G_X'].astype(np.float32); n1 = GX.shape[-1]; r = len(x0) // n1
assert r * n1 == len(x0)
xi = x0 / x1[1]; i0 = np.floor(xi).astype(int) % n1; w = (xi - np.floor(xi)).astype(np.float32)
D['G_X'] = ((1 - w) * GX[..., i0] + w * GX[..., (i0 + 1) % n1]).astype(np.float32)
m = json.loads(str(D['meta'])); m['hybrid'] = dict(G_X_from=gxp, G_X_meta=mg, note='Y kernels production; X from diagnosis product, linear x upsampling')
D['meta'] = json.dumps(m)
np.savez(out, **D)
print('wrote', out, D['G_X'].shape, 'gap_frac', mg.get('esl_gap_frac'), 'X total prompt', float(D['G_X'][:, 0].sum(0).mean()))
