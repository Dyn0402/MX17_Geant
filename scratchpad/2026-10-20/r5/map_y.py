"""Does the Y view carry its own smooth asymmetry field? (2026-10-20) Same 8 x 8 cells as map_fit.py (position from the
matched X column and Y row), a1Y = alphaY(+1) - alphaY(-1). Reports the constant-model chi2 of the cell map for X and
Y side by side (X structure vs Y structure), and the 4 x 4 Y map. Stripes along y can only make X structure; a vector
misregistration (via/drill or inner-layer offsets) would make both. Usage: map_y.py"""
import json, os, numpy as np
import r5lib as R
OUT = {}
def cellchi(cx, cy, v, nc):
    i = cx * nc // 512; j = cy * nc // 512; M = []
    for a in range(nc):
        for b in range(nc):
            m = (i == a) & (j == b)
            if m.sum() >= 25: M.append((v[m].mean(), v[m].std() / np.sqrt(m.sum())))
    M = np.array(M); w = 1 / M[:, 1] ** 2; m0 = (M[:, 0] * w).sum() / w.sum()
    return float((w * (M[:, 0] - m0) ** 2).sum()), len(M) - 1, float(np.sqrt(np.average((M[:, 0] - m0) ** 2, weights=w)))
for b in ('det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1'):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    ax = E['alx'][:, 5] - E['alx'][:, 3]; ay = E['aly'][:, 5] - E['aly'][:, 3]
    cX, nX, rX = cellchi(cx, cy, ax, 8); cY, nY, rY = cellchi(cx, cy, ay, 8)
    OUT[b] = dict(X=(cX, nX, rX), Y=(cY, nY, rY), meanX=float(ax.mean()), meanY=float(ay.mean()))
    print(f'{b:13s} X: mean {ax.mean():+.3f}, cell scatter rms {rX:.3f}, chi2 {cX:6.1f}/{nX} | Y: mean {ay.mean():+.3f}, rms {rY:.3f}, chi2 {cY:6.1f}/{nY}')
    G = [[ay[(cx // 128 == i) & (cy // 128 == j)].mean() for j in range(4)] for i in range(4)]
    print('      Y 4x4: ' + ' | '.join(' '.join(f'{v:+.2f}' for v in r) for r in G))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'map_y.json'), 'w'), indent=1)
