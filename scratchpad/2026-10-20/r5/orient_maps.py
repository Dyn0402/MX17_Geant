"""Pin the field to the resistive strips (2026-10-20 evening). det0 (4-1 He/ethane HV scan) and det1 (Jan) have their
strips along x (X-slow, strip_orient2.py); every other board along y. Strip-induced => the asymmetry field lives in
the view that measures ACROSS the strips: Y on det0/det1, X elsewhere. Board/readout-induced => X on all.
Per board: 8 x 8 cell maps (matched X true column x Y true row) of a1X and a1Y, constant-model chi2, weighted cell rms.
Usage: orient_maps.py"""
import os, glob, json, numpy as np
import r5lib as R, asym_more as A
ORIENT = '/media/dylan/data/claude_tmp/tmp/nb/raw/orient'
def from_pairs(pairs):
    E = {k: [] for k in ('px', 'py', 'ax', 'ay', 'shx', 'shy')}
    for fx, fy in pairs:
        try: dx, dy, ix, iy = A.matched(fx, fy)
        except Exception as e: print('skip', fx, e); continue
        if len(ix) < 20: continue
        E['px'].append(dx['p'][ix] % 512); E['py'].append(dy['p'][iy] % 512)
        E['ax'].append(dx['a1'][ix]); E['ay'].append(dy['a1'][iy])
        E['shx'].append((dx['al'][1] + dx['al'][-1])[ix]); E['shy'].append((dy['al'][1] + dy['al'][-1])[iy])
    return {k: np.concatenate(v) for k, v in E.items()} if E['px'] else None
def from_build(b):
    return from_pairs([(rx['f'], ry['f']) for rx, ry in A.pairs(R.rows(), lambda r: r['build'] == b)])
def cmap(E, key, nc=8):
    cx = 64 * (E['px'] // 64) + 63 - E['px'] % 64; cy = 64 * (E['py'] // 64) + 63 - E['py'] % 64
    i = cx * nc // 512; j = cy * nc // 512; M = []; G = np.full((nc, nc), np.nan)
    for a in range(nc):
        for b in range(nc):
            m = (i == a) & (j == b)
            if m.sum() >= 25: M.append((E[key][m].mean(), E[key][m].std() / np.sqrt(m.sum()))); G[a, b] = M[-1][0]
    M = np.array(M); w = 1 / M[:, 1] ** 2; m0 = (M[:, 0] * w).sum() / w.sum()
    return float((w * (M[:, 0] - m0) ** 2).sum()), len(M) - 1, float(np.sqrt(np.average((M[:, 0] - m0) ** 2, weights=w))), G
def main():
    S = {}
    pairs = sorted(glob.glob(f'{ORIENT}/d0all_*_x.npz'))
    S['det0 (strips along x)'] = from_pairs([(f, f[:-6] + '_y.npz') for f in pairs if os.path.exists(f[:-6] + '_y.npz')])
    S['det1 (strips along x)'] = from_build('det1')
    S['det1 all 5 files (strips along x)'] = from_pairs([(f'{ORIENT}/d1all_x.npz', f'{ORIENT}/d1all_y.npz')])
    q = sorted(glob.glob(f'{ORIENT}/d0ar/d0ar_*_x.npz'))
    S['det0 Ar/iso 4-7 (strips along x)'] = from_pairs([(f, f[:-6] + '_y.npz') for f in q if os.path.exists(f[:-6] + '_y.npz')])
    if os.path.exists(f'{ORIENT}/a24_x.npz'): S['Apr-24 bulk ("horizontal", broken)'] = from_pairs([(f'{ORIENT}/a24_x.npz', f'{ORIENT}/a24_y.npz')])
    for b in ('det3 June-15', 'det2', 'det3 May', 'det4 May'):
        if b == 'det4 May': S[b + ' (along y)'] = from_pairs([(f'{ORIENT}/{t}_x.npz', f'{ORIENT}/{t}_y.npz') for t in ('d4m_490', 'd4m_510')])
        else: S[b + ' (along y)'] = from_build(b)
    OUT = {}
    print(f'{"board":38s} {"n":>6s} | X: chi2/ndf   rms   mean | Y: chi2/ndf   rms   mean')
    for nm, E in S.items():
        if E is None: print(f'{nm:38s} too few matched events'); continue
        cX, nX, rX, GX = cmap(E, 'ax'); cY, nY, rY, GY = cmap(E, 'ay')
        OUT[nm] = dict(n=int(len(E['ax'])), X=[cX, nX, rX, float(E['ax'].mean())], Y=[cY, nY, rY, float(E['ay'].mean())], GX=GX.tolist(), GY=GY.tolist())
        print(f'{nm:38s} {len(E["ax"]):6d} | X: {cX:6.0f}/{nX:<3d} {rX:.3f} {E["ax"].mean():+.3f} | Y: {cY:6.0f}/{nY:<3d} {rY:.3f} {E["ay"].mean():+.3f}')
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'orient_maps.json'), 'w'), indent=1)
if __name__ == '__main__':
    main()
