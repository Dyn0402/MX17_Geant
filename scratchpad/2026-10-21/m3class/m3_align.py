"""Step 1 of the M3 pad-class test (2026-10-21): align M3 to each view and measure the pointing resolution.
Tagged dumps (dump_tagged.py; position-ordered +-5 windows, M3 x, y, ax, ay at the detector plane).
Per view: centroid of the 3 central channels (mm, from the dump's position map), affine fit
pos = a + b*m3x + c*m3y + d*ax (robust, 3 passes of 4-sigma clipping), residual core sigma.
Usage: m3_align.py [TAGGED_DIR]"""
import sys, os, json, numpy as np
T = sys.argv[1] if len(sys.argv) > 1 else '/media/dylan/data/claude_tmp/tmp/nb/tagged'
TAGS = ('o22', 'o23', 'g_det3_wknd', 'g_det2', 'g_det4', 'g_det6_long', 'g_det7_long')
def load(tag):
    z = np.load(f'{T}/{tag}.npz', allow_pickle=True); o = {}
    for v in 'xy':
        W = np.nan_to_num(z[f'{v}_wf']); ch = z[f'{v}_ch']; pos = z[f'{v}_pos']; A = W.max(2)
        a3 = np.clip(A[:, 4:7], 0, None); p3 = np.where(ch[:, 4:7] >= 0, pos[np.maximum(ch[:, 4:7], 0)], np.nan)
        cen = (a3 * p3).sum(1) / a3.sum(1)
        o[v] = dict(ev=z[f'{v}_event'], W=W, A=A, ch=ch, pos=pos, cen=cen, pk=pos[z[f'{v}_peak_ch']], m3=z[f'{v}_m3'])
    return o
def fit(m3, y, coord):
    X = np.column_stack([np.ones(len(y)), m3[:, 0], m3[:, 1], m3[:, 2 + coord]])
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    for _ in range(4):
        c, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None); r = y - X @ c
        s = 1.4826 * np.nanmedian(np.abs(r[ok] - np.nanmedian(r[ok]))); ok = ok & (np.abs(r) < 4 * max(s, .05))
    return c, r, s, ok
if __name__ == '__main__':
    OUT = {}
    for tag in TAGS:
        if not os.path.exists(f'{T}/{tag}.npz'): continue
        d = load(tag); OUT[tag] = {}
        for i, v in enumerate('xy'):
            e = d[v]; good = (e['A'][:, 5] > 200) & (e['A'][:, 5] < 3500)
            c, r, s, ok = fit(e['m3'][good], e['cen'][good], i)
            core = r[ok]; OUT[tag][v] = dict(coef=c.tolist(), sigma_mad=float(s), rms_core=float(core.std()), n=int(ok.sum()), frac_out=float(1 - ok.mean()))
            print(f'{tag:12s} {v}: pos = {c[0]:+8.2f} {c[1]:+.4f}*x {c[2]:+.4f}*y {c[3]:+.4f}*a   sigma(MAD) {s*1000:5.0f} um  rms {core.std()*1000:5.0f} um  n {ok.sum()}  out {1-ok.mean():.2%}')
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'm3_align.json'), 'w'), indent=1)
