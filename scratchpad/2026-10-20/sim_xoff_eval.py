"""Evaluate the X-view-offset sim arms (condor 4505153: v5a recipe, th00, --x-view-offset-um 0/-100/-200/-300) through
det3's analysis, against the det3 June-15 6-27/28 1000 V session (the v5a target: FEU 7/8, 490/1000 V) and all det3
June-15 data. Inputs: raw_asym.py npz of each arm's X (FEU 7) and Y (FEU 8), in ~/x17/nb_asym_raw/simxoff/.
A plain copy a (0..0.2) is also injected into each arm's raw X windows. Absolute Y-share bins (the full sim has the
resistive spreading and shaping). Usage: sim_xoff_eval.py"""
import os, sys, json, glob, re, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import asym_lib as L, make_asym_deck as M, asym_more as A
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9]); SRC = os.path.expanduser('~/x17/nb_asym_raw/simxoff')
def classes(Wx, px0, py, ysh):
    A_ = Wx.max(2); sh = A_[:, 3:6].argmax(1) - 1; c0 = 4 + sh; ii = np.arange(len(Wx))
    c = Wx[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
    d = (np.einsum('ij,ij->i', Wx[ii, c0 + 1], c) - np.einsum('ij,ij->i', Wx[ii, c0 - 1], c)) / n2
    px = px0 + sh; par = (px + py) % 2; ok = (px % 64 >= 4) & (px % 64 <= 59)
    yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean()); b = np.digitize(ysh, ED) - 1; r = []
    for k in range(len(ED) - 1):
        for cl in (True, False):
            m = ok & (b == k) & (yp == cl); r.append((d[m].mean() if m.sum() > 20 else np.nan, d[m].std() / np.sqrt(max(m.sum(), 1))))
        mb = ok & (b == k); r.append((yp[mb].mean() if mb.sum() > 20 else np.nan, 0.02))
    return np.array(r), float(np.median(ysh)), float(d[ok].mean()), int(ok.sum())
def from_pair(fx, fy):
    dx, dy = L.load(fx), L.load(fy)
    kx = dx['fi'].astype(np.int64) * 10 ** 9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10 ** 9 + dy['ev']
    _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
    return dx['W'][ix].astype(float), dx['p'][ix] % 512, dy['p'][iy] % 512, (dy['al'][1] + dy['al'][-1])[iy]
def data_sel(sel):
    rows = M.survey_rows(); out = [[], [], [], []]
    for rx, ry in A.pairs(rows, sel):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        for o, v in zip(out, (dx['W'][ix].astype(float), dx['p'][ix] % 512, dy['p'][iy] % 512, (dy['al'][1] + dy['al'][-1])[iy])): o.append(v)
    return [np.concatenate(o) for o in out]
REF = {'6-27/28 1000 V': data_sel(lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] == '1000'),
       'all det3 June-15': data_sel(lambda r: r['build'] == 'det3 June-15')}
if os.environ.get('XOFF_SESSIONS'):   # the two narrow-cluster sessions on their own (fwd_session_fit.py SESS)
    REF['6-22/23 1000 V'] = data_sel(lambda r: r['build'] == 'det3 June-15' and 'det2_det3_overnight_6-22' in r['run'])
    REF['6-23/24 600 V'] = data_sel(lambda r: r['build'] == 'det3 June-15' and '6-23-26' in r['run'])
REFS = {k: classes(*v) for k, v in REF.items()}
fmt = lambda t, j: ' '.join('  -  ' if not np.isfinite(v) else f'{v:+.2f}' for v in t[j::3, 0])
for k, (t, med, a, n) in REFS.items():
    print(f'DATA {k:18s} n={n:5d} median Ysh {med:.2f} a1 {a:+.3f} | Ypad [{fmt(t,0)}] Xpad [{fmt(t,1)}] fY [{fmt(t,2)}]')
OUT = dict(data={k: dict(t=v[0].tolist(), med=v[1], a1=v[2], n=v[3]) for k, v in REFS.items()}, sim=[])
# arms: xo<off>_<ang>_s<seed> (run_xoff_job.sh), pooled per offset over angles and seeds; old single-angle xoff* kept apart
tags = sorted({re.sub(r'_[xy]\.npz$', '', os.path.basename(f)) for f in glob.glob(f'{SRC}/*_x.npz')})
groups = {}
for t in tags:
    if os.path.exists(f'{SRC}/{t}_y.npz'): groups.setdefault(t.split('_')[0] if t.startswith('xo') and not t.startswith('xoff') else t, []).append(t)
for arm, members in sorted(groups.items()):
    # steep arms (thx/thy 15-30 in the measured view) keep no isolated cluster and are written empty: skip them
    members = [t for t in members if re.search(os.environ.get('XOFF_ARMS', '.'), t) and all(np.load(f'{SRC}/{t}_{v}.npz')['wf'].ndim == 3 for v in 'xy')]
    if not members: continue
    parts = [from_pair(f'{SRC}/{t}_x.npz', f'{SRC}/{t}_y.npz') for t in members]
    Wx, px, py, ysh = [np.concatenate([q[i] for q in parts]) for i in range(4)]
    for a in [float(c) for c in os.environ.get('XOFF_COPIES', '0,.05,.1,.15,.2').split(',')]:
        W = Wx.copy()
        if a: W[:, 1:] += a * Wx[:, :-1]
        t, med, a1, n = classes(W, px, py, ysh)
        chi = {}
        for k, (td, *_r) in REFS.items():
            sel = np.r_[[3 * b for b in range(7)], [3 * b + 1 for b in range(7)]]
            ok = np.isfinite(t[sel, 0]) & np.isfinite(td[sel, 0])
            chi[k] = float((((t[sel][ok, 0] - td[sel][ok, 0]) ** 2) / (t[sel][ok, 1] ** 2 + td[sel][ok, 1] ** 2 + .02 ** 2)).sum())
        OUT['sim'].append(dict(arm=arm, copy=a, t=t.tolist(), med=med, a1=a1, n=n, chi=chi))
        print(f'SIM {arm:10s} ({len(members)} files) +copy {a:.2f} n={n:5d} median Ysh {med:.2f} a1 {a1:+.3f} | Ypad [{fmt(t,0)}] Xpad [{fmt(t,1)}] fY [{fmt(t,2)}]'
              ' | chi2 ' + ', '.join(f'{k} {v:.0f}' for k, v in chi.items()))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.environ.get('XOFF_OUT', 'sim_xoff_eval.json')), 'w'), default=float)
