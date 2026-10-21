"""Is det3 June-15's prompt-wide, time-flat X sharing an instantaneous symmetric pad coupling? (2026-10-21 night; FOLLOWUP §8)
A coupling a per neighbour (neighbour = a x central, at every instant) predicts:
  (1) r(f) = (W+1 + W-1)/W_c flat down to the onset of the central pulse (5-10 % of its rise); physics (diffusion,
      spreading) makes it grow through the pulse, because early electrons come from near the mesh and diffuse least;
  (2) a floor under s1 = alpha+1 + alpha-1 at ~2a: no event shares less;
  (3) the same floor in every session (gain, drift field, time), being a property of the board, not of the clouds.
Also the Y view with the same cut (does the Y view of the same board show it?).
Near-vertical cut as prompt_share.py SEL=vert (|dt_nb| < DTN_MAX ns, FWHM >= its FW_PCT percentile).
Usage: coupling_checks.py [out.json]   (NB_RAW defaults to ~/x17/nb_asym_raw)"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NB_RAW', os.path.expanduser('~/x17/nb_asym_raw'))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import prompt_share as PS
import asym_lib as L
RAW = os.environ['NB_RAW']; HW = 4; DT = 60.
DTN_MAX = float(os.environ.get('DTN_MAX', 80.)); FW_PCT = float(os.environ.get('FW_PCT', 25.))
FR = (.05, .1, .15, .25, .5, .75); BINS = ((.2, .45), (.45, .7))

def events(files):
    W, s1 = [], []
    for f in files:
        d = L.load(f); p = d['p'] % 512; m = (p % 64 >= 4) & (p % 64 <= 59)
        W.append(d['W'][m]); s1.append(d['al'][1][m] + d['al'][-1][m])
    W = np.concatenate(W).astype(float); s1 = np.concatenate(s1); c = W[:, HW]; pk = c.argmax(1)
    ts = np.arange(W.shape[2]); pos = lambda w: np.clip(w, 0, None)
    tc = lambda w: (pos(w) * ts).sum(1) / np.maximum(pos(w).sum(1), 1e-9)
    dtn = (tc(W[:, HW + 1]) - tc(W[:, HW - 1])) * DT
    fwhm = (PS.cross(c, pk, .5, False) - PS.cross(c, pk, .5)) * DT
    vert = (np.abs(dtn) < DTN_MAX) & (fwhm >= np.nanpercentile(fwhm, FW_PCT))
    return W[vert], s1[vert], int(len(s1))

def profile(W, s1):
    c = W[:, HW]; nb = W[:, HW - 1] + W[:, HW + 1]; pk = c.argmax(1)
    pts = {f'r{int(f * 100):02d}': PS.cross(c, pk, f) for f in FR}; pts['peak'] = pk.astype(float)
    pts['f75'] = PS.cross(c, pk, .75, False); pts['f50'] = PS.cross(c, pk, .5, False)
    out = {}
    for lo, hi in BINS:
        m = (s1 >= lo) & (s1 < hi); row = {'n': int(m.sum())}
        for k, t in pts.items():
            nu, de = PS.at(nb, t), PS.at(c, t); g = m & np.isfinite(nu) & np.isfinite(de)
            if g.sum() > 20:
                ii = np.flatnonzero(g); bs = PS.RNG.integers(0, len(ii), (200, len(ii)))
                row[k] = [float(nu[g].sum() / de[g].sum()), float(np.std(nu[ii][bs].sum(1) / de[ii][bs].sum(1)))]
        out[f'{lo}-{hi}'] = row
    return out

def floor(s1):
    q = np.percentile(s1, [1, 5, 10, 25, 50]); return dict(p1=q[0], p5=q[1], p10=q[2], p25=q[3], p50=q[4],
                                                          f_lt02=float((s1 < .2).mean()), f_lt03=float((s1 < .3).mean()))

def show(name, n, nsel, fl, pr):
    print(f'\n{name}  ({nsel} of {n} near-vertical)  s1 p1/p5/p10/p25/p50 = '
          f'{fl["p1"]:.2f}/{fl["p5"]:.2f}/{fl["p10"]:.2f}/{fl["p25"]:.2f}/{fl["p50"]:.2f}  f(s1<0.2) {fl["f_lt02"]:.3f}  f(<0.3) {fl["f_lt03"]:.3f}')
    keys = [f'r{int(f * 100):02d}' for f in FR] + ['peak', 'f75', 'f50']
    print('   s1 bin      n   ' + '  '.join(f'{k:>11s}' for k in keys))
    for b, row in pr.items():
        print(f'   {b:9s} {row["n"]:5d}  ' + '  '.join(f'{row[k][0]:5.3f}±{row[k][1]:.3f}' if k in row else '     -     ' for k in keys))

R = {'views': {}, 'sessions': {}}
SETS = {'det3 June-15': 'd3jun', 'det3 May': 'd3may', 'det2': 'd2'}
for view in ('x', 'y'):
    for name, tag in SETS.items():
        W, s1, n = events([f'{RAW}/{tag}_{view}.npz']); fl = floor(s1); pr = profile(W, s1)
        R['views'][f'{name} {view.upper()}'] = dict(n=n, n_sel=len(s1), floor=fl, profile=pr); show(f'{name} {view.upper()}', n, len(s1), fl, pr)
    W, s1, n = events(sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_{view}.npz'))) if view == 'x' else (None, None, 0)
    if W is not None:
        fl = floor(s1); pr = profile(W, s1); R['views'][f'sim glg0 {view.upper()}'] = dict(n=n, n_sel=len(s1), floor=fl, profile=pr)
        show(f'sim glg0 {view.upper()}', n, len(s1), fl, pr)

# per session (survey), X view: s1 floor and the prompt ratio at 10 % and 25 % in the 0.2-0.45 bin
import r5lib
print('\nPer session, X view (survey), near-vertical: s1 p5 / p25 / p50, r10 and r25 in s1 0.2-0.45')
for r in r5lib.rows():
    if r['build'] not in SETS or r['view'] != 'X': continue
    try: W, s1, n = events([r['f']])
    except Exception as e: continue
    if len(s1) < 60: continue
    fl = floor(s1); pr = profile(W, s1)['0.2-0.45']
    key = f"{r['build']} | {r['run'][-14:]} {r['sub'][:24]} | mesh {r['mesh']} drift {r['drift']}"
    R['sessions'][key] = dict(n=len(s1), floor=fl, r10=pr.get('r10'), r25=pr.get('r25'), a=r['a'])
    g = lambda k: f'{pr[k][0]:.3f}±{pr[k][1]:.3f}' if k in pr else '     -     '
    print(f'  {key:80s} n {len(s1):4d}  s1 {fl["p5"]:.2f}/{fl["p25"]:.2f}/{fl["p50"]:.2f}  r10 {g("r10")}  r25 {g("r25")}  a {r["a"]:+.3f}')
json.dump(R, open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'coupling_checks.json'), 'w'), indent=1, default=float)
