"""Round-5 helpers (2026-10-20 day): matched X/Y events per build with pad class, from the raw survey.
ev(build) -> dict of per-event arrays: Wx, Wy (9x32 windows), px, py, ypad (bool), al_x/al_y (in-phase content
k=-4..4 as arrays (n, 9)), qx/qy (integrals per channel over samples, (n, 9)), ysh, sess (run name)."""
import os, sys, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
os.environ.setdefault('NB_RAW', '/media/dylan/data/claude_tmp/tmp/nb/raw')
import make_asym_deck as M, asym_lib as L, asym_more as A
_rows = None
def rows():
    global _rows
    if _rows is None: _rows = M.survey_rows()
    return _rows
def ev(sel, edge=True):
    E = {k: [] for k in ('Wx', 'Wy', 'px', 'py', 'par', 'alx', 'aly', 'sess')}
    for rx, ry in A.pairs(rows(), sel):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        E['Wx'].append(dx['W'][ix].astype(np.float32)); E['Wy'].append(dy['W'][iy].astype(np.float32))
        E['px'].append(dx['p'][ix]); E['py'].append(dy['p'][iy]); E['par'].append((dx['p'][ix] + dy['p'][iy]) % 2)
        E['alx'].append(np.stack([dx['al'][k][ix] for k in range(-4, 5)], 1))
        E['aly'].append(np.stack([dy['al'][k][iy] for k in range(-4, 5)], 1))
        E['sess'].append(np.array([rx['run']] * len(ix)))
    E = {k: np.concatenate(v) for k, v in E.items()}
    E['ysh'] = E['aly'][:, 3] + E['aly'][:, 5]
    yp = int(E['ysh'][E['par'] == 1].mean() < E['ysh'][E['par'] == 0].mean()); E['ypad'] = E['par'] == yp
    if edge:
        m = (E['px'] % 64 >= 4) & (E['px'] % 64 <= 59) & (E['py'] % 64 >= 4) & (E['py'] % 64 <= 59)
        E = {k: (v[m] if hasattr(v, '__len__') and len(v) == len(m) else v) for k, v in E.items()}
    E['qx'] = np.clip(E['Wx'], 0, None).sum(2); E['qy'] = np.clip(E['Wy'], 0, None).sum(2)
    return E
BUILDS = {'det3 June-15': lambda r: r['build'] == 'det3 June-15', 'det3 May': lambda r: r['build'] == 'det3 May',
          'det2': lambda r: r['build'] == 'det2', 'det1': lambda r: r['build'] == 'det1', 'det4': lambda r: r['build'] == 'det4',
          'det6': lambda r: r['build'] == 'det6', 'det7': lambda r: r['build'] == 'det7'}
