"""Second-round observables for the X-asymmetry deck (2026-10-19 evening), cached to asym_more.json.

All from the raw survey (make_asym_deck.survey_rows) and the EXTRA EOS runs, X and Y matched per sub-run by (file index,
eventId). Parity (p_x + p_y) % 2 picks the pad type under the charge; which class is "over a Y pad" is read from the data
(the class with the smaller Y +-1 share and the smaller A_x/A_y), so a board's own channel convention cannot flip it.
  boards   : per board/build, copy a1 in the X-pad and Y-pad classes (+ A_x/A_y, Y share)
  chips    : det3 June-15, Y-pad copy per X chip and per Y chip, per session, with the local pulse width
  subruns  : det3 June-15 X, per sub-run: centre-chip copy, Y-pad copy, excess/Y, Y share, pulse width w30, median A
Usage: asym_more.py   (writes asym_more.json next to this file)"""
import os, re, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import make_asym_deck as M, asym_lib as L

def fidx(files):
    return np.array([int(re.search(r'_(\d{3})_\d\d\.root$', str(f)).group(1)) for f in files])

def matched(fx, fy):
    dx, dy = L.load(fx), L.load(fy)
    kx = fidx(dx['z']['files'])[dx['fi']].astype(np.int64) * 10**9 + dx['ev']
    ky = fidx(dy['z']['files'])[dy['fi']].astype(np.int64) * 10**9 + dy['ev']
    _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
    return dx, dy, ix, iy

def w30(c):
    m = (c / c.max(1, keepdims=True)).mean(0); a = np.flatnonzero(m > 0.3 * m.max()); return float((a[-1] - a[0] + 1) * 60)

def mse(x):
    x = np.asarray(x); return (float(x.mean()), float(x.std() / np.sqrt(len(x))), int(len(x))) if len(x) > 10 else (np.nan, np.nan, int(len(x)))

def pairs(rows, sel):
    g = {}
    for r in rows:
        if sel(r): g.setdefault((r['run'], r['sub']), {})[r['view']] = r
    return [(v['X'], v['Y']) for v in g.values() if 'X' in v and 'Y' in v]

def pooled(rows, sel):
    acc = {k: [] for k in ('par', 'a1', 'axay', 'ysh', 'px', 'py', 'wy', 'exc')}
    for rx, ry in pairs(rows, sel):
        dx, dy, ix, iy = matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        Wx, Wy = dx['W'][ix], dy['W'][iy]; sy = Wy[:, 4]
        acc['par'].append((dx['p'][ix] + dy['p'][iy]) % 2); acc['a1'].append(dx['a1'][ix])
        acc['axay'].append(dx['A'][ix][:, 4] / dy['A'][iy][:, 4]); acc['ysh'].append((dy['al'][1] + dy['al'][-1])[iy])
        acc['px'].append(dx['p'][ix] % 512); acc['py'].append(dy['p'][iy] % 512); acc['wy'].append(sy)
        acc['exc'].append(np.einsum('ij,ij->i', Wx[:, 5] - Wx[:, 3], sy) / np.einsum('ij,ij->i', sy, sy))
    if not acc['par']: return None
    a = {k: np.concatenate(v) for k, v in acc.items()}
    # the Y-pad class: smaller Y share
    ypar = int(a['ysh'][a['par'] == 1].mean() < a['ysh'][a['par'] == 0].mean())
    a['ypad'] = a['par'] == ypar
    return a

def main():
    rows = M.survey_rows(); out = {}
    # per board
    out['boards'] = []
    for b, bu in sorted({(r['board'], r['build']) for r in rows}):
        a = pooled(rows, lambda r, b=b, bu=bu: r['board'] == b and r['build'] == bu)
        if a is None: continue
        out['boards'].append(dict(board=b, build=bu, n=int(len(a['a1'])), ypad=mse(a['a1'][a['ypad']]), xpad=mse(a['a1'][~a['ypad']]),
                                  axay_y=float(np.median(a['axay'][a['ypad']])), axay_x=float(np.median(a['axay'][~a['ypad']])),
                                  ysh_y=float(a['ysh'][a['ypad']].mean()), ysh_x=float(a['ysh'][~a['ypad']].mean())))
    # per chip, det3 June-15 sessions
    sess = {'6-22/23, 1000 V, FEU 3': lambda r: r['build'] == 'det3 June-15' and 'det2_det3_overnight_6-22' in r['run'],
            '6-23/24, 600 V, FEU 3': lambda r: r['build'] == 'det3 June-15' and '6-23-26' in r['run'],
            '6-27/28, 1000 V, FEU 7': lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] == '1000'}
    out['chips'] = {}
    for lab, sel in sess.items():
        a = pooled(rows, sel); d = {}
        for ax, key in (('xchip', 'px'), ('ychip', 'py')):
            ch = a[key] // 64
            d[ax] = [dict(copy=mse(a['a1'][a['ypad'] & (ch == c)]), w30=w30(a['wy'][ch == c]) if (ch == c).sum() > 80 else None,
                          n=int((ch == c).sum())) for c in range(8)]
        out['chips'][lab] = d
    # per sub-run, det3 June-15 X (+ EXTRA, X only)
    out['subruns'] = []
    for rx, ry in pairs(rows, lambda r: r['build'] == 'det3 June-15'):
        dx, dy, ix, iy = matched(rx['f'], ry['f'])
        if len(ix) < 100 or not rx['t']: continue
        Wx, Wy = dx['W'][ix], dy['W'][iy]; sy = Wy[:, 4]; par = (dx['p'][ix] + dy['p'][iy]) % 2
        ysh = (dy['al'][1] + dy['al'][-1])[iy]; ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean()); yp = par == ypar
        exc = np.einsum('ij,ij->i', Wx[:, 5] - Wx[:, 3], sy) / np.einsum('ij,ij->i', sy, sy)
        out['subruns'].append(dict(t=rx['t'], run=rx['run'], sub=rx['sub'], feu=rx['feu'], drift=rx['drift'], mesh=rx['mesh'],
                                   ac=rx['ac'], eac=rx['eac'], ypad=mse(dx['a1'][ix][yp]), exc=mse(exc[yp]), ysh=float(ysh[yp].mean()),
                                   w30=w30(dx['W'][:, 4]), medA=float(np.median(dx['A'][:, 4]))))
    for f, run, sub, feu, t, drift in M.EXTRA:
        p = f'{L.RAW}/{f}'
        if not os.path.exists(p): continue
        d = L.load(p); ac, eac, nc = M.central(d)
        out['subruns'].append(dict(t=t, run=run, sub=sub, feu=feu, drift=drift, mesh='500', ac=ac, eac=eac, ypad=None, exc=None,
                                   ysh=None, w30=w30(d['W'][:, 4]), medA=float(np.median(d['A'][:, 4]))))
    json.dump(out, open(os.path.join(HERE, 'asym_more.json'), 'w'), indent=0, default=float)
    print('wrote asym_more.json:', {k: len(v) for k, v in out.items()})

if __name__ == '__main__':
    main()
