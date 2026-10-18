#!/usr/bin/env python3
"""Hit-time ladder compression, sim (gun truth) vs det3 data (M3 angles) — 2026-10-18, item 3 of RECO_TRUTH.

RECONSTRUCTION_BASIS.md (nTof_x17) says per-strip hit times compress the drift-time ladder by 20-30 %
because each strip carries delayed copies of its neighbours. Here the ladder is measured the same way on
both legs and the compression is C = (fitted dt/dx) / (true dt/dx), with true dt/dx = 1 / (v tan theta):
  sim:  theta = gun angle (pencil beam), v = the sim's own drift velocity (38.0 um/ns, Stage B log)
  data: theta = M3 angle in that view,   v = 36.6 um/ns (wft fit, = wet Magboltz)
C = 1 means an unbiased ladder; C < 1 means compressed (the time span is too short for the strip span).

Per event, tilted view only: strips >= 5 sigma and >= 100 ADC, unsaturated (peak < 3500), >= 4 of them.
Hit time per strip, two estimators: leading-edge 50 % of the strip's own peak (linear interpolation, "cfd"),
and the parabolic peak time ("pk"). Ladder slope = amplitude-weighted least squares of t on strip position.
The other view must be within 3 deg of vertical (the sim's other view is exactly 0).

Usage: ladder.py OUTDIR SIMROOT PT [PT ...]    (SIMROOT/<PT>_<ang>/*_07.root, *_08.root)
"""
import os, sys, glob, json
import numpy as np

NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'):
    sys.path.insert(0, p)
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
from wft.io import FeuReader, strip_position_map

cfg = get_config('sat_det3')
PM = strip_position_map(cfg)
DEC = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
SAT, DT = 3500.0, 60.0
V_SIM, V_DATA = 38.0, 36.6
ANGS = {'thx05': ('x', 5), 'thx10': ('x', 10), 'thx15': ('x', 15), 'thx20': ('x', 20), 'thx30': ('x', 30),
        'thy05': ('y', 5), 'thy10': ('y', 10), 'thy15': ('y', 15), 'thy20': ('y', 20), 'thy30': ('y', 30)}


def strip_times(w):
    """(t_cfd, t_pk) per strip row of w [n, ns], in ns from the window start."""
    k = w.argmax(axis=1)
    n, ns = w.shape
    pk = w[np.arange(n), k]
    t_pk = np.full(n, np.nan); t_cfd = np.full(n, np.nan)
    for i in range(n):
        ki = k[i]
        if 0 < ki < ns - 1:
            a, b, c = w[i, ki - 1], w[i, ki], w[i, ki + 1]
            den = a - 2 * b + c
            t_pk[i] = (ki + (0.5 * (a - c) / den if den < 0 else 0.0)) * DT
        half = 0.5 * pk[i]
        j = ki
        while j > 0 and w[i, j - 1] >= half:
            j -= 1
        if j > 0:
            y0, y1 = w[i, j - 1], w[i, j]
            t_cfd[i] = (j - 1 + (half - y0) / (y1 - y0)) * DT
    return t_cfd, t_pk


def ladder_events(fn, feu, want=None):
    """{eid: dict(slope_cfd, slope_pk, n, span)} for one decoded file."""
    pos = PM[feu]
    r = FeuReader(fn)
    noise = np.where(r.noise > 0, r.noise, np.inf)
    ids = set(int(x) for x in r.event_ids) if want is None else (set(int(x) for x in r.event_ids) & want)
    out = {}
    for e, _t, w in r.iter_events(ids):
        amp = w.max(axis=1)
        sel = (amp / noise >= 5) & (amp >= 100) & np.isfinite(pos)
        if amp[sel].max(initial=0) >= SAT:       # any saturated strip: skip the event (same on both legs)
            continue
        if sel.sum() < 4:
            continue
        idx = np.nonzero(sel)[0]
        tc, tp = strip_times(w[idx])
        p, a = pos[idx], amp[idx]
        res = {'n': int(len(idx)), 'span': float(p.max() - p.min())}
        for nm, t in (('cfd', tc), ('pk', tp)):
            g = np.isfinite(t)
            if g.sum() < 4:
                res[nm] = res[nm + '_xt'] = np.nan; continue
            ww = a[g]; x = p[g]; y = t[g]
            xm = (ww * x).sum() / ww.sum(); ym = (ww * y).sum() / ww.sum()
            res[nm] = float((ww * (x - xm) * (y - ym)).sum() / (ww * (x - xm) ** 2).sum())
            # the reverse regression (position on time, as a track fit x(z) does), stored as dt/dx
            b = (ww * (x - xm) * (y - ym)).sum() / (ww * (y - ym) ** 2).sum()
            res[nm + '_xt'] = float(1.0 / b) if b != 0 else np.nan
        out[int(e)] = res
    return out


def collect(d, view):
    feu = 7 if view == 'x' else 8
    out = {}
    for f in sorted(glob.glob(f'{d}/*_{feu:02d}.root')):
        out.update({(f, k): v for k, v in ladder_events(f, feu).items()})
    return list(out.values())


def summarise(rows, tan_true, v):
    """Median of the per-event C = |fitted dt/dx| * v * tan(theta). For data, tan(theta) is the event's own
    M3 angle (r['th']), so the +-3 deg window does not bias C through the 1/tan of the true slope."""
    o = {'n': len(rows)}
    for nm in ('cfd', 'pk', 'cfd_xt', 'pk_xt'):
        c = np.array([abs(r[nm]) * v * 1e-3 * (abs(np.tan(np.radians(r['th']))) if 'th' in r else tan_true)
                      for r in rows], float)
        c = c[np.isfinite(c)]
        if len(c) < 10:
            o[nm] = None; continue
        rng = np.random.default_rng(1)
        bs = [np.median(rng.choice(c, len(c))) for _ in range(200)]
        q = np.percentile(c, [25, 75])
        o[nm] = {'C': float(np.median(c)), 'C_err': float(np.std(bs)), 'q25': float(q[0]), 'q75': float(q[1])}
    return o


def main():
    outdir, simroot, pts = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(outdir, exist_ok=True)
    res = {'sim': {}, 'data': {}}
    # ---- data, M3 angles (same recipe as scratchpad/2026-10-17/xy_extent_m3.py)
    from M3RefTracking import M3RefTracking
    rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
    x0, y0, ev = map(np.asarray, rays.get_xy_positions(702.0))
    x1, y1, _ = map(np.asarray, rays.get_xy_positions(802.0))
    ax = dict(zip(ev.astype(int), np.degrees(np.arctan((x1 - x0) / 100))))
    ay = dict(zip(ev.astype(int), np.degrees(np.arctan((y1 - y0) / 100))))
    th = {'x': {k: ay[k] for k in ay}, 'y': {k: -ax[k] for k in ax}}   # X view <- M3 y, Y view <- -M3 x
    for view, feu in (('x', 7), ('y', 8)):
        other = 'y' if view == 'x' else 'x'
        rows = {}
        for f in sorted(glob.glob(f'{DEC}/*_{feu:02d}.root')):
            for e, r in ladder_events(f, feu, want=set(th[view])).items():
                if e in th[other] and np.isfinite(th[view][e]) and np.isfinite(th[other][e]):
                    r['th'] = th[view][e]; r['tho'] = th[other][e]
                    rows[(f, e)] = r
        rows = list(rows.values())
        for a in (5, 10, 15, 20, 30):
            sel = [r for r in rows if abs(abs(r['th']) - a) < (2 if a == 5 else 3) and abs(r['tho']) < 3]
            res['data'][f'th{view}{a:02d}'] = summarise(sel, np.tan(np.radians(a)), V_DATA)
            print('data', view, a, json.dumps(res['data'][f'th{view}{a:02d}']), flush=True)
    # ---- sim
    for pt in pts:
        res['sim'][pt] = {}
        for ang, (view, a) in ANGS.items():
            d = f'{simroot}/{pt}_{ang}'
            if not glob.glob(f'{d}/*_07.root'):
                continue
            rows = collect(d, view)
            res['sim'][pt][ang] = summarise(rows, np.tan(np.radians(a)), V_SIM)
            print(pt, ang, json.dumps(res['sim'][pt][ang]), flush=True)
    with open(f'{outdir}/ladder.json', 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
