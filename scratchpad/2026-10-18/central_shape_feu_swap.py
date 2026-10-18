"""Is the X/Y difference of the central-strip pulse shape (FRESH_LOOK §21.8: data X keeps a slow late tail, data Y
falls faster after the peak and undershoots ~-0.12) a property of the VIEW (chamber) or of the FEU (electronics)?
det3 ran with X/Y on FEU 7/8 (6-27 saturday scan, sat_det3) and on FEU 3/4 (6-22 long run, g_det3).
Same selection on both: unsaturated peak strip (< 3500 ADC), baseline = mean of the first 3 samples, <= NMAX strips >= 5 sigma in each view (near-vertical
proxy, no M3 needed), events with both views. Event-aligned average of the peak strip, normalised to its peak,
aligned at the leading 50 % crossing (as xtalk_fit.py). Optional: a sim point dir as the 3rd leg.
Usage: central_shape_feu_swap.py [SIM_TH00_DIR]

STATUS 2026-10-18: NOT USABLE. The multiplicity proxy (<= 5 or 8 strips) and the first-3-sample baseline give late
plateaus of 0.2-0.4 of the peak in the data, an artefact (the M3-selected sample, dump_pm3_m3.py, is clean). A real FEU
swap test needs the 6-22 run's own pedestals and an M3 vertical selection (its m3_tracking_root exists on yoga)."""
import os, sys, glob
import numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis', f'{NT}/mx17_sim_wft'): sys.path.insert(0, p)
from qa_config import get_config
from wft.io import FeuReader, strip_position_map
NMAX, MAXEV, DT = 8, 6000, 60.0
grid = np.arange(-4, 20, 0.25)

def cross(w, frac, side):
    # as 2026-10-12/wf_shape.py (which reads a file at import, so it is not imported)
    ipk = int(np.argmax(w)); lvl = frac * w[ipk]
    if side < 0:
        b = np.flatnonzero(w[:ipk + 1] < lvl)
        if not len(b): return np.nan
        i = b[-1]; j = i + 1
    else:
        b = np.flatnonzero(w[ipk:] < lvl)
        if not len(b): return np.nan
        j = ipk + b[0]; i = j - 1
    return i + (lvl - w[i]) / (w[j] - w[i]) if w[j] != w[i] else float(i)

def per_view(d, feu, pm):
    pos = pm[feu]; valid = ~np.isnan(pos); res = {}
    for path in sorted(glob.glob(f'{d}/*_{feu:02d}.root')):
        r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
        ids = list(int(e) for e in r.event_ids)[:MAXEV]
        for e, _t, w in r.iter_events(set(ids)):
            # own baseline per channel (mean of the first 3 samples): the 6-22 files keep a pedestal
            w = w - w[:, :3].mean(axis=1, keepdims=True)
            amp = w.max(axis=1); over = valid & (amp / noise >= 5)
            n = int(over.sum())
            if n == 0 or n > NMAX: continue
            c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])])
            if amp[p] >= 3500: continue
            res[int(e)] = w[p] / amp[p]
        if len(res) >= MAXEV: break
    return res

def avg(ws):
    C = []
    for w in ws:
        l = cross(w, .5, -1)
        if np.isnan(l): continue
        C.append(np.interp(grid + l, np.arange(len(w)), w))
    return np.mean(C, 0), len(C)

legs = []
for key in ('sat_det3', 'g_det3'):
    cfg = get_config(key); pm = strip_position_map(cfg)
    legs.append((f'{key} (X/Y FEU {cfg.MX17_FEU_X}/{cfg.MX17_FEU_Y})', f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root',
                 cfg.MX17_FEU_X, cfg.MX17_FEU_Y, pm))
if len(sys.argv) > 1:
    cfg = get_config('sat_det3'); legs.append((f'sim {os.path.basename(sys.argv[1].rstrip("/"))}', sys.argv[1], 7, 8,
                                                 strip_position_map(cfg)))
print(f'near-vertical proxy: <= {NMAX} strips >= 5 sigma per view, unsaturated, both views; t = 0 at the leading 50 %')
print('t[ns]                          ' + ' '.join(f'{g*DT:6.0f}' for g in grid[::4][:24]))
for name, d, fx, fy, pm in legs:
    X = per_view(d, fx, pm); Y = per_view(d, fy, pm); both = X.keys() & Y.keys()
    for v, S in (('X', X), ('Y', Y)):
        c, n = avg([S[e] for e in both])
        print(f'{name[:24]:24s} {v} n={n:5d} ' + ' '.join(f'{x:6.3f}' for x in c[::4][:24]))
