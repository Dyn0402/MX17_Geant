"""Last-ditch tests (2026-10-23): per-event table from the wft bench caches (raw 32-sample windows + M3 truth).
Per view the central channel is the max-amplitude channel of the window; alpha_k = in-phase content of the channel at
position index c+k (projection on the central pulse, as asym_lib), so a1 = alpha_+1 - alpha_-1 is in POSITION order
(+ = larger pos_mm), s1 = alpha_+1 + alpha_-1. Cuts as asym_lib.good() (200 < A_c < 3500, peak sample 4..22, flat
pre-samples) and the central channel at least 1 from the window edge. Positions: pos_mm of the central channel, the
3-channel amplitude centroid, M3 ref_x/ref_y/tan_x/tan_y (raw strip frame, pos_mm = 0.78 x column).
Usage: extract.py   (nTof_x17 venv) -> tables/<tag>.npz"""
import os, pickle, numpy as np
A = '/media/dylan/data/x17/cosmic_bench/Analysis'
RUNS = {'d3_622': ('det3 June-15', 'mx17_det2_det3_overnight_6-22-26/long_run/mx17_3'),
        'd3_627s': ('det3 June-15', 'mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check/mx17_3'),
        'd3_627': ('det3 June-15', 'mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V/mx17_3'),
        'd2_622': ('det2', 'mx17_det2_det3_overnight_6-22-26/long_run/mx17_2'),
        'd2_622l': ('det2', 'mx17_det2_det3_overnight_6-22-26/longer_run/mx17_2'),
        'd4_624': ('det4', 'mx17_det4_day_6-24-26/long_run/mx17_4'),
        'd6_626': ('det6', 'mx17_det6_det7_overnight_6-26-26/long_run/mx17_6'),
        'd7_626': ('det7', 'mx17_det6_det7_overnight_6-26-26/long_run/mx17_7')}
HERE = os.path.dirname(os.path.abspath(__file__)); os.makedirs(f'{HERE}/tables', exist_ok=True)

def view(w):
    W = np.asarray(w['W'], float); pos = np.asarray(w['pos'], float); o = np.argsort(pos); W, pos = W[o], pos[o]
    if len(pos) < 3 or np.any(np.abs(np.diff(pos) - 0.78) > 0.01): return None
    Am = W.max(1); c = int(Am.argmax())
    if c < 1 or c > len(pos) - 2: return None
    w0 = W[c]; A0 = w0.max(); pk = int(w0.argmax())
    if not (200 < A0 < 3500 and 4 <= pk <= 22 and abs(w0[:3].mean()) < 0.1 * A0): return None
    n2 = w0 @ w0; am, ap = W[c - 1] @ w0 / n2, W[c + 1] @ w0 / n2
    q = np.clip(Am[c - 1:c + 2], 0, None); cen = (q * pos[c - 1:c + 2]).sum() / q.sum()
    return pos[c], cen, A0, ap - am, ap + am, int(np.asarray(w['ch'])[o][c])

for tag, (build, rel) in RUNS.items():
    f = f'{A}/{rel}/wft/bench_cache.pkl'
    if not os.path.exists(f): print(tag, 'missing'); continue
    E = pickle.load(open(f, 'rb'))['events']; rows = []
    for k, e in E.items():
        t = e.get('truth')
        if not t or not e['wins'].get('x') or not e['wins'].get('y') or e.get('spark'): continue
        vx, vy = view(e['wins']['x'][0]), view(e['wins']['y'][0])
        if vx is None or vy is None: continue
        rows.append((k, t['ref_x'], t['ref_y'], t['tan_x'], t['tan_y'], *vx, *vy))
    r = np.array(rows, float)
    names = 'ev ref_x ref_y tan_x tan_y xc xcen Ax a1x s1x chx yc ycen Ay a1y s1y chy'.split()
    np.savez_compressed(f'{HERE}/tables/{tag}.npz', build=build, **{n: r[:, i] for i, n in enumerate(names)})
    print(f'{tag:8s} {build:13s} events {len(E):6d} kept {len(r):6d}  <a1x> {r[:, 8].mean():+.3f}  <a1y> {r[:, 14].mean():+.3f}', flush=True)
