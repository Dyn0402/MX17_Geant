"""What makes the data's partial-drift-column pulses? (2026-10-22, FOLLOWUP §12.) Per session and view: central-pad pulse
metrics (onset = 5 % rise, f50 = 50 % fall, FWHM, |dt_nb| = neighbour centroid-time difference, as cloud/short_pulses.py)
joined event by event to the M3 reference track (m3lib: Chi2X, Chi2Y < 5) at the detector plane.
Classes among |dt_nb| < 80 ns: short FWHM < 350 ns, mid 350-450, long >= 450. Per class: number of M3 tracks (0 / 1 good /
>1), the true angles of the single good track along the view's measured coordinate (th_meas) and along the pads (th_along),
the track's position (edge distance), and the pulse timing. Axis (checked from chip medians, diag.py): view X <-> M3 y in May and June
(the May run_config's rot 0 notwithstanding). Usage: partial_m3.py [SESSION ...]  (-> partial_m3.json, small)"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-19', 'asym_deck')); sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import asym_lib as L, prompt_share as PS, m3lib
RAW = os.path.expanduser('~/x17/nb_asym_raw'); HW = 4; DT = 60.
#            name: (npz stem, det z mm, rotated 90?, drift V)
SESS = {'d3may': ('d3may', 232., True, 900), 'd3jun': ('d3jun', 232., True, 600),
        'drift100': ('survey/s027', 702., True, 100), 'drift300': ('survey/s029', 702., True, 300),
        'drift500': ('survey/s030', 702., True, 500), 'drift700': ('survey/s031', 702., True, 700),
        'drift900': ('survey/s032', 702., True, 900), 'drift1100': ('survey/s028', 702., True, 1100)}
CLS = (('short', 0, 350), ('mid', 350, 450), ('long', 450, 5000))

def metrics(W):
    c = W[:, HW]; pk = c.argmax(1); ts = np.arange(W.shape[2]); pos = lambda w: np.clip(w, 0, None)
    tc = lambda w: (pos(w) * ts).sum(1) / np.maximum(pos(w).sum(1), 1e-9)
    on, f50, r50 = PS.cross(c, pk, .05), PS.cross(c, pk, .5, False), PS.cross(c, pk, .5)
    s1 = (np.einsum('ij,ij->i', W[:, HW + 1] + W[:, HW - 1], c) / np.einsum('ij,ij->i', c, c))
    return dict(on=on, f50=f50, fw=(f50 - r50) * DT, dtn=(tc(W[:, HW + 1]) - tc(W[:, HW - 1])) * DT, amp=c[np.arange(len(c)), pk], s1=s1, pk=pk)

def session(name, view):
    stem, z, rot, hv = SESS[name]; f = f'{RAW}/{stem}_{view}.npz'
    if not os.path.exists(f): return None
    d = L.load(f); files = d['z']['files']; T = {i: m3lib.m3_table(fn, z) for i, fn in enumerate(files)}
    if any(t is None for t in T.values()): return None
    M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])])
    m = metrics(d['W']); meas = (view == 'x') != rot          # measured coordinate is M3 x?
    thm, tha = (M[:, 2], M[:, 3]) if meas else (M[:, 3], M[:, 2]); um, ua = (M[:, 4], M[:, 5]) if meas else (M[:, 5], M[:, 4])
    v = np.abs(m['dtn']) < 80; res = dict(n=int(len(v)), nv=int(v.sum()), drift=hv, classes={})
    for cl, lo, hi in CLS:
        k = v & (m['fw'] >= lo) & (m['fw'] < hi); g = k & np.isfinite(thm)
        if k.sum() < 10: continue
        q = lambda x, s=g: [float(u) for u in np.nanpercentile(x[s], [25, 50, 75])] if s.sum() else None
        res['classes'][cl] = dict(n=int(k.sum()), frac=float(k.sum() / v.sum()), m3_none=float((M[k, 0] == 0).mean()), m3_unmatched=float((M[k, 0] < 0).mean()),
            m3_1good=float(g.sum() / k.sum()), m3_multi=float((M[k, 0] > 1).mean()), abs_th_meas=q(np.abs(thm)), abs_th_along=q(np.abs(tha)),
            th_meas=q(thm), th_along=q(tha), u_meas=q(um), u_along=q(ua), onset=q(m['on'], k), f50=q(m['f50'], k), amp=float(np.median(m['amp'][k])), s1=float(np.median(m['s1'][k])))
    res['_ev'] = dict(M=M, thm=thm, tha=tha, um=um, ua=ua, **m)
    return res

def show(name, view, r):
    print(f"\n{name} {view.upper()}  drift {r['drift']} V   n {r['n']}  |dt_nb|<80: {r['nv']}")
    print('   class   n   frac  M3:0trk  1good  >1trk   |th_meas| p50 (p25-p75)   |th_along| p50    onset p50  f50 p50   amp   s1')
    for cl, c in r['classes'].items():
        a, b = c['abs_th_meas'] or [np.nan] * 3, c['abs_th_along'] or [np.nan] * 3
        print(f"   {cl:6s} {c['n']:4d} {c['frac']:5.2f}   {c['m3_none']:5.2f}  {c['m3_1good']:5.2f}  {c['m3_multi']:5.2f}    {a[1]:5.1f} ({a[0]:4.1f}-{a[2]:4.1f})        "
              f"{b[1]:5.1f} ({b[0]:4.1f}-{b[2]:4.1f})   {c['onset'][1]:5.1f}   {c['f50'][1]:5.1f}   {c['amp']:5.0f}  {c['s1']:.2f}")

if __name__ == '__main__':
    out = {}
    for name in (sys.argv[1:] or SESS):
        for view in ('x', 'y'):
            r = session(name, view)
            if r is None: print('skip', name, view); continue
            show(name, view, r); r.pop('_ev'); out[f'{name}_{view}'] = r
    json.dump(out, open(os.path.join(HERE, 'partial_m3.json'), 'w'), indent=1)
