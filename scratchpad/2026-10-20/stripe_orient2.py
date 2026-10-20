"""Stripe direction per board, every build we have, X vs Y side by side (2026-10-20, follow-up to stripe_orient.py).
Charge landing on a resistive stripe spreads ALONG it, so the view measuring the coordinate along the stripes sees
slow, late neighbour charge; the other view's neighbours are prompt. Same metric as stripe_orient.py (mean +-1 / +-2
waveforms over the peak channel, aligned on the central peak; 'late' = 600-960 ns after it), now per (build, view)
including the boards the survey lacked: det0 (4-1 He/ethane HV scan), det4 May build (5-7, Ar/CF4), and the
April-24 build ('Strip resist horizontal', EOS det_1_old 4-25, M3-filtered arrays, raw_asym.py --array).
The Apr-24 board needs a 60 ADC floor instead of 200 (THR).
Verdict per build: Y-slow (stripes along y, as det3 May 'vertical') or X-slow (stripes along x).
Usage: NB_RAW=<raw dir> stripe_orient2.py"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L
OR = f'{L.RAW}/orient'
EXTRA = {'det0 (4-1 He)': ['d0_505', 'd0_520'], 'det4 May (5-6 bulk)': ['d4m_490', 'd4m_510'],
         'Apr-24 bulk "horizontal"': ['a24']}
acc = {}
for r in M.survey_rows():
    acc.setdefault((r['build'], r['view']), []).append(r['f'])
for b, tags in EXTRA.items():
    for t in tags:
        for v in 'xy':
            acc.setdefault((b, v.upper()), []).append(f'{OR}/{t}_{v}.npz')
THR = {'Apr-24 bulk "horizontal"': 60}   # that board's cosmics sit mostly below 200 ADC (noise sigma ~10 ADC)
def load_thr(f, thr):
    # asym_lib.load with the 200 ADC floor replaced by thr; hot-channel veto as in asym_lib (the Apr-24 board has a
    # fixed pulsing spot at X 487-494 / Y 48-52 that it removes)
    z = np.load(f); W = z['wf']; c = W[:, L.HW]; A = c.max(1); pk = c.argmax(1); pre = c[:, :3].mean(1)
    g = (A > thr) & (A < 3500) & (pk >= 4) & (pk <= 22) & (np.abs(pre) < 0.1 * A)
    g &= ~np.isin(z['p'], L.hot(z['p'][g])); return dict(W=W[g], p=z['p'][g].astype(int))
def metric(files, thr=None):
    Ws = []
    for f in files:
        d = load_thr(f, thr) if thr else L.load(f); W = d['W']; ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59); Ws.append(W[ok])
    W = np.concatenate(Ws)
    if len(W) < 30: return None
    c = W[:, 4]; A = c.max(1); pk = c.argmax(1); Z = 24; al = np.zeros((len(W), 9, 64))
    for i in range(len(W)):
        s = Z - pk[i]; al[i, :, s:s + 32] = W[i] / A[i]
    m = al.mean(0); n1 = (m[3] + m[5]) / 2; n2 = (m[2] + m[6]) / 2; late = slice(Z + 10, Z + 16)
    # bootstrap error on the +-2 late content
    rng = np.random.default_rng(1); bs = []
    for _ in range(100):
        j = rng.integers(0, len(W), len(W)); mm = al[j].mean(0); bs.append(((mm[2] + mm[6]) / 2)[late].mean())
    return dict(n=len(W), n1=float(n1.max()), n1_delay=float((np.argmax(n1) - Z) * 60), n2_late=float(n2[late].mean()),
                n2_late_err=float(np.std(bs)), n2_delay=float((np.argmax(n2) - Z) * 60))
OUT = {}
builds = sorted({b for b, _ in acc}, key=lambda b: ('det' not in b, b))
print(f'{"build":26s} | {"n X/Y":>11s} | ±1 share X/Y | ±1 delay X/Y ns | ±2 delay X/Y ns | ±2 late X / Y (x1000)   | slow view')
for b in builds:
    mx, my = (metric(acc[(b, v)], THR.get(b)) if (b, v) in acc else None for v in 'XY')
    if not mx or not my: continue
    dl = my['n2_late'] - mx['n2_late']; e = np.hypot(mx['n2_late_err'], my['n2_late_err'])
    verdict = ('Y' if dl > 0 else 'X') + f'  ({dl / e:+.1f} σ)'
    OUT[b] = dict(X=mx, Y=my, slow=verdict)
    print(f'{b:26s} | {mx["n"]:5d}/{my["n"]:5d} |  {mx["n1"]:.2f} / {my["n1"]:.2f} |  {mx["n1_delay"]:+5.0f} / {my["n1_delay"]:+5.0f} |'
          f'  {mx["n2_delay"]:+5.0f} / {my["n2_delay"]:+5.0f} | {1e3*mx["n2_late"]:5.1f}±{1e3*mx["n2_late_err"]:.1f} / {1e3*my["n2_late"]:5.1f}±{1e3*my["n2_late_err"]:.1f} | {verdict}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stripe_orient2.json'), 'w'), indent=1)
