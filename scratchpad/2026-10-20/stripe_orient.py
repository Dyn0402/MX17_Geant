"""Which way do the resistive stripes run on each board? (2026-10-20)
Charge on a stripe spreads along it, so the view measuring the coordinate ALONG the stripes gets slow, wide sharing
(stripes along y -> Y view slow). Per board/build and view, from the raw survey: mean +-1 and +-2 neighbour waveforms
over the central (Y-share-matched is not needed for a direction check), early (prompt) and late (> 600 ns after the
central peak) content, and the late/early ratio. run_config: det3 May 'Strip resist vertical', det1 April 'horizontal'.
Usage: stripe_orient.py"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L
rows = M.survey_rows(); OUT = {}
acc = {}
for r in rows:
    d = L.load(r['f']); W = d['W']; ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
    W = W[ok]
    if len(W) < 30: continue
    acc.setdefault((r['build'], r['view']), []).append(W)
print(f'{"build":13s} view   n   | ±1 peak/central  ±1 late/peak  ±2 late/c  | ±1 peak delay (ns)')
for (bu, v), Ws in sorted(acc.items()):
    W = np.concatenate(Ws); c = W[:, 4]; A = c.max(1, keepdims=True); pk = c.argmax(1)
    # align on the central peak sample
    T = 32; Z = 24; al = np.zeros((len(W), 9, 64)); 
    for i in range(len(W)):
        s = Z - pk[i]; al[i, :, s:s + T] = W[i] / A[i]
    m = al.mean(0); c0 = m[4]; n1 = (m[3] + m[5]) / 2; n2 = (m[2] + m[6]) / 2
    late = slice(Z + 10, Z + 16)          # 600-960 ns after the central peak
    o = dict(n=len(W), n1_peak=float(n1.max()), n1_late=float(n1[late].mean() / n1.max()), n2_late=float(n2[late].mean()),
             c_late=float(c0[late].mean()), n1_delay=float((np.argmax(n1) - Z) * 60), n2_delay=float((np.argmax(n2) - Z) * 60))
    OUT[f'{bu}|{v}'] = o
    print(f'{bu:13s} {v}  {o["n"]:6d} |   {o["n1_peak"]:.3f}          {o["n1_late"]:.3f}        {o["n2_late"]:.3f}   (c late {o["c_late"]:.3f}) | ±1 {o["n1_delay"]:+.0f}  ±2 {o["n2_delay"]:+.0f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stripe_orient.json'), 'w'))
