"""Audit (2026-10-21): does the session-to-session a1 of det3 June-15 follow the cloud width, once the pad labels
are dropped? Per session: a1 (all events), median Y share, median s1X, median central pulse width (samples above half
max, X view), and a1 within the session vs Y share (all events). -> sess_width.json.
Usage: NB_RAW=~/x17/nb_asym_raw python sess_width.py"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import r5lib as R
E = R.ev(R.BUILDS['det3 June-15'])
a1 = E['alx'][:, 5] - E['alx'][:, 3]; ysh = E['ysh']; sx = E['alx'][:, 5] + E['alx'][:, 3]
c = E['Wx'][:, 4]; wid = (c >= 0.5 * c.max(1, keepdims=True)).sum(1)
cy = E['Wy'][:, 4]; widy = (cy >= 0.5 * cy.max(1, keepdims=True)).sum(1)
OUT = {}
for s in np.unique(E['sess']):
    m = E['sess'] == s
    if m.sum() < 300: continue
    q = np.quantile(wid[m], [0.25, 0.5, 0.75])
    by_w = [float(a1[m & (wid <= q[0])].mean()), float(a1[m & (wid > q[2])].mean())]
    OUT[s] = dict(n=int(m.sum()), a1=float(a1[m].mean()), a1e=float(a1[m].std() / np.sqrt(m.sum())), ysh_med=float(np.median(ysh[m])),
                  s1x_med=float(np.median(sx[m])), wid_med=float(np.median(wid[m])), widy_med=float(np.median(widy[m])), a1_short_long=by_w)
    print(f"{s:45s} n {m.sum():5d} a1 {OUT[s]['a1']:+.3f}±{OUT[s]['a1e']:.3f}  Ysh med {OUT[s]['ysh_med']:.2f}  s1X med {OUT[s]['s1x_med']:.2f}  pulse width X/Y {OUT[s]['wid_med']:.0f}/{OUT[s]['widy_med']:.0f} smp  a1 short/long-pulse quartile {by_w[0]:+.3f}/{by_w[1]:+.3f}")
json.dump(OUT, open(os.path.join(HERE, 'sess_width.json'), 'w'), indent=1)
