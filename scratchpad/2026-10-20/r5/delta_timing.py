"""Is the sim's stripe-registration asymmetry prompt, like the data's? (2026-10-21)
The production kernel's registration asymmetry builds up over 200-300 ns of sheet transport (kernel_delta.py), while
the data's excess peaks with the central (timing.py, route 21). Same X-only estimator on both: per event, waveforms
normalised by the central amplitude and aligned on the central's peak sample; mean excess e = W(n+1) - W(n-1) and
mean symmetric neighbour s = (W(n+1) + W(n-1)) / 2. Sim: |delta| >= 150 um, sign-folded so the excess is positive
(delta from sim_delta.py, 'col' orientation: a1 > 0 at delta > 0). Data: det3 June-15 X (all sub-runs), det2 X.
Reports peak time and 50 % rise of e, s and the central (ns, relative to the central's peak sample).
Usage: delta_timing.py [extra sim arm prefix and dir, e.g. fu1 ~/x17/nb_asym_raw/simfunnel] (data are skipped then)"""
import os, sys, glob, re, json, numpy as np
import r5lib as R
import timing as T
L = R.L; HW = L.HW; Z = 8; NT = 40
def prof(W, sgn):
    A = W[:, HW].max(1); pk = W[:, HW].argmax(1)
    X = T._align(W / A[:, None, None], pk, Z, NT)
    c = X[:, HW].mean(0); e = (sgn[:, None] * (X[:, HW + 1] - X[:, HW - 1])).mean(0); s = (.5 * (X[:, HW + 1] + X[:, HW - 1])).mean(0)
    f = lambda w: dict(peak=T._pk(w) - Z * 60, t50=T._t50(w) - Z * 60, amp=float(w.max()))
    return dict(n=int(len(W)), c=f(c), e=f(e), s=f(s), w=dict(c=c.tolist(), e=e.tolist(), s=s.tolist()))
OUT = {}
PRE = sys.argv[1] if len(sys.argv) > 1 else 'xo0'
SIM = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else '~/x17/nb_asym_raw/simxoff'); P = 780.; ORG = -(512 - 1) / 2 * P
Ws, sg, dl = [], [], []
for f in sorted(glob.glob(f'{SIM}/{PRE}_*_x.npz')):
    if not re.search('_th00_|_thy', f) or np.load(f)['wf'].ndim != 3: continue
    d = L.load(f); p = d['p'] % 512; ok = (p % 64 >= 4) & (p % 64 <= 59)
    col = 64 * (p // 64) + 63 - p % 64; delta = ((ORG + col * P + 400.) % 800.) - 400.
    m = ok & (np.abs(delta) >= 150); Ws.append(d['W'][m].astype(np.float32)); sg.append(np.sign(delta[m])); dl.append(delta[m])
Ws, sg, dl = np.concatenate(Ws), np.concatenate(sg), np.concatenate(dl)
OUT[f'sim {PRE} |delta|>=150'] = prof(Ws, sg)
m = np.abs(dl) >= 250; OUT[f'sim {PRE} |delta|>=250'] = prof(Ws[m], sg[m])
for b in (('det3 June-15', 'det2', 'det3 May') if PRE == 'xo0' else ()):
    W = [L.load(r['f'])['W'].astype(np.float32) for r in R.rows() if r['build'] == b and r['view'] == 'X']
    W = np.concatenate(W); OUT[f'data {b} X'] = prof(W, np.ones(len(W)))
print(f'{"sample":22s} {"n":>7s} | central pk/t50 | excess amp pk/t50 | sym nb amp pk/t50  (ns)')
for k, o in OUT.items():
    print(f'{k:22s} {o["n"]:7d} | {o["c"]["peak"]:+5.0f}/{o["c"]["t50"]:+5.0f} | {o["e"]["amp"]:.3f} {o["e"]["peak"]:+5.0f}/{o["e"]["t50"]:+5.0f} '
          f'| {o["s"]["amp"]:.3f} {o["s"]["peak"]:+5.0f}/{o["s"]["t50"]:+5.0f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'delta_timing{"" if PRE == "xo0" else "_" + PRE}.json'), 'w'), indent=1)
