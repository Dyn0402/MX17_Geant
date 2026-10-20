"""Inverse question with the fast model (2026-10-20): line L gets extra sensitivity w over a region R on its bus side
(+x, toward column L+1 = toward the data's n-1 side of an event whose central is L+1), taken from a donor.
Indicator regions at the resistive layer, folded with the same clouds as fwd_xpos.py:
  gap    x - L P in [+340, +440] (the inter-column gap where line L's bus runs)      donor: Y (no X change)
  ypad   the Y pads of column L+1                                                     donor: Y
  ynear  the half of those Y pads nearest line L                                      donor: Y
  xskip  the X pads of column L+2                                                     donor: X line L+2
  xnext  the X pads of line L+1 (rows +-P in the +x column)                           donor: X line L+1
  col    all of column L+1 (its X pads and Y pads)                                    donor: half X line L+1, half Y
Donor Y means the Y view loses the same charge (Y share and peak re-found with it). Scored like fwd_xpos.py.
Usage: fwd_region.py"""
import numpy as np, json, os
import fwd_xpos as F
P, K, n = F.P, F.K, F.n; px, py, W = F.px, F.py, F.W
def ind(kind, L):
    # line L's frame, read off the field map (fW): X_L's pad at (L P + 390, L P); its column alternates X_L / Y pads;
    # the +x column holds Y_L's pad in row L P and X_{L+1}'s pads at rows L P +- P. y is periodic with 2P.
    x = px - P * L - 390.; y = ((py - P * L + P) % (2 * P)) - P
    inpad = lambda xc, yc: (np.abs(x - xc) <= 340) & (np.minimum(np.abs(y - yc), 2 * P - np.abs(y - yc)) <= 340)
    if kind == 'gap':     r = (x >= 340) & (x <= 440)
    elif kind == 'ypad':  r = inpad(P, 0.)                 # Y_L's pad in the +x column, same row
    elif kind == 'ynear': r = inpad(P, 0.) & (x <= P)       # its half nearest line L
    elif kind == 'xskip': r = inpad(2 * P, 0.)             # X_{L+2}'s pad, same row
    elif kind == 'xnext': r = inpad(P, P)                  # X_{L+1}'s pads (rows +- P) in the +x column
    elif kind == 'col':   r = (np.abs(x - P) <= 340)
    return (r * W).sum((1, 2))
def model(kind, w):
    S = F.SX(0.).copy(); SYm = F.SY.copy()
    gain = w * 0.83 * np.stack([ind(kind, L) for L in K], 1)      # 0.83 = a full pad's direct weighting potential
    S = S + gain
    if kind in ('ypad', 'ynear', 'gap'): SYm = SYm - gain           # donor: Y line L (its pad in the +x column)
    elif kind == 'xskip': S[:, :-2] -= gain[:, 2:] * 0 ; S[:, 2:] -= gain[:, :-2]
    elif kind == 'xnext': S[:, 1:] -= gain[:, :-1]                  # donor: X line L+1
    elif kind == 'col':   S[:, 1:] -= 0.5 * gain[:, :-1]; SYm = SYm - 0.5 * gain
    return S, SYm
def main():
    base = F.analyse(F.SX(0.)); fb = F.fingerprint(base)
    f = lambda a: ' '.join(f'{x:+.2f}' for x in a)
    print(f'{"DATA":16s} AY [{f(F.DAT["AY"])}] AX [{f(F.DAT["AX"])}] dfY [{f(F.DAT["dfY"])}] dvu [{f(F.DAT["dvu"])}]')
    OUT = {}
    for kind in ('gap', 'ypad', 'ynear', 'xnext', 'xskip', 'col'):
        for w in (.1, .2, .35, .5):
            S, SYm = model(kind, w); SY0 = F.SY; F.SY = SYm; e = F.analyse(S); F.SY = SY0; fm = F.fingerprint(e)
            r = dict(AY=fm['AY'], AX=fm['AX'], dfY=fm['fY'] - fb['fY'], dvu=fm['vu'] - F.fingerprint(base, qu=np.array(fm['uq']))['vu'])
            chi = {k: float(np.sum((r[k] - F.DAT[k]) ** 2) / .03 ** 2) for k in F.DAT}
            OUT[f'{kind} {w}'] = dict(a1=fm['a1'], chi=chi, **{k: v.tolist() for k, v in r.items()})
            print(f'{kind+" "+str(w):16s} a1 {fm["a1"]:+.3f} | AY [{f(r["AY"])}] AX [{f(r["AX"])}] dfY [{f(r["dfY"])}] dvu [{f(r["dvu"])}] | ' +
                  ' '.join(f'{k} {v:4.0f}' for k, v in chi.items()) + f' tot {sum(chi.values()):5.0f}')
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fwd_region.json'), 'w'), indent=1)
if __name__ == '__main__':
    main()
