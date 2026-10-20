"""Can stripe registration make the asymmetry? (2026-10-20 evening) The sim's ESL stripes (800 um, centred on x = 0)
over 780 um pads give each pad column c a stripe offset delta(c) = signed distance from the column centre
(PAD_ORIGIN + c P) to the nearest stripe centre; it changes by only 20 um per column, so locally each delta is a
pitch-locked board. X-only events of every xo0 arm (th00 + thy*, seeds 5..16): a1 and coin per 50 um delta bin, both
column-mapping orientations (the DAQ->x direction of 'col' is a convention), and the narrow/wide halves by X share.
Usage: sim_delta.py"""
import os, glob, re, json, numpy as np
import r5lib as R
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff'); P = 780.; ORG = -(512 - 1) / 2 * P
p, a1, A, s1 = [], [], [], []
for f in sorted(glob.glob(f'{SIM}/xo0_*_x.npz')):
    if not re.search('_th00_|_thy', f) or np.load(f)['wf'].ndim != 3: continue
    d = R.L.load(f); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
    p.append(d['p'][ok] % 512); a1.append(d['a1'][ok]); A.append(d['A'][ok]); s1.append((d['al'][1] + d['al'][-1])[ok])
p, a1, A, s1 = map(np.concatenate, (p, a1, A, s1)); col = 64 * (p // 64) + 63 - p % 64
print(f'X-only events {len(p)}, distinct columns {len(np.unique(col))}')
OUT = {}
for orient, cc in (('col', col), ('511-col', 511 - col)):
    x = ORG + cc * P; delta = ((x + 400.) % 800.) - 400.
    ed = np.arange(-400, 401, 50); b = np.digitize(delta, ed) - 1; med = np.median(s1); rows = []
    for k in range(len(ed) - 1):
        m = b == k
        if m.sum() < 40: rows.append(None); continue
        rows.append(dict(d=float((ed[k] + ed[k + 1]) / 2), n=int(m.sum()), a1=float(a1[m].mean()), e=float(a1[m].std() / np.sqrt(m.sum())),
                         coin=float((A[m, 5] > A[m, 3]).mean()), narrow=float(a1[m & (s1 < med)].mean()), wide=float(a1[m & (s1 >= med)].mean()),
                         s1=float(s1[m].mean())))
    OUT[orient] = rows
    print(f'-- delta from column {orient}: delta(um)  n  a1  coin  a1(narrow half)  a1(wide half)  s1')
    for r in rows:
        if r: print(f'   {r["d"]:+5.0f} {r["n"]:5d} {r["a1"]:+.3f}±{r["e"]:.3f} {r["coin"]:.3f} {r["narrow"]:+.3f} {r["wide"]:+.3f}  {r["s1"]:.3f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_delta.json'), 'w'), indent=1)
