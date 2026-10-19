"""Copy survey over many dumps (2026-10-19). For each dump and view: the one-sided copy from
 (1) the per-event in-phase content: a_ev = mean(alpha_-1 - alpha_+1) over interior peaks (±1 on the same chip),
 (2) the ±2 difference: mean(alpha_-2 - alpha_+2) / mean(alpha_+1)  (copy model: n_-2 - n_+2 ~ a s1),
plus the clean-side share alpha_+1 and the mean central amplitude. alpha_k = <w_k, w_c> / <w_c, w_c> per event.
Peak strip 200..3500 ADC. Errors are standard errors over events.
Usage: nb_copy_survey.py LABEL=NPZ ...   (prints a table; --csv OUT writes it)"""
import sys, numpy as np
HW = 5
rows = []; args = [a for a in sys.argv[1:] if not a.startswith('--')]
csv = sys.argv[sys.argv.index('--csv') + 1] if '--csv' in sys.argv else None
if csv: args = [a for a in args if a != csv]
for item in args:
    lab, f = item.split('=', 1); z = np.load(f)
    for v in 'xy':
        W = z[f'{v}_wf']; ch = z[f'{v}_ch'].astype(int)
        if len(W) < 50: continue
        amp = np.nanmax(W, 2); c = amp[:, HW]
        m = (c > 200) & (c < 3500) & np.isfinite(W[:, HW - 2:HW + 3, :20]).all((1, 2))
        W, ch = np.nan_to_num(W[m]), ch[m]; wc = W[:, HW]; n2 = np.einsum('ij,ij->i', wc, wc)
        al = {k: np.einsum('ij,ij->i', W[:, HW + k], wc) / n2 for k in (-2, -1, 1, 2)}
        inner = (ch[:, HW - 1] // 64 == ch[:, HW] // 64) & (ch[:, HW + 1] // 64 == ch[:, HW] // 64)
        d1 = (al[-1] - al[1])[inner]; d2 = al[-2] - al[2]
        a1, e1 = d1.mean(), d1.std() / np.sqrt(len(d1))
        a2 = d2.mean() / al[1].mean(); e2 = d2.std() / np.sqrt(len(d2)) / al[1].mean()
        rows.append((lab, v, int(m.sum()), a1, e1, a2, e2, al[1].mean(), al[-1].mean(), float(np.median(c[m]))))
print(f'{"dataset":38s} view     n   copy a (±1)        copy a (±2/s1)     alpha+1  alpha-1  median c')
for r in rows:
    print(f'{r[0]:38s} {r[1].upper():>3s} {r[2]:6d}  {r[3]:+.3f} ± {r[4]:.3f}     {r[5]:+.3f} ± {r[6]:.3f}     {r[7]:.3f}    {r[8]:.3f}   {r[9]:6.0f}')
if csv:
    with open(csv, 'w') as fh:
        fh.write('dataset,view,n,a1,e1,a2,e2,alpha_p1,alpha_m1,median_c\n')
        for r in rows: fh.write(','.join(str(x) for x in r) + '\n')
