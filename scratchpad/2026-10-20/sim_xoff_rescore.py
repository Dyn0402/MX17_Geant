"""Rescore sim_xoff_eval*.json (2026-10-20): the eval's chi2 uses only bins where the sim has events in that pad class,
so arms that make no narrow X-pad events (plain copy, no offset) never pay for the data's narrow X-pad sign flip.
Here the Y-pad fraction fY per Y-share bin is added (binomial error on the data count, 0.02 floor), which does
penalise that. Usage: sim_xoff_rescore.py [sim_xoff_eval_sym.json]"""
import sys, json, numpy as np
J = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'sim_xoff_eval_sym.json'))
for ref, dd in J['data'].items():
    td = np.array(dd['t']); n = dd['n']
    print(f'== vs {ref}')
    print(f'{"arm":8s} copy | chi2 asym (nbins)  chi2 fY  total | narrow X-pad sim (data {td[1,0]:+.2f} {td[4,0]:+.2f} {td[7,0]:+.2f})')
    for s in J['sim']:
        t = np.array(s['t']); ca, na = 0., 0
        for j in range(7):
            for c in (0, 1):
                i = 3 * j + c
                if np.isfinite(t[i, 0]) and np.isfinite(td[i, 0]):
                    ca += (t[i, 0] - td[i, 0]) ** 2 / (t[i, 1] ** 2 + td[i, 1] ** 2 + .02 ** 2); na += 1
        f, fd = t[2::3, 0], td[2::3, 0]; ok = np.isfinite(f) & np.isfinite(fd)
        ef = np.sqrt(fd * (1 - fd) / max(n / 7, 1) + f * (1 - f) / max(s['n'] / 7, 1) + .02 ** 2)
        cf = float((((f - fd) / ef) ** 2)[ok].sum())
        nx = ' '.join('  -  ' if not np.isfinite(v) else f'{v:+.2f}' for v in t[[1, 4, 7], 0])
        print(f'{s["arm"]:8s} {s["copy"]:.2f} | {ca:6.0f} ({na:2d})        {cf:6.0f}  {ca + cf:6.0f} | {nx}')
