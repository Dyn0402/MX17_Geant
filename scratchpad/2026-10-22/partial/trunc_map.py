"""Where and when are the 6-23 vertical-track columns truncated? (2026-10-22, §12.) Truncated = 9-channel-sum end
(last 30 % crossing) < 14 samples. Fraction vs M3 track position at det3 (4 x 4 board bins), vs |th_meas| / |th_along|,
vs event order (time in the run), vs M3 track multiplicity (all tracks, any chi2), and vs peak sample.
Usage: trunc_map.py TAG [VIEW]"""
import sys, numpy as np
import session_profile as SP, event_ends as EE
tag = sys.argv[1]; S = {s['tag']: s for s in SP.sessions()}; s = S[tag]
for view in (sys.argv[2:] or 'xy'):
    d = SP.L.load(s['npz'][view]); files = d['z']['files']; T = {i: SP.m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
    M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
    spread = [np.std([np.median(M[g & (p // 64 == c), k]) for c in range(8) if (g & (p // 64 == c)).sum() > 15]) for k in (4, 5)]
    km = 4 if spread[0] > spread[1] else 5; thm, tha = M[:, km - 2], M[:, 3 if km == 4 else 2]
    W = d['W'].astype(float); S9 = W.sum(1); A9 = S9.max(1); e9 = EE.last_cross(S9, .3); tr = e9 < 14
    print(f'\n{tag} {view} {s["sess"]}: all events n {len(tr)}, truncated {tr.mean():.2f};  M3 1-good {g.mean():.2f}')
    for nm, x, bins in (('|th_meas|', np.abs(thm), (0, 2, 5, 10, 20, 60)), ('|th_along|', np.abs(tha), (0, 2, 5, 10, 20, 60)),
                        ('M3 x', M[:, 4], np.percentile(M[g, 4], [0, 25, 50, 75, 100])), ('M3 y', M[:, 5], np.percentile(M[g, 5], [0, 25, 50, 75, 100])),
                        ('event order (file, ev)', d['fi'] * 1e6 + d['ev'], np.percentile(d['fi'] * 1e6 + d['ev'], [0, 25, 50, 75, 100])),
                        ('A9', A9, np.percentile(A9, [0, 25, 50, 75, 100]))):
        row = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (x >= lo) & (x < hi) & (g if 'M3' in nm or 'th' in nm else True)
            if m.sum() >= 10: row.append(f'[{lo:.4g},{hi:.4g}) {tr[m].mean():.2f} (n {m.sum()})')
        print(f'   {nm:24s}', '  '.join(row))
    for nt in (0, 1, 2):
        m = (M[:, 0] == nt) if nt < 2 else (M[:, 0] >= 2)
        if m.sum() >= 10: print(f'   M3 tracks {nt}{"+" if nt == 2 else ""}: truncated {tr[m].mean():.2f} (n {m.sum()})')
    pk = S9.argmax(1); print('   peak sample p10/50/90: truncated', np.percentile(pk[tr], [10, 50, 90]), ' full', np.percentile(pk[~tr], [10, 50, 90]))
