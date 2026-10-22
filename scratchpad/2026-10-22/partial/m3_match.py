"""Does the M3 track point at the det3 hit? (2026-10-22, §12.) Per session/view: chip centres in the measured M3 coordinate
from the chip medians of 'full-column' events; residual = M3 coord - centre of the hit's chip (a 64-channel chip spans
~50 mm, so a matched track has |res| <~ 25-30 mm, an unrelated one is spread over the board). Fractions matched for
truncated (9-ch end < 14) vs full events. Usage: m3_match.py TAG ..."""
import sys, numpy as np
import session_profile as SP, event_ends as EE
S = {s['tag']: s for s in SP.sessions()}
for tag in sys.argv[1:]:
    for view in 'xy':
        s = S[tag]; d = SP.L.load(s['npz'][view]); files = d['z']['files']; T = {i: SP.m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
        M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
        W = d['W'].astype(float); S9 = W.sum(1); tr = EE.last_cross(S9, .3) < 14; ch = p // 64
        spread = [np.std([np.median(M[g & ~tr & (ch == c), k]) for c in range(8) if (g & ~tr & (ch == c)).sum() > 15]) for k in (4, 5)]
        u = M[:, 4 if spread[0] > spread[1] else 5]
        cen = {c: np.median(u[g & ~tr & (ch == c)]) for c in range(8) if (g & ~tr & (ch == c)).sum() > 15}
        ok = g & np.isin(ch, list(cen)); res = np.full(len(u), np.nan); res[ok] = u[ok] - np.array([cen[c] for c in ch[ok]])
        mt = lambda m: (np.abs(res[m & ok]) < 30).mean()
        print(f"{tag} {view} {s['sess'][-38:]:38s} trunc {tr.mean():.2f}  matched |res|<30mm: full {mt(~tr):.2f} (n {(~tr & ok).sum()})  trunc {mt(tr):.2f} (n {(tr & ok).sum()})"
              f"   |res| p50 full {np.nanmedian(np.abs(res[~tr])):.0f} trunc {np.nanmedian(np.abs(res[tr])):.0f} mm")
