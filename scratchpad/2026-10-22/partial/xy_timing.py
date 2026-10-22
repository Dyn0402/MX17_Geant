"""det3 X vs Y (FEU 3 vs 4, same avalanches) per-event timing (2026-10-22, §12): joined by (decoded file number, eventId).
If the random onsets of 6-23 are physical (charge arriving late in det3) X and Y agree event by event; an FEU timing
problem would decorrelate them. Usage: xy_timing.py TAG ..."""
import re, sys, numpy as np
import session_profile as SP, event_ends as EE
for tag in sys.argv[1:]:
    s = {x['tag']: x for x in SP.sessions()}[tag]; R = {}
    for v in 'xy':
        d = SP.L.load(s['npz'][v]); fn = [re.search(r'_(\d{3})_\d\d\.root$', str(f)).group(1) for f in d['z']['files']]
        S9 = d['W'].astype(float).sum(1); A = S9.max(1); on = np.argmax(S9 >= .3 * A[:, None], 1); e = EE.last_cross(S9, .3)
        pk = S9.argmax(1); R[v] = {(fn[fi], int(ev)): (on[i], e[i], pk[i], A[i]) for i, (fi, ev) in enumerate(zip(d['fi'], d['ev']))}
    k = sorted(set(R['x']) & set(R['y'])); X = np.array([R['x'][q] for q in k], float); Y = np.array([R['y'][q] for q in k], float)
    print(f"{tag} {s['sess'][-40:]:40s} joined {len(k)}  corr onset {np.corrcoef(X[:, 0], Y[:, 0])[0, 1]:.2f}  end {np.corrcoef(X[:, 1], Y[:, 1])[0, 1]:.2f}  "
          f"peak {np.corrcoef(X[:, 2], Y[:, 2])[0, 1]:.2f}  amp {np.corrcoef(np.log(X[:, 3]), np.log(Y[:, 3]))[0, 1]:.2f};  |dpeak|<=1: {np.mean(np.abs(X[:, 2] - Y[:, 2]) <= 1):.2f}  "
          f"peak p10/50/90 X {np.percentile(X[:, 2], [10, 50, 90])} Y {np.percentile(Y[:, 2], [10, 50, 90])}")
