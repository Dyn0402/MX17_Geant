"""Is the 6-23 det3 column a few separate bursts? (2026-10-22, §12.) Per event (9-ch sum, both views): number of separated
local maxima above 30 % of the event peak (>= 3 samples apart, dip between them below 70 % of the smaller one), and the
'fill' = integral / (peak x 30 %-width) (a box -> ~1, isolated spikes -> low). M3-unselected, all events.
Usage: sparse.py TAG ..."""
import sys, numpy as np
import session_profile as SP

def bursts(w):
    pk = w.max(); thr = .3 * pk; m = [i for i in range(1, len(w) - 1) if w[i] >= thr and w[i] >= w[i - 1] and w[i] >= w[i + 1]]
    keep = []
    for i in m:
        if keep and i - keep[-1] < 3: 
            if w[i] > w[keep[-1]]: keep[-1] = i
            continue
        if keep and w[keep[-1]:i + 1].min() > .7 * min(w[keep[-1]], w[i]):
            if w[i] > w[keep[-1]]: keep[-1] = i
            continue
        keep.append(i)
    return len(keep)

for tag in sys.argv[1:]:
    s = {x['tag']: x for x in SP.sessions()}[tag]
    S9 = np.concatenate([SP.L.load(s['npz'][v])['W'] for v in 'xy']).astype(float).sum(1)
    nb = np.array([bursts(w) for w in S9]); A = S9.max(1); on = np.argmax(S9 >= .3 * A[:, None], 1)
    end = 31 - np.argmax((S9 >= .3 * A[:, None])[:, ::-1], 1); fill = np.clip(S9, 0, None).sum(1) / (A * np.maximum(end - on + 1, 1))
    print(f"{tag} {s['sess'][-40:]:40s} n {len(nb)}  bursts 1/2/3+: {np.mean(nb == 1):.2f} {np.mean(nb == 2):.2f} {np.mean(nb >= 3):.2f}   fill p25/50/75 {np.percentile(fill, [25, 50, 75]).round(2)}   A9 p50 {np.median(A):.0f}")
