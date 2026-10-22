"""Vertical-track (M3 |th_meas| < 2 deg) mean waveforms per session, by channel distance (2026-10-22, §12): is the late
(cathode-side) charge missing (attachment / no drift) or spread sideways (diffusion)? Prints the mean central, +-1, +-2..4
and 9-channel-sum profiles (ADC, not normalised) per sample, and the charge fraction outside the central pad per time slice.
Usage: vert_profiles.py TAG [TAG ...]   (survey tags s0NN or d3may; -> vert_profiles.json)"""
import os, sys, json, numpy as np
import session_profile as SP
S = {s['tag']: s for s in SP.sessions()}; out = {}
for tag in sys.argv[1:]:
    for view in 'xy':
        s = S[tag]; d = SP.L.load(s['npz'][view]); files = d['z']['files']; T = {i: SP.m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
        M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
        spread = [np.std([np.median(M[g & (p // 64 == c), k]) for c in range(8) if (g & (p // 64 == c)).sum() > 15]) for k in (4, 5)]
        thm = M[:, 2 if spread[0] > spread[1] else 3]; v = g & (np.abs(thm) < 2); W = d['W'][v].astype(float)
        c, n1, far = W[:, 4].mean(0), (W[:, 3] + W[:, 5]).mean(0), (W[:, :3].sum(1) + W[:, 6:].sum(1)).mean(0)
        out[f'{tag}_{view}'] = dict(sess=s['sess'], n=int(v.sum()), c=c.round(1).tolist(), n1=n1.round(1).tolist(), far=far.round(1).tolist())
        print(f"\n{tag} {view} {s['sess']}  n {v.sum()}\n  s : " + ' '.join(f'{i:5d}' for i in range(2, 30, 2)))
        for nm, a in (('c', c), ('n1', n1), ('far', far), ('sum', c + n1 + far)): print(f'  {nm:3s}' + ' '.join(f'{a[i]:5.0f}' for i in range(2, 30, 2)))
        print('  out/sum' + ''.join(f'{(n1[i] + far[i]) / max(c[i] + n1[i] + far[i], 1):5.2f} ' for i in range(2, 30, 2)))
json.dump(out, open('vert_profiles.json', 'w'))
