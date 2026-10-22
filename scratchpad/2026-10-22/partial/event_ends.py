"""Per-event drift-column length of M3-vertical tracks (2026-10-22, §12): 9-channel sum and central pad, onset (10 %) and
end (last crossing of 30 % of the event's peak, walking back from the window end), per session; plus the end vs amplitude.
A per-event truncation (column cut at a random depth) shows up as a broad, early end distribution in the 9-ch sum itself.
Usage: event_ends.py TAG ...   (-> event_ends.json)"""
import sys, json, numpy as np
import session_profile as SP
S = {s['tag']: s for s in SP.sessions()}; out = {}

def last_cross(w, f):
    pk = w.max(1, keepdims=True); above = w >= f * pk; ns = w.shape[1]
    return ns - 1 - np.argmax(above[:, ::-1], 1)

if __name__ == '__main__':
  print('tag  view  n   sum9: onset p50   end p10/p50/p90 (samples)   central end p10/p50/p90   corr(end, A9)   A9 p50   frac end9 < 14')
  for tag in sys.argv[1:]:
      for view in 'xy':
          s = S[tag]; d = SP.L.load(s['npz'][view]); files = d['z']['files']; T = {i: SP.m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
          M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
          spread = [np.std([np.median(M[g & (p // 64 == c), k]) for c in range(8) if (g & (p // 64 == c)).sum() > 15]) for k in (4, 5)]
          thm = M[:, 2 if spread[0] > spread[1] else 3]; v = g & (np.abs(thm) < 2); W = d['W'][v].astype(float)
          S9 = W.sum(1); c = W[:, 4]; A9 = S9.max(1)
          on = np.argmax(S9 >= .1 * A9[:, None], 1); e9 = last_cross(S9, .3); ec = last_cross(c, .3)
          q = lambda x: np.percentile(x, [10, 50, 90])
          out[f'{tag}_{view}'] = dict(sess=s['sess'], n=int(v.sum()), on9=float(np.median(on)), end9=q(e9).tolist(), endc=q(ec).tolist(),
                                      hist_end9=np.bincount(e9, minlength=32).tolist(), r=float(np.corrcoef(e9, A9)[0, 1]), A9=float(np.median(A9)))
          print(f"{tag:5s} {view}  {v.sum():4d}        {np.median(on):4.1f}        {q(e9)[0]:4.1f} {q(e9)[1]:4.1f} {q(e9)[2]:4.1f}                  "
                f"{q(ec)[0]:4.1f} {q(ec)[1]:4.1f} {q(ec)[2]:4.1f}         {np.corrcoef(e9, A9)[0, 1]:+.2f}       {np.median(A9):5.0f}   {np.mean(e9 < 14):.2f}   {s['sess'][-40:]}")
  json.dump(out, open('event_ends.json', 'w'))
