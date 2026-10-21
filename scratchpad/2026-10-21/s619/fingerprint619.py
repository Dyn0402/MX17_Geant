"""Which board is 6-19/20's X data? (2026-10-21). Per-channel noise sigma (pedestal rms, 512 ch) of each sub-run is a
board+cable fingerprint (dead/hot/noisy channels, per-connector levels). Correlate the 6-19 extractions
(extra/d619_longer_x, extra/d619_longer_feu3_x) and 6-17 (extra/d617_long_x) with every survey sub-run's X and Y sig,
after normalising each 64-channel block (removes FEU/gain level), and print the best matches. Also a1 per chip.
Usage: NB_RAW=~/x17/nb_asym_raw fingerprint619.py"""
import os, sys, csv, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L
RAW = L.RAW; HERE = os.path.dirname(os.path.abspath(__file__))
meta = {(r['run'], r['sub'], r['det'], r['view']): r for r in csv.DictReader(open(f'{HERE}/../../2026-10-19/copy_survey.csv'))}
def norm(s):
    s = np.asarray(s, float).copy(); s[~np.isfinite(s)] = 0
    b = s.reshape(8, 64); med = np.median(b, 1, keepdims=True); med[med == 0] = 1
    return np.log(np.clip(b / med, .05, 20)).ravel()
ref = []
for l in open(f'{RAW}/survey/index.tsv'):
    f, rd, sub, det, v, feu = l.rstrip('\n').split('\t'); f = f'{RAW}/survey/{os.path.basename(f)}'
    m = meta.get((rd.split('/')[-1], sub, det, v.upper()))
    if not m or not os.path.exists(f): continue
    z = np.load(f); ref.append((f'{m["phys"]}{"(May)" if m["phys"]=="det3" and ("5-" in rd or "6-6" in rd) else ""} {v} FEU{feu} {rd.split("/")[-1][:30]} {sub[:14]}', norm(z['sig'])))
if __name__ == "__main__":
    for q in ('d617_long_x', 'd619_longer_x', 'd619_longer_feu3_x'):
      z = np.load(f'{RAW}/extra/{q}.npz'); s = norm(z['sig'])
      cc = sorted(((float(np.corrcoef(s, r)[0, 1]), name) for name, r in ref), reverse=True)
      print(f'== {q}: best noise-pattern matches')
      for c, name in cc[:6]: print(f'   {c:.3f}  {name}')
      # best per board/view
      best = {}
      for c, name in cc:
          k = ' '.join(name.split()[:2]); best.setdefault(k, (c, name))
      print('   best per board/view: ' + ', '.join(f'{k} {v[0]:.2f}' for k, v in sorted(best.items(), key=lambda t: -t[1][0])))
      d = L.load(f'{RAW}/extra/{q}.npz'); p = d['p'] % 512; ch = p // 64
      print('   a1 per chip (DAQ 64-ch blocks): ' + ' '.join(f'{d["a1"][ch==k].mean():+.3f}({(ch==k).sum()})' if (ch==k).sum() > 50 else '  -  ' for k in range(8)))
