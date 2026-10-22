"""Column-collection health of every bench detector-session with M3 tracks (2026-10-22, §13). For M3-vertical tracks
(|th_meas| < 2 deg; the measured M3 axis found from chip medians): 9-pad-sum 30 % onset and end percentiles (60 ns samples),
fraction of columns ending before sample 14, the 9-pad peak amplitude A9, central FWHM < 350 ns fraction, alpha+1-like
sharing s1, plus stack slot (z), HV channels and drift voltage from run_config. Long/longer runs only unless TAGS given.
Usage: [TAGS=s001,s004] slot_table.py   (-> slot_table.json)"""
import os, sys, json, numpy as np
import session_profile as SP, event_ends as EE, m3lib
RAW = SP.RAW

def sessions():
    rows = [l.rstrip('\n').split('\t') for l in open(f'{RAW}/survey/index.tsv')]; by = {}
    for f, sd, sub, det, view, feu in rows:
        by.setdefault((sd, sub, det), {})[view] = f'{RAW}/survey/' + os.path.basename(f)
    out = []
    for (sd, sub, det), npz in by.items():
        cfg = json.load(open(f'{sd}/run_config.json')); d = [x for x in cfg['detectors'] if x['name'] == det]
        if not d or 'x' not in npz or not os.path.exists(npz['x']): continue
        d = d[0]; card, ch = d['hv_channels']['drift']; sr = [s for s in cfg.get('sub_runs', []) if s['sub_run_name'] == sub]
        drift = sr[0]['hvs'].get(str(card), {}).get(str(ch)) if sr else None
        out.append(dict(tag=os.path.basename(npz['x'])[:4], npz=npz, run=os.path.basename(sd), sub=sub, det=det, z=float(d['det_center_coords']['z']),
                        hv=d['hv_channels'], drift=drift, desc=d.get('description', '')[:45]))
    return out

def health(s, view='x'):
    d = SP.L.load(s['npz'][view]); files = d['z']['files']; T = {i: m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
    if any(t is None for t in T.values()): return None
    M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
    sp = [np.std([np.median(M[g & (p // 64 == c), k]) for c in range(8) if (g & (p // 64 == c)).sum() > 15] or [0]) for k in (4, 5)]
    thm = M[:, 2 if sp[0] > sp[1] else 3]; v = g & (np.abs(thm) < 2)
    if v.sum() < 15: return dict(n_vert=int(v.sum()))
    W = d['W'][v].astype(float); S9 = W.sum(1); A9 = S9.max(1); on = np.argmax(S9 >= .3 * A9[:, None], 1); e9 = EE.last_cross(S9, .3)
    c = W[:, 4]; fw = (SP.xcross(c, .5, False) - SP.xcross(c, .5, True)) * 60.; s1 = np.einsum('ij,ij->i', W[:, 3] + W[:, 5], c) / np.einsum('ij,ij->i', c, c)
    q = lambda x: [float(u) for u in np.percentile(x, [10, 50, 90])]
    return dict(n_vert=int(v.sum()), axis='xy'[0 if sp[0] > sp[1] else 1], on=q(on), end=q(e9), early_end=float(np.mean(e9 < 14)), A9=float(np.median(A9)),
                short=float(np.nanmean(fw < 350)), s1=float(np.median(s1)))

if __name__ == '__main__':
    want = os.environ.get('TAGS'); out = []
    print('tag  date/run                          sub                  det     z    drift  HV(drift,resist)   nV   onset p10/50/90   end p10/50/90   end<14  short   A9     s1')
    for s in sessions():
        if want and s['tag'] not in want.split(','): continue
        if not want and not s['sub'].startswith('long'): continue
        h = health(s)
        if not h or h['n_vert'] < 15: continue
        out.append(dict(**{k: v for k, v in s.items() if k != 'npz'}, **h))
        print(f"{s['tag']} {s['run'][-32:]:32s} {s['sub'][:20]:20s} {s['det']} {s['z']:5.0f} {str(s['drift']):>6s}  {str(s['hv']['drift']):7s}{str(s['hv']['resist']):8s}  {h['n_vert']:4d}   "
              f"{'/'.join(f'{x:.0f}' for x in h['on']):11s}      {'/'.join(f'{x:.0f}' for x in h['end']):11s}   {h['early_end']:.2f}   {h['short']:.2f}  {h['A9']:5.0f}  {h['s1']:.2f}", flush=True)
    json.dump(out, open('slot_table.json', 'w'), indent=1)
