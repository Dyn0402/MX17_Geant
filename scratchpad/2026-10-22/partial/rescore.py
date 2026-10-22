"""Step 2 of §12: det3 June X timing re-scored on CLEAN sessions at the sim's drift field (2026-10-22).
The §8-10 det3 June reference (d3jun = 6-23 long_run, 600 V) is a degraded session (§12: columns collected as a ~0.5 us
burst at a random depth, 9-ch amplitude 4-5x low and flat vs resist HV, alpha+1 0.13). Here det3 June = 6-27 long_run
490 V / 1000 V (s051) + 6-27 overnight sanity run 490 V / 1000 V (s025); det3 May (900 V) and 6-23 for reference.
Neighbour/central ratio r = (W+1 + W-1) / Wc at the central's 5 % and 25 % rise, peak and 50 % fall (ratio of sums),
in absolute central-FWHM x s1 cells (cloud/fwhm_bins.py), two selections:
  dtn : |dt_nb| < 80 ns (the §8-11 cut)         vs sim glg0 all angles (same cut)
  m3v : M3 |th_meas| < 2 deg (true near-vertical) vs sim glg0 th00 only (no cut)
Usage: rescore.py   (-> rescore.json)"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import prompt_share as PS, session_profile as SP
RAW = SP.RAW; HW = 4; DT = 60.; FWB = (0, 350, 450, 550, 650, 2000); S1B = ((.2, .45), (.45, .7))
SESS = {s['tag']: s for s in SP.sessions()}

def prep(W):
    W = W.astype(float); c = W[:, HW]; pk = c.argmax(1); n2 = np.einsum('ij,ij->i', c, c)
    s1 = np.einsum('ij,ij->i', W[:, HW + 1] + W[:, HW - 1], c) / n2
    ts = np.arange(W.shape[2]); pos = np.clip(W, 0, None); tc = (pos * ts).sum(2) / np.maximum(pos.sum(2), 1e-9)
    dtn = (tc[:, HW + 1] - tc[:, HW - 1]) * DT; fw = (SP.xcross(c, .5, False) - SP.xcross(c, .5, True)) * DT
    pts = dict(r05=SP.xcross(c, .05, True), r25=SP.xcross(c, .25, True), peak=pk.astype(float), f50=SP.xcross(c, .5, False))
    return dict(W=W, s1=s1, dtn=dtn, fw=fw, pts=pts)

def data_set(tags, view='x'):
    out = []
    for t in tags:
        s = SESS[t]; d = SP.L.load(s['npz'][view]); p = d['p'] % 512; keep = (p % 64 >= 4) & (p % 64 <= 59)
        files = d['z']['files']; T = {i: SP.m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
        M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); g = np.isfinite(M[:, 2]); pp = d['p']
        spread = [np.std([np.median(M[g & (pp // 64 == c), k]) for c in range(8) if (g & (pp // 64 == c)).sum() > 15]) for k in (4, 5)]
        thm = M[:, 2 if spread[0] > spread[1] else 3]
        r = prep(d['W'][keep]); r['m3v'] = (np.abs(thm[keep]) < 2); out.append(r)
    return {k: (np.concatenate([o[k] for o in out]) if k != 'pts' else {q: np.concatenate([o['pts'][q] for o in out]) for q in out[0]['pts']}) for k in out[0]}

def sim_set(pat):
    W = []
    for f in sorted(glob.glob(pat)):
        d = SP.L.load(f); p = d['p'] % 512; W.append(d['W'][(p % 64 >= 4) & (p % 64 <= 59)])
    r = prep(np.concatenate(W)); r['m3v'] = np.ones(len(r['s1']), bool); return r

def table(name, D, sel):
    S = D['W'][:, HW - 1] + D['W'][:, HW + 1]; c = D['W'][:, HW]; base = (np.abs(D['dtn']) < 80) if sel == 'dtn' else D['m3v']
    rows = []; print(f'\n{name} [{sel}]  n {base.sum()}\n   FWHM ns     s1 bin       n    r05          r25          peak         f50          early/pk late/pk')
    for lo, hi in zip(FWB[:-1], FWB[1:]):
        for slo, shi in S1B:
            m = base & (D['fw'] >= lo) & (D['fw'] < hi) & (D['s1'] >= slo) & (D['s1'] < shi)
            if m.sum() < 20: continue
            v = {}
            for k, t in D['pts'].items():
                nu, de = PS.at(S, t), PS.at(c, t); gg = m & np.isfinite(nu) & np.isfinite(de); ii = np.flatnonzero(gg)
                bs = PS.RNG.integers(0, len(ii), (200, len(ii))); v[k] = (float(nu[gg].sum() / de[gg].sum()), float(np.std(nu[ii][bs].sum(1) / de[ii][bs].sum(1))))
            rows.append(dict(fw=[lo, hi], s1=[slo, shi], n=int(m.sum()), **v))
            print(f'   {lo:4d}-{hi:<5d}  {slo:.2f}-{shi:.2f}  {m.sum():5d}  ' + '  '.join(f'{v[k][0]:.3f}±{v[k][1]:.3f}' for k in D['pts'])
                  + f'   {v["r05"][0] / v["peak"][0]:.2f}  {v["f50"][0] / v["peak"][0]:.2f}')
    fwv = D['fw'][base]; print(f'   FWHM p25/50/75 {np.nanpercentile(fwv, [25, 50, 75]).round(0)}   fraction < 350 ns {np.mean(fwv < 350):.2f}')
    return dict(rows=rows, fw_q=np.nanpercentile(fwv, [25, 50, 75]).tolist(), short=float(np.mean(fwv < 350)), n=int(base.sum()))

if __name__ == '__main__':
    R = {}
    sets = {'det3 June 6-27 1000 V (s051+s025)': data_set(['s051', 's025']), 'det3 June 6-23 600 V (d3jun)': data_set(['s052']),
            'det3 May 900 V': data_set(['d3may'])}
    sims = {'sim glg0 (all angles)': sim_set(f'{RAW}/simgapleak/glg0_th*_x.npz'), 'sim glg0 th00': sim_set(f'{RAW}/simgapleak/glg0_th00_*_x.npz')}
    for n, D in sets.items():
        for sel in ('dtn', 'm3v'): R[f'{n} [{sel}]'] = table(n, D, sel)
    R['sim glg0 [dtn]'] = table('sim glg0 (all angles)', sims['sim glg0 (all angles)'], 'dtn')
    R['sim glg0 th00 [m3v]'] = table('sim glg0 th00', sims['sim glg0 th00'], 'm3v')
    json.dump(R, open(os.path.join(HERE, 'rescore.json'), 'w'), indent=1)
