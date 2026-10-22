"""Drift-column arrival profile per det3 session from M3-vertical tracks (2026-10-22, FOLLOWUP §12).
With a fixed (scintillator) trigger, the 9-channel sum waveform of a track vertical in the measured coordinate
(|th_meas| < 2 deg from M3) is the column's arrival-time distribution x the shaping: onset = electrons from the mesh side,
end = electrons from the cathode (30 mm). Per session: mean normalised sum profile, its 5 % onset / 50 % fall / FWHM, the
late/early charge ratio, the peak amplitude, the fraction of SHORT central pulses (FWHM < 350 ns) among all |dt_nb| < 80
events and among M3-vertical ones, drift voltage and start time. The measured M3 axis is found per session from the chip
medians of the track position (FEU chip = 64 channels = ~50 mm). Sessions: det3 May (d3may) + every det3 June survey npz.
Usage: session_profile.py  (-> session_profile.json, small)"""
import os, re, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L, m3lib
RAW = os.path.expanduser('~/x17/nb_asym_raw'); HW = 4; DT = 60.

def xcross(w, f, rising):
    """Vectorised fractional crossing of f * peak on each row of w (first crossing walking away from the peak)."""
    n, ns = w.shape; pk = w.argmax(1); thr = f * w[np.arange(n), pk]; s = np.arange(ns)
    if rising:
        ok = (s[None, :-1] < pk[:, None]) & (w[:, :-1] < thr[:, None]) & (w[:, 1:] >= thr[:, None])
        i = np.where(ok.any(1), ns - 2 - np.argmax(ok[:, ::-1], 1), -1)
    else:
        ok = (s[None, :-1] >= pk[:, None]) & (w[:, :-1] >= thr[:, None]) & (w[:, 1:] < thr[:, None])
        i = np.where(ok.any(1), np.argmax(ok, 1), -1)
    j = np.clip(i, 0, ns - 2); a, b = w[np.arange(n), j], w[np.arange(n), j + 1]
    return np.where(i >= 0, j + (thr - a) / np.where(b != a, b - a, 1), np.nan)

def sessions():
    S = [dict(tag='d3may', npz={v: f'{RAW}/d3may_{v}.npz' for v in 'xy'}, sess='det3 May 5-5 resist_470V_drift_900V', z=232., drift=900)]
    rows = [l.rstrip('\n').split('\t') for l in open(f'{RAW}/survey/index.tsv')]
    by = {}
    for f, sd, sub, det, view, feu in rows:
        if det == 'mx17_3': by.setdefault((sd, sub), {})[view] = f'{RAW}/survey/' + os.path.basename(f)
    for (sd, sub), npz in by.items():
        cfg = json.load(open(f'{sd}/run_config.json')); d = [x for x in cfg['detectors'] if x['name'] == 'mx17_3'][0]
        card, ch = d['hv_channels']['drift']; sr = [s for s in cfg['sub_runs'] if s['sub_run_name'] == sub][0]
        S.append(dict(tag=os.path.basename(npz['x'])[:4], npz=npz, sess=f'{os.path.basename(sd)} {sub}', z=float(d['det_center_coords']['z']),
                      drift=sr['hvs'][str(card)].get(str(ch))))
    return S

def one(s, view):
    d = L.load(s['npz'][view]); files = d['z']['files']; T = {i: m3lib.m3_table(fn, s['z']) for i, fn in enumerate(files)}
    t0 = re.search(r'datrun_(\d{6})_(\d\d)H(\d\d)', str(files[0])).groups()
    M = np.array([T[fi].get(int(e), [-1, -1] + [np.nan] * 6) for fi, e in zip(d['fi'], d['ev'])]); p = d['p']; g = np.isfinite(M[:, 2])
    spread = [np.std([np.median(M[g & (p // 64 == c), k]) for c in range(8) if (g & (p // 64 == c)).sum() > 15]) for k in (4, 5)]
    km = 4 if spread[0] > spread[1] else 5; thm = M[:, km - 2]
    W = d['W'].astype(float); c = W[:, HW]; fw = (xcross(c, .5, False) - xcross(c, .5, True)) * DT
    ts = np.arange(W.shape[2]); pos = np.clip(W, 0, None); tc = (pos * ts).sum(2) / np.maximum(pos.sum(2), 1e-9)
    dtn = (tc[:, HW + 1] - tc[:, HW - 1]) * DT; v = np.abs(dtn) < 80; vert = g & (np.abs(thm) < 2)
    S9 = W.sum(1); A = S9.max(1); prof = (S9[vert] / A[vert, None]).mean(0) if vert.sum() >= 5 else None
    r = dict(view=view, axis='xy'[km - 4], spread=[float(x) for x in spread], start=f'20{t0[0][:2]}-{t0[0][2:4]}-{t0[0][4:]} {t0[1]}:{t0[2]}',
             n=int(len(c)), n_vert=int(vert.sum()), short_dtn80=float(np.mean(fw[v] < 350)) if v.sum() else None,
             short_vert=float(np.mean(fw[vert] < 350)) if vert.sum() else None, fwhm_vert=float(np.nanmedian(fw[vert])) if vert.sum() else None,
             amp9_vert=float(np.median(A[vert])) if vert.sum() else None)
    if prof is not None:
        on, f50, r50 = (float(xcross(prof[None], f, rr)[0]) for f, rr in ((.05, True), (.5, False), (.5, True)))
        r.update(prof=[round(float(x), 4) for x in prof], onset=on, f50=f50, fwhm_prof=(f50 - r50) * DT)
    return r

if __name__ == '__main__':
    out = []
    for s in sessions():
        rec = dict(tag=s['tag'], sess=s['sess'], drift=s['drift'], z=s['z'])
        for view in 'xy':
            if view not in s['npz'] or not os.path.exists(s['npz'][view]): continue
            try: rec[view] = one(s, view)
            except Exception as e: rec[view] = dict(err=str(e))
        out.append(rec); x = rec.get('x', {})
        print(f"{rec['tag']:5s} {x.get('start', ''):16s} drift {str(rec['drift']):>5s}  {rec['sess'][-45:]:45s} " + '  '.join(
            f"{v.upper()}[{rec[v].get('axis', '?')}] nV {rec[v].get('n_vert', 0):4d} shortV {rec[v].get('short_vert') or np.nan:4.2f} short80 {rec[v].get('short_dtn80') or np.nan:4.2f} "
            f"fwV {rec[v].get('fwhm_vert') or np.nan:4.0f} on {rec[v].get('onset', np.nan):4.1f} f50 {rec[v].get('f50', np.nan):4.1f} A9 {rec[v].get('amp9_vert') or np.nan:5.0f}"
            for v in 'xy' if v in rec), flush=True)
    json.dump(out, open(os.path.join(HERE, 'session_profile.json'), 'w'))
