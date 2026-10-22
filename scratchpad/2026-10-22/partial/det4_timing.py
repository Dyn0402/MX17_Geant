"""6-23: does det4 (same events, z = 702 mm, FEU 6 = X) see the same late/random onsets as det3? (2026-10-22, §12.)
det4 per event: pedestal-subtracted (median of the first 1500 events) raw ADC, peak channel, sum of +-4 channels around it,
30 % onset and end. Joined to det3's raw npz (X view, file 000) by eventId. Same for the 6-22 det2+det3 overnight? (det2 FEU)
Usage: det4_timing.py   (prints; small)"""
import sys, numpy as np, uproot, awkward as ak
import session_profile as SP, event_ends as EE
RUN = '/home/dylan/x17/cosmic_bench/det3_det4/mx17_det3_det4_overnight_6-23-26/long_run/decoded_root/MX17_long_run_datrun_260624_03H16_000_{:02d}.root'

def other_det(fn, evs):
    t = uproot.open(fn)['nt']; ped = None; out = {}
    for arr in t.iterate(['eventId', 'amplitude'], step_size=1500, library='ak'):
        n = np.asarray(ak.num(arr['amplitude'])); keep = n == 32 * 512
        x = ak.to_numpy(arr['amplitude'][keep]).reshape(-1, 32, 512).astype(np.float32); eid = np.asarray(arr['eventId'])[keep]
        if ped is None: ped = np.median(x, axis=(0, 1))
        r = x - ped; sel = np.isin(eid, evs)
        for i in np.flatnonzero(sel):
            A = r[i].max(0); p = int(np.clip(A.argmax(), 4, 507)); s9 = r[i][:, p - 4:p + 5].sum(1)
            out[int(eid[i])] = (A[p], s9)
    return out

s = {x['tag']: x for x in SP.sessions()}['s052']; d = SP.L.load(s['npz']['x']); f0 = d['fi'] == 0
W = d['W'][f0].astype(float); S9 = W.sum(1); A9 = S9.max(1); on3 = np.argmax(S9 >= .3 * A9[:, None], 1); e3 = EE.last_cross(S9, .3)
o4 = other_det(RUN.format(6), d['ev'][f0]); ev = d['ev'][f0]
has = np.array([e in o4 and o4[e][0] > 150 for e in ev]); S4 = np.array([o4[e][1] if e in o4 else np.zeros(32) for e in ev])
A4 = S4.max(1); on4 = np.argmax(S4 >= .3 * A4[:, None], 1); e4 = EE.last_cross(S4, .3)
print(f'det3 X file 000 events {len(ev)}, det4 X peak > 150 ADC: {has.mean():.2f}')
print(f'det4 onset p10/50/90 {np.percentile(on4[has], [10, 50, 90])}  width p10/50/90 {np.percentile((e4 - on4)[has], [10, 50, 90])}')
print(f'det3 onset p10/50/90 {np.percentile(on3[has], [10, 50, 90])}  width p10/50/90 {np.percentile((e3 - on3)[has], [10, 50, 90])}')
print('corr(det3 onset, det4 onset) %.2f' % np.corrcoef(on3[has], on4[has])[0, 1])
for lo, hi in ((0, 7), (7, 10), (10, 14), (14, 32)):
    m = has & (on3 >= lo) & (on3 < hi)
    if m.sum() > 10: print(f'   det3 onset {lo:2d}-{hi:2d} (n {m.sum():4d}): det4 onset p50 {np.median(on4[m]):.1f}  p10/90 {np.percentile(on4[m], [10, 90])}')
