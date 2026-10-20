"""Raw-ADC neighbour extractor for the X-asymmetry deck (2026-10-19). Independent of wft/FeuReader on purpose:
reads the DREAM decoded_root 'amplitude' (raw 12-bit ADC, the run had no ZS, no on-board pedestal or common-noise
subtraction), subtracts a per-channel pedestal (median over the first 1000 events, all samples) and nothing else.
No strip map, no clustering, no common-mode subtraction: everything is in DAQ channel numbers.
Per event: A_ch = max over samples; channels over 5 sigma (sigma = per-channel MAD noise); keep events with 1..12 channels
over threshold, all within p-3..p+3 of the peak channel p, and p at 4..59 inside its 64-channel DREAM chip.
Stores the 9x32 raw waveforms of channels p-4..p+4 (channel order, NOT position order).
Noise: on quiet events (no channel over 5 sigma) accumulates the channel-distance correlation corr(x_n, x_n+d), d = 1..6,
within each chip, raw and after a per-chip-per-sample median (common-mode) subtraction, plus the same with a copy
x_n+1 += 0.15 x_n injected, to show what an electronics-side copy would do to the noise.
--cm (beam data): also subtract a per-chip, per-sample median (common mode) before thresholding, which removes the
gamma-flash ringing, and cut each 32-sample window around the peak channel's maximum instead of samples 0-31;
channels noisier than 3x the median are masked and the isolation cut is applied within the peak's own chip.
--edge (2026-10-20): keep ONLY peaks within 4 channels of a 64-channel chip/connector edge, with the isolation cut and
the peak search done in TRUE pad-column order (each 64-channel connector reversed); stores the column-order window
(wf, physical neighbours) and the DAQ-order window around the same channel (wf_daq). At a connector edge the two
differ: DAQ neighbours sit on another connector/chip 127 columns away, physical ones are 127 DAQ channels away.
--array (2026-10-20): read decoded_array / filtered_root files (amplitude uint16[512][32], channel-major), e.g. the
April-24 build on EOS, which has no decoded_root. Filtered files hold M3-triggered events only; the first-1000-event
pedestal is then a median over events that mostly carry a track in some chip, which the median tolerates.
--mix-with DIR2 --mix-w W (2026-10-20): sim only. Reads the same events from a second arm that differs only in the
X response (digitizer --x-view-offset-um; same seed, so identical noise and Y) and uses (1-W)*A + W*A2 sample by
sample: a fraction W of each avalanche's X induction read at the shifted position (linear induction, linear DAQ).
Usage: raw_asym.py DECODED_DIR FEU OUT.npz [--files 3] [--skip 0] [--cm] [--edge] [--array] [--mix-with DIR2 --mix-w W]"""
import sys, glob, argparse, numpy as np, uproot, awkward as ak
ap = argparse.ArgumentParser(); ap.add_argument('dec'); ap.add_argument('feu', type=int); ap.add_argument('out')
ap.add_argument('--files', type=int, default=3); ap.add_argument('--skip', type=int, default=0)
ap.add_argument('--cm', action='store_true'); ap.add_argument('--edge', action='store_true')
ap.add_argument('--array', action='store_true')
ap.add_argument('--mix-with', default=None); ap.add_argument('--mix-w', type=float, default=0.0)
a = ap.parse_args()
if a.array:   # decoded_array / filtered_root files: amplitude is a fixed uint16[512][32] (channel, sample) per event
    files = sorted(glob.glob(f'{a.dec}/*_{a.feu:02d}_decoded_array*.root'))
else:
    files = [f for f in sorted(glob.glob(f'{a.dec}/*_{a.feu:02d}.root')) if 'pedestal' not in f and 'decoded_array' not in f]
files = files[a.skip:a.skip + a.files]
if a.mix_with:   # same files (same events, same noise) of another sim arm, matched by sorted order
    files2 = [f for f in sorted(glob.glob(f'{a.mix_with}/*_{a.feu:02d}.root')) if 'pedestal' not in f and 'decoded_array' not in f]
    files2 = files2[a.skip:a.skip + a.files]; assert len(files2) == len(files)
HW, D, AINJ = 4, 6, 0.15
WF, P, EV, FI, DISP, WFD = [], [], [], [], [], []
COL = np.array([64 * (c // 64) + 63 - c % 64 for c in range(512)])   # DAQ ch -> column; its own inverse
cs = {k: np.zeros(D + 1) for k in ('raw', 'cm', 'inj')}; cn = {k: np.zeros(D + 1) for k in cs}; cq = {k: np.zeros(D + 1) for k in cs}
nquiet = 0
def corr_acc(x, key):
    # x: (n_ev, ns, 8, 64) ; accumulate sum x_n x_n+d and normalisers over chips/samples/events
    for d in range(D + 1):
        u, v = x[..., :64 - d], x[..., d:]
        cs[key][d] += (u * v).sum(); cn[key][d] += (u * u).sum(); cq[key][d] += (v * v).sum()
for fi, f in enumerate(files):
    t = uproot.open(f)['nt']; ped = None
    it2 = uproot.open(files2[fi])['nt'].iterate(['eventId', 'amplitude'], step_size=1500, library='ak') if a.mix_with else None
    for arr in t.iterate(['eventId', 'amplitude'], step_size=1500, library='ak'):
        if it2 is not None:
            arr2 = next(it2); assert np.array_equal(np.asarray(arr['eventId']), np.asarray(arr2['eventId']))
            amp = (1 - a.mix_w) * arr['amplitude'] + a.mix_w * arr2['amplitude']
        elif a.array:
            amp = ak.Array(ak.to_numpy(arr['amplitude']).transpose(0, 2, 1).reshape(len(arr), -1))
        else:
            amp = arr['amplitude']
        n = ak.num(amp); ns_mode = int(np.bincount(np.asarray(n) // 512).argmax())
        keep = np.asarray(n) == ns_mode * 512
        x = ak.to_numpy(amp[keep]).reshape(-1, ns_mode, 512).astype(np.float32); eid = np.asarray(arr['eventId'])[keep]
        if ped is None:
            ped = np.median(x, axis=(0, 1)); r0 = x - ped
            if a.cm:
                r0 = r0 - np.repeat(np.median(r0.reshape(len(r0), ns_mode, 8, 64), axis=3), 64, axis=2)
            sig = 1.4826 * np.median(np.abs(r0 - np.median(r0, axis=(0, 1))), axis=(0, 1)); sig[sig <= 0] = np.inf
            if a.cm:
                sig[sig > 3 * np.median(sig[np.isfinite(sig)])] = np.inf
        r = x - ped                                        # (ev, ns, 512) raw - pedestal
        if a.cm:                                           # (beam) minus the per-chip, per-sample median
            r = r - np.repeat(np.median(r.reshape(len(r), ns_mode, 8, 64), axis=3), 64, axis=2)
        A = r.max(axis=1); over = A / sig > 5; nov = over.sum(1)
        q = nov == 0
        if q.any() and nquiet < 6000:
            xq = r[q][:, :, :].reshape(-1, ns_mode, 8, 64); nquiet += len(xq)
            corr_acc(xq, 'raw')
            xc = xq - np.median(xq, axis=3, keepdims=True); corr_acc(xc, 'cm')
            xi = xc.copy(); xi[..., 1:] += AINJ * xc[..., :-1]; corr_acc(xi, 'inj')
        cand = np.flatnonzero(nov >= 1) if a.cm else np.flatnonzero((nov >= 1) & (nov <= 12))
        if a.edge:
            for i in cand:
                cc = COL[np.flatnonzero(over[i])]; Ac = A[i, COL]               # A in column order
                pc = int(cc[np.argmax(Ac[cc])])
                if not (4 <= pc <= 507) or (3 < pc % 64 < 60) or cc.min() < pc - 3 or cc.max() > pc + 3: continue
                pd = int(COL[pc])
                if not (4 <= pd <= 507): continue
                WF.append(r[i, :32][:, COL[pc - HW:pc + HW + 1]].T.astype(np.float32)); WFD.append(r[i, :32, pd - HW:pd + HW + 1].T.astype(np.float32))
                P.append(pd); EV.append(int(eid[i])); FI.append(fi)
            continue
        for i in cand:
            c = np.flatnonzero(over[i]); p = int(c[np.argmax(A[i, c])])
            if a.cm:
                c = c[c // 64 == p // 64]
                if len(c) > 12: continue
            if not (4 <= p % 64 <= 59) or c.min() < p - 3 or c.max() > p + 3: continue
            s0 = min(max(int(r[i, :, p].argmax()) - 8, 0), ns_mode - 32) if a.cm else 0
            blk = r[i, s0:s0 + 32, p - HW:p + HW + 1].T
            if blk.shape[1] < 32: continue
            WF.append(blk.astype(np.float32)); P.append(p); EV.append(int(eid[i])); FI.append(fi)
            if len(DISP) < 40 and 600 < A[i, p] < 2500: DISP.append(r[i, :32, (p // 64) * 64:(p // 64) * 64 + 64].T.astype(np.float32))
    print(f, len(WF), 'kept', nquiet, 'quiet', flush=True)
corr = {k: cs[k] / np.sqrt(cn[k] * cq[k]) for k in cs}
np.savez_compressed(a.out, wf=np.array(WF, np.float32), wf_daq=np.array(WFD, np.float32), p=np.array(P, np.int16), ev=np.array(EV, np.int64), fi=np.array(FI, np.int8),
                    disp=np.array(DISP, np.float32), ped=ped, sig=sig, files=np.array(files), nquiet=nquiet,
                    corr_raw=corr['raw'], corr_cm=corr['cm'], corr_inj=corr['inj'])
print('DONE', a.out, len(WF), {k: np.round(v[1:4], 3).tolist() for k, v in corr.items()})
