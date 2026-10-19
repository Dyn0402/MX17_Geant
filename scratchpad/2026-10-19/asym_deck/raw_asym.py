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
Usage: raw_asym.py DECODED_DIR FEU OUT.npz [--files 3] [--skip 0]"""
import sys, glob, argparse, numpy as np, uproot, awkward as ak
ap = argparse.ArgumentParser(); ap.add_argument('dec'); ap.add_argument('feu', type=int); ap.add_argument('out')
ap.add_argument('--files', type=int, default=3); ap.add_argument('--skip', type=int, default=0)
a = ap.parse_args()
files = [f for f in sorted(glob.glob(f'{a.dec}/*_{a.feu:02d}.root')) if 'pedestal' not in f and 'decoded_array' not in f]
files = files[a.skip:a.skip + a.files]
HW, D, AINJ = 4, 6, 0.15
WF, P, EV, FI, DISP = [], [], [], [], []
cs = {k: np.zeros(D + 1) for k in ('raw', 'cm', 'inj')}; cn = {k: np.zeros(D + 1) for k in cs}; cq = {k: np.zeros(D + 1) for k in cs}
nquiet = 0
def corr_acc(x, key):
    # x: (n_ev, ns, 8, 64) ; accumulate sum x_n x_n+d and normalisers over chips/samples/events
    for d in range(D + 1):
        u, v = x[..., :64 - d], x[..., d:]
        cs[key][d] += (u * v).sum(); cn[key][d] += (u * u).sum(); cq[key][d] += (v * v).sum()
for fi, f in enumerate(files):
    t = uproot.open(f)['nt']; ped = None
    for arr in t.iterate(['eventId', 'amplitude'], step_size=1500, library='ak'):
        amp = arr['amplitude']; n = ak.num(amp); ns_mode = int(np.bincount(np.asarray(n) // 512).argmax())
        keep = np.asarray(n) == ns_mode * 512
        x = ak.to_numpy(amp[keep]).reshape(-1, ns_mode, 512).astype(np.float32); eid = np.asarray(arr['eventId'])[keep]
        if ped is None:
            ped = np.median(x, axis=(0, 1)); r0 = x - ped
            sig = 1.4826 * np.median(np.abs(r0 - np.median(r0, axis=(0, 1))), axis=(0, 1)); sig[sig <= 0] = np.inf
        r = x - ped                                        # (ev, ns, 512) raw - pedestal
        A = r.max(axis=1); over = A / sig > 5; nov = over.sum(1)
        q = nov == 0
        if q.any() and nquiet < 6000:
            xq = r[q][:, :, :].reshape(-1, ns_mode, 8, 64); nquiet += len(xq)
            corr_acc(xq, 'raw')
            xc = xq - np.median(xq, axis=3, keepdims=True); corr_acc(xc, 'cm')
            xi = xc.copy(); xi[..., 1:] += AINJ * xc[..., :-1]; corr_acc(xi, 'inj')
        for i in np.flatnonzero((nov >= 1) & (nov <= 12)):
            c = np.flatnonzero(over[i]); p = int(c[np.argmax(A[i, c])])
            if not (4 <= p % 64 <= 59) or c.min() < p - 3 or c.max() > p + 3: continue
            blk = r[i, :32, p - HW:p + HW + 1].T
            if blk.shape[1] < 32: continue
            WF.append(blk.astype(np.float32)); P.append(p); EV.append(int(eid[i])); FI.append(fi)
            if len(DISP) < 40 and 600 < A[i, p] < 2500: DISP.append(r[i, :32, (p // 64) * 64:(p // 64) * 64 + 64].T.astype(np.float32))
    print(f, len(WF), 'kept', nquiet, 'quiet', flush=True)
corr = {k: cs[k] / np.sqrt(cn[k] * cq[k]) for k in cs}
np.savez_compressed(a.out, wf=np.array(WF, np.float32), p=np.array(P, np.int16), ev=np.array(EV, np.int64), fi=np.array(FI, np.int8),
                    disp=np.array(DISP, np.float32), ped=ped, sig=sig, files=np.array(files), nquiet=nquiet,
                    corr_raw=corr['raw'], corr_cm=corr['cm'], corr_inj=corr['inj'])
print('DONE', a.out, len(WF), {k: np.round(v[1:4], 3).tolist() for k, v in corr.items()})
