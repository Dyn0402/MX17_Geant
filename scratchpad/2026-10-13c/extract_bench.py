#!/usr/bin/env python3
"""Generic bench extractor for the pre-June chambers (det0 April, old det3 May).

Per subrun and view, for every M3-good ray (chi2<1, the June recipe), pointing
anywhere: the threshold-free peak-strip max (same observable as nTof_x17
mx17_sim_wft/hv_slope/extract.py), the +-1 neighbour max, the ray position at
the detector plane and the ray slopes (dx/dz, dy/dz from positions at z and
z+100 mm), plus the peak strip's +-2 neighbourhood waveforms (5 x 32, pedestal
and per-64 common-mode subtracted) for shape studies.

Reads decoded_root (flat sample-major vectors) or filtered_root
(uint16[512][32], channel-major: verified identical to decoded_root transposed,
2026-10-13). Pedestal = per-channel median over the first 300 events, CNS =
per-64-channel median per sample: the wft.io.FeuReader recipe.

Optional track window (needed on det0, whose X FEU has hot channels that win a
threshold-free max): ALIGN = 'x:OFF:SIGN:REF,y:OFF:SIGN:REF' maps a view's strip
position to the M3 ray, pos = OFF + SIGN*REF at the detector plane (fitted from
the data, see README). The peak is then searched only among strips within
10 mm + |tan|*30 mm of that prediction, and channels whose noise exceeds 3x the
view's median noise are masked.

  python extract_bench.py BASE RUN SUBRUN DET_Z OUTPREFIX [ALIGN]
"""
import glob, os, re, sys
import numpy as np, pandas as pd, uproot
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'):
    sys.path.insert(0, p)
from qa_config import _Config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
from wft.io import strip_position_map
from M3RefTracking import M3RefTracking

N_PED, BLK, CHUNK = 300, 64, 400


def reader(path):
    t = uproot.open(path)['nt']
    filt = 'filtered' in os.path.basename(path)
    def get(lo, hi):
        a = t.arrays(['eventId', 'amplitude', 'ftst'], entry_start=lo, entry_stop=hi, library='np')
        if filt:
            W = np.stack(a['amplitude']).astype(np.float32)               # (n, 512, 32)
        else:
            W = np.stack([x.reshape(-1, 512).T for x in a['amplitude']]).astype(np.float32)
        return a['eventId'].astype(np.int64), W
    n = t.num_entries
    _, P = get(0, min(N_PED, n))
    ped = np.median(P, axis=(0, 2))                                       # (512,)
    S = P - ped[None, :, None]
    S -= np.repeat(np.median(S.reshape(len(S), 512 // BLK, BLK, -1), axis=2), BLK, axis=1)
    reader.noise = 1.4826 * np.median(np.abs(S), axis=(0, 2))
    for lo in range(0, n, CHUNK):
        ids, W = get(lo, min(lo + CHUNK, n))
        W -= ped[None, :, None]
        cm = np.median(W.reshape(len(W), 512 // BLK, BLK, -1), axis=2)    # (n, 8, ns)
        W -= np.repeat(cm, BLK, axis=1)
        yield ids, W


def main():
    base, run, sub, z, outp = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
    align = {}
    if len(sys.argv) > 6:
        for tok in sys.argv[6].split(','):
            v, off, sg, ref = tok.split(':'); align[v] = (float(off), float(sg), ref)
    volt = int(re.search(r'resist_(\d+)V', sub).group(1)) if 'resist_' in sub else -1
    cfg = _Config('x', run, sub, feus=[3, 4], det_z=z, base_path=base, det_name='mx17_1')
    pos = strip_position_map(cfg)
    rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
    x0, y0, ev = map(np.asarray, rays.get_xy_positions(z))
    x1, y1, ev1 = map(np.asarray, rays.get_xy_positions(z + 100.0))
    assert np.array_equal(ev, ev1)
    ok = np.isfinite(x0) & np.isfinite(y0)
    R = {int(e): (a, b, (c - a) / 100.0, (d - b) / 100.0) for e, a, b, c, d in zip(ev[ok], x0[ok], y0[ok], x1[ok], y1[ok])}
    print(f'[{sub}] {volt} V  {len(R):,} M3-good rays', flush=True)
    sub_dir = os.path.join(base, run, sub)
    files = glob.glob(f'{sub_dir}/decoded_root/*.root') or glob.glob(f'{sub_dir}/filtered_root/*.root')
    rows, wfs = [], []
    for view, feu in (('x', 3), ('y', 4)):
        p = pos[feu]; valid = np.isfinite(p)
        chs = np.flatnonzero(valid)[np.argsort(p[valid])]
        rank = np.full(512, -1); rank[chs] = np.arange(len(chs))
        for f in sorted(x for x in files if re.search(rf'_0{feu}(_decoded_array_filtered)?\.root$', x)):
            try:
                ok = os.path.getsize(f) > 0 and 'nt' in [k.split(';')[0] for k in uproot.open(f).keys()]
            except Exception:
                ok = False
            if not ok:
                print(f'   [skip unreadable] {os.path.basename(f)}', flush=True)
                continue
            gen = reader(f)
            for ids, W in gen:
                use = valid.copy()
                if view in align:
                    nz = reader.noise
                    use &= nz < 3 * np.median(nz[valid])
                for k in np.flatnonzero(np.isin(ids, list(R))):
                    w = W[k]; amp = w.max(axis=1)
                    rx, ry, sx, sy = R[int(ids[k])]
                    cand = use
                    if view in align:
                        off, sg, ref = align[view]
                        pred = off + sg * (rx if ref == 'ref_x' else ry)
                        half = 10.0 + abs(sx if ref == 'ref_x' else sy) * 30.0
                        cand = use & (np.abs(np.nan_to_num(p, nan=1e9) - pred) < half)
                        if not cand.any():
                            continue
                    pk = int(np.argmax(np.where(cand, amp, -np.inf))); r = rank[pk]
                    nbr = [chs[r + d] if 0 <= r + d < len(chs) else None for d in (-2, -1, 0, 1, 2)]
                    blk = np.array([w[c] if c is not None else np.full(w.shape[1], np.nan) for c in nbr], np.float32)
                    rows.append(dict(run=run, subrun=sub, volt=volt, view=view, event_id=int(ids[k]),
                                     peak_ch=pk, peak_pos=float(p[pk]), peak_amp=float(amp[pk]),
                                     nb_amp=float(np.nansum(blk[[1, 3]].max(axis=1))),
                                     peak_sample=int(np.argmax(w[pk])),
                                     ref_x=rx, ref_y=ry, tan_x=sx, tan_y=sy))
                    wfs.append(blk[:, :32] if blk.shape[1] >= 32 else np.pad(blk, ((0, 0), (0, 32 - blk.shape[1])), constant_values=np.nan))
        print(f'   {view}: {sum(r["view"] == view for r in rows):,} events', flush=True)
    pd.DataFrame(rows).to_parquet(outp + '.parquet')
    np.save(outp + '_wf.npy', np.array(wfs, np.float16))
    print('wrote', outp)


if __name__ == '__main__':
    main()
