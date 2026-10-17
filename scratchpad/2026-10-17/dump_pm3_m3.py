"""Peak-strip +-3 neighbour waveforms (position-ordered), vertical tracks, for the neighbour-delay fit
(2026-10-13/xtalk_fit.py). Desktop version of 2026-10-12/laptop_dump_waveforms.py: the data leg is
selected by the TRUE M3 angle (|ax|, |ay| < 3 deg, as ledger_m3.py) instead of the T14 reco window;
the sim leg is every event of a th00 point. Same seeding (5 sigma) and peak-strip choice as t14.
Usage: dump_pm3_m3.py SIM_TH00_DIR OUT.npz   (SIM dir holds *_07.root / *_08.root)"""
import os, sys, glob
import numpy as np, pandas as pd
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis', f'{NT}/mx17_sim_wft'): sys.path.insert(0, p)
from qa_config import get_config
from wft.io import FeuReader, strip_position_map

cfg = get_config('sat_det3'); PM = strip_position_map(cfg)
DEC = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
pk = pd.read_parquet(os.path.expanduser('~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet'))
pk = pk[(pk.ref_x > -200) & (pk.ref_x < 120) & (pk.ref_y.abs() < 170)]
vert = set(pk[(pk.ax.abs() < 3) & (pk.ay.abs() < 3)].event_id.astype(int))

def dump(d, feu, wanted):
    pos = PM[feu]; valid = ~np.isnan(pos)
    order = np.flatnonzero(valid)[np.argsort(pos[valid])]; rank = {c: i for i, c in enumerate(order)}
    W = []
    for path in sorted(glob.glob(f'{d}/*_{feu:02d}.root')):
        r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
        ids = set(int(e) for e in r.event_ids)
        for _e, _t, w in r.iter_events(ids if wanted is None else ids & wanted):
            amp = w.max(axis=1); over = valid & (amp / noise >= 5)
            if not over.any(): continue
            c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])])
            if p not in rank: continue
            i = rank[p]; blk = np.full((7, w.shape[1]), np.nan)
            for k in range(-3, 4):
                if 0 <= i + k < len(order): blk[k + 3] = w[order[i + k]]
            W.append(blk)
    return np.array(W, np.float32)

out = {}
for view, feu in (('x', 7), ('y', 8)):
    out[f'data_{view}_wf'] = dump(DEC, feu, vert)
    out[f'sim_{view}_wf'] = dump(sys.argv[1], feu, None)
    print(view, len(out[f'data_{view}_wf']), len(out[f'sim_{view}_wf']), flush=True)
np.savez_compressed(sys.argv[2], **out)
