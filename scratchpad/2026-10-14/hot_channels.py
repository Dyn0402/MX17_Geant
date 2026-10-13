#!/usr/bin/env python3
"""Sim noise model: which channels cross 5 sigma far (>60 strips) from the track, and how often.
yoga: hot_channels.py <sim dir> <data dir>"""
import sys, os, glob, numpy as np, collections
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17')]
from wft.io import FeuReader
for d in sys.argv[1:]:
  for feu in (7, 8):
    f = sorted(glob.glob(f'{d}/decoded_root/*_{feu:02d}.root'))[0]
    rdr = FeuReader(f); noise = np.where(rdr.noise > 0, rdr.noise, np.inf)
    cnt = collections.Counter(); n = 0; mx_all = []
    for eid, _f, w in rdr.iter_events(set(list(map(int, rdr.event_ids))[:600])):
        amp = w.max(axis=1); sig = amp / noise; pk = np.argmax(amp)
        far = np.flatnonzero((sig >= 5) & (np.abs(np.arange(len(amp)) - pk) > 60))
        cnt.update(far.tolist()); n += 1; mx_all.append(sig)
    S = np.array(mx_all)
    print(os.path.basename(d.rstrip('/')), feu, 'events', n, 'n_samples', w.shape, 'noise pctl', np.percentile(rdr.noise, [1, 10, 50, 90, 99]).round(2),
          'top far channels', cnt.most_common(8), 'distinct', len(cnt))
    print('   median over events of max/σ per channel, percentiles over channels:', np.percentile(np.median(S, 0), [1, 10, 50, 90, 99]).round(2))
