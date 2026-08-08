#!/usr/bin/env python3
"""
truth.py — the per-event truth sidecar (plan §7 step 7, audit C9).

The decoded files are deliberately indistinguishable from data: ADC counts,
pedestal, noise, and nothing else. That is the whole point of T13, and it also
means a simulated run cannot answer "how well did the reconstruction do", only
"did it run". This writes the other half — what actually happened — as a
sidecar keyed by eventId, so any reconstruction output can be joined against it.

WHY NPZ AND NOT PARQUET. The fix order allows either. The response chain runs
on the system python3, which has numpy and uproot but NO pyarrow or pandas
(checked 2026-08-08), so a parquet sidecar could not be written by the process
that produces it. npz needs nothing extra and the S1 products are already npz.
The per-channel arrays are jagged, so they are stored CSR-style — one flat
array plus an offsets array — which is exact and needs no object dtype.

WHAT IS AND IS NOT TRUTH HERE. The charges recorded are the induced charge per
channel BEFORE the DREAM shaper. That distinction used to be cosmetic and is
not any more: with the P6 residual-pole model int(h) = (1 - beta) int(h_nom),
so shaping no longer preserves the integral and "true charge" measured after it
would silently carry the shaper's beta.

    from .truth import TruthWriter
    tw = TruthWriter()
    tw.add(event_id, ...)
    tw.write("sim_truth.npz")
"""

from __future__ import annotations

import numpy as np

SCHEMA = "mx17_truth/1"


class TruthWriter:
    """Accumulates per-event truth, then writes one npz."""

    def __init__(self):
        self.ev = []            # scalars, one row per event
        self.ch_view = []       # jagged, flattened
        self.ch_idx = []
        self.ch_q = []
        self.ch_off = [0]

    def add(self, event_id, *, n_electron_in, n_seed, q_seed_total,
            x_true_mm, y_true_mm, t0_ns, off_x_mm, off_y_mm, channels=None):
        """
        One event.

        `channels` is {(view, index): charge_in_electrons}; None or empty is
        valid and is exactly the case C10 exists to preserve — an event that
        produced no charge is still a trigger the DAQ recorded.
        """
        self.ev.append((int(event_id), int(n_electron_in), int(n_seed),
                        float(q_seed_total), float(x_true_mm),
                        float(y_true_mm), float(t0_ns),
                        float(off_x_mm), float(off_y_mm)))
        for (view, idx), q in sorted((channels or {}).items()):
            self.ch_view.append(0 if view == "X" else 1)
            self.ch_idx.append(int(idx))
            self.ch_q.append(float(q))
        self.ch_off.append(len(self.ch_q))

    def write(self, path, meta=None):
        import json
        a = np.array(self.ev, dtype=np.float64) if self.ev else np.zeros((0, 9))
        np.savez_compressed(
            path,
            schema=SCHEMA,
            event_id=a[:, 0].astype(np.uint64),
            n_electron_in=a[:, 1].astype(np.int64),
            n_seed=a[:, 2].astype(np.int64),
            q_seed_total=a[:, 3],
            x_true_mm=a[:, 4], y_true_mm=a[:, 5], t0_ns=a[:, 6],
            offset_x_mm=a[:, 7], offset_y_mm=a[:, 8],
            # CSR-style jagged per-channel truth: channels of event i are
            # chan_*[chan_offset[i]:chan_offset[i+1]].
            chan_offset=np.asarray(self.ch_off, dtype=np.int64),
            chan_view=np.asarray(self.ch_view, dtype=np.uint8),   # 0=X, 1=Y
            chan_index=np.asarray(self.ch_idx, dtype=np.uint16),
            chan_q_electrons=np.asarray(self.ch_q, dtype=np.float64),
            meta=json.dumps(meta or {}),
        )
        return path


def load(path):
    """Read a truth sidecar back as a dict of arrays (jagged stays CSR)."""
    with np.load(path, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def channels_of(t, i):
    """Per-channel truth of row i as [(view, index, charge), ...]."""
    lo, hi = int(t["chan_offset"][i]), int(t["chan_offset"][i + 1])
    return [("X" if v == 0 else "Y", int(c), float(q))
            for v, c, q in zip(t["chan_view"][lo:hi], t["chan_index"][lo:hi],
                               t["chan_q_electrons"][lo:hi])]
