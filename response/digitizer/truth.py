#!/usr/bin/env python3
"""
truth.py — the per-event truth sidecar (plan §7 step 7, audit C9).

The decoded files are deliberately indistinguishable from data: ADC counts,
pedestal, noise, and nothing else. That is the whole point of T13, and it also
means a simulated run cannot answer "how well did the reconstruction do", only
"did it run". This writes the other half — what actually happened — as a
sidecar keyed by eventId, so any reconstruction output can be joined against it.

FORMAT: PARQUET, with npz as an explicit fallback. The fix order asked for
parquet and the rest of the chain reads parquet (wft writes events.parquet), so
that is the default. It needs pyarrow, which the nTof_x17 venv has —
`~/PycharmProjects/nTof_x17/.venv/bin/python` is the interpreter to run this
chain with. Where pyarrow is genuinely absent (a bare condor worker, say) the
writer falls back to npz and SAYS SO on stdout rather than switching format
silently, because a sidecar that is sometimes one format and sometimes another
without telling you is worse than either.

The per-channel arrays are jagged. In parquet they are list columns; in npz
they are stored CSR-style — one flat array plus an offsets array — which is
exact and needs no object dtype. `load()` hides the difference.

WHAT IS AND IS NOT TRUTH HERE. The charges recorded are the induced charge per
channel BEFORE the DREAM shaper. That distinction used to be cosmetic and is
not any more: with the P6 residual-pole model int(h) = (1 - beta) int(h_nom),
so shaping no longer preserves the integral and "true charge" measured after it
would silently carry the shaper's beta.

THERE ARE TWO TRUTH IMPACT POINTS AND THEY ARE NOT THE SAME NUMBER. Stage A's
EventTree records the primary VERTEX (where the particle entered, audit C14);
`x_true_mm` / `y_true_mm` here are the CHARGE-WEIGHTED MEAN ARRIVAL POINT at
the ESL. They differ by track inclination over the 30 mm gap and, much more
violently, by delta rays dragging the charge centroid away from the track.
Measured on 40 normal-incidence muons: mean difference 0.16-0.19 mm, sd
0.7-1.3 mm, worst case 8.2 mm. Neither is wrong — a position reconstruction
should be scored against the arrival point, a tracking study against the
vertex — but comparing one to the other and calling the difference an error
would be a mistake.

    from .truth import TruthWriter
    tw = TruthWriter()
    tw.add(event_id, ...)
    tw.write("sim_truth.parquet")
"""

from __future__ import annotations

import numpy as np

SCHEMA = "mx17_truth/1"


class TruthWriter:
    """Accumulates per-event truth, then writes one parquet (or npz) file."""

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
        """
        Write the sidecar. `path` may end in .parquet or .npz; anything else
        gets the extension of whichever backend is actually used.
        """
        import json
        base = path[:-8] if path.endswith(".parquet") else (
            path[:-4] if path.endswith(".npz") else path)
        try:
            import pyarrow            # noqa: F401
        except ImportError:
            print(f"  [truth] pyarrow not available — writing npz instead of "
                  f"parquet. Use ~/PycharmProjects/nTof_x17/.venv/bin/python "
                  f"for the parquet sidecar the rest of the chain reads.")
            return self._write_npz(base + ".npz", meta)
        return self._write_parquet(base + ".parquet", meta)

    def _write_parquet(self, path, meta):
        import json
        import pandas as pd
        cols = ("event_id", "n_electron_in", "n_seed", "q_seed_total",
                "x_true_mm", "y_true_mm", "t0_ns", "offset_x_mm", "offset_y_mm")
        df = pd.DataFrame(self.ev, columns=cols) if self.ev else \
            pd.DataFrame({c: [] for c in cols})
        df["event_id"] = df["event_id"].astype("int64")
        for c in ("n_electron_in", "n_seed"):
            df[c] = df[c].astype("int64")
        # Per-channel truth as list columns, one row per event.
        off = self.ch_off
        df["chan_view"] = [["X" if v == 0 else "Y"
                            for v in self.ch_view[off[i]:off[i + 1]]]
                           for i in range(len(self.ev))]
        df["chan_index"] = [list(self.ch_idx[off[i]:off[i + 1]])
                            for i in range(len(self.ev))]
        df["chan_q_electrons"] = [list(self.ch_q[off[i]:off[i + 1]])
                                  for i in range(len(self.ev))]
        df.to_parquet(path, index=False)
        with open(path.replace(".parquet", ".meta.json"), "w") as f:
            json.dump({"schema": SCHEMA, **(meta or {})}, f, indent=1)
        return path

    def _write_npz(self, path, meta=None):
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
    """
    Read a truth sidecar written in either format.

    Returns a dict of arrays. Jagged per-channel data stays in whatever shape
    its backend used; use `channels_of` rather than touching it directly.
    """
    if path.endswith(".parquet"):
        import pandas as pd
        return pd.read_parquet(path)
    with np.load(path, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def channels_of(t, i):
    """Per-channel truth of row i as [(view, index, charge), ...]."""
    if hasattr(t, "iloc"):                       # a pandas DataFrame
        r = t.iloc[i]
        return [(str(v), int(c), float(q)) for v, c, q in
                zip(r["chan_view"], r["chan_index"], r["chan_q_electrons"])]
    lo, hi = int(t["chan_offset"][i]), int(t["chan_offset"][i + 1])
    return [("X" if v == 0 else "Y", int(c), float(q))
            for v, c, q in zip(t["chan_view"][lo:hi], t["chan_index"][lo:hi],
                               t["chan_q_electrons"][lo:hi])]
