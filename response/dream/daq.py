#!/usr/bin/env python3
"""
daq.py — FEU sampling, ADC conversion, ZS, and the decoded_root writer (T12/T13).

This is the last stage before the simulation becomes indistinguishable from
data on disk: it turns the 1 ns shaped channel currents into exactly the file
`wft/io.py` opens, so wft runs on simulation UNCHANGED (plan §2's
`sim_decoded.root` contract, T13).

THE OUTPUT CONTRACT, read off wft/io.py rather than assumed:

    tree `nt`, branches `eventId`, `ftst`, `amplitude`
    amplitude is FLAT, row-major [sample][channel], length n_sample * 512

wft reshapes with `.reshape(-1, 512)`, takes its pedestal as the per-channel
median over the first events, and subtracts a per-64-channel-block median as
common-mode. So the simulation has to supply pedestals and coherent noise for
wft to have anything to subtract — see response/dream/noise.py.

DENSE, NOT ZERO-SUPPRESSED, IS THE DEFAULT, and that is not laziness. wft's
reader requires a dense n_sample x 512 block; a zero-suppressed file has
variable-length events (the det3 ArIso run's modal amplitude length is 0) and
`.reshape(-1, 512)` cannot consume it. The det3 reference runs used for the
noise model are themselves dense (`Feu_RunCtrl_ZS 0`), so dense is also the
like-for-like comparison. `zero_suppress()` is provided for modelling ZS runs
and is off by default.

THE ADC SCALE IS DERIVED, NOT FITTED. From the manual (see shaper.py): the
selected gain code is 2, giving a 200 fC input range into a fixed 2 V p-p
output, and the FEU digitises that with a 12-bit ADC. So

    ADC per fC = 4096 counts / 200 fC = 20.48

with one stated assumption: that the ADC spans the DREAM's full 2 V output.
The det3 pedestal sits at 341 ADC, leaving ~3750 counts of headroom, i.e. an
effective range of ~183 fC before clipping — and clipping is modelled, because
the data reaches 4095.

TIMING PHASE IS RANDOM PER EVENT. The 60 ns sampling clock has no relationship
to the arrival time of a cosmic, so a fixed phase would produce a sampling
artefact that no real run has: every pulse would be sampled at the same point
on its rise. Each event draws a uniform phase over one sample period.
"""

from __future__ import annotations

import json

import numpy as np

N_CHAN = 512
CNS_BLOCK = 64
ADC_BITS, ADC_MAX = 12, 4095

# Manual: gain code 2 -> 200 fC full range -> fixed 2 V p-p output.
INPUT_RANGE_FC = 200.0
ADC_PER_FC = (1 << ADC_BITS) / INPUT_RANGE_FC      # 20.48

# UNITS, and the trap in them. digitize.induce() returns current in ELEMENTARY
# CHARGES per second, and shaper.apply() multiplies by dt and convolves with a
# PEAK-NORMALISED h — so the shaped waveform is in units of e, and its peak
# equals the input charge for instantaneous delivery. Feeding that straight
# into an fC scale over-scales by 1/1.6e-4 ~ 6200x, which saturated the 12-bit
# ADC on every event and was caught by the median simulated MIP pegging at
# 3888 of 3748 available counts. Convert explicitly.
FC_PER_E = 1.602176634e-4                          # 1 e = 1.602e-19 C = 1.6e-4 fC


def charges_to_fc(q_e):
    """Shaped waveform in elementary charges -> fC, the unit to_adc expects."""
    return np.asarray(q_e, dtype=np.float32) * FC_PER_E

# CosmicTb_MX17.cfg (det3): sample period 60 ns, 32 samples, ZsTyp 1 (tpc),
# ZsChkSmp 1, CmOffset 256, "Sys PedRun Threshold 5" sigma.
SAMPLE_PERIOD_NS = 60.0
N_SAMPLE = 32
CM_OFFSET = 256
ZS_NSIGMA = 5.0
ZS_CHK_SMP = 1


class Daq:
    """FEU sampling + ADC + optional ZS, with a noise model attached."""

    def __init__(self, noise_spec, n_sample=N_SAMPLE,
                 sample_period_ns=SAMPLE_PERIOD_NS, seed=1234):
        from .noise import NoiseModel
        self.noise = NoiseModel(noise_spec)
        self.n_sample = n_sample
        self.dt_ns = sample_period_ns
        self.rng = np.random.default_rng(seed)

    # ── sampling ─────────────────────────────────────────────────────────────

    def sample(self, wf_1ns, phase_ns=None, t0_ns=0.0):
        """
        Sample 1 ns waveforms onto the 60 ns DAQ grid.

        wf_1ns: (n_chan, n_t) in the shaper's output units.

        The FEU samples an analogue level, so this is a point sample of the
        shaped waveform, NOT an average over the 60 ns bin: the shaper has
        already band-limited the signal and averaging again would smear it a
        second time. With a 180 ns peaking time against 60 ns sampling there
        are ~3 samples on the rise, which is what the data shows.
        """
        if phase_ns is None:
            phase_ns = self.rng.uniform(0.0, self.dt_ns)
        idx = np.rint((t0_ns + phase_ns
                       + np.arange(self.n_sample) * self.dt_ns)).astype(int)
        n_t = wf_1ns.shape[1]
        out = np.zeros((wf_1ns.shape[0], self.n_sample), dtype=np.float32)
        ok = (idx >= 0) & (idx < n_t)
        out[:, ok] = wf_1ns[:, idx[ok]]
        return out, float(phase_ns)

    # ── ADC ──────────────────────────────────────────────────────────────────

    def to_adc(self, q_fc, n_ev=1):
        """
        (n_ev, n_sample, n_chan) signal in fC -> quantised ADC with pedestal
        and noise. Signal, pedestal and noise are summed in float and
        quantised ONCE, because quantising twice would add a second, fictitious
        LSB of dither.
        """
        sig = np.asarray(q_fc, dtype=np.float32) * ADC_PER_FC
        bkg = self.noise.sample(self.rng, n_ev, self.n_sample, N_CHAN)
        return self.noise.quantise(sig + bkg)

    # ── zero suppression ─────────────────────────────────────────────────────

    def zero_suppress(self, adc, ped, sigma, nsigma=ZS_NSIGMA,
                      chk_smp=ZS_CHK_SMP):
        """
        Firmware ZS, ZsTyp=1 ("tpc"), as documented in nTof_x17's
        26_zs_sim_extract.py: a channel-sample crosses when it exceeds its own
        pedestal by nsigma * sigma_ch, and the readout keeps the crossing
        samples plus `chk_smp` further samples per crossing run.

        Returns a boolean keep-mask of the same shape; writing a suppressed
        file is a separate concern from deciding what survives.
        """
        over = adc > (ped[None, None, :] + nsigma * sigma[None, None, :])
        keep = over.copy()
        for k in range(1, chk_smp + 1):
            keep[:, k:, :] |= over[:, :-k, :]
        return keep


# ── decoded_root writer ──────────────────────────────────────────────────────

def write_decoded(path, adc, event_ids=None, ftst=None):
    """
    Write `sim_decoded.root` in the exact schema wft/io.py reads.

    adc: (n_ev, n_sample, n_chan) uint16.

    The flattening order is [sample][channel] because wft does
    `.reshape(-1, 512)`, which reads the fastest-varying axis as the channel.
    Getting this transposed produces a file that loads without error and is
    silently scrambled, so it is asserted on read-back by `verify()`.
    """
    import awkward as ak
    import uproot

    n_ev, n_sample, n_chan = adc.shape
    n_hit = n_sample * n_chan
    flat = adc.reshape(n_ev * n_hit).astype(np.uint16)
    counts = np.full(n_ev, n_hit, dtype=np.int64)

    if event_ids is None:
        event_ids = np.arange(n_ev, dtype=np.uint64)
    if ftst is None:
        ftst = np.zeros(n_ev, dtype=np.uint16)

    # The index branches, matching the real file: channel is the FAST axis and
    # sample the slow one (verified against det3 -- channel[:512] == 0..511).
    smp1 = np.repeat(np.arange(n_sample, dtype=np.uint16), n_chan)
    chn1 = np.tile(np.arange(n_chan, dtype=np.uint16), n_sample)

    def jag(x1, n):
        """One event's vector, tiled to all events, as vector<uint16_t>."""
        return ak.unflatten(np.tile(x1, n).astype(np.uint16), counts)

    # TTree, NOT RNTuple, and the branch must be genuinely jagged. Two traps,
    # both of which produce a file that opens fine and that wft cannot read:
    #
    #   * `f["nt"] = {...}` with awkward input writes a ROOT::RNTuple under
    #     uproot 5.7. The data is a TTree, so mktree is used explicitly.
    #   * handing uproot a 2-D array writes a fixed-length RegularArray;
    #     uproot's `library='np'` then cannot convert it and wft dies in
    #     ak.to_numpy. The real branches are std::vector<uint16_t>, so the
    #     payloads are unflattened into var-length lists.
    var_u16 = ak.types.from_datashape("var * uint16", highlevel=False)
    with uproot.recreate(path) as f:
        f.mktree("nt", {"eventId": np.uint64, "timestamp": np.uint64,
                        "ftst": np.uint16, "sample": var_u16,
                        "channel": var_u16, "amplitude": var_u16})
        f["nt"].extend({
            "eventId": np.asarray(event_ids, dtype=np.uint64),
            "timestamp": np.zeros(n_ev, dtype=np.uint64),
            "ftst": np.asarray(ftst, dtype=np.uint16),
            "sample": jag(smp1, n_ev),
            "channel": jag(chn1, n_ev),
            "amplitude": ak.unflatten(flat, counts),
        })
    return path


def verify(path, adc):
    """
    Read the file back THE WAY WFT DOES and require it to reproduce `adc`.

    This is the T13 acceptance in one function: not "the file opened" but "the
    consumer's own access pattern returns the array we wrote". A transposed
    flatten passes every other check and fails this one.
    """
    import uproot
    t = uproot.open(path)["nt"]
    a = t.arrays(["amplitude"], library="np")["amplitude"]
    n_sample, n_chan = adc.shape[1], adc.shape[2]
    bad = 0
    for i, row in enumerate(a):
        got = np.asarray(row).reshape(-1, n_chan)
        if got.shape != (n_sample, n_chan) or not np.array_equal(got, adc[i]):
            bad += 1
    return bad


def load_noise(path):
    with open(path) as f:
        return json.load(f)
