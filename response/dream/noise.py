#!/usr/bin/env python3
"""
noise.py — DREAM front-end noise, characterised from det3 pedestals (plan T12).

WHY THIS IS AN INPUT AND NOT A TUNING KNOB. Electronics noise is a property of
the apparatus, measured from its own pedestals, in the same category as the gas
mixture or the mesh voltage. It is not one of §9's response observables
(sharing fractions, peak ratios, peak-time shifts), which stay blind until T14.

WHAT THE DATA ACTUALLY LOOKS LIKE, and why the obvious model is wrong. Measured
on det3 `mx17_det3_long_run_5-6-26`, 32 samples x 512 channels at 60 ns:

    pedestal                    341 ADC mean, 34 ADC spread across channels
    common mode per 64-block    83 ADC rms          <-- dominant
    per-channel after CNS       10.4 ADC median (3.4 spread)
    per-channel before CNS      19.3 ADC
    sample-lag autocorrelation  1.00 0.99 0.92 0.81 0.66 0.51

Two things follow, and both bite if ignored:

  1. **The noise is dominated by COHERENT common mode**, 83 ADC against 10 ADC
     incoherent. wft removes it with a per-64-channel-block median, so a
     simulation that emits only per-channel white noise gives wft nothing to
     subtract and lands in a regime the reconstruction never sees. The ZS
     decision is also made on the raw, pre-CNS waveform in the DAQ, so the
     common mode is what sets the occupancy.

  2. **The noise is strongly correlated sample to sample** — the
     autocorrelation is still 0.5 at lag 5 (300 ns), which is the 180 ns
     shaping doing exactly what it is supposed to. White noise of the same rms
     would produce a completely different threshold-crossing rate, so it is
     generated here from the measured POWER SPECTRUM rather than from an rms.

Generation is spectral: draw white noise in frequency, scale by the measured
amplitude spectrum, transform back. That reproduces the measured
autocorrelation by construction instead of approximating it with an AR fit, and
`selftest()` closes the loop by characterising the generated noise with the
same estimator used on the data.

    python3 -m response.dream.noise --characterise <decoded.root> --out noise.json
    python3 -m response.dream.noise --selftest noise.json
"""

from __future__ import annotations

import argparse
import glob
import json

import numpy as np

CNS_BLOCK = 64          # must match wft/io.py
N_CHAN = 512
ADC_MAX = 4095


def _stack(path, n_events):
    """Dense (n_ev, n_sample, 512) ADC block from a decoded FEU file."""
    import uproot
    t = uproot.open(path)["nt"]
    a = t.arrays(["amplitude"], entry_stop=n_events, library="np")["amplitude"]
    lens = np.array([len(x) // N_CHAN for x in a])
    ns = int(np.bincount(lens).argmax())
    if ns == 0:
        raise SystemExit(
            f"{path}: modal amplitude length is 0 — this file is "
            "zero-suppressed and carries no pedestal traces. Use a dense run.")
    return np.stack([x.reshape(ns, N_CHAN)
                     for x, l in zip(a, lens) if l == ns]).astype(np.float32), ns


def _cns(sub):
    """Subtract the per-block, per-sample median. Same operation as wft."""
    nblk = N_CHAN // CNS_BLOCK
    cms = np.median(sub.reshape(sub.shape[0], sub.shape[1], nblk, CNS_BLOCK),
                    axis=3)
    return sub - np.repeat(cms, CNS_BLOCK, axis=2), cms


def characterise(path, n_events=1500, quiet_k=3.0):
    """
    Measure the noise structure of a dense decoded run.

    QUIET SELECTION. These are cosmic data runs, not dedicated pedestal runs,
    so real tracks sit in them. Any spectrum taken over everything would carry
    signal, and the model would then generate noise shaped like muons. Traces
    whose residual excursion exceeds quiet_k times their own robust sigma are
    dropped, and the fraction kept is reported so the cut is visible rather
    than implicit.
    """
    S, ns = _stack(path, n_events)
    ped = np.median(S, axis=(0, 1))
    sub = S - ped
    res, cms = _cns(sub)

    sig = 1.4826 * np.median(np.abs(res), axis=(0, 1))       # per channel
    quiet = np.abs(res).max(axis=1) < quiet_k * sig[None, :]  # (n_ev, n_chan)
    frac = float(quiet.mean())

    # --- incoherent spectrum, quiet traces only -----------------------------
    # The DC bin is KEPT. Per-trace mean subtraction before the FFT would zero
    # f=0, but the pedestal has already been removed above, so what sits at DC
    # is genuine event-to-event baseline wander — real noise that the DAQ sees
    # and that wft's CNS is there to remove. Dropping it makes the generated
    # noise quieter than the data by exactly that term, which is most of the
    # common mode.
    q = res.transpose(0, 2, 1)[quiet]              # (n_quiet, n_sample)
    amp_chan = np.sqrt((np.abs(np.fft.rfft(q, axis=1)) ** 2).mean(axis=0))

    # --- common-mode spectrum ----------------------------------------------
    c = cms.transpose(0, 2, 1).reshape(-1, ns)
    amp_cm = np.sqrt((np.abs(np.fft.rfft(c, axis=1)) ** 2).mean(axis=0))

    # Autocorrelation is a SHAPE, so it is reported per-trace mean-subtracted;
    # selftest applies the identical estimator to the generated noise.
    lag = _autocorr(q - q.mean(axis=1, keepdims=True))
    return {
        "schema": "dream_noise/1",
        "source": path.split("/")[-1],
        "n_events_used": int(len(S)), "n_sample": int(ns),
        "quiet_fraction": frac,
        "pedestal_mean_adc": float(ped.mean()),
        "pedestal_spread_adc": float(ped.std()),
        "pedestal_per_channel": ped.astype(float).tolist(),
        "sigma_chan_median_adc": float(np.median(sig)),
        "sigma_chan_spread_adc": float(sig.std()),
        "sigma_chan": sig.astype(float).tolist(),
        "sigma_cm_adc": float(cms.std()),
        "amp_spectrum_chan": amp_chan.astype(float).tolist(),
        "amp_spectrum_cm": amp_cm.astype(float).tolist(),
        "autocorr_lag": lag,
    }


def _autocorr(q, nlag=6):
    q = q - q.mean()
    denom = float((q * q).mean())
    n = q.shape[1]
    return [float((q[:, :n - k] * q[:, k:]).mean() / denom) for k in range(nlag)]


def _parseval_rms(amp, ns):
    """
    Rms of a real trace whose one-sided amplitude spectrum is `amp`.

    The same identity `_coloured` is built on, read the other way round:

        var = (1/ns^2) [ A_0^2 + A_nyq^2 + 2 sum_mid A_k^2 ]

    Used to normalise, so it must follow numpy's rfft binning exactly — the
    Nyquist bin exists only for even ns and, like DC, is NOT doubled.
    """
    a = np.asarray(amp, dtype=np.float64)
    tot = a[0] ** 2
    if ns % 2 == 0:
        tot += a[-1] ** 2
        mid = a[1:-1]
    else:
        mid = a[1:]
    return float(np.sqrt(tot + 2.0 * np.sum(mid ** 2)) / ns)


class NoiseModel:
    """Generates ADC noise with the measured spectrum, per 64-channel block."""

    def __init__(self, spec):
        self.spec = spec
        self.ns = int(spec["n_sample"])
        self.ped = np.asarray(spec["pedestal_per_channel"], dtype=np.float32)
        self.sigma_chan = np.asarray(spec["sigma_chan"], dtype=np.float32)
        self.amp_chan = np.asarray(spec["amp_spectrum_chan"], dtype=np.float64)
        self.amp_cm = np.asarray(spec["amp_spectrum_cm"], dtype=np.float64)

    def _coloured(self, rng, shape, amp, ns):
        """
        Noise with a prescribed amplitude spectrum.

        For real white noise x ~ N(0,1) of length ns, numpy's rfft gives
        E|X_k|^2 = ns in every bin. Setting Y_k = X_k * A_k / sqrt(ns)
        therefore gives E|Y_k|^2 = A_k^2, i.e. exactly the measured spectrum,
        and Parseval then fixes the variance with no free constant:

            var(y) = (1/ns^2) [ A_0^2 + A_nyq^2 + 2 sum_mid A_k^2 ]

        which is the same identity that relates the data's own spectrum to its
        variance. numpy's irfft already carries the 1/ns, so nothing further is
        applied. An earlier version multiplied by sqrt(ns) here and generated
        noise 5.7x too loud — the autocorrelation still matched, because a
        pure scale error leaves the normalised autocorrelation untouched. That
        is why selftest checks amplitude and shape separately.
        """
        w = np.fft.rfft(rng.standard_normal(shape + (ns,)), axis=-1)
        w *= (amp / np.sqrt(ns))[(None,) * len(shape)]
        return np.fft.irfft(w, ns, axis=-1)

    def sample(self, rng, n_ev, n_sample=None, n_chan=N_CHAN):
        """(n_ev, n_sample, n_chan) of pedestal + common mode + per-channel."""
        ns = n_sample or self.ns
        nblk = n_chan // CNS_BLOCK

        cm = self._coloured(rng, (n_ev, nblk), self.amp_cm, ns)   # (ev,blk,ns)
        cm = np.repeat(cm, CNS_BLOCK, axis=1).transpose(0, 2, 1)  # (ev,ns,chan)

        ch = self._coloured(rng, (n_ev, n_chan), self.amp_chan, ns)
        ch = ch.transpose(0, 2, 1)
        # Restore the per-channel sigma spread: the spectrum is an average over
        # channels, so scale each to its own measured sigma.
        #
        # THE DENOMINATOR IS THE TRACE'S OWN RMS, NOT median(sigma_chan)
        # (fix 2026-08-07, audit C1). `_coloured` produces a trace whose rms is
        # the Parseval rms of `amp_chan` — and `amp_chan` is the MEAN spectrum
        # over channels, so that rms is sqrt(mean variance). Dividing by the
        # MEDIAN sigma instead is only correct if the two coincide, and the
        # sigma distribution is right-skewed, so mean-square > median^2 and
        # every channel came out ~3.7 % loud. That is exactly the ~3.5 %
        # residual `selftest` has always reported and always tolerated.
        scale = self.sigma_chan[:n_chan] / max(
            _parseval_rms(self.amp_chan, ns), 1e-9)
        ch = ch * scale[None, None, :]

        return cm + ch + self.ped[None, None, :n_chan]

    def quantise(self, x):
        """ADC is an unsigned 12-bit integer; clipping is real, so model it."""
        return np.clip(np.rint(x), 0, ADC_MAX).astype(np.uint16)


def selftest(spec, n_ev=1500, seed=3):
    """
    Round-trip: generate, then characterise with the SAME estimator used on the
    data. This is the internal check that the model reproduces its input — it
    does not look at any §9 observable.
    """
    m = NoiseModel(spec)
    rng = np.random.default_rng(seed)
    S = m.sample(rng, n_ev).astype(np.float32)
    ped = np.median(S, axis=(0, 1))
    res, cms = _cns(S - ped)
    sig = 1.4826 * np.median(np.abs(res), axis=(0, 1))
    q = res.transpose(0, 2, 1).reshape(-1, m.ns)
    lag = _autocorr(q - q.mean(axis=1, keepdims=True))

    # Tolerances tightened 2026-08-07 with the C1 scale fix: the per-channel
    # residual was a systematic +3.5 % (the median-vs-Parseval bug) and is now
    # 0.2 %, so a 10 % bar no longer tests anything. 3 % leaves room for the
    # MAD estimator's own sampling error at this n_ev.
    rows = [
        ("sigma per channel [ADC]", spec["sigma_chan_median_adc"],
         float(np.median(sig)), 0.03),
        ("sigma common mode [ADC]", spec["sigma_cm_adc"], float(cms.std()), 0.05),
    ]
    ok = True
    print(f"  {'quantity':28s} {'data':>9s} {'model':>9s}  {'rel':>7s}")
    for name, want, got, tol in rows:
        rel = abs(got - want) / max(abs(want), 1e-9)
        good = rel <= tol
        ok &= good
        print(f"  {name:28s} {want:9.2f} {got:9.2f}  {rel:6.1%}  "
              f"{'OK' if good else 'FAIL'}")

    print(f"\n  {'autocorr lag':28s} " +
          "".join(f"{k:>7d}" for k in range(len(lag))))
    print(f"  {'  data':28s} " +
          "".join(f"{v:7.3f}" for v in spec["autocorr_lag"]))
    print(f"  {'  model':28s} " + "".join(f"{v:7.3f}" for v in lag))
    dmax = max(abs(a - b) for a, b in zip(spec["autocorr_lag"], lag))
    good = dmax <= 0.05
    ok &= good
    print(f"  max |diff| {dmax:.3f}  {'OK' if good else 'FAIL'}")
    print("\n" + ("PASS" if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--characterise")
    ap.add_argument("--n-events", type=int, default=1500)
    ap.add_argument("--out")
    ap.add_argument("--selftest")
    a = ap.parse_args()

    if a.characterise:
        p = a.characterise
        if "*" in p:
            p = sorted(glob.glob(p))[0]
        spec = characterise(p, a.n_events)
        print(f"  {spec['source']}: {spec['n_events_used']} events x "
              f"{spec['n_sample']} samples, quiet fraction "
              f"{spec['quiet_fraction']:.3f}")
        print(f"  pedestal   {spec['pedestal_mean_adc']:.1f} ADC "
              f"(spread {spec['pedestal_spread_adc']:.1f})")
        print(f"  sigma      per-channel {spec['sigma_chan_median_adc']:.2f} "
              f"+- {spec['sigma_chan_spread_adc']:.2f}, "
              f"common mode {spec['sigma_cm_adc']:.2f} ADC")
        print(f"  autocorr   " +
              " ".join(f"{v:.3f}" for v in spec["autocorr_lag"]))
        if a.out:
            with open(a.out, "w") as f:
                json.dump(spec, f)
            print(f"  wrote {a.out}")
        return 0

    if a.selftest:
        with open(a.selftest) as f:
            spec = json.load(f)
        return 0 if selftest(spec) else 1

    ap.error("need --characterise or --selftest")


if __name__ == "__main__":
    raise SystemExit(main())
