#!/usr/bin/env python3
"""
shaper.py — the DREAM analogue chain, from the manual (plan §8 step 1).

Everything here is read out of `~/x17/Documents/dream/DREAM_User Manual_prod_v3.pdf`
and the run config `CosmicTb_MX17.cfg`, not assumed.

WHAT THE MANUAL SAYS
--------------------
"The analog filter consists in a Pole Zero Cancellation stage followed by a
2 complex pole Sallen-Key low pass filter. The peaking time of the global filter
is selectable among sixteen values in the range of 50 ns to 900 ns."

Peaking times are Table 9, quoted as **5 % -> 100 %**, not 0 -> peak. That
distinction is not pedantic: it is a ~10 % difference in the shaping constant,
and taking 180 ns as 0->peak would make every simulated pulse too slow.

Gain is the CSA feedback capacitor: input ranges 50 / 100 / 200 / 600 fC into a
fixed 2 V p-p output.

WHAT THE RUN CONFIG SAYS
------------------------
    Feu * Dream *  1 0x081F 0xD023 0x0000 0x0000
    Feu * Dream *  6 0xAAAA 0xAAAA 0xAAAA 0xAAAA
    Feu * Dream *  7 0xAAAA 0xAAAA 0xAAAA 0xAAAA

Peaking time lives in state1<4:7>, so (0xD023 >> 4) & 0xF = 2 -> **180 ns**,
which confirms the number the nTof_x17 notes carry.

The GAIN does not live in register 1 at all. The manual is explicit that it is
"done channel by channel with two 64-bit slow control registers state6<63:0> &
state7<63:0> (2 bits per channel)". Those are the 0xAAAA words = 1010..., i.e.
code 2 for every channel -> the **200 fC** range, hence 2 V / 200 fC =
**10 mV/fC**. Anyone looking for the gain in register 1 will not find it.

THE MODEL, AND ITS ONE APPROXIMATION
------------------------------------
A PZC stage followed by a 2-pole low pass is a semi-Gaussian shaper. This
implements it as CR-RC^2, whose impulse response is (t/tau)^2 exp(-t/tau) —
the standard form, with the shaping constant calibrated so the 5 %->100 % rise
matches the manual's table rather than by assuming tau = t_peak / 2.

The approximation: a Sallen-Key section has COMPLEX poles, while CR-RC^2 has a
real double pole. The two agree closely near the peak and differ mildly on the
tail, where the complex-pole response undershoots slightly. Since the measured
undershoot is a §9 validation target (-4 to -6 %), do NOT read this model's
undershoot as a prediction — it is a property of the CR-RC^2 stand-in. Fixing
it means implementing the actual Sallen-Key biquad, which is the natural next
step if the undershoot comparison matters.
"""

from __future__ import annotations

import numpy as np

# Table 9 of the manual, 5 % -> 100 % peaking time [ns], indexed by the 4-bit code.
PEAKING_TIME_NS = (76, 123, 180, 228, 283, 328, 388, 433,
                   578, 618, 675, 717, 781, 818, 880, 919)

# Table 1: selectable CSA input charge ranges [fC], indexed by the 2-bit
# per-channel gain code in state6/state7.
INPUT_RANGE_FC = (50.0, 100.0, 200.0, 600.0)
OUTPUT_SWING_V = 2.0

# Bench run config, CosmicTb_MX17.cfg.
BENCH_REG1 = 0xD023
BENCH_REG6 = 0xAAAA


def peaking_code(reg1=BENCH_REG1):
    """state1<4:7> — the 4-bit peaking-time selector."""
    return (reg1 >> 4) & 0xF


def gain_code(reg6=BENCH_REG6, channel=0):
    """2 bits per channel out of state6/state7. 0xAAAA -> code 2 everywhere."""
    return (reg6 >> (2 * (channel % 8))) & 0x3


def mV_per_fC(code=None):
    code = gain_code() if code is None else code
    return OUTPUT_SWING_V * 1e3 / INPUT_RANGE_FC[code]


def _tau_for_peaking(t_peak_ns, n=2):
    """
    Shaping constant tau such that CR-RC^n rises 5 % -> 100 % in t_peak_ns.

    h(u) = u^n e^-u peaks at u = n. The 5 % point is the SMALL root of
    u^n e^-u = 0.05 n^n e^-n, and the manual's quoted time is the gap between
    them. Solving it rather than assuming tau = t_peak/n matters: for n=2 the
    two differ by ~9 %.
    """
    peak = (n ** n) * np.exp(-n)
    target = 0.05 * peak
    # Plain bisection on [0, n], where u^n e^-u rises monotonically. Fifty
    # halvings take the bracket below 1e-15; not worth a scipy dependency for
    # a scalar root, and response/ is meant to stay numpy-only for the core.
    lo, hi = 0.0, float(n)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (mid ** n) * np.exp(-mid) < target:
            lo = mid
        else:
            hi = mid
    return t_peak_ns / (n - 0.5 * (lo + hi))


class DreamShaper:
    """CR-RC^2 semi-Gaussian shaper matched to the DREAM manual."""

    def __init__(self, reg1=BENCH_REG1, reg6=BENCH_REG6, n=2, dt_ns=1.0,
                 t_max_ns=None):
        self.code = peaking_code(reg1)
        self.t_peak_ns = float(PEAKING_TIME_NS[self.code])
        self.n = n
        self.tau_ns = _tau_for_peaking(self.t_peak_ns, n)
        self.gain_code = gain_code(reg6)
        self.mV_per_fC = mV_per_fC(self.gain_code)
        self.dt_ns = dt_ns
        # Long enough that the tail is dead: 12 tau is <1e-4 of the peak.
        t_max_ns = t_max_ns or 12 * self.tau_ns
        t = np.arange(0.0, t_max_ns, dt_ns)
        h = (t / self.tau_ns) ** n * np.exp(-t / self.tau_ns)
        # Normalise to UNIT PEAK, so the shaped amplitude of a delta of charge
        # Q is exactly Q. Normalising to unit area instead would make the
        # amplitude depend on the peaking time, which is not what a
        # charge-sensitive front end does.
        self.h = h / h.max()
        self.t = t

    def apply(self, current, dt_ns=None):
        """
        Shape an induced-current waveform. Returns the shaped output in units
        of charge (multiply by mV_per_fC for volts).

        Convolve with the CURRENT, not with its integral. h(t) = (t/tau)^n
        e^(-t/tau) is by definition the response to a DELTA OF CHARGE, so the
        CSA integration is already inside it — the chain is CSA (integrate) ->
        CR (differentiate) -> RC^n, and the first two cancel. Integrating the
        current first and then applying h therefore counts the integration
        twice: the self-check below caught exactly that, reporting a peak at
        1186 ns and a 1105 ns rise for a shaper whose manual value is 180 ns,
        because a step convolved with h just climbs to a plateau.
        """
        dt = dt_ns or self.dt_ns
        i = np.asarray(current, dtype=float) * (dt * 1e-9)
        return np.convolve(i, self.h, mode="full")[:len(i)]

    def describe(self):
        return {
            "peaking_code": self.code,
            "t_peak_ns_5to100": self.t_peak_ns,
            "tau_ns": self.tau_ns,
            "order_n": self.n,
            "gain_code": self.gain_code,
            "input_range_fC": INPUT_RANGE_FC[self.gain_code],
            "mV_per_fC": self.mV_per_fC,
            "filter_model": "CR-RC^2 (real double pole) standing in for "
                            "PZC + 2-complex-pole Sallen-Key",
            "caveat": "undershoot is a property of the CR-RC^2 stand-in, "
                      "NOT a prediction — see module docstring",
        }


def sample(waveform, dt_ns, n_samp=32, period_ns=60.0, rng=None, phase_ns=None):
    """
    SCA sampling: `n_samp` samples every `period_ns`, at a random trigger phase.

    The phase is uniform over one sampling period because the trigger is
    asynchronous with the sampling clock — that jitter is a real contribution
    to the measured peak-amplitude spread and must not be left out.
    """
    rng = rng or np.random.default_rng()
    if phase_ns is None:
        phase_ns = rng.uniform(0.0, period_ns)
    idx = np.rint((phase_ns + np.arange(n_samp) * period_ns) / dt_ns).astype(int)
    idx = np.clip(idx, 0, len(waveform) - 1)
    return waveform[idx], phase_ns


if __name__ == "__main__":
    import pprint
    sh = DreamShaper()
    pprint.pprint(sh.describe())
    print(f"\n  tau = {sh.tau_ns:.1f} ns   (naive t_peak/n would give "
          f"{sh.t_peak_ns/sh.n:.1f} ns, "
          f"{100*(sh.tau_ns/(sh.t_peak_ns/sh.n)-1):+.0f} %)")
    imp = np.zeros(2000); imp[0] = 1.0 / 1e-9    # a delta of unit charge
    w = sh.apply(imp)
    pk = int(np.argmax(w))
    above = np.flatnonzero(w >= 0.05 * w.max())
    print(f"  check: peak at {pk*sh.dt_ns:.0f} ns, "
          f"5%->100% = {(pk-above[0])*sh.dt_ns:.0f} ns "
          f"(manual says {sh.t_peak_ns:.0f} ns)")
