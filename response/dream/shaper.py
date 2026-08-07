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

On the SK filter (§2.1.4): "The filter damping factor is 0.75, so that the
global filter response exhibits a 1% only undershoot." So the biquad's own
undershoot is ~1 % BY DESIGN — the measured −4 to −6 % cannot come from the
Sallen-Key poles, and implementing "the actual biquad" (the original P6 plan)
would not have produced it.

On the CSA (§2.1.2): the DC feedback emulates a high-value resistor Rf, and
"the Rf.Cf time constant can be set to 5 µs ... or 50 µs independently of the
Gain" via state1<31>. On the PZC (§2.1.3): it exists "to avoid long duration
undershoot at the shaped output due to the CSA main pole" — i.e. the CSA pole
IS the undershoot mechanism, and the PZC's job is to cancel it. Any residual
(mismatched) fraction of that pole appears as a high-pass with recovery time
Rf·Cf, giving a shallow, µs-persistent negative tail. That matches both
run_71 features: the −4 to −6 % level still present at the 3.84 µs window end,
and the −30 % end-of-drift-ladder lobe at 700 V (a compact unipolar h could
never go negative after the current stops; a quasi-differentiating one must).

Peaking times are Table 9, quoted as **5 % -> 100 %**, not 0 -> peak. That
distinction is not pedantic: it is a ~10 % difference in the shaping constant,
and taking 180 ns as 0->peak would make every simulated pulse too slow.

Gain is the CSA feedback capacitor: input ranges 50 / 100 / 200 / 600 fC into a
fixed 2 V p-p output.

WHAT THE RUN CONFIG SAYS
------------------------
    Feu * Dream *  1 0x081F 0xD023 0x0000 0x0000
    Feu 1 Dream *  1 0x881F 0xD043 0x0000 0x0000
    Feu * Dream *  6 0xAAAA 0xAAAA 0xAAAA 0xAAAA
    Feu * Dream *  7 0xAAAA 0xAAAA 0xAAAA 0xAAAA

The FEU 1 (trigger-plane) override fixes the word order: its first word flips
the top bit and its second moves the peaking code, and only the reading
"first word = state1<31:16>, second = state1<15:0>" makes both sensible
(trigger FEU: 50 µs integration + 283 ns peaking; data FEUs: 5 µs + 180 ns).
So for the detector FEUs: state1 = 0x081FD023 →

  state1<7:4>  = 2  -> 180 ns peaking (Table 9), as the nTof_x17 notes carry.
  state1<31>   = 0  -> Rf·Cf = 5 µs   (manual: "1 -> 50 µs, 0 -> 5 µs").
  state1<10>   = 0  -> normal PZC mode (not the 450 ns integrator mode).

The GAIN does not live in register 1 at all. The manual is explicit that it is
"done channel by channel with two 64-bit slow control registers state6<63:0> &
state7<63:0> (2 bits per channel)". Those are the 0xAAAA words = 1010..., i.e.
code 2 for every channel -> the **200 fC** range, hence 2 V / 200 fC =
**10 mV/fC**. Anyone looking for the gain in register 1 will not find it.

THE MODEL
---------
Nominal filter, per the manual's topology: a real pole (the PZC replacement
pole Rs·Cs, same shaping time as the SK per §2.1.3) followed by the 2-complex-
pole Sallen-Key at damping ζ = 0.75:

    H_nom(s) = ω0³ / [ (s + ω0)(s² + 2·0.75·ω0·s + ω0²) ]

with ω0 calibrated so the 5 %->100 % rise matches Table 9. Validation that
this is the right reconstruction: its impulse response undershoots by −1.1 %,
against the manual's stated "1 % only".

CSA residual: the un-cancelled fraction β of the 5 µs pole, applied as

    H(s) = H_nom(s) · [ 1 − β / (1 + s·τ_f) ] ,   τ_f = Rf·Cf = 5 µs.

β = 0 is a perfect PZC (the design intent); β = 1 is full AC coupling (zero
at DC). β is a hardware property of the PZC trim that the manual does not
quote — treat it as a SCAN/NUISANCE parameter like ρ_s, not a constant.
Indicative values from run_71 (2026-08-07 prototype, scratch): the −4 to −6 %
undershoot depth corresponds to β ≈ 0.6–0.9 and the −30 % end-of-ladder lobe
at 700 V independently to β ≈ 0.85–0.9; the default 0.75 sits in that band
but MUST be treated as constrained tuning at T14, where the undershoot target
selects β and the remaining shape observables become the cross-checks.

Side effects to be aware of: the residual high-pass shaves 0.6–2.3 % off the
peak (β 0.25→1) and makes the long-window area integral β-dependent — the
LTI area-invariance argument in the plan holds per filter, but ∫h is now
(1−β)·∫h_nom, so area budgets on windows ≳ τ_f are no longer β-blind.

The legacy CR-RC² stand-in ((t/τ)² e^(−t/τ), real double pole) is kept as
`filter_model="crrc2"` for A/B tests. Note it has EXACTLY ZERO undershoot —
its impulse response is nonnegative — so any undershoot measured on it is
numerics, full stop.
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

# SK filter damping factor, manual §2.1.4.
SK_ZETA = 0.75

# Bench run config, CosmicTb_MX17.cfg. REG1 is state1<15:0>, REG1_HI is
# state1<31:16> (word order pinned by the FEU 1 trigger-plane override).
BENCH_REG1 = 0xD023
BENCH_REG1_HI = 0x081F
BENCH_REG6 = 0xAAAA


def peaking_code(reg1=BENCH_REG1):
    """state1<4:7> — the 4-bit peaking-time selector."""
    return (reg1 >> 4) & 0xF


def tau_csa_us(reg1_hi=BENCH_REG1_HI):
    """state1<31>: CSA integration time constant Rf·Cf. 0 -> 5 µs, 1 -> 50 µs."""
    return 50.0 if (reg1_hi >> 15) & 1 else 5.0


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


def _h_sk(t_ns, w0, zeta=SK_ZETA):
    """
    Impulse response of ω0³/[(s+ω0)(s²+2ζω0s+ω0²)] by residues, on t_ns [ns]
    with w0 in rad/ns. Exact closed form — no time stepping.
    """
    wd = w0 * np.sqrt(1.0 - zeta * zeta)
    a_real = w0 / (2.0 - 2.0 * zeta)                    # residue at s = -w0
    sp = complex(-zeta * w0, wd)
    r_cplx = w0 ** 3 / ((sp + w0) * (2j * wd))          # residue at s_+
    h = a_real * np.exp(-w0 * t_ns) + 2.0 * np.real(r_cplx * np.exp(sp * t_ns))
    return np.where(t_ns < 0, 0.0, h)


def _rise_5_100_ns(t_ns, h):
    pk = int(np.argmax(h))
    above = np.flatnonzero(h >= 0.05 * h[pk])
    return t_ns[pk] - t_ns[above[0]]


def _w0_for_peaking(t_peak_ns, zeta=SK_ZETA):
    """ω0 [rad/ns] such that the SK model rises 5 % -> 100 % in t_peak_ns."""
    t = np.arange(0.0, 30.0 * t_peak_ns, t_peak_ns / 400.0)
    lo, hi = 1e-2 / t_peak_ns, 1e2 / t_peak_ns
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if _rise_5_100_ns(t, _h_sk(t, mid, zeta)) > t_peak_ns:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


class DreamShaper:
    """DREAM shaper matched to the manual: PZC pole + ζ=0.75 Sallen-Key,
    times a residual CSA-pole high-pass (the undershoot mechanism)."""

    def __init__(self, reg1=BENCH_REG1, reg6=BENCH_REG6, reg1_hi=BENCH_REG1_HI,
                 filter_model="sk", pzc_residual=0.75, n=2, dt_ns=1.0,
                 t_max_ns=None):
        self.code = peaking_code(reg1)
        self.t_peak_ns = float(PEAKING_TIME_NS[self.code])
        self.filter_model = filter_model
        self.n = n
        self.gain_code = gain_code(reg6)
        self.mV_per_fC = mV_per_fC(self.gain_code)
        self.dt_ns = dt_ns
        self.tau_csa_ns = tau_csa_us(reg1_hi) * 1e3
        self.beta = float(pzc_residual)
        if filter_model == "sk":
            self.w0 = _w0_for_peaking(self.t_peak_ns)
            self.tau_ns = 1.0 / self.w0                 # bookkeeping only
        elif filter_model == "crrc2":
            self.tau_ns = _tau_for_peaking(self.t_peak_ns, n)
        else:
            raise ValueError(f"unknown filter_model {filter_model!r}")
        # The fast lobe is dead after ~16 tau; the undershoot tail recovers
        # with tau_csa, so when beta > 0 the kernel must extend a few tau_csa
        # or the truncated negative area silently un-conserves the response.
        if t_max_ns is None:
            t_max_ns = 16 * self.tau_ns
            if self.beta > 0:
                t_max_ns = max(t_max_ns, 4 * self.tau_csa_ns)
        self.t = t = np.arange(0.0, t_max_ns, dt_ns)
        if filter_model == "sk":
            h = _h_sk(t, self.w0)
        else:
            h = (t / self.tau_ns) ** n * np.exp(-t / self.tau_ns)
        # Normalise the NOMINAL response to UNIT PEAK, so the shaped amplitude
        # of a delta of charge Q is exactly Q, matching the derived ADC scale.
        # The residual high-pass is applied AFTER normalisation: its ~1-2 %
        # peak loss is hardware behaviour, not a normalisation choice.
        h = h / h.max()
        if self.beta > 0:
            h = h - self.beta * self._csa_pole_conv(h)
        self.h = h

    def _csa_pole_conv(self, h):
        """(k ⊛ h)(t) with k = (1/τ_f)e^(−t/τ_f), exactly, via the scaled
        cumulative sum — O(N), no quadrature error beyond the 1 ns grid."""
        r = self.dt_ns / self.tau_csa_ns
        g = np.exp(-self.t / self.tau_csa_ns)
        # exp(+t/τ_f) stays ≤ e^(t_max/τ_f); with t_max = 4 τ_f that is e^4.
        return r * g * np.cumsum(h / g)

    def apply(self, current, dt_ns=None):
        """
        Shape an induced-current waveform. Returns the shaped output in units
        of charge (multiply by mV_per_fC for volts).

        Convolve with the CURRENT, not with its integral. h(t) is by
        definition the response to a DELTA OF CHARGE, so the CSA integration
        is already inside it — the chain is CSA (integrate) -> PZC
        (differentiate) -> low pass, and the first two cancel. Integrating the
        current first and then applying h therefore counts the integration
        twice: the self-check below caught exactly that, reporting a peak at
        1186 ns and a 1105 ns rise for a shaper whose manual value is 180 ns,
        because a step convolved with h just climbs to a plateau.
        """
        dt = dt_ns or self.dt_ns
        if abs(dt - self.dt_ns) > 1e-12:
            raise ValueError(f"shaper kernel is on a {self.dt_ns} ns grid but "
                             f"apply() was called with dt_ns={dt}; rebuild the "
                             f"shaper with the matching dt_ns")
        i = np.asarray(current, dtype=float) * (dt * 1e-9)
        # Kernel samples beyond the output window cannot contribute; truncate
        # before the FFT so the transform length tracks the waveform.
        h = self.h[:len(i)]
        nfft = 1
        while nfft < len(i) + len(h) - 1:
            nfft *= 2
        y = np.fft.irfft(np.fft.rfft(i, nfft) * np.fft.rfft(h, nfft), nfft)
        return y[:len(i)]

    def describe(self):
        d = {
            "peaking_code": self.code,
            "t_peak_ns_5to100": self.t_peak_ns,
            "gain_code": self.gain_code,
            "input_range_fC": INPUT_RANGE_FC[self.gain_code],
            "mV_per_fC": self.mV_per_fC,
            "filter_model": self.filter_model,
            "tau_csa_ns": self.tau_csa_ns,
            "pzc_residual_beta": self.beta,
            "undershoot_pred_pct": round(100 * self.h.min() / self.h.max(), 2),
        }
        if self.filter_model == "sk":
            d["sk_zeta"] = SK_ZETA
            d["w0_rad_per_ns"] = self.w0
            d["provenance"] = ("topology+tau_csa from DREAM manual and "
                               "CosmicTb_MX17.cfg state1=0x081FD023; beta is a "
                               "scan/nuisance parameter (PZC residual, not in "
                               "the datasheet)")
        else:
            d["tau_ns"] = self.tau_ns
            d["caveat"] = ("CR-RC^2 legacy stand-in: impulse response is "
                           "nonnegative, undershoot is exactly zero by "
                           "construction")
        return d


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
    for kwargs in ({"filter_model": "sk", "pzc_residual": 0.0},
                   {"filter_model": "sk", "pzc_residual": 0.75},
                   {"filter_model": "crrc2", "pzc_residual": 0.0}):
        sh = DreamShaper(**kwargs)
        pprint.pprint(sh.describe())
        imp = np.zeros(20000); imp[0] = 1.0 / 1e-9   # a delta of unit charge
        w = sh.apply(imp)
        pk = int(np.argmax(w))
        above = np.flatnonzero(w >= 0.05 * w.max())
        print(f"  peak {w.max():+.4f} at {pk*sh.dt_ns:.0f} ns, "
              f"5%->100% = {(pk-above[0])*sh.dt_ns:.0f} ns "
              f"(manual {sh.t_peak_ns:.0f} ns), "
              f"min {100*w.min()/w.max():+.2f} % at "
              f"{np.argmin(w)*sh.dt_ns:.0f} ns\n")
