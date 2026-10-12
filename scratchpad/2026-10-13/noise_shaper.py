#!/usr/bin/env python3
"""
noise_shaper.py -- what time scale does the DREAM shaper ACTUALLY have, from the noise?

The register says 180 ns (FRESH_LOOK §7.6), but chip spread / PZC trim are not in any
register, and at 180 ns the inclined-track data need the ion signal ~x5 suppressed.
The FEU pulser would measure the real response; this asks the data's own noise first.

Noise at the CSA input is a mix of parallel (white current) and series (white voltage
x detector capacitance -> current ~ f) noise. The output autocorrelation is therefore

    R(tau) = A (h * h)(tau) + B (h' * h')(tau) + C

with h the shaper's response to a delta of charge and C the slow baseline wander that
survives pedestal + common-mode subtraction. h is the chain's own DreamShaper (SK filter
+ residual CSA-pole high-pass), stretched in time by s; s = 1 is the nominal 180 ns.
A, B, C are linear and fitted at every s; chi2(s) says which time scale the noise wants.

Data: det3 27 June long run (the T14 target), decoded_root (raw, pedestal and CNS not
applied), quiet traces only, pedestal = per-channel median, CNS = per-64-block median.
The unbiased lag estimator averages every in-trace pair; no per-trace mean is removed.

Usage: python3 noise_shaper.py <decoded_07.root> [<decoded_08.root> ...]   (from the repo root)
"""
import sys
import numpy as np
import uproot

sys.path.insert(0, ".")
from response.dream.shaper import DreamShaper

N_CHAN, BLK, DT = 512, 64, 60.0


def stack(path, n_ev=4000):
    t = uproot.open(path)["nt"]
    a = t.arrays(["amplitude"], entry_stop=n_ev, library="np")["amplitude"]
    lens = np.array([len(x) // N_CHAN for x in a])
    ns = int(np.bincount(lens).argmax())
    return np.stack([x.reshape(ns, N_CHAN) for x, l in zip(a, lens) if l == ns]).astype(np.float32), ns


def acf_data(path, quiet_k=3.0):
    S, ns = stack(path)
    ped = np.median(S, axis=(0, 1))
    sub = S - ped
    cms = np.median(sub.reshape(sub.shape[0], ns, N_CHAN // BLK, BLK), axis=3)
    res = sub - np.repeat(cms, BLK, axis=2)
    sig = 1.4826 * np.median(np.abs(res), axis=(0, 1))
    quiet = np.abs(res).max(axis=1) < quiet_k * sig[None, :]
    q = res.transpose(0, 2, 1)[quiet]                 # (n_trace, ns)
    q = q / sig[np.nonzero(quiet)[1]][:, None]        # each channel in its own sigma units
    R = np.array([(q[:, :ns - k] * q[:, k:]).mean() for k in range(ns)])
    return R / R[0], float(quiet.mean()), ns


def model_terms(s, ns, filter_model="sk"):
    dt = 1.0
    sh = DreamShaper(dt_ns=dt, filter_model=filter_model)
    t = sh.t
    h = np.interp(t / s, t, sh.h, right=0.0)          # stretched response
    hp = np.gradient(h, dt)
    def ac(x):
        n = len(x); f = np.fft.rfft(x, 2 * n)
        r = np.fft.irfft(f * np.conj(f))[:n]
        return r
    Rp, Rs = ac(h), ac(hp)
    lags = (np.arange(ns) * DT).astype(int)
    return Rp[lags] / Rp[0], Rs[lags] / Rs[0]


def fit(R, ns, filter_model):
    out = []
    w = np.ones(ns); w[0] = 0.0                       # lag 0 is the normalisation
    for s in np.arange(0.40, 1.81, 0.02):
        Rp, Rs = model_terms(s, ns, filter_model)
        M = np.stack([Rp, Rs, np.ones(ns)], 1)
        # A, B >= 0 enforced by a coarse scan on the ratio, C free
        best = None
        for fs in np.linspace(0, 1, 41):
            Mm = np.stack([(1 - fs) * Rp + fs * Rs, np.ones(ns)], 1)
            coef, *_ = np.linalg.lstsq(Mm * w[:, None], R * w, rcond=None)
            if coef[0] < 0:
                continue
            r = (R - Mm @ coef) * w
            c2 = float(np.sum(r ** 2))
            if best is None or c2 < best[0]:
                best = (c2, fs, coef)
        out.append((s, *best))
    return out


def main():
    for path in sys.argv[1:]:
        R, qf, ns = acf_data(path)
        print(f"\n{path.split('/')[-1]}: quiet fraction {qf:.3f}, {ns} samples x {DT:.0f} ns")
        print("data R(k):", " ".join(f"{x:.3f}" for x in R[:12]))
        for fm in ("sk", "crrc2"):
            res = fit(R, ns, fm)
            c2 = np.array([r[1] for r in res]); i = int(np.argmin(c2))
            s, c2b, fs, coef = res[i]
            ok = np.array([r[0] for r in res])[c2 <= c2b * 1.5]
            Rp, Rs = model_terms(s, ns, fm)
            mod = coef[0] * ((1 - fs) * Rp + fs * Rs) + coef[1]
            r1 = [r for r in res if abs(r[0] - 1.0) < 1e-9][0]
            print(f"  {fm:5s}: best stretch s = {s:.2f} (peaking ~{180*s:.0f} ns), series fraction {fs:.2f}, "
                  f"const {coef[1]:+.3f}; chi2 best {c2b:.4f} vs s=1: {r1[1]:.4f}; s within x1.5 chi2: "
                  f"{ok.min():.2f}-{ok.max():.2f}")
            print("        model R(k):", " ".join(f"{x:.3f}" for x in mod[:12]))


if __name__ == "__main__":
    main()
