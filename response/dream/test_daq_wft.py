#!/usr/bin/env python3
"""
test_daq_wft.py — T13 acceptance: wft reads the simulation UNCHANGED.

Not "the file opened" but "the real consumer's own access pattern returns what
we wrote, with the right numbers in the right places". Four checks:

  1. the file is a TTree with the data's branch types (a RegularArray or an
     RNTuple opens fine and is unreadable to wft — both were hit in
     development, see daq.write_decoded);
  2. wft's FeuReader recovers the pedestal and the noise we generated;
  3. a known injected charge comes back at the expected ADC amplitude, which
     is the only end-to-end check of the derived 20.48 ADC/fC scale;
  4. channel/sample orientation is not transposed — the signal must appear on
     the channel it was injected into, at the sample it was injected at. A
     transposed flatten survives every other check.

wft imports scipy via its package __init__, which the system python lacks, so
wft/io.py is loaded directly by path. If nTof_x17 is absent the wft-dependent
checks are SKIPPED and reported as skipped, never silently passed.

    python3 -m response.dream.test_daq_wft
"""

from __future__ import annotations

import importlib.util
import json
import os

import numpy as np

from .daq import ADC_PER_FC, Daq, N_CHAN, verify, write_decoded

WFT_IO = os.path.expanduser("~/PycharmProjects/nTof_x17/wft/io.py")
NOISE = os.path.expanduser("~/x17/response_sim/dream/noise_det3.json")

INJ_CHAN, INJ_SAMPLE, INJ_FC = 100, 11, 5.0


def _load_wft_io():
    """Import wft/io.py without triggering the wft package __init__."""
    if not os.path.exists(WFT_IO):
        return None
    spec = importlib.util.spec_from_file_location("wft_io_standalone", WFT_IO)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:                       # scipy/uproot missing etc.
        print(f"  wft/io.py did not import ({type(e).__name__}) — SKIP")
        return None
    return mod


def main(n_ev=40, seed=7):
    if not os.path.exists(NOISE):
        raise SystemExit(f"missing {NOISE}; run response.dream.noise "
                         "--characterise first")
    with open(NOISE) as f:
        spec = json.load(f)

    d = Daq(spec, seed=seed)
    q = np.zeros((n_ev, d.n_sample, N_CHAN), dtype=np.float32)
    q[:, INJ_SAMPLE, INJ_CHAN] = INJ_FC
    adc = d.to_adc(q, n_ev=n_ev)

    path = os.path.join(os.path.dirname(NOISE), "sim_decoded_acceptance.root")
    write_decoded(path, adc)

    ok = True

    # --- 1. schema ----------------------------------------------------------
    import uproot
    f = uproot.open(path)
    cls = list(f.classnames().values())[0]
    good = cls == "TTree"
    ok &= good
    print(f"  container            {cls:34s} {'OK' if good else 'FAIL'}")
    # What governs readability is the INTERPRETATION, not the typename string.
    # A known and accepted difference from the data: uproot's TTree writer
    # emits jagged branches as counter-based variable-length arrays
    # ("uint16_t[]", AsJagged with no header) while the DAQ writes true STL
    # vectors ("std::vector<uint16_t>", AsJagged with header_bytes=10). uproot
    # reads both as AsJagged and wft is uproot-based, so this does not affect
    # T13 -- but a ROOT/C++ consumer expecting std::vector would see the
    # difference, so it is asserted as AsJagged and stated rather than hidden.
    want = {"eventId": "AsDtype", "timestamp": "AsDtype", "ftst": "AsDtype",
            "sample": "AsJagged", "channel": "AsJagged",
            "amplitude": "AsJagged"}
    for b, kind in want.items():
        got = type(f["nt"][b].interpretation).__name__
        good = got == kind
        ok &= good
        print(f"  branch {b:12s}  {f['nt'][b].typename:20s} {got:12s} "
              f"{'OK' if good else 'FAIL'}")

    # --- 2. our own round-trip ---------------------------------------------
    bad = verify(path, adc)
    ok &= bad == 0
    print(f"  self round-trip      {bad} mismatched events"
          f"{'':17s}{'OK' if bad == 0 else 'FAIL'}")

    # --- 3/4. wft ----------------------------------------------------------
    io = _load_wft_io()
    if io is None:
        print("\n  wft checks SKIPPED (nTof_x17 not importable)")
        print("\n" + ("PASS (partial)" if ok else "FAIL"))
        return 0 if ok else 1

    r = io.FeuReader(path, n_ped=n_ev)
    checks = [
        ("n_sample", r.n_sample, d.n_sample, 0),
        ("pedestal median [ADC]", float(np.median(r.ped)),
         spec["pedestal_mean_adc"], 0.05),
        ("noise median [ADC]", float(np.median(r.noise)),
         spec["sigma_chan_median_adc"], 0.15),
    ]
    for name, got, wnt, tol in checks:
        rel = 0.0 if wnt == 0 else abs(got - wnt) / abs(wnt)
        good = (got == wnt) if tol == 0 else rel <= tol
        ok &= good
        print(f"  {name:20s} {got:10.2f}  want {wnt:8.2f}  "
              f"{'OK' if good else 'FAIL'}")

    W = None
    for _ev, _ftst, w in r.iter_events():
        W = w
        break
    good = W.shape == (N_CHAN, d.n_sample)
    ok &= good
    print(f"  W shape              {str(W.shape):20s} want "
          f"{(N_CHAN, d.n_sample)}  {'OK' if good else 'FAIL'}")

    # The injected charge, after wft's pedestal + CNS. CNS subtracts a
    # 64-channel block median, which for a single hot channel removes almost
    # nothing, so the peak should sit close to q*ADC_PER_FC.
    # Averaged over events: a single event carries ~10 ADC of noise on a
    # ~100 ADC signal, so a one-event check is noise-dominated and would pass
    # or fail on the seed rather than on the ADC scale.
    want_adc = INJ_FC * ADC_PER_FC
    amps = []
    for _ev, _ftst, w in r.iter_events():
        amps.append(float(w[INJ_CHAN, INJ_SAMPLE]))
    got_adc = float(np.mean(amps))
    rel = abs(got_adc - want_adc) / want_adc
    good = rel <= 0.10
    ok &= good
    print(f"  injected {INJ_FC:.0f} fC       {got_adc:10.1f}  want "
          f"{want_adc:8.1f} ADC  ({rel:.1%}, mean of {len(amps)})  "
          f"{'OK' if good else 'FAIL'}")

    # Orientation: the peak must be at the injected (channel, sample).
    pk = np.unravel_index(np.argmax(W), W.shape)
    good = pk == (INJ_CHAN, INJ_SAMPLE)
    ok &= good
    print(f"  peak location        {str(pk):20s} want "
          f"{(INJ_CHAN, INJ_SAMPLE)}  {'OK' if good else 'FAIL'}")

    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
