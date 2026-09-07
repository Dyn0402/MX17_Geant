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
     transposed flatten survives every other check;
  5. ftst survives the round trip with the semantics the real FEUs use, and
     the X/Y pair carries the fixed board offset (Fix 3);
  6. the decoded file NAMES are ones wft can both discover and pair (Fix 8).

THE PEDESTAL EVENTS ARE SIGNAL-FREE, and that is check 3's whole validity.
The pre-fix version injected into EVERY event and then let FeuReader build its
pedestal from those same events, so the injection rode into the pedestal
through the common mode and 5 fC came back as 98.0 ADC against 102.4 predicted
— a 4 % "discrepancy" that was purely the test contaminating itself (audit B2;
with clean pedestal events it closes at 102.0/102.4). The first `n_ped` events
are now empty, and the tolerance is tightened accordingly: at 10 % this check
would have passed a +13 % gain error.

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

from .daq import (ADC_PER_FC, Daq, N_CHAN, TriggerPhase, verify,
                  write_decoded)

WFT_IO = os.path.expanduser("~/PycharmProjects/nTof_x17/wft/io.py")
NOISE = os.path.expanduser("~/x17/response_sim/dream/noise_det3.json")

INJ_CHAN, INJ_SAMPLE, INJ_FC = 100, 11, 5.0
# Tolerance on the end-to-end ADC scale. 10 % (the pre-fix bar) would pass a
# +13 % gain error; with a clean pedestal the check closes to well under 3 %.
ADC_SCALE_TOL = 0.03


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


def main(n_ev=40, n_ped=40, seed=7):
    if not os.path.exists(NOISE):
        raise SystemExit(f"missing {NOISE}; run response.dream.noise "
                         "--characterise first")
    with open(NOISE) as f:
        spec = json.load(f)

    d = Daq(spec, seed=seed)
    # SIGNAL-FREE pedestal events FIRST (FeuReader takes its pedestal from the
    # first n_ped entries), then the injected ones. See the header.
    n_tot = n_ped + n_ev
    q = np.zeros((n_tot, d.n_sample, N_CHAN), dtype=np.float32)
    q[n_ped:, INJ_SAMPLE, INJ_CHAN] = INJ_FC
    adc = d.to_adc(q, n_ev=n_tot)

    # One trigger phase per event, both views, real ftst semantics (Fix 3).
    trig = TriggerPhase(sample_period_ns=d.dt_ns)
    rng = np.random.default_rng(seed + 3)
    ph = [trig.draw(rng) for _ in range(n_tot)]
    ftst = {v: np.array([p[v][1] for p in ph], dtype=np.uint16)
            for v in ("X", "Y")}

    path = os.path.join(os.path.dirname(NOISE), "sim_decoded_acceptance.root")
    write_decoded(path, adc, ftst=ftst["X"])

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

    r = io.FeuReader(path, n_ped=n_ped)
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
    got_ftst = []
    for ev, ft, w in r.iter_events():
        got_ftst.append((int(ev), int(ft)))
        if W is None and ev >= n_ped:
            W = w                                # a SIGNAL event, not a pedestal one
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
    for ev, _ftst, w in r.iter_events():
        if ev < n_ped:                  # pedestal events carry no injection
            continue
        amps.append(float(w[INJ_CHAN, INJ_SAMPLE]))
    got_adc = float(np.mean(amps))
    rel = abs(got_adc - want_adc) / want_adc
    good = rel <= ADC_SCALE_TOL
    ok &= good
    print(f"  injected {INJ_FC:.0f} fC       {got_adc:10.1f}  want "
          f"{want_adc:8.1f} ADC  ({rel:.1%} vs {ADC_SCALE_TOL:.0%}, mean of "
          f"{len(amps)})  {'OK' if good else 'FAIL'}")

    # Orientation: the peak must be at the injected (channel, sample).
    pk = np.unravel_index(np.argmax(W), W.shape)
    good = pk == (INJ_CHAN, INJ_SAMPLE)
    ok &= good
    print(f"  peak location        {str(pk):20s} want "
          f"{(INJ_CHAN, INJ_SAMPLE)}  {'OK' if good else 'FAIL'}")

    # --- 5. ftst round trip (Fix 3) ----------------------------------------
    got = np.array([f for _e, f in sorted(got_ftst)], dtype=np.uint16)
    good = np.array_equal(got, ftst["X"])
    ok &= good
    print(f"  ftst round trip      {len(got):5d} events matched"
          f"{'':10s}{'OK' if good else 'FAIL'}")
    # Data semantics: uniform over the 6 states, and the X/Y pair separated by
    # the fixed board offset on EVERY event (never independently drawn).
    n_states = len(np.unique(ftst["X"]))
    good = n_states == trig.n_ftst
    ok &= good
    print(f"  ftst states          {n_states:5d} of {trig.n_ftst}"
          f"{'':21s}{'OK' if good else 'FAIL'}")
    dif = np.unique((ftst["Y"].astype(int) - ftst["X"].astype(int))
                    % trig.n_ftst)
    good = len(dif) == 1 and int(dif[0]) == trig.offset_ticks % trig.n_ftst
    ok &= good
    print(f"  ftst X->Y offset     {str(dif):20s} want "
          f"[{trig.offset_ticks % trig.n_ftst}]      "
          f"{'OK' if good else 'FAIL'}")

    # --- 6. wft can discover AND pair the decoded names (Fix 8) ------------
    from ..digitizer.run import decoded_path
    fx = decoded_path("/tmp/mx17/run/sub/decoded_root/MX17_sim", "sim_000", 3)
    fy = decoded_path("/tmp/mx17/run/sub/decoded_root/MX17_sim", "sim_000", 4)
    good = fx.endswith("_03.root") and fy.endswith("_04.root")
    ok &= good
    print(f"  feu suffixes         {os.path.basename(fx)} / "
          f"{os.path.basename(fy)}  {'OK' if good else 'FAIL'}")
    # The pairing tag must MATCH across the two views or reco's by_tag dict
    # silently never assembles a pair (wft/reco.py:449).
    tx, ty = io.file_tag(fx), io.file_tag(fy)
    good = tx == ty
    ok &= good
    print(f"  wft pairing tag      {tx!r} == {ty!r}  "
          f"{'OK' if good else 'FAIL'}")

    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
