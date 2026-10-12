# DREAM front-end settings — read back from the archived run cfgs

**Settled 2026-10-12. Do not re-open without a new cfg or a bench measurement.**

The DAQ saves the cfg it actually loaded next to the raw data
(`<area>/Run|runs/<run>/<subrun>/raw_daq_data/*.cfg` on EOS
`/eos/experiment/ntof/data/x17/`). Decode any of them with

```
python ~/PycharmProjects/nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py <file.cfg>          # full decode
python ~/PycharmProjects/nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py --sweep <dir> ...   # distinct settings
```

Per-run table for every archived cfg (cosmic bench June, n_TOF July, SPS P2 July;
3574 files, pulled 2026-10-12): `dream_cfg/dream_settings_by_run.csv`, made by
`scratchpad/2026-10-12/dream_settings_by_run.py`.

## The T14 target run

`cosmic_bench/june_tests/Run/mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V`
(copy: `dream_cfg/satdet3_long_run_resist_490V_drift_1000V_CosmicTb_MX17.cfg`).
All 25 subruns of the scan, and their `dream_run/` copies, have identical
DREAM register lines. August's audit said this run "archived no .cfg"; it did.

| setting | data FEUs 7, 8 (det3 X, Y) | FEU 1 (M3 trigger) | source |
|---|---|---|---|
| register 1 | `0x081F 0xD023` → 0x081FD023 | `0x881F 0xD043` → 0x881FD043 | |
| **shaper peaking time** (5→100 %, Table 9) | **180 ns** (code 2) | 283 ns (code 4) | reg1<7:4> |
| CSA Rf·Cf | 5 µs | 50 µs | reg1<31> |
| input polarity / VICM | negative (CSA current out) / 1.35 V | same | reg1<11>, <13:12> |
| external input | none (CSA in path) | none | reg1<29:28> |
| integrator mode | off (normal PZC) | off | reg1<10> |
| CSA input current | 1.2 mA | 1.2 mA | reg1<1:0> |
| test capacitor / test mode | 100 fF / nothing | same | reg1<3:2>, <9:8> |
| discriminator | 17.5 % range, ThDAC −127, HIT enabled, 4 µs veto | same | reg1<30>, <21:14>, <27>, <22> |
| **gain** | code 2 on all 64 ch → **200 fC range, 10 mV/fC** | same | reg 6/7, 2 bits/ch |
| discriminator mask | reg 8/9 = 0xFFFFFFFF (all masked; ext. trigger) | same | |
| register 2 | 0x00000008 (meaning not decoded; manual not to hand) | same | |
| pulser (reg 3/4), FEU pulser | 0 / disabled | same | |
| latency | 45 samples = 2.7 µs | same | reg 12 |
| samples × period | 32 × **60 ns** (WrClk_Div 6 on 100 MHz; RdClk_Div 6) | same | DrmClk, Main_Conf_Samples |
| trigger / mode | external (TCM), raw, no ZS, no CM subtraction | | Sys DaqRun, Feu_RunCtrl |

## Word order of the 4-word register lines (why 180 and not 123 ns)

`Feu f Dream d R w0 w1 w2 w3` writes register R as 16-bit words **most
significant first, left-justified**: reg 1 = (w0 << 16) | w1. Read the
other way (0xD023081F), the same line would give 123 ns. That reading is ruled
out by the register itself and by the other cfgs:

1. **Field coherence.** Low-word-first gives a front end that cannot work:
   CSA bypassed (external input = SK filter), positive input polarity with
   a negative-polarity VICM (1.25 V), 5 % discriminator range, HIT off.
   High-word-first gives a coherent Micromegas setup: negative polarity +
   matching VICM 1.35 V, CSA in, 17.5 %, ThDAC −127.
2. **Deliberate overrides decode to deliberate changes.** The trigger-FEU
   override `0x881F 0xD043` changes only CSA τ (5 → 50 µs) and peaking
   (180 → 283 ns) under high-word-first. Under low-word-first its peaking
   would be unchanged, even though bits 4–7 of the second word are exactly
   what the cfg comment ("Peaking time : 0:76 1:123 2:180 …") and the
   commented-out variants (`0xD043`, `0xD063`) edit. The SiPM overrides in
   `dream_configs/dream_config_saclay_sipm_cosmics` (`0x181F 0xE823`,
   `0x281F 0xE823`) flip three fields together: positive polarity, the
   positive VICM pair (1.55 V), CSA bypassed (SK / gain-2 input). That is one
   coherent "SiPM mode", and it only reads that way high word first.
3. **The hardware followed the high-word-first encoding.** The July
   self-trigger scan (`make_selftrig_cfgs.py`) and the SPS P2 self-trigger
   runs set ThDAC (bits 14–20, straddling the two words) by encoding
   high word first, and produced working threshold scans (SPS P2 `run_1`
   Trg FEUs: ThDAC −19/−22/−29/−60).

Not yet used as a witness: the CEA FEU software that parses these lines
(`~/Feu/Firmware5/Distribution/Sources/Software/` on the DAQ PC, not
reachable from the desktop on 2026-10-12). It would be a fourth confirmation,
but the three above are already conclusive.

## Every campaign (data FEUs)

| campaign | peaking | CSA τ | range (gain) | samples × period | latency |
|---|---|---|---|---|---|
| cosmic bench June (10 runs) | 180 ns | 5 µs | 200 fC (10 mV/fC) | 32 × 60 ns | 45 |
| n_TOF July (207 runs) | 180 ns | **50 µs** | **600 fC (3.3 mV/fC)** | 32 × **20 ns** (runs 1–12, 14, 15, 20–22, 25–28, 38, 41, 44, 46; runs 30 and 32 mixed 20/60) or 32 × 60 ns (the rest) | 2–100, run-dependent |
| SPS P2 July (73 runs) | 180 ns | 5 µs | 200 fC | **16** × 60 ns | 24–40 |

All data FEUs in all three campaigns: 180 ns, negative polarity, CSA in.
Per-run detail is in the CSV. The July cfgs also carry `Sys NbOfSamples`
(20 / 400) that differs from the FEU's `Main_Conf_Samples` (32). The table
reports the FEU value; which one the readout honours was not checked here.

## What the cfg cannot tell you

- **The actual analog response.** 180 ns is the manual's nominal value for
  code 2 (Table 9). Chip spread and the PZC trim (β in
  `response/dream/shaper.py`) are not in any register. One FEU pulser shot
  through the test capacitor would measure them.
- Register 2's bit meaning (0x00000008 here; per the manual's Table 3,
  bit 26 is HIT_multi level 4, as used in `make_selftrig_cfgs.py`).
