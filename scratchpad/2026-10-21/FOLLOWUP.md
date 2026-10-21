# det3 X asymmetry: follow-up to the audit (2026-10-21 evening, dylan-MS-7C84)

Continues `audit/AUDIT.md` §4.3. Scripts and JSON in the subdirectories named below; data = raw survey (NB_RAW=~/x17/nb_asym_raw),
sim arms copied from yoga `~/x17/nb_asym_raw/{simgapleak,simfunnel,simxoff}`.

## 1. The even term retires the whole attractor family (`even/`)
- **Every sim arm** (`sim_even.py` → `sim_even.json`, 25 arms: gap leak, funnel, strip width 410/690, rho_s 0.5/2/5, anisotropy
  3/10, rigid offsets): E/O = 1.2–2.5, slope of s1 vs |a1| across delta bins +2.2 … +4.5. No lever changes the ratio.
- **Data, same PCB** (`same_pcb.py`): det3 May vs det3 June-15 (one PCB bulked twice), X sharing at matched Y share (all events,
  no labels): equal to ±0.02 over 7 bins (0.568/0.550, 0.653/0.654, 0.708/0.728, 0.770/0.794, 0.845/0.824, 0.884/0.868,
  0.846/0.851) while a1 goes 0.00 → 0.14. The attractor family predicts +0.3 (sim) or ≥ +0.1 (toy).
- **Data, within det2** (`det2_profile.py`, 32 cells, a1 −0.30 … +0.26, Y share and Y ±2 as width controls): per unit a1 the
  X content changes by −0.54 / +0.47 at ∓1, −0.08 / +0.04 at ∓2, ~0 at ∓3 (odd), and by −0.09 … −0.10 at ±1, −0.07 at ±2 per
  unit |a1| (even, slightly **narrower**). det3 June-15 has the same shape (−0.32 / +0.19, ∓0.02).
- **Shape of the transfer** (`antisym.py`, May → June-15 at matched width): the ±1 change is antisymmetric, but neither linear
  operator fits the rest. A nearest-neighbour antisymmetric coupling puts ~40 % of the ±1 odd part at ±2 (data ~12 %); a
  one-sided copy gives n−1 no loss. The ±2/±3 even drop between May and June-15 (−0.05/−0.02) is mostly a bulk difference
  (det2's within-board even slope predicts only −0.01 at a1 0.14), e.g. a different photon-feedback halo.
- **Reading:** a direct n−1 → n+1 transfer, local to where the charge sits, that conserves the sum and leaves the cluster no
  wider (slightly narrower). Any explanation must have Δs1 ≤ 0 with |a1|. The registration/attractor sims fail by ~20×.

## 2. Time, gain and the session jumps (`time/time_a1.py` → `time_a1_det3.json`)
- On the direct a1 (not the copy estimator, which misreads low-amplitude runs: the 6-23 early runs "0.035/0.057" are
  +0.06 ± 0.05 / +0.16 ± 0.02 on a1, consistent with the night's 0.19), every det3 June-15 session falls on one curve in the
  median Y share: 0.19–0.20 at 0.5, 0.16 at 0.56–0.65, 0.14–0.15 at 0.75, 0.07–0.11 at 1.0.
- Flat with gain at fixed width: 6-23 resist 525 → 465 V (pulse ×6 smaller), a1 0.15–0.21 throughout; 6-22 likewise.
- Flat over hours: 6-22's 8 h scan, 6-27 afternoon → overnight. No charging transient at the survey's resolution (≤ 2 files
  per sub-run; within-run splits agree).
- Charging-up would scale with deposited charge; the data scale with cloud width only.

## 3. 6-19/20 was det2, mislabelled (`s619/`)
- The "det3 FEU 7" 6-19 data (a1 0.03, all 8 chips lit) has det2's readout-chain noise pattern (r = 0.957 with det2 X on
  6-22) and not det3's on the same FEU 7 slot on 6-17 (`fingerprint619.py`). The pattern follows the chain + board, not the
  FEU slot (`fp_validate.py`: same board, other FEU, r ≈ 0; different boards median 0.04).
- Board properties agree (`board_props.py`): FEU 7 lights all 8 chips and reads a1 +0.026 like det2; det3's X readout on
  6-17 lit only chips 0–3, and on 6-19 that 4-chip readout is the FEU 3 file (372 events, a1 +0.077 ± ~0.03, s1 0.75).
- So the one falsifier of "fixed map × cloud width" is a labelling error. Real det3 on 6-19 is low-statistics and ~2σ low.
  The raw_daq_data HV scan + long run could still be decoded for more FEU 3 events.

## What this leaves
- The fingerprint (audit §4.1) gains: (2) the transfer is n−1 → n+1 direct, ±2 odd ≈ 12 % of ±1, Δs1 ≤ 0 (−0.19 ± 0.06 per |a1|);
  (5) gain- and time-independent at fixed width; no session outlier.
- Excluded on top of the audit: every attractor arm simulated this week, a fixed linear channel operator (copy or
  antisymmetric coupling), charging-up.
- Open: what physical, position-level, sum-conserving mechanism moves charge from the n−1 pad to the n+1 pad without
  widening the cluster, X-only (across the strips), with a smooth two-signed geometry per board, June bulks only, and an
  amplitude set by the cloud width. The skewed-induction family (audit §4.2: an offset, partly conducting strip steering the
  induced signal, or tilted ion drift in a lateral field between strips at different potentials) is the remaining
  untested class; it should now be scored on the odd shape (±1 vs ±2) and the even term before anything else.
