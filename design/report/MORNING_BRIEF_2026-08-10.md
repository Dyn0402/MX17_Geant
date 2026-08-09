# Morning brief — overnight of 2026-08-09 → 08-10

Five-minute read. Full working record with every table:
`design/report/OVERNIGHT_2026-08-10.md`.

---

## Two things need you

**1. The desktop needs a Tailscale SSH re-auth.** Every new session is refused
pending an interactive check-in; the host itself is fine (ping 16 ms, port 22
open). The link is regenerated per attempt, so just run `ssh desktop true` and
visit whatever URL it prints.

Nothing was lost to it: the T7 slope-hunt chain is running there unattended and
self-merges, and the one job I needed off that machine (the 08-08→09 HV-scan
verification) turned out to be doable from EOS instead.

**2. When the slope hunt lands, it is one command:**

```bash
python3 -m response.validation.slopehunt_verdict \
    --calib  ~/x17/response_sim/avalanche/aval_calib_slopehunt.json \
    --slopes ~/x17/response_sim/hv_slope/slopes.json
```

The decision rule is pre-registered in that script's header, written before the
data existed. See "the decision that is yours" below.

---

## The three threads, and where each now stands

### Thread 1 — rise time. The ion term is the whole discrepancy. ✅ measured

The rise mismatch used to rest on one statistic (p5 of a vertical sample). It
now rests on the whole distribution at two inclinations.

At 10° and 20°, the ion term acts as a **near-rigid ~95 ns delay** across every
quantile, and the data needs about **75 of those 95 ns removed**. Removing the
ion term entirely slightly *overshoots* — the sim becomes 10 ns (10°) to 24 ns
(20°) too fast everywhere.

> **The data's rise is reproduced at f_eff ≈ 0–0.25, on six quantiles at two
> independent inclinations, against a defended f_ion = 0.9056.**

Nothing was fitted. The legs were already on disk; they were read at matched
quantiles. **This is the same factor-4–5 contradiction as before, but it is now
much harder to dislodge.**

Track inclination — billed as the last live candidate — is **eliminated**, and
along the way the committed claim that "the sim barely responds to inclination"
turned out to be an artifact of a 200 ns threshold sitting below the sim's own
rise floor. At 240 ns the sim goes 0.045 → 0.749 across 0–20°. Corrected in the
nTof_x17 note.

### Thread 2 — amplitude. One candidate left, and two numbers point at it

The deficit is best quoted **in charge, not peak**: `q_sum` is *invariant* under
f_ion (0.626–0.642 across the entire range) because f_ion redistributes charge
in time and conserves it. So the charge deficit **cannot be double-counted**
against thread 1 — however f_eff resolves, ×0.63 survives.

| | |
|---|---|
| charge deficit | **×0.63** (f_ion-independent) |
| gain it demands at 490 V | **~4 × 10⁴** against the sim's 24 094 (**×1.6**) |
| independent HV-slope error | sim **0.3106 ± 0.0033** vs data **0.4487 ± 0.0093** per 10 V (**×1.44**) |

Every other row of the ledger is now dead, controlled, or separated — including
the two that owed numbers last night:

* **primary ionisation — measured and correct.** 91.2 e⁻/cm (lit 90–100),
  implied W = 25.97 eV (lit ~26), 3.83 e⁻/cluster (lit ~3.5). Cannot carry ×1.6.
* **diffusion dilution — closed by construction.** An integral cannot be diluted
  by spreading, and the deficit is quoted in an integral.

Also dead: ADC scale (1.2 %, a control row), electronics gain range (44/44
archived cfgs), ρ_s (×1.15 over a *factor 10*), β (0.6–2.3 %), mesh
transparency, W2 prompt capture (already applied).

> **The avalanche gain is the only survivor, and a single α(E) / Penning error
> at the operating point would produce both the ×1.6 and the ×1.44.** Neither
> number was derived from the other.

### Thread 3 — X/Y asymmetry. Real, sim-side, and modest

The detector is X/Y symmetric; the simulation is not. Both terms are now
established on 738 pooled paired events:

| | value | significance |
|---|---|---|
| peak Y/X, sim ÷ data | **0.8190 ± 0.0171** | 10.6 σ |
| charge-partition term | 0.8823 ± 0.0156 | 7.6 σ |
| extra-spreading term | 1.0777 ± 0.0254 | 3.1 σ |

**An 18 % effect made of two ~10 % pieces** — the sim puts ~10 % less charge
into Y and spreads it ~10 % more. Charge is *not* lost: both legs show Y
carrying more charge with a lower peak, because Y spreads over more channels.
The detector does this too; the sim overdoes it.

Excluded: the data (symmetric once de-biased, converging to 1.0 by two
independent corrections), selection (the sim's Y/X is immune to pairing and
de-saturation), and the S1 electrostatics (a strips-vs-uniform A/B moves Y/X by
0.02–0.03, and the wrong way). **It lives in Stage B/C.** Note `kY` is probably
*not* the culprit — it is applied identically to both legs, so it cannot create
a sim-only asymmetry on its own.

---

## The standing open question — now fully cornered

**What suppresses the ion-induction term at the readout by ×4–5?**

Every mechanism proposed since 2026-08-08 is eliminated by measurement or by
structure: the charge split (0.9056 through the real mesh, and the shift is the
*wrong way*), the template shape (validated to 3.8 % at every quantile), ion
species (Blanc: 4 %, wrong direction), T10 lateral factorisation (3.7 ns), β
(4 ns across its full range), any missing high-pass (falsified *as a class* by
the rise/undershoot exchange rate), the peaking register (code 2, 44/44 cfgs),
the amplification-gap geometry (150 µm, from the pillar gerber), track
inclination (closed last night), and — re-derived independently last night —
resistive-sheet screening.

That last one was the soft spot, so it got a second, deliberately different
derivation. The result **hardens** the retirement rather than merely repeating
it:

> The weighting potential factorises as `Ψ_sheet(k,τ) × cont(k,z)` with the
> continuation factor **time-independent**, because the gas is source-free and
> bounded by the grounded mesh. An in-gas source therefore sees *exactly* the
> same temporal sheet response as an on-sheet source. **The sheet cannot
> distinguish induced charge from deposited charge at all.** The ρ_s-dependent
> part is 0.34 % across a factor 10 in ρ_s, where ×4.5 was needed.

So the elimination is not just complete, it is **structural**. Which sharpens
what remains: the suppression must live either **outside the current chain
decomposition**, or in an **assumption every defended piece inherits**. The
remaining assumption surface, named explicitly rather than left vague:

* the **weighting-field family** — every piece above is computed against the
  same electrostatic model of the stack;
* the **DAQ frame / t0 definition** shared by both legs;
* the **mapping of the rise metric onto the model** (10–90 % on a shaped
  waveform, against a model quantity that is not obviously the same thing).

⚠️ One disagreement left open on purpose: my re-derivation is flat in ρ_s, the
original swings 26 %. **Both agree the effect cannot deliver ×4.5**, so the
retirement is robust either way, but they cannot both be right about the
ρ-dependence. My reading is that the original double-counts the sheet — I have
not proven it, and I am not calling their number wrong on the strength of a
derivation that itself needed a bug fix last night. Daylight item.

---

## The decision that is yours

**When the slope hunt lands**, the pre-registered rule reads:

| outcome | condition | consequence |
|---|---|---|
| **A** | best rP fixes the slope **and** gain ≥ ×1.4 | **The ledger closes on a single defect.** Adopt rP\*, re-run Stage B/C, amplitude explained. |
| **B** | slope fixed, gain < ×1.2 | **No candidate survives** — the chain decomposition itself is suspect. The interesting outcome. Do *not* patch it with a fitted gain factor. |
| **C** | no rP reaches the data slope | Penning is not the knob; α(E) is wrong in the cross sections or the field map. |
| **D** | slope fixed, gain > ×2.2 | The two threads are not one defect. Report that; don't split the difference. |

Note the ×1.6 (gain demanded) and ×1.44 (slope error) sit either side of
outcome A's ×1.4 gate, which is deliberate.

**Standing rule, same as β and undershoot:** rP\* is *fitted* to the slope, so
slope agreement is not evidence. The gain is the out-of-sample prediction and
the only part that can confirm anything.

**Two other things are yours, neither urgent:**

* **The f_eff question above** — whether to open the assumption surface, and
  which part first.
* **Whether the 18 % X/Y term is worth chasing.** It is real and now
  well-measured, but it is 18 %, not a factor.

---

## Also landed, briefly

* **08-08→09 HV-scan chain verified and re-certified** without the desktop (the
  merged product was already on EOS). All six certs pass on both gases —
  including "one distinct field map per voltage" at 8/8 and 7/7, the direct
  guard against the T7 voltage-label incident recurring. Cross-checks at 490 V:
  gain 24 172 vs the pooled 24 094, survival 0.9518 vs T6's independent 0.955.
  *Tomorrow's first chore is already done.*
  One line worth keeping: f_ion reads **0.9006 exactly** from the template's own
  integrals, confirming the shipped calib still uses the parallel-plate ψ rather
  than the through-mesh 0.9056 — a known ≤0.005 understatement, not a new bug.

* **det3 contaminant family — prefers O₂-like attachment.** Across the
  30-mixture water grid, **29 give η = 0 exactly**: water does not attach at any
  fraction, so *any* finite decay excludes attachment-free transport, and det3
  decays at 5.6–11.1 mm at every drift field. Agrees independently with the
  freeze queue's observation that 0.8 % H₂O + ~1 % air best reproduces the bench
  v_drift. **Caveat: the λ(E) *shape* test does not close** — every Magboltz
  curve rises with field while det3's falls — so the preference rests on the
  decay existing, not on its shape. Family-constraint inference, not a
  measurement; no concentration may be read off it.

* **Wet gain bracket submitted** (condor 16705137, 32 jobs, running). Two
  findings while setting it up: the roadmap's "~2 new Magboltz jobs" was wrong —
  **zero are needed**, the tables already exist at amp range (the `.gas` files
  store *E/p*, which reads as a drift table until multiplied by 745.83 Torr) —
  and the bracket as specified would have confounded water with a Penning-model
  change, so all arms now run at rP = 0.40 with dry-at-auto carried separately.

---

## Corrections made to the committed record

Four, all dated and loud, per project habit — including two of my own:

| document | what was withdrawn |
|---|---|
| `nTof_x17 ANGLED_LADDER_2026-08-09.md` §4 | "the sim barely responds" to inclination — threshold artifact |
| same | open thread #4 (vertical broadening = one-view θ window) — **tested, answered NO** |
| `GAS_AND_DRIFT_CAGE_ROADMAP_2026-08-08.md` §2 | "~2 new Magboltz jobs" — zero needed; plus the Penning confound |
| `OVERNIGHT_2026-08-10.md` §8 → §10 | **my own** factor-0.74 X/Y claim — it is 18 %, and the comparison behind 0.74 was invalid |

Two further self-caught errors are recorded in the working log rather than
corrected in the record, because they were caught before anything rested on
them: the sheet-screening re-derivation's pre-registered falsifier fired on its
*author's* bug first, and an avalanche certification read FAIL on my own wrong
sign convention rather than on the data.
