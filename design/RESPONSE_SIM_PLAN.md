# MX17 Full Detector-Response Simulation — Implementation Plan

**Status:** approved plan, ready for implementation. **Date:** 2026-08-06.
**Authority:** this document. Where it conflicts with older notes, this wins. Update it as decisions land ("living document").
**Audience:** an implementing agent/developer who has NOT read the research behind this plan. Every stage gives the equation or the reference, the file contract, the host to run on, and acceptance criteria. When something is genuinely open, it is listed in §12 with the default to proceed with — do not block on §12 items.

**Coordination rule (active now):** a parallel effort is updating the Geant4 geometry (`src/DetectorConstruction.cc`, `shared/MX17ModuleGeometry.hh`). Until that merges, do NOT edit existing C++ files in this repo. All new work in this plan lives in NEW files/directories, and Stage A (§6) starts only after the geometry work is committed.

---

## 0. Goal and principles

Simulate the complete MX17 response to MIPs from first principles: Geant4 ionization → drift/diffusion → mesh transparency → avalanche gain → **charge spreading on the resistive strips + signal induction on the X/Y readout** → DREAM electronics → digitized waveforms that the existing `wft` reconstruction (in `~/PycharmProjects/nTof_x17`) can process **unchanged**.

Principles:
1. **First-principles first.** No parameter is tuned to MX17 detector data on the first pass. Measured MX17 quantities (sharing fractions, peak-time shifts, response library — §9) are *blind validation targets*. Only after the blind comparison do we iterate.
2. **Staged, file-joined pipeline** (adopted from `~/CLionProjects/p2_geant/docs/SIM_CAMPAIGN_PLAN.md` §2): expensive stages run once, cheap stages re-run per parameter point. Every stage reads/writes documented files; any stage can be re-run in isolation.
3. **Two-tier fidelity.** A rigorous slow path (Garfield++ delayed-signal formalism) validates a fast path (precomputed response-function convolution). The fast path is the production digitizer. Same physics, different caching.
4. **Everything scanned that is unknown.** Resistivity, coverlay thickness, amplification gap are scan axes, not guesses.

Pipeline overview:

```
S1  Weighting solver   Ψ_n(x,y,z,t)  dynamic weighting potentials      [new python]
S2  Mesh field solver  E(x,y,z) unit cell, transparency, funneling     [Garfield++ neBEM / Elmer]
S3  Avalanche calib    gain(V), Polya θ, σ0, ion current shape         [Garfield++ microscopic]
A   Geant4             ClusterTree upgrade: time + local coords        [this repo, C++]
B   Digitizer          clusters → per-strip induced current waveforms  [new python package]
C   DREAM electronics  currents → ADC samples, noise, ZS               [new python package]
V   Validation         wft reconstruction on sim output vs data        [nTof_x17 tooling]
```

---

## 1. Detector description and parameter table

Anode stack, top (gas) to bottom. z=0 at the top surface of the ESL resistive layer, +z toward the mesh.

| # | Layer | Value | Source | Confidence |
|---|---|---|---|---|
| 1 | Micromesh (grounded ref. for amp field; at −HV_mesh in reality) | woven SS, wire 2×19 µm crossing, fill factor 0.223 | `shared/MX17ModuleGeometry.hh` | good |
| 2 | Amplification gap | **150 µm — CONFIRMED (user, 2026-08-06). Not 128. Do not scan.** | user; supersedes the "unverified" flag in `design/GEOMETRY_IMPLEMENTATION_NOTES.md` | fixed |
| 2b | Pillars (inside the amp gap) | **Dynamask** dry-film (photoimageable solder mask; user 2026-08-06), Ø 0.60 mm on a 4.68 mm grid (3571 pillars), height = gap 150 µm; area coverage ≈ 1.3% | material: user; geometry: `3498A_bulk.gbr` | good |
| 3 | **ESL resistive strips** | **550 µm wide, 250 µm gap, 800 µm pitch, running along y ("vertical")**; terminated on copper bus strips at both y-ends only (user 2026-08-07); print thickness ~10 µm (user estimate, unconfirmed — does NOT enter the model: S1 treats the ESL as a zero-thickness sheet, thickness is absorbed into ρ_s) | user/fab knowledge; NOT in gerbers (screen-printed after fab) | geometry good; ρ_s unknown → scan **{0.5, 1, 2, 5} MΩ/sq** |
| 4 | Insulator between pad Cu and ESL | **kapton, thickness unknown** (NOT Dynamask — that is the pillar material only; user 2026-08-06). Default 75 µm, ε_r = 3.5 | user (material class); thickness assumed | scan **{50, 75, 125} µm** |
| 5 | Pad plane (Cu) | 512×512 pads, **0.68 mm square on 0.78 mm pitch**, active area 399.36 mm | `design/gerbers/readout_pcb/DFS3498A_L2-pads.gbr` | exact |
| 6 | Buried interconnect | pads bussed by vias to 512 Y-strip traces (L3-TrackY) and 512 X-strip traces (L4-TrackX), 0.1 mm traces on 0.78 mm pitch | gerbers | exact (in-plane); layer z-spacing unknown |

Critical geometric facts:
- **Pitch mismatch / beat:** resistive pitch 800 µm vs readout pitch 780 µm. The ESL-to-pad registration phase advances 20 µm per pitch and repeats every **LCM = 31.2 mm** (39 resistive strips = 40 readout pitches). The response kernel is therefore **position-dependent with a 31.2 mm superperiod** in x. The solvers and the digitizer must carry the absolute x position, not just position-within-one-strip. This beat is itself a physics prediction to look for in data (position-dependent sharing/residuals with 31.2 mm period).
- **Orientation and sharing anisotropy:** resistive strips run along y. Resistive transport moves charge along y → spreads signal across the *y-measuring* channels → the Y view has stronger, slower sharing (data: τ_Y ≈ 410 ns vs τ_X ≈ 230 ns; kY ≈ 1.8–2.9). Across x, gaps block DC transport; X-view sharing is diffusion + induction + weak inter-strip capacitance. The simulation must reproduce this asymmetry *from geometry alone* — it is a headline validation target.
- **Strip biasing/grounding — A1 CONFIRMED as hardware, τ_g REINTERPRETED (2026-08-07):** the ESL strips contact copper bus strips at both y-ends of the active area and **nothing in between** (user, 2026-08-07 — this is now a fact, not an assumption). The implied global drain is L²/(π²D) = 4–41 ms across the ρ_s scan. The quantity previously mislabeled here as "measured global drain constant τ_g = 5.3–7.3 µs" (nTof_x17 `rc_line_step2.py`) is **not a drain measurement**: it is fitted as T_Y = T_X ⊛ [δ + d/dt(Gaussian diffusion × e^(−t/τ_g))] on a ≤1.4 µs window with the measured X template as reference, so (i) any decay common to both views — electronics, a true global drain, HV sag — cancels by construction, and (ii) an ms-scale drain is a factor 0.9998 over the window, invisible. τ_g is the *extra decay of the Y view relative to X beyond the Gaussian toy model*, i.e. a kernel-shape observable. A drain-free S1 kernel of the true comb channel, pushed through the same toy fit, reproduces an apparent τ_g ≈ 2.2–3.3 µs for on-strip deposits and ∞ for gap deposits (charge level, ρ_s 1–2 MΩ/sq; `response/validation/tau_g_reinterpretation.py`); a strip/gap mixture lands in the measured 5.3–7.3 µs band, and the mechanism is position-flat along the strip and fleet-consistent, both as observed. Consequence for the solver: **`tau_drain_s=None` is the production default** (§3). The digitizer-level closure (run the actual rc_line fit on simulated waveforms) is task T13b.
- Gas: Ar/iso 95/5 (+~1% H2O on the bench — matters, measured v_drift 36.6 µm/ns is far below dry Magboltz). Drift gap 30 mm at 1000 V typical bench.

---

## 2. Repository layout and file contracts

New code lives in this repo under `response/` (Python ≥3.10, numpy/scipy only for the core; ROOT via `uproot` for I/O):

```
response/
  solver/          S1 semi-spectral weighting solver
  meshcell/        S2 unit-cell scripts (Garfield++/neBEM driver + collectors)
  avalanche/       S3 calibration campaign scripts + frozen JSON outputs
  digitizer/       B  cluster → strip-current pipeline
  dream/           C  electronics + ZS emulation
  validation/      V  closure scripts, blind-comparison figures
  common/          geometry constants (single source: parse shared/MX17ModuleGeometry.hh
                   values into python once, assert against the header at import)
  params/          *.yaml parameter sets (one file per scan point; git-tracked)
```

Large products (grids, libraries, waveform files) never go in git. **Canonical bulk store is EOS: `/eos/experiment/ntof/data/x17/response_sim/`** (user, 2026-08-07 — effectively unlimited under the nTOF allocation), laid out as `s1/ s2/ avalanche/raw/ clusters/ currents/`. `~/x17/response_sim/` on the laptop is a *working copy*, not the archive, and the desktop holds nothing permanently (20 GB free there is the tightest disk in the fleet).

Transfer routing. **All three hosts can write EOS** (desktop enabled 2026-08-07). Two details that are not obvious and will waste an afternoon if forgotten:

- There is **no `/eos` mount and no `eos`/`xrdcp` client on the desktop.** Its access is via `ssh lxplus` with GSSAPI credential delegation — `~/.ssh/config` there already sets `GSSAPIAuthentication`/`GSSAPIDelegateCredentials yes`, so a plain `rsync … lxplus:/eos/…` works. There is no local `/eos/...` path to write to, so any script that assumes one fails with "No such file or directory".
- That route depends on a **Kerberos ticket, which expires** (24 h; `klist` on the desktop to check). Before the ticket existed, desktop→lxplus could read EOS but writes failed with "Operation not permitted" and even `.bashrc` was unreadable — so a *silent* AFS/EOS permission failure in an unattended desktop job means the ticket lapsed, not that the code is wrong. Long campaigns should either renew it or fall back to desktop → laptop → EOS.

Anything produced on the desktop must be drained on a `.done` marker written *after* the producer returns, never on file existence: rsync of a file numpy is still writing produces a silently truncated archive (cost us one bad 1 GB transfer). Every product file embeds: git hash of `response/`, the parameter YAML content, and a UTC timestamp. **No un-manifested runs** — one CSV manifest row per production run (copy the discipline from p2 `SIM_CAMPAIGN_PLAN.md` §"bookkeeping").

Key file contracts (formats frozen here; extend, don't break):

| Product | Producer | Format |
|---|---|---|
| `wpot_<ch>_<params>.npz` | S1 | Ψ_n on grid: axes `x` (over 31.2 mm superperiod, ≤10 µm step near strip edges), `y` (relative, to ±25 mm), `z` {0, 8–16 slices to mesh}, `t` (60 log-spaced, 0.1 ns–10 µs); arrays `psi[t,z,y,x]`, prompt slice `psi_p[z,y,x]` |
| `greens_<params>.npz` | S1 post | G_n(x0,y0,t): induced charge on channel n (X and Y sets) for unit point charge landing at (x0,y0); same axes; this is Ψ on the z=0 plane by reciprocity |
| `meshfield.root/.txt` | S2 | E-field map of one weave unit cell, amp gap + last 200 µm of drift |
| `aval_calib.json` | S3 | per (gas, HV, gap): mean gain ḡ, Polya θ, transverse avalanche σ0, funneling offset map, ion current shape i_ion(t) (normalized), ion fraction to mesh, transparency ε(E_d/E_a) |
| `clusters.root` | A | ClusterTree (schema §6) |
| `currents_<ev>.npz` / batched | B | per-event dict: `{channel_id: i(t)}` on a 1 ns grid, plus truth block |
| `sim_decoded.root` | C | **identical schema to data** `decoded_root/*_<feu>.root` tree `nt` (`eventId`, `amplitude[n_samp×512]`, `ftst`) so `wft/io.py` reads it unchanged |

---

## 3. S1 — semi-spectral dynamic weighting-potential solver

**Physics.** Extended Ramo–Shockley for resistive elements (Riegler, JINST 11 (2016) P11002, arXiv:1602.07949; Janssens et al., arXiv:2304.01883). To get the signal induced on readout channel n: apply V_w·Θ(t) to channel n's electrode (all its pads), everything else (other pads, mesh) grounded; solve the **quasi-static** relaxation of the potential in the stack where the ESL layer is a thin sheet with patterned surface conductivity. The time-dependent solution Ψ_n(x,y,z,t) is the dynamic weighting potential. Induced current from a charge q at position x_q(t):

```
i_n(t) = -(q/V_w) ∇ψ_p·ẋ_q(t)  -  (q/V_w) ∫₀ᵗ H_d[x_q(t'), t-t']·ẋ_q(t') dt'
H_d(x,t) = -∇ ∂ψ_d(x,t)/∂t ,   ψ_d = Ψ - ψ_p ,  ψ_p = Ψ(t→0⁺)
```

At t=0⁺ the ESL acts as a dielectric (prompt/static solution); as t→∞ it acts as a grounded conductor. By reciprocity, Ψ_n(x0,y0,0⁻ surface, t)/V_w is exactly the Green's function G_n: the charge induced on channel n by a unit point charge *sitting* on the ESL at (x0,y0) since t=0.

**Geometry model W1 (baseline).** Layers in z: grounded mesh plane at z=+g (treat as solid — justified: weave-scale field ripple decays as e^(−2πz/p_weave), negligible at the ESL for p_weave ≪ g; verify in V5) / gas ε=1 thickness g / ESL sheet at z=0 with σ_s(x) = 1/ρ_s on strips, 0 in gaps (periodic in x, period 800 µm, uniform in y) / kapton ε_r=3.5 thickness d_k / segmented pad plane at z=−d_k. (Pillars — Dynamask, 1.3% coverage — are ignored in the weighting solve; they matter only as dead/perturbed spots, checked in validation, and in S2 if included in the unit cell.) Buried trace layers are screened by the pad plane and ignored in W1 (W2 check in V6: expose 100 µm inter-pad gaps and re-run — expected percent-level).

**How many distinct kernels there are — RESOLVED 2026-08-07 (T2b).** The checkerboard plus the 31.2 mm superperiod collapses 1024 channels to an exact, small set, and this is what makes the comb solve affordable:

- **2 distinct Y kernels.** Moving a Y channel by one row shifts its comb by one pad pitch in x and flips which columns it owns; the row's own y position is a pure translation of a y-uniform sheet. Only the column *parity* survives.
- **40 distinct X kernels**, one per column phase — the column's x position relative to the ESL strips does not drop out. **So the 31.2 mm beat lives in the X view.** The Y view averages over all phases within its own comb and sees the beat only through where the *deposit* lands. That sharpens the §9 prediction: look for the 31.2 mm period in X-view sharing/residuals, not in Y.
- An X channel's drive is exactly 1.56 mm periodic in y, so its solver box can *be* 1.56 mm — the periodic images are the rest of the comb, not an artifact. The 40 X solves are therefore cheap; the Y kernel, localised in y, is the expensive one.

**Sum-rule validation (new, and the strongest test in S1).** Driving every pad is a uniform Dirichlet plane, which excites only k=0, where the sheet has no in-plane gradient and never relaxes. So with pitch-sized (tiling) pads the total induced charge must equal S(0)/C(0) = ε_r/d ÷ (ε_r/d + 1/g) = 0.875 at **every** time. Reproduced to 4×10⁻⁷, with the channels partitioning the plane to 2×10⁻¹². It exercises the comb assembly, the Bloch blocks, all 40 X phases and the row sum against a closed form, and it caught three bugs that would each have silently biased every §9 sharing number (row-sum roll direction; X pads on the wrong row parity; hard-pixel pad masks, which made the total swing 0.76→0.57 with resolution — replaced by exact fractional-area coverage).

**First predictions from the comb kernels (charge level, ρ_s 1 MΩ/sq, d_k 75 µm):**
- **X/Y charge balance is exactly 0.500/0.500 at t=0, necessarily** — at t=0 the sheet is a plain dielectric and the checkerboard is 90°-symmetric, so any measured departure from 0.5 (data: 0.49/0.51) is generated by the *time-dependent* part plus electronics, not by geometry.
- ~76 % of a deposit's image charge lands on the channels, the rest on the mesh and on the grounded 100 µm inter-pad copper.
- **Point-charge prompt sharing between adjacent channels is near zero.** In the checkerboard a d=±1 neighbour's nearest pad is 0.78 mm away in x *and* 0.78 mm in y, i.e. 1.10 mm diagonally, while the prompt kernel is localised on the ~0.2 mm scale of the stack. So essentially all of the measured ±1 sharing (§9: c1 = 0.23–0.28) has to come from the avalanche's own transverse size, electron diffusion, and resistive spreading — none of it from prompt induction geometry. That is a strong, falsifiable statement and it is a headline check for T9.

Production point (ρ_s 1 MΩ/sq, d_k 75 µm, nx 3120, ny 1024, 61 log times, drain-free), deposit **exactly on a pad centre**, shares within d = −3..+3:

| deposit | view | prompt share d=−3..+3 | view total | late (10 µs) share | view total |
|---|---|---|---|---|---|
| on-strip (pad owned by X) | X | 0 0 0 **1.000** 0 0 0 | 0.870 | 0 0 0.006 **0.990** 0.004 0 0 | 0.366 |
| | Y | 0 0 0.258 0.485 0.258 0 0 | **0.001** | 0.241 0.003 0.255 0.001 0.255 0.003 0.241 | 0.143 |
| in-gap (pad owned by Y) | X | 0 0 0.243 0.515 0.243 0 0 | **0.001** | 0 0 0.004 **0.991** 0.005 0 0 | 0.365 |
| | Y | 0 0 0 **1.000** 0 0 0 | 0.870 | 0.002 0.323 0.005 0.338 0.005 0.323 0.002 | 0.111 |

Read three things off it. (i) **Prompt charge division is essentially total**: the view owning the pad under the deposit takes 0.870, the other view takes 0.001 — a factor ~870. The checkerboard does not "share" between views at t=0, it assigns. (ii) **The other view is fed only by sheet transport**, rising from 0.001 to 0.11–0.14 by 10 µs, and it arrives spread over ±3 channels rather than concentrated — so cross-view charge is a *late, diffuse* component. (iii) **The owning view's own total falls** (0.870 → 0.366) as charge migrates outside the ±3 window; the drain-free sheet conserves charge globally (the all-channel total is time-independent at 0.665) but not locally.

**ρ_s scan, 2026-08-07 (d_k 75 µm, four points, desktop).** The position-independent quantities are exactly ρ_s-independent, as they must be: the all-channel total (0.6650), the 0.500/0.500 prompt balance, and the sum rule (passes at 4.6–5.1×10⁻⁷ at every point) are all fixed by the k=0 mode and by symmetry, not by transport. The ρ_s dependence is entirely in the *timing*, and it is clean:

⚠️ **A first version of this section quoted a `τ(1/e)` per view and concluded that the measured τ sits "at the top of the scanned range, nearest 5 MΩ/sq", recommending the scan be extended above 5. That was wrong and is withdrawn.** The estimator — time for the d=+1 neighbour to cover 1−1/e of the way from its prompt value to its maximum — is unsound whenever `tau_drain_s=None`, because a charge-conserving sheet never saturates: the "maximum" is just the value at the 10 µs window edge, so the normalisation is set by where the window was cut. It moved τ_X from 42 ns to 1169 ns on a d_k change that can shift real timescales by at most 1.44. See `response/solver/spread.py`.

**The sound measurement.** Charge on the sheet spreads diffusively, so the induced-charge profile obeys σ²(t) = σ_p² + 2Dt, with σ_p the prompt electrostatic width before the sheet conducts at all. Both terms come from moments of the kernel — no normalisation, no saturation, no window choice:

| ρ_s [MΩ/sq] | σ_p (d_k 75) | D fit | D exact | t(σ = one pad pitch), d_k 75 | d_k 50 |
|---|---|---|---|---|---|
| 0.5 | 222.6 µm | 3.712 | 4.235 | **75 ns** | 105 ns |
| 1.0 | 222.4 µm | 1.856 | 2.118 | **151 ns** | 211 ns |
| 2.0 | 221.8 µm | 0.928 | 1.059 | **301 ns** | 422 ns |
| 5.0 | 220.4 µm | 0.372 | 0.424 | **753 ns** | 1055 ns |

Two free validations fall out. σ_p is flat to <1 % across a factor 10 in ρ_s (it is a prompt geometric width, so it must be) while correctly *depending* on d_k (222 µm at 75 µm kapton, 212 µm at 50 µm — thicker kapton lets the image spread further). And D_fit/D_exact = 0.89 ± 0.02, flat to 1.7 % over that same factor 10 — the solver reproduces the 1/ρ_s law without being told it. The ~10 % deficit is expected rather than an error: D_exact uses c′ = C(k→0), while a profile of finite width samples k > 0 where C(k) > c′, hence slower effective diffusion.

**What this says about ρ_s.** Taking the measured τ_X ≈ 230 ns / τ_Y ≈ 410 ns as the time to spread of order one pad pitch puts ρ_s at roughly **1.5–2.7 MΩ/sq** (d_k 75) or 1.1–2.0 (d_k 50) — the **middle of the existing scan**. The scanned range 0.5–5 MΩ/sq already brackets the data comfortably and does **not** need extending. Still indicative rather than a measurement, for the same reason as before: the data number is a fit to shaped waveforms with 180 ns peaking, and the legitimate comparison is T13b/T14.

The X/Y asymmetry cannot be read off this table at all — σ_y here measures spreading *along* the strips, which is the Y-sharing direction only. Getting the §9 X-vs-Y asymmetry requires the x-direction profile, where transport is blocked by the inter-strip gaps; that is a T9 observable.

Caveat, and it matters for reading any of the above against §9: these are **point-charge, charge-level** kernels. A real avalanche has transverse extent and the electrons arrive with diffusion spread, both comparable to the pad pitch, so the "all on one channel" prompt result will be substantially smeared before it reaches a waveform. The τ(1/e) values the run prints (7–42 ns) are likewise *not* comparable to the measured τ_X ≈ 230 ns / τ_Y ≈ 410 ns, which come from fits to shaped waveforms with a 180 ns peaking time. Both comparisons become legitimate only at T9/T11.

**Electrode definition — RESOLVED 2026-08-07 (T2, gerber-derived).** The checkerboard is confirmed, from the 0.1 mm connector stubs (the vias do NOT answer it — the plated through hole rings every layer on all 512×512 pads). A Y channel = one row's pads at (col+row) even; an X channel = one column's pads at (col+row) odd. **A channel is therefore a comb of 256 pads on a 1.56 mm pitch, not a solid line of 512 pads on 0.78 mm.** Extraction and proof: `response/common/channel_map.py` (262030/262144 pads assigned; the 114 stragglers are a 0.04% edge effect). All channel-level kernels must be built with comb drive patterns — the solid-line `pad_pattern`/`strip_pattern_y` in `wpot.py` are per-pad/limiting cases only (task T2b). The x/y sign convention vs `Mx17StripMap.py` connector numbering is still unverified.

**Method.** Expand in lateral Fourier modes. In y (uniform sheet direction) modes decouple: continuous wavenumber k_y (discretize, ~200 log+linear points to k_y·L=π·400). In x the periodic conductivity couples k_x ↔ k_x + m·2π/p_ESL (Bloch): for each (k_y, k_x∈[0, 2π/p)) truncate to |m| ≤ M (start M=12, converge-test). Each dielectric layer gives an algebraic transfer relation per mode; the sheet gives the junction condition

```
n·(J_above − J_below) = −∇_T · [ σ_s(x) ∇_T V ]|_{z=0}   (+ ε0 ∂/∂t displacement terms in the layers)
```

Result: per (k_y, Bloch block) a linear ODE system `dV/dt = A·V + b·Θ(t)` of dimension ~(2M+1); solve by eigendecomposition of A (exact exponentials — no time-stepping error). Assemble Ψ on the output grid at the 60 requested times. Runtime target: minutes per channel type on the laptop. Pure numpy.

**Boundary/drain — AMENDED 2026-08-07.** Production default is **`tau_drain_s=None`** (strict charge conservation on the sheet). The real drain — copper bus at the two y-ends, §1 — is ms-scale, a factor ≥0.999 over any waveform window, and indistinguishable from no drain in every §9 observable; implementing it as a boundary condition buys nothing. The solver keeps the uniform-leak term (`tau_drain_s` explicit) strictly as a nuisance/systematic knob. Do NOT pin it to "the measured τ_g" — that number is not a drain (§1). Superseded: ~~implement A1 as a lumped leak 1/τ_g on the k_y→0 mode and check the tail against the measured 5–7 µs~~.

**Validation (all must pass before S1 output is used):**
- V1: set σ_s uniform (no gaps) → matches Riegler's closed-form uniform-layer solution (JINST 2016 eqns; also the Gaussian limit σ²(t)=2t/(ρ_s c′), c′=ε0(ε_r/d_k + 1/g)).
- V2: single isolated strip, k_y only → matches 1D telegraph solution (Galan arXiv:1110.6640).
- V3: charge conservation: Σ_n G_n(x0,y0,t) + charge remaining on sheet + mesh charge = 1 at all t; each G_n → its t→∞ electrostatic value.
- V4 — REDEFINED 2026-08-07 (the original "drain tail within factor ~2 of τ_g, else revisit A1" compared incommensurable observables and correctly "failed"; that failure is documented in `design/report/s1_validation.txt` and resolved in §1): the drain-free comb-channel kernel, fitted with the rc_line toy model (Gaussian line diffusion × exp(−t/τ_g), ≤1.4 µs window), must yield an apparent τ_g of the measured order for on-strip deposits and a much slower/no decay for gap deposits. Passed at charge level (2.2–3.3 µs on-strip, ∞ in-gap; `response/validation/tau_g_reinterpretation.py`). Final closure at waveform level is T13b.
- V5: mesh-as-plane check: perturbative estimate of weave ripple amplitude at z=0 < 1%.
- V6: W2 (exposed inter-pad gaps) shifts G_n by < few %.

Numbers to expect (sanity): c′ ≈ 5×10⁻⁷ F/m²; sheet diffusivity D = 1/(ρ_s c′) ≈ 2.0 m²/s at 1 MΩ/sq → relaxation of a 130 µm feature in ~8 ns, of one 800 µm pitch in ~0.3 µs.

**Host:** laptop (light). **Scan matrix:** ρ_s {0.5,1,2,5} MΩ/sq × d_k {50,75,125} µm = 12 points (gap fixed at 150 µm); each point = 2 channel types → 24 solver runs, laptop-scale; if slow, desktop.

---

## 4. S2 — realistic mesh field (unit cell)

Purpose: (a) electron transparency ε(E_drift/E_amp) from geometry instead of an assumed constant; (b) funneling map (entry (x,y) in the weave cell → avalanche seed position/spread below); (c) field non-uniformity feeding S3; (d) ion drift endpoints (fraction terminating on mesh vs escaping to drift — sets the ion-tail shape split).

Method: one weave unit cell (period from `shared/MX17ModuleGeometry.hh` wire diameter + fill factor — **resolve the weave pitch from these two numbers and cross-check against the bulk-MM standard 400 lpi; if inconsistent, flag in §12 and use the header**), woven wire geometry (two orthogonal sinusoid-ish wires, standard Garfield++ neBEM wire/primitive representation or an Elmer tetra mesh), periodic lateral BCs, plates: drift cathode far above (apply E_drift), anode plane at ESL surface potential below. Solve electrostatics; export field map.

Tools: **first choice Garfield++ neBEM** (built on all three hosts at the pinned commit, §5a; no meshing needed). Fallback: install `gmsh`+`Elmer` (apt/pip, ~30 min) if neBEM struggles with the woven geometry — note `ComponentElmer2d` gained delayed weighting-potential support in the pinned version. **Host: desktop** (one-off heavy; hours). Deliverable includes a transparency curve figure vs field ratio compared to the generic bulk-MM curve from literature.

**Start from the upstream example, not from scratch:** `Examples/ResistiveMicromegas/` in the pinned Garfield (added 2026-05-21 from the DRD1 GDSimS 2026 tutorial) is a working resistive-MM chain that ships COMSOL maps for a **woven mesh** and a **dynamic weighting potential** with strip electrodes, and drives them with `SetDynamicWeightingPotential` + `CopyWeightingPotential` (one solved map translated onto many electrodes — the answer to "512 channels"), `EnableDelayedSignal`, `AvalancheMicroscopic::GetIons()` → `AvalancheMC` for the ion tail, and a `Shaper`. Its maps are also an independent cross-check for S1's V1–V3.

## 5. S3 — avalanche calibration campaign

Garfield++ `AvalancheMicroscopic` + `MediumMagboltz` in the S2 field map (or uniform-field fallback for first pass): 10³–10⁴ single electrons per (gas, HV) point.

Extract per point into `aval_calib.json`: mean gain ḡ, Polya θ (fit P(g) ∝ (g/ḡ)^θ e^(−(1+θ)g/ḡ)), transverse avalanche spread σ0 at the ESL, longitudinal α(z) profile, normalized ion-induced current shape i_ion(t) in the gap (uniform-field ion mobility from Magboltz tables; **flag: ion mobility is the single softest parameter** — carry ±30% as a systematic), fraction of ions to mesh.

Gases: Ar/iso 95/5 dry AND +1% H2O (tables partly exist in `~/PycharmProjects/nTof_x17/garfield_sim/results/` and on EOS — reuse; the condor workflow in that repo is the template). HV: 480–540 V mesh in 10 V steps (bench operating 490 V, SPS up to 625 V different gas — add Ar/CO2/iso 95/3/2 and Ar/CF4/iso 88/10/2 later).

**Host: lxplus condor** (systematic campaign; reuse `garfield_sim/mm_condor_*` submission machinery). Quick single-point smoke tests: laptop.

**Result, first pass (uniform field, Ar/iC₄H₁₀ 95/5 dry, 150 µm gap), 2026-08-07.** `aval_calib.json`, figure `design/figures/response/s3_avalanche_calib.png`:

| V_mesh | mean gain | Polya θ | σ₀ at the ESL | nev |
|---|---|---|---|---|
| 470 | 23 107 | 1.42 | 34.3 µm | 1600 |
| 480 | 31 387 | 1.54 | 34.0 µm | 1360 |
| 490 | 44 472 | 1.64 | 33.7 µm | 1120 |
| 500 | 60 309 | 1.84 | 33.4 µm | 880 |
| 510 | 81 646 | 1.66 | 33.2 µm | 640 |
| 520 | 112 447 | 1.95 | 32.6 µm | 480 |
| 530 | 156 731 | 1.66 | 32.4 µm | 320 |

Gain rises by ×6.8 over 60 V — an e-folding every **31 V**, which is the normal slope for a 150 µm bulk MM. At the 490 V bench operating point the gain is ≈4.4×10⁴. θ ≈ 1.4–1.9 with no real trend; the non-monotonicity is statistics, since nev falls 1600 → 320 as the avalanche size cap bites at high V, and θ here is a moment estimator. If θ(V) is ever wanted as a physical trend rather than a per-point constant, the high-V points need more events.

**σ₀ ≈ 33 µm is the number to carry forward, and it is a strong constraint on T9.** The avalanche footprint at the ESL is **~24× smaller than the 780 µm pad pitch** and barely moves with voltage. Combined with the T2b finding that point-charge prompt sharing to d=±1 is near zero (§3), this means the avalanche's own size contributes essentially *nothing* to the measured c1 = 0.23–0.28. The sharing has to come from transverse diffusion over the 30 mm drift gap plus resistive spreading — nothing else is left. That is a sharp, testable corner for T9: if the digitizer cannot reach c1 ≈ 0.25 from diffusion + sheet transport alone, the model is missing physics, and per §9 the answer is to find it rather than to widen σ₀.

## 5a. Garfield++ version — PINNED (decided 2026-08-06)

**Pin: garfieldpp master `927e5c21`.** Built from source on all three hosts; all pass upstream `ctest` 22/22.

| Host | ROOT | Garfield install (full path) | source clone |
|---|---|---|---|
| laptop | 6.36.06 | `/home/dylan/garfield/install` | `/home/dylan/garfield` |
| desktop | 6.30.02 (`/home/dylan/Software/root_6_30`) | `/home/dylan/Software/garfield/install` | `/home/dylan/Software/garfield` |
| lxplus | 6.38.00 (LCG_109 view) | `/afs/cern.ch/user/d/dneff/work/garfield_install/lcg109-927e5c21` | `/afs/cern.ch/user/d/dneff/work/git/garfieldpp` |

**lxplus paths in full** (these are the ones to paste into job scripts):

```
# the install
/afs/cern.ch/user/d/dneff/work/garfield_install/lcg109-927e5c21
# the tarball shipped to condor workers (6.7 MB)
/afs/cern.ch/user/d/dneff/work/garfield_install/garfield-927e5c21.tar.gz
# the source clone it was built from, pinned at 927e5c21
/afs/cern.ch/user/d/dneff/work/git/garfieldpp
# the LCG view supplying ROOT/python/compiler (its Garfield is NOT used)
/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc14-opt
```

The install sets `GARFIELD_INSTALL`, `LD_LIBRARY_PATH`, `PYTHONPATH`
(`lib64/python3.13/site-packages`) and `HEED_DATABASE`. Upstream examples,
including `Examples/ResistiveMicromegas` and its COMSOL maps, live under the
source clone.

Single entry point on every host: `source nTof_x17/garfield_sim/setup_garfield.sh`. It is the only file that names a Garfield or LCG path; nothing else in the toolchain hard-codes one. It resolves the paths above automatically (and on a condor worker unpacks the tarball instead), so prefer it over pasting a path.

**Do not use the CVMFS Garfield.** LCG_108 ships `6fb94b35` (2025-07-07, 664 commits behind the pin) and LCG_109 ships `78fe1bd3` (2026-02-02, 281 behind). The APIs this plan names in §7 exist in all of them, so this is not about being blocked — it is that everything aimed at *this* problem landed between March and August 2026: the `ResistiveMicromegas` example (§4), `AvalancheMicroscopic::GetIons()` for the ion component of §7 step 5, the neBEM OpenMP race fix in the SVD inversion (S2 correctness), interface-crossing checks (electrons no longer tunnel through mesh wires — that *is* the S2 transparency observable), the FFT-convolution fix and arbitrary-PSD noise generators (§8), and the regression test suite itself.

**Magboltz needs no separate upgrade** and the existing gas tables stay valid: Magboltz is vendored inside Garfield at version 11.19 (January 2024), and between the LCG_108 Garfield and the pin `Magboltz/magboltz.f` changes only by the fixed-form continuation marker in column 6 (`/` → `&`, 354 lines) plus one missing comma in a `FORMAT` *print* statement. No cross-section or transport change. Garfield's built-in Penning table is likewise unchanged — re-probed with `garfield_sim/probe_penning.py`, every rP in `mm_config.py` reproduces exactly.

Re-run `probe_penning.py` and reconcile it against `mm_config.py` whenever the pin moves.

## 6. Stage A — Geant4 upgrades (AFTER geometry work merges)

Port the p2_geant ClusterTree schema (`p2_geant/docs/OUTPUT_FORMAT.md`) into this repo's `SteppingAction`/`EventData`/`RunAction`:
1. Add **`time`** (globalTime, ns) per ionization cluster. (Blocking for everything downstream.)
2. Add module-local coordinates or store the world→active-area transform in the file header. Active-area frame: origin at active-area center, x/y per strip-map convention (`nTof_x17/common/Mx17StripMap.py`), z=0 at ESL surface. Beware: PCB plates are offset (+15,+15) mm from the active-area axis — use the active-area axis, not the plate center.
3. Optional (cheap, valuable): the p2 provenance block (creator process, origin volume, ancestor).
4. Do NOT model strips/pads/coverlay as Geant4 volumes — material budget is unchanged at the level that matters and the response chain owns that geometry.

Physics settings: keep EM opt4; evaluate PAI model in a gas G4Region (p2 `TOOLCHAIN_NOTES.md` argues default condensed-history straggling is inadequate in thin gas — for our 30 mm drift gap it matters less than for p2's 3 mm, but PAI in the gas region is cheap: turn it on, compare cluster statistics, keep it).
Acceptance: a 10⁴-event muon run whose ClusterTree loads in `response/digitizer` and produces sensible (x,y,z,t,nPrimary) distributions.

## 7. Stage B — digitizer

Per event, per ionization cluster (vectorize over clusters):
1. Electrons: n = round(edep/W) with Fano-factor binomial correction (F≈0.2 Ar); or later Heed re-ionization mode.
2. Drift each electron packet: arrival time t = t_cluster + z/v_d + Gauss(σ_L√z /v_d); transverse Gauss(σ_T√z); attachment survival e^(−z/λ). Parameters interpolated from Magboltz tables (wet gas!). Use the packet approximation (per-cluster, not per-electron) until profiling says otherwise.
3. Mesh transparency ε (S2) — binomial thin.
4. Per surviving electron: gain g ~ Polya(ḡ,θ) (S3); avalanche lands at (x+funnel offset, y) with spread σ0.
5. **Induction — fast path (production):** for an avalanche of total charge Q=g·e at (x0,y0) at time t0: per readout channel n,
   `i_n(t) = Q·[ f_e·δ_fast(t−t0) ⊛ (−∂G̃_n/∂t) + f_ion·i_ion ⊛ ... ]`
   concretely implemented as precomputed **response templates** `R_n(x0 mod 31.2 mm, t)`: the full current on channel n for a standard avalanche at (x0,y0), built once by running the slow path (below) over a grid of source positions and caching. Ion component included (ions drift up; their induction uses Ψ at moving z — this is why templates come from the slow path, not from surface G alone).
6. **Induction — slow path (validation, and template generation):** Garfield++ `Sensor` with `ComponentGrid` loading the S1 Ψ time slices (`LoadWeightingField(file, fmt, t_k, true)` per slice; `SetDelayedSignalTimes`; `EnableDelayedSignal(true)`; `SetWeightingFieldOffset` to place channels), drift/avalanche trajectories parameterized (electron spike + ion line current), NOT microscopic per event. Run on O(100) events per parameter point to certify the fast path (<2% waveform residual target).
7. Sum currents per channel on a 1 ns grid; write `currents` product + truth block (true (x,y), t0, per-channel true charge).

**Host:** laptop for development/small runs; **lxplus condor** for productions (pure python + npz — trivially portable); desktop for medium one-offs.

## 8. Stage C — DREAM electronics

1. Shaper: build the DREAM transfer function **from the manual** (`~/x17/Documents/dream/DREAM_User Manual_prod_v3.pdf`) at the register settings in the run config (`CosmicTb_MX17.cfg`, local copies in `~/x17/cosmic_bench/det_3/*/raw_daq_data/`; peaking-time code `(0xd023>>4)&0xF = 2` → 180 ns class per nTof_x17 notes). Convolve channel currents.
2. Sampling: 60 ns (bench, 32 samp) / 60 ns 64 samp (SPS config) with uniform-random trigger phase; ftst semantics as in data.
3. Gain/ADC: charge→ADC scale from the manual's mV/fC + ADC full scale; leave one global scale factor free-but-recorded (this is the one place absolute calibration enters; it does not affect shapes/sharing).
4. Noise: per-channel Gaussian σ from the det3 pedestal runs (raw fdf/decoded root at `~/x17/cosmic_bench/det3/mx17_det3_saturday_scan_6-27-26/*/raw_daq_data/MX17_pedestals_pedthr_260627_16H35_*` and the standalone 6-22 run) **plus** the common-mode component per 64-ch block (measure covariance from pedestal data; inject correlated noise; the analysis CNS step then removes most of it, as in data). Optional pink/coherent extras only if pedestal PSD demands.
5. Saturation at 3550 ADC (clip, plus reproduce the repeated-constant pathology only if needed later). ZS: port `nTof_x17/mx_july_beam_qa/26_zs_sim_extract.py` (DREAM firmware ZS: `ZsTyp=1`, `ZsChkSmp=4`, N·σ thresholds); RAW mode = no ZS.
6. Output in the exact `decoded_root` schema (§2) so `wft/io.py` reads simulation as if it were data.

**Host:** laptop.

## 9. Validation & closure (blind targets — do not tune to these)

Run the **unmodified** wft chain (`nTof_x17/wft/`) + the SPS-style kernel analyses on simulated det3-bench cosmics and simulated normal-incidence tracks. Compare:

| Observable | Measured value | Source |
|---|---|---|
| Dispersed ±1 share c1 | 0.23–0.28 (gain/gas/drift-invariant) | `sps_beam_test_26/analysis/M70V_FLAT_ANALYSIS.md`, `FLAT_CF4_RUN63.md` |
| ±1 median peak-time shift | +54–61 ns | `RAW_RUN71_REANALYSIS_2026-08-04.md` |
| Charge budget d=0/±1/±2/±3 (area) | 1.00 / 0.71–0.77 / 0.40–0.48 / 0.15–0.18 | same (trim20, clean) |
| Peak-amplitude ratios | 1.00 / 0.16–0.19 / 0.06–0.08 / 0.03 | same |
| Full W_d(t) library, 3 drift fields | npz archive | `staging/run_71/reanalysis_2026-08-04/` (data disk) |
| X vs Y sharing asymmetry | τ 230 vs 410 ns; kY 1.8–2.9 | `mx_june_wft/ANALYSIS_STATE_2026-07-31.md` |
| Apparent τ_g from the rc_line fit (differential Y-vs-X kernel decay, NOT a drain — §1) | 5.3–7.3 µs; predicted two-component (gap charge ~no decay) | `rc_line_step2.py` results; closure = T13b |
| X/Y charge balance | 0.49/0.51 (det3) | `bench_constants.py` |
| Undershoot | −4 to −6% | run_71 reanalysis |
| Angular/position resolution | σ_θ 1.08–1.11°, core σ|r| 0.46 mm | `ANALYSIS_STATE_2026-07-31.md` |
| Prompt diffusion onto ±1 | 0.19–0.21 | M70V analysis (checks steps B.2–B.4 alone) |
| **Predicted, look for in data:** 31.2 mm beat in sharing/residuals | — | this plan §1 |

Procedure: predictions FIRST for the full ρ_s × d_k grid, as a band; then overlay data; identify which scan point matches; only then permit tuning, restricted to the physical parameter set {ρ_s, d_k, ion mobility, absolute gain, ENC scale} — the p2 "don't tune past these" firewall applies: if *shapes* disagree beyond these knobs, the model is missing physics; find it, don't fudge it.

## 10. Compute distribution policy

| Host | Hardware | Use for | Don't use for |
|---|---|---|---|
| **laptop** (this machine) | i7-8550U 4c/8t, 16 GB, GTX 1050 4 GB | development, S1 solver, Stage B/C on ≤10⁴ events, analysis/plots, single-event Garfield checks | anything >1 h wall or >8 GB RAM |
| **desktop** (`ssh desktop`) | Ryzen 7 5800X 8c/16t, 62 GB, RTX 3060 Ti 8 GB. Garfield++ ready at the pin (§10a); ⚠ home disk 20 GB free | one-off heavy: S2 neBEM/Elmer solves, medium Garfield campaigns, big single Geant4 runs, S1 if it outgrows laptop | systematic multi-point campaigns (no batch system); storing bulk output (ship results back to `~/x17/response_sim/`) |
| **lxplus** (`ssh lxplus`) | HTCondor + CVMFS; LCG_109 view for the runtime, our own Garfield (§5a) | ALL systematic campaigns: S3 avalanche grid, gas tables, Stage A productions, Stage B parameter sweeps | interactive iteration |

Existing lxplus workflow to reuse: `nTof_x17/garfield_sim/mm_condor_submit.py` and wrappers (AFS work dir `/afs/cern.ch/user/d/dneff/work/git/...`, EOS for tables). Copy the pattern, don't reinvent.

### 10a. Desktop (`ssh desktop`, `dylan-MS-7C84`) — state as of 2026-08-06

Correcting this plan's original claim that the desktop had "no Geant4/ROOT/Garfield installed yet": ROOT and Garfield were already there. A bare `ssh desktop <cmd>` runs a *non-login* shell that does not source `.bashrc`, so `which root` finds nothing and the box looks emptier than it is — use `ssh desktop 'bash -lc "..."'` when probing.

**Present and working:**

| | |
|---|---|
| OS / toolchain | Ubuntu 22.04.5, gcc 11.4.0, cmake 3.22.1, system python 3.10.12, gfortran |
| ROOT | **6.30.02** at `~/Software/root_6_30`, sourced from `.bashrc` (line 121) |
| Garfield++ | `~/Software/garfield`, **at the pin `927e5c21`** — rebuilt 2026-08-06 (72 s, `make -j12`), `CMAKE_BUILD_TYPE=Release`, installed to `~/Software/garfield/install`. Upstream `ctest`: **22/22 pass**. PyROOT smoke test passes (`import Garfield`, `GetIons`, `SetDynamicWeightingPotential`, neBEM all present) |
| Repos | `~/PycharmProjects/nTof_x17`, `~/CLionProjects/MX17_Geant` |

Use it via `source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh` — it detects the desktop and wires up ROOT + the pinned Garfield. Nothing else needs configuring for S2/S3/Stage-B work.

**Still needed (the remainder of T0):**

1. ~~**Python stack**~~ — **done (2026-08-06).** No conda/mamba, system python 3.10 only, but the analysis stack is present and verified importable: `numpy 2.2.6`, `scipy 1.15.3`, `matplotlib 3.10.8`, `uproot 5.7.5`, `awkward 2.12.0`, `PyYAML 5.4.1`, plus PyROOT 6.30/02. (`uproot`/`awkward` installed via `pip --user`; the rest were already present.) Conda was deliberately *not* installed: PyROOT here is hard-bound to python 3.10 (`libcppyy3_10.so` → system `libpython3.10.so`), so a conda python 3.11+ env would have no PyROOT and could not use the pinned Garfield without a full rebuild. Same reason the note below matters — the pinned Garfield's python module lives under `lib/python3.10/site-packages`, so a *different* python would not see it without re-pointing `PYTHONPATH`.
2. **Geant4** — genuinely absent. Only needed if big single Geant4 runs are actually wanted here; lxplus covers Stage A productions.
3. **Disk** — 20 GB free on `/` (plus 28 GB on `/media/ucla`, 7 GB on `/media/dylan/data`). The garfieldpp clone alone is 639 MB because the `ResistiveMicromegas` COMSOL maps are ~100 MB each. Keep bulk products off this box, per the table above.
4. ~~**Cleanup**~~ — **done (2026-08-06).** `~/Software/garfield/build_lcg110/` and `install_lcg110/` (a `3627927` build; the name was misleading — it used local ROOT 6.30 and `/usr/bin/c++`, not an LCG_110 view) have been deleted, reclaiming 267 MB. **Note for whoever wrote that step:** they were *not* safe to delete as written — `.bashrc` sourced `install_lcg110/share/Garfield/setupGarfield.sh`, so deleting first would have broken every login shell. `.bashrc` now sources `install/` (the `927e5c21` pin) and was fixed before the deletion. A login shell and `setup_garfield.sh` now agree on one Garfield path.

**Constraints to remember:** there is **no CVMFS** on this machine (`/cvmfs` does not exist, the client is not installed) — so no LCG views, and the Garfield build must go against the local ROOT. `sudo` works but **requires a password**, so nothing here can be provisioned non-interactively.

## 11. Task graph

| ID | Task | Depends | Host | Acceptance |
|---|---|---|---|---|
| T0 | Desktop env setup — **Garfield++ part DONE** (2026-08-06: pinned build, ctest 22/22, smoke test passes). Remaining: python stack (numpy/scipy/uproot) and, if wanted there, Geant4 | — | desktop | `python -c "import numpy"`; ~~garfield smoke test~~ ✅ |
| T1 | `response/` package skeleton + `common/` geometry constants (parse/assert vs `MX17ModuleGeometry.hh`) + params YAML schema — **DONE 2026-08-07** | — | laptop | unit tests pass |
| T2 | Pad↔X/Y channel mapping from gerbers — **DONE 2026-08-07** (checkerboard confirmed from connector stubs, §3; `response/common/channel_map.py`) | T1 | laptop | map figure; agrees with `Mx17StripMap.py` channel count |
| **T2b** | **Comb channel kernels** — **DONE 2026-08-07** (`response/solver/kernels.py`). 2 distinct Y kernels + 40 distinct X kernels (§3); validated against the closed-form sum rule to 4e-7 with an exact plane partition; first charge-level predictions in §3. Remaining: the x/y sign convention vs `Mx17StripMap.py` connector numbering is still unverified | T2, T5 | laptop/desktop | comb G_n(t) exported for both channel types ✅; sum rule passes ✅ |
| T3 | **S1 solver core** (uniform sheet first: V1) — **DONE 2026-08-07** (V1 at machine precision) | T1 | laptop | V1 passes to <1% |
| T4 | S1 Bloch patterning (strips) + V2, V3 — **DONE 2026-08-07** (superperiod-commensurate box, no truncation parameter) | T3 | laptop | V2, V3 pass |
| T5 | S1 boundary/drain + full grid export; V4, V5, V6 — solver + V4 (as redefined, charge level) **DONE 2026-08-07**; production grid export NOT started (export-strategy decision pending, several GB/point vs local-window vs on-demand); V5, V6 not run | T4 | laptop | all V pass; 48-run scan produced |
| T6 | S2 mesh unit cell | T0 | desktop | transparency curve; field map export |
| T7 | S3 avalanche campaign — **DONE 2026-08-07**. 56/56 slices merged into `aval_calib.json` (7 voltage points, 470–530 V; result table below). Raw seed JSONs moved AFS → EOS `response_sim/avalanche/raw/`. Still open: verify the 6 Magboltz-table jobs' `.gas` outputs landed (they run through the nTof_x17 `garfield_sim` EOS workflow, NOT the mx17_response tree), and re-run once T6 supplies a real field map — this pass is uniform-field, with S2-dependent fields as explicit nulls carrying `uniform_field` provenance | T6 (or uniform-field first pass without T6) | lxplus | `aval_calib.json` grid ✅ |
| T8 | Stage A schema upgrade — **DONE 2026-08-07** (`time`, `creator`, Meta tree with world→active-area transform) | geometry merge | laptop | §6 acceptance |
| T9 | Stage B fast path (templates from S1 surface G + analytic ion term, first pass) | T5, T7 | laptop | digitizes 10³ events; energy/charge bookkeeping closes |
| T10 | Stage B slow path (Garfield++ ComponentGrid delayed signals) + certify fast path | T5, T9 | desktop | <2% waveform residual fast vs slow |
| T11 | Stage C DREAM (shaper from manual, sampling, gain) | T9 | laptop | shaped single-avalanche pulse; rise/peak vs measured template compared (report, don't tune) |
| T12 | Stage C noise from det3 pedestals + ZS port | T11 | laptop | simulated pedestal PSD/σ matches data pedestals |
| T13 | End-to-end: sim decoded_root through wft unchanged | T8–T12 | laptop | wft runs, produces events.parquet |
| **T13b** | **τ_g closure at waveform level**: run the *unmodified* nTof_x17 `rc_line_step1/2.py` fit machinery on simulated waveforms produced with `tau_drain_s=None`; the apparent τ_g must land in the measured 5.3–7.3 µs and be position-flat along the strip. Also look for the predicted two-component tail (gap-deposit charge barely decays — kept ~0.70 at 1.4 µs deep in a gap). Charge-level version already passes: `response/validation/tau_g_reinterpretation.py` | T13 | laptop | apparent τ_g within the measured band with NO drain in the model |
| **Td** | *(nTof_x17 repo, data-side, can run any time)* Two-component refit of the measured rc_line templates: replace the single exp(−t/τ_g) with a strip+gap pair (amplitude ratio ~ the 550/250 area split as a starting point). A resolved second component is direct evidence for the kernel-shape reading of τ_g | — | laptop | refit result recorded either way |
| T14 | Blind comparison (§9) + report | T13 | laptop/lxplus | prediction-band figures vs data |
| T15 | Iterate: constrained tuning, systematics, feed kernel back to wft as physical model | T14 | — | documented parameter posterior |

~~Parallelizable now (before geometry merge): T0–T7. T3–T5 is the critical path and the highest-skill task.~~ **State 2026-08-07:** geometry merged; T1–T5 (solver), T2, T7 (jobs), T8 done. The critical path is now **T2b (comb kernels) → S1 export decision → T9 (digitizer)**, with T6 (desktop mesh cell) and the T7 collect step parallel to it. Td can run immediately in nTof_x17.

**State 2026-08-07, later:** T2b done and validated (§3). Running: the nominal comb-kernel production point on the **desktop** (ρ_s 1 MΩ/sq, d_k 75 µm, nx 3120 = 10 µm, ny 1024 = 48.75 µm, 61 log times; ~2.5 min per Y kernel), and the T7 `collect.py` merge on **lxplus** — run *there*, not locally, because the 56 raw seed JSONs are 19 GB and only the merged `aval_calib.json` needs to travel. Critical path is now **S1 export decision → T9**. The export decision has hard numbers at last, measured not estimated: **one grid point is 2.3 GB compressed and takes 6 minutes of desktop wall clock** (2 Y kernels at 139 s each, 40 X kernels in 67 s together, at nx 3120 / ny 1024 / nt 61). So the 12-point scan is ~28 GB and ~75 minutes of compute. Compute is a non-issue, and as of 2026-08-07 **storage is no longer a constraint either**: EOS `/eos/experiment/ntof/data/x17/response_sim/s1/` is effectively unlimited (§2). The earlier decision to defer the full grid until T9 was made purely on the desktop's 20 GB, so it no longer applies — **the full 12-point ρ_s × d_k grid is being produced** (~28 GB, ~75 min), with the desktop throttled to at most 3 resident points by the drain. The drain is routed desktop → laptop → EOS, which also leaves the laptop a complete working copy for analysis; the desktop can now push to EOS directly, so future *bulk-only* products (Stage B waveforms, S2 maps) should go straight there and skip the hop. Trimming the Y kernel's y range from ±25 mm to ±12.5 mm remains available if a *resident working set* is ever wanted small, but it is no longer needed for the archive.

## 12. Outstanding questions — REFERENCE ONLY

**These are NOT being sent to anyone and will likely never be answered. They exist to record exactly which inputs are assumed rather than known. Never block on, or wait for, any item here — the listed default is the answer for all purposes.**

**Fab-side unknowns (Saclay/CEA would know):**
1. ESL resistive paste: measured surface resistivity (Ω/sq) and paste type/batch for our modules. *Default: scan 0.5–5 MΩ/sq. The T2b spread measurement (§3) puts the data at roughly 1.5–2.7 MΩ/sq, i.e. comfortably inside that range — **the scan does not need extending**. (An earlier note here said to extend above 5 MΩ/sq; that came from a broken τ estimator and is withdrawn — see §3.)* (2026-08-07: the τ_g reinterpretation in §1 removes the apparent pressure toward absurdly low ρ_s; the scan stands. Print thickness ~10 µm per user, unconfirmed — irrelevant to the thin-sheet model, useful only for sanity-checking ρ_s against paste volume resistivity.)
2. Kapton thickness between pad Cu and ESL (material IS kapton — user 2026-08-06; pillars are Dynamask, a separate fact). *Default: 75 µm, ε_r 3.5, scan 50–125.*
3. Screen-print registration: nominal alignment/tolerance of the 800 µm ESL pattern vs the 780 µm pad pattern (and vs active-area center). *Default: nominal aligned at center; the beat makes absolute phase measurable from data later.*
4. ~~ESL strip termination: how strips connect to the HV/ground bus~~ **RESOLVED (user, 2026-08-07): copper bus strips at both y-ends of the active area, no connection anywhere in between. A1 confirmed as hardware; see §1 for why this never conflicted with the data.**
5. ~~Amplification gap as built: 128 or 150 µm~~ **RESOLVED: 150 µm (user, 2026-08-06).**
6. PCB internal stackup: layer z-spacings (pads→L3→L4). *Default: pads-only model W1; spacing irrelevant in W1.*
7. ~~Confirm pad↔X/Y bussing pattern (checkerboard?)~~ **RESOLVED (T2, 2026-08-07): exact checkerboard, extracted from the connector stubs (not the vias — those ring every pad on both layers). A channel is a 256-pad comb on 1.56 mm pitch. See §3 and `response/common/channel_map.py`.**

**DAQ-side unknowns:**
8. Confirm `CosmicTb_MX17.cfg` on `daq:/mnt/cosmic_data/MX17/dream_config/` is byte-identical to the May copies we have locally (for the 6-27 det3 runs). *Default: trust local copies.*
9. DREAM ADC full-scale and mV/fC at our register settings if not unambiguous from the manual. *Default: manual values + one recorded global scale factor.*

**Internal (resolve by doing):**
10. ~~Mesh weave pitch from fill factor 0.223 + 19 µm wire vs standard 400 lpi — reconcile in T6.~~ **RESOLVED 2026-08-07, and there was never a conflict.** The header does not leave the pitch to be inferred from the fill factor — it states the weave directly (`meshWire_um = 19.0`, `meshOpen_um = 48.0`, so pitch = 67 µm = 379 lpi) and *derives* the fill as the areal-mass scale for a 2d-thick slab, `fill = πd/(4·pitch)` = 0.2227, which is the quoted 0.223. Against the bulk-MM standard 400 lpi / 18 µm wire (pitch 63.5 µm):

    | | header 19/48 | standard 400 lpi / 18 µm | difference |
    |---|---|---|---|
    | d / pitch | 0.2836 | 0.2835 | **+0.04 %** |
    | optical open fraction | 0.5133 | 0.5134 | −0.03 % |
    | areal-mass fill | 0.2227 | 0.2226 | +0.04 % |
    | pitch | 67.0 µm | 63.5 µm | +5.5 % |

    The two are the **same weave geometry up to a uniform 5.5 % scale** — identical d/pitch, hence identical optical transparency, identical areal mass, and an identical *shape* of the electrostatic field pattern. The only physical difference is the absolute scale against the 150 µm gap: pitch/gap 0.447 vs 0.423. So per plan §4, use the header, and carry the 5.5 % pitch as a **systematic** in T6 rather than treating it as an inconsistency to be resolved. Note the header still flags the 19 µm wire as a "P2-like weave, placeholder" — the *scale* is the assumption, the weave family is not.
11. Whether per-electron (vs per-cluster-packet) treatment changes closure observables — profile in T9.

## 13. References

Dixit & Rankin, NIM A 566 (2006) 281 (physics/0605121) — dispersion model · Riegler, JINST 11 (2016) P11002 (arXiv:1602.07949) — resistive-layer weighting theory, THE math reference for S1 · Janssens et al., arXiv:2304.01883 — 2D resistive-strip bulk MM with delayed weighting potentials, the rigor benchmark · Galan et al., arXiv:1110.6640, 1304.2057 — strip telegraph line · Alexopoulos et al., arXiv:2409.19297 — NSW strip spreading solutions · T2K ERAM, arXiv:2303.04481 — data-constrained RC workflow · Garfield++ User Guide 2025.1 §7.2 (delayed signals), `ComponentGrid::LoadWeightingField`, `Sensor::EnableDelayedSignal` · p2_geant `docs/SIM_CAMPAIGN_PLAN.md`, `docs/TESTBEAM_PLAN.md` — architecture and tuning-firewall templates · nTof_x17 `RAW_RUN71_REANALYSIS_2026-08-04.md`, `ANALYSIS_STATE_2026-07-31.md` — validation targets.
