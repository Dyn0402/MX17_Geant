# Penning transfer in gaseous-detector simulation: physics, implementation, and why a fitted r_P is still standard

Notes compiled 2026-10-15. The primary sources read in full were:
- the full text of Şahin, Tapan, Özmutlu, Veenhof, *Penning transfer in argon-based gas mixtures*, JINST 5 (2010) P05002 ("Şahin 2010");
- the full text of H. Schindler's PhD thesis, CERN-THESIS-2012-208 ("Schindler 2012"), which contains the microscopic de-excitation model that is in Garfield++;
- the current Garfield++ master source (`MediumGas.cc`, `MediumMagboltz.cc` and the `Data/Deexcitation/*.txt` files), downloaded from gitlab.cern.ch on 2026-10-15;
- the full text of Andronic et al. 2004 (arXiv:physics/0402044).

For the later Şahin/Kowalski papers (2014 NIM A, 2016, 2017, 2020, 2021 JINST) only the abstracts could be accessed, because IOP and Elsevier blocked the downloads. Their fitted parameters are taken from the Garfield++ source, which implements them with DOIs in comments.

Labels used below: **[fact]** = established or measured; **[model]** = a modelling assumption; **[inference]** = my own reasoning or speculation.

---

## 1. Mechanism of Penning transfer in Ar (and Ne) based detector gases

### Takeaway
"Penning transfer" in detector physics is shorthand for *every* channel by which electron-impact excitation energy of the noble gas ends up as an extra free electron. The channels are:
- collisional Penning ionisation, A* + B → A + B⁺ + e⁻;
- heteronuclear associative ionisation, giving (AB)⁺;
- photo-ionisation of the quencher by resonance VUV photons;
- Hornbeck–Molnar (homonuclear associative) ionisation, Ar** + Ar → Ar₂⁺ + e⁻.

These compete with radiative decay to lower levels, non-ionising quenching (into vibrational, rotational or dissociative channels of B), collisional mixing between Ar levels, and three-body excimer formation (Ar₂*). Each channel scales differently with total pressure p and admixture fraction c, and Şahin et al. use exactly that scaling to tell the channels apart. Which Ar levels can contribute depends on the quencher's ionisation potential:

| Quencher | Ionisation potential | Ar levels that can transfer |
|---|---|---|
| CO₂ | 13.77 eV | only 3p⁵3d / 5s and higher (≥13.85 eV) |
| CH₄ | 12.65 eV | 4p and higher |
| C₂H₆, C₃H₈, iC₄H₁₀ | 11.52, 10.95, 10.67 eV | all levels including the 4s manifold |

### Cited Findings

**Ar level structure and lifetimes**
- **[fact]** The 3p⁵4s manifold (Paschen 1s₅, 1s₄, 1s₃, 1s₂) lies at 11.55, 11.62, 11.72 and 11.83 eV. 1s₅ (³P₂, J=2) and 1s₃ (³P₀, J=0) are metastable, with lifetimes of seconds. The J=1 levels mix. The lower resonant level (1s₄) decays in 8.6 ± 0.4 ns by emitting a VUV photon; an older measurement of 21 ± 2 ns is attributed to radiation trapping. 1s₂ decays in 2.2 ± 0.2 ns. — [Şahin 2010, §2.2](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** The ten 3p⁵4p levels (2p₁₀–2p₁, 12.91–13.48 eV) decay mostly to 4s levels, emitting red or IR light at 697–912 nm, with lifetimes of 21.7–40.5 ns. Their radiative decay to the ground state is forbidden, so photo-ionisation from 4p is suppressed. — [Şahin 2010, §1.2, §2.2](https://doi.org/10.1088/1748-0221/5/05/P05002); [Schindler 2012, ch. 5](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[fact]** The 3p⁵3d levels have a 13.85 eV threshold and are produced in amounts comparable to 4p. Their lifetimes cluster around 50 ns, except the J=1 levels, which live about 3.5 ns. Higher levels live longer: 200–360 ns for 4d, 100–210 ns for 5p, about 75 ns for 6s. — [Şahin 2010, §2.2](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Magboltz has resolved 44 separate Ar excitation levels since version 8.6; before that they were lumped into three terms. Above the individual ns and nd levels up to 8s and 6d, a lumped "HIGH" term at 15.66 eV remains. — [Schindler 2012, §5.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Microscopic transfer mechanisms**
- **[fact/model]** Collisional transfer proceeds by one of two mechanisms:
  - Electron exchange (Hotop–Niehaus), which is Auger-like, needs close approach of the two particles, and works even for radiatively forbidden states.
  - Direct dipole–dipole transfer: a pair of optical transitions with 1/r³ coupling, which can act at longer range.

  Both are cited in [Şahin 2010, §1.2](https://doi.org/10.1088/1748-0221/5/05/P05002).
- **[fact]** Hornbeck–Molnar ionisation in Ar has a threshold of 14.710 ± 0.009 eV, near the 3p⁵5p levels. It is dominated by highly excited, non-metastable states and needs two-body binding. Contributions from 3d and 5s are improbable because the gap of 0.5–1 eV is much larger than kT ≈ 0.04 eV. — [Şahin 2010, §1.2, §2.5](https://doi.org/10.1088/1748-0221/5/05/P05002)

**Excimers**
- **[fact]** In argon at atmospheric pressure, excimer formation Ar* + 2Ar → Ar₂* + Ar is a three-body process. The excimer emission bands are:
  - the first continuum at about 11.3 eV (110 nm), which is pressure-suppressed;
  - the second continuum at 9.8 eV (126 nm), which is dominant;
  - the third continuum at 5.4–6.9 eV.

  The 9.8 eV photons cannot ionise CO₂, CH₄, C₂H₆, C₃H₈ or iC₄H₁₀. Only the first continuum can reach the ionisation potentials of C₃H₈ and iC₄H₁₀. The main effect of excimers is therefore *removing* excited states, not ionising. — [Şahin 2010, §2.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Three-body excimer rate constants used in Garfield++ (from Kolts & Setser 1978) are 1.1×10⁻⁴¹ cm⁶ ns⁻¹ for 1s₅ and 0.83×10⁻⁴¹ cm⁶ ns⁻¹ for 1s₃. Two-body 1s₅→1s₄ and 1s₃→1s₄ mixing rates are 2.1×10⁻²⁴ and 5.3×10⁻²⁴ cm³ ns⁻¹; Tachibana 1986 gives similar values. — [Garfield++ `Data/Deexcitation/RateConstants_Ar_Ar.txt`](https://gitlab.cern.ch/garfield/garfieldpp/-/tree/master/Data/Deexcitation); [Schindler 2012, §5.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Photo-absorption and the Şahin dissipation model**
- **[fact]** Argon is transparent below about 10.5 eV. It has discrete absorption lines between 11 eV and its ionisation potential, peaking at 580 Mb at the upper 4s J=1 line, and a 25–40 Mb continuum above the ionisation potential. CO₂ absorbs about 30 Mb at the 4s lines and about 60 Mb at the 3d level, so decay photons are absorbed inside the avalanche. CO₂ is transparent below about 10 eV, so excimer photons escape. Methane and heavier alkanes have higher absorption; essentially no isobutane photo-absorption data were found. — [Şahin 2010, §2.6, §4.2](https://doi.org/10.1088/1748-0221/5/05/P05002); [Schindler 2012, §5.4](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[model]** Şahin et al. follow an excited atom A* through time steps, with competing per-step loss terms:

  | Process | Term per step |
  |---|---|
  | ionising collision with B | p·c·f_B⁺·dt/τ_A*B |
  | homonuclear associative ionisation | p(1−c)·f_A⁺·dt/τ_A*A |
  | excimer formation | p²(1−c)²·f_A2·dt/τ_AAA |
  | non-ionising loss on B | p·c·f_B̄·dt/τ_A*B |
  | non-ionising loss on A | p(1−c)·f_Ā·dt/τ_A*A |
  | radiative decay | dt/τ_A* |

  This gives r(p,c) = k·τ_P, with
  - k = pc·f_B⁺/τ_A*B + p(1−c)·f_A⁺/τ_A*A + p²(1−c)²·f_A2·f_exc/τ_AAA + f_rad/τ_A*
  - 1/τ_P = pc·(f_B⁺ + f_B̄)/τ_A*B + p(1−c)·(f_A⁺ + f_Ā)/τ_A*A + p²(1−c)²·f_A2/τ_AAA + 1/τ_A*
  - f_rad = f_γ·η_B(γ), where η_B is the photo-ionisation yield of the quencher. f_rad is independent of p and c as long as the photon is absorbed in the gas. The approximation fails as c → 0, when photons reach the cathode (exp(−d·N·c·σ_pa) must be ≪ 1; for example c ≫ 0.3% for d = 0.5 cm and 25 Mb).

  — [Şahin 2010, §1.3, eqs. 1.10–1.15](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[model]** Dropping the excimer and Ar*–Ar terms reduces this to the rational form r(c) = (a1·c + a3)/(c + a2) used for most mixtures. Here:
  - a1 = f_B⁺/(f_B⁺ + f_B̄) is the asymptotic collisional branching ratio at high c;
  - a2 is set by the ratio of collision time to lifetime (τ_A*B/τ_A*, divided by p);
  - a3/a2 carries the radiative (photo-ionisation) and Hornbeck–Molnar contributions.

  The p² terms (a5, a6 in the Garfield++ Ar/CO₂ and Ne/CO₂ formulae) represent excimer and three-body physics. — [Şahin 2010, eqs. 4.x](https://doi.org/10.1088/1748-0221/5/05/P05002); [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)

**Collision-time scale and its uncertainty**
- **[fact]** Typical collision times at 1 atm with c = 1:

  | Pair | τ_A*B |
  |---|---|
  | Ar 3d – CO₂ | 77 ps |
  | Ar – C₂H₂ | 150 ps |
  | Ar – C₃H₈ | 130 ps |
  | Ar 4p – Xe | 130 ps |
  | Ar 3d – Xe | 88 ps |

  These are much shorter than the 10–50 ns radiative lifetimes, so at 10% admixture an excited Ar atom suffers of order 10–100 quencher collisions before it decays radiatively. The collision diameters are "the most problematic ingredient": covalent, kinetic, Deutsch–Märk and Lennard-Jones definitions differ by factors of 2–3. Excited Ar is 2.2–6.3 times larger than ground-state Ar. — [Şahin 2010, §2.1, Table 1, §4](https://doi.org/10.1088/1748-0221/5/05/P05002)

**Pressure and concentration dependence**
- **[fact]** In Ar–CH₄ (Ye et al. data, p up to 21–28 atm) the transfer probability clearly rises with pressure. Şahin et al. attribute this to Hornbeck–Molnar ionisation at high pressure and high CH₄ fraction. In Ar–C₂H₂, r is roughly constant in c because radiation trapping lengthens the effective lifetime of the 4s states. The resonant 4s states average about 5.4 ns naturally, which cannot explain the observed rates; trapped radiation is stored in the nearby metastables. — [Şahin 2010, §4.3–4.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Conclusion of Şahin 2010: the main gain-producing processes besides direct ionisation are photo-ionisation and collisional Ar*–quencher transfer. There are hints of Hornbeck–Molnar ionisation. "Loss of excited states through excimer formation is not observed" in the mixtures and pressures studied. — [Şahin 2010, §5](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** In Ne–CO₂ and Ne–N₂, transfer rates *decrease* at high molecular-gas concentration. The author takes this as evidence of an unusual transfer mechanism for excited Ne. — [Şahin, JINST 16 (2021) P03026, abstract](https://iopscience.iop.org/article/10.1088/1748-0221/16/03/P03026)
- **[fact]** Adding N₂ to Ar–CO₂ (85:15) lowers the fitted Penning fraction from 37% with no N₂ to 20% at 10% N₂ and 8% at 20% N₂ (Magboltz 2, single tuned parameter). N₂ quenches Ar* without ionising, since its highest level in that calculation is 13.0 eV. — [Andronic et al., NIM A 523 (2004) 302, arXiv:physics/0402044](https://arxiv.org/abs/physics/0402044)

### Inferences
- **[inference]** For the most common mixtures, Ar/CO₂ and Ar/iC₄H₁₀, the physics is quite different:
  - In Ar/CO₂ only 3d, 5s and higher levels count. These are a modest fraction of excitations, which explains why r_P is about 0.4–0.6 but the gain boost is "moderate".
  - In Ar/iC₄H₁₀ the large 4s population, including the metastables, is eligible, and radiation trapping and metastable kinetics matter most.
- **[inference]** Wet gases and ternary mixtures add non-ionising quenching channels (H₂O, N₂) that lower r_P, as the Andronic N₂ result shows. No dedicated r_P fit for water-contaminated Ar/iC₄H₁₀ or Ar/CO₂ was found.

### Gaps
- No full-text access to the 2014, 2017 and 2021 papers, so their mechanism attributions beyond the abstracts are not given here.
- No source found quantifying Holstein escape factors for Ar 1s₄/1s₂ in detector-relevant geometries (μm–mm cells). Garfield++ implements trapping explicitly by photon transport with Voigt line profiles instead (section 2).
- Branching ratios of metastable Ar quenched by C₂H₂: Şahin 2010 cites 0.61 ± 0.14 and 0.74 ± 0.07. A search snippet attributed to an older Ar(³P₂) study said ionisation channels for NO and C₂H₂ were weak, with efficiencies of about 0.2 and 0.1. I could not verify that snippet; treat it as an unverified conflict.

---

## 2. What Magboltz, Degrad and Garfield++ implement today

### Takeaway
Magboltz does **not** apply Penning transfer. It outputs level-resolved excitation and ionisation collision frequencies, with 44 Ar levels resolved since v8.6. Garfield++ then offers three options:
1. **Macroscopic (gas table):** the Townsend coefficient is rescaled using the excitation frequencies of eligible levels times r_P.
2. **Microscopic "simplified Penning":** at each eligible excitation, with probability r_P an extra electron is created on the spot, or within a sphere of radius λ.
3. **Detailed de-excitation model** (`EnableDeexcitation`, Ar only, plus optional `EnableRadiationTrapping`): explicit level kinetics with radiative cascades, photon transport with Voigt lines, quenching rate constants and Ar–Ar kinetics. It is mutually exclusive with the simple Penning mode.

A default parameter table, `EnablePenningTransfer()` with no arguments, covers only *binary* mixtures.

### Cited Findings

**API (current master, downloaded 2026-10-15)**
- **[fact]** The calls are:
  - `EnablePenningTransfer()`, which uses the built-in parameterisations and returns an error unless the mixture has exactly 2 components;
  - `EnablePenningTransfer(r, lambda)`, global;
  - `EnablePenningTransfer(r, lambda, gasname)`, per component;
  - `DisablePenningTransfer([gas])`, `GetPenningTransfer(...)`, `GetNumberOfPenningTransfers()`.

  MediumMagboltz also stores r and λ per excitation *level* (`m_rPenning[nMaxLevels]`) and exposes `GetPenningTransfer(i, r, lambda)` per level. — [Garfield++ MediumGas.hh / MediumMagboltz.hh](https://gitlab.cern.ch/garfield/garfieldpp/-/tree/master/Include/Garfield)

**Built-in default table**

Formulas are copied from `MediumGas::EnablePenningTransfer()` (c = admixture fraction, p in atm). I evaluated the r_P values with the coded constants.

| Mixture | Source in code | Form | r_P examples |
|---|---|---|---|
| Ar/CO₂, p = 1 atm | NIM A 768 (2014) 104 | (0.6643c + 0.0028)/(c + 0.0518) | 3%: 0.28; 5%: 0.35; 7%: 0.41; 10%: 0.46; 20%: 0.54; 30%: 0.57 |
| Ar/CO₂, p ≠ 1 | JINST 12 (2017) C01035 | (a5·(p(1−c))² + a1·pc + a4·c + a3)/(a6·(p(1−c))² + pc + a2), with a1 = 0.627898, a2 = 0.041394, a3 = 0.004716, a4 = 0.001562, a5 = 0.002422, a6 = 0.027115 | 10%: 0.36 at 0.4 atm, 0.43 at 1 atm, 0.43 at 1.8 atm |
| Ar/CH₄ | JINST 5 (2010) P05002 | (b4·pc + b1·p(1−c) + b2·c + b5)/(b4·pc + p(1−c) + b3), with b1 = 0.1956, b2 = 16.38, b3 = 22.12, b4 = 3.842, b5 = 2.992 | 10%: 0.22 at 1 atm, 0.26 at 5 atm, 0.33 at 20 atm |
| Ar/C₂H₆ | 2010 | constant 0.31 (10%, 1 atm only) | |
| Ar/C₃H₈ | 2010 | 0.4536c/(c + 0.0035), 1 atm | 5%: 0.42; 10%: 0.44 |
| Ar/iC₄H₁₀ | 2010 | **constant 0.40** (only 10%, 1 atm measured; prints a warning otherwise) | |
| Ar/C₂H₂ | 2010 | constant 0.72 (cylindrical-chamber value) | |
| Ar/Xe | 2010 | (1.248c + 0.008)/(c + 0.039) | exceeds 1 above about 15% Xe |
| Ne/CO₂ | JINST 16 (2021) P03026 | 7-parameter rational form in pc, c², (p(1−c))² | 10%: 0.58 at 1 atm; 30%: 0.52 at 1 atm |
| Ne/N₂ | JINST 16 (2021) P03026 | similar form | |
| Xe/TMA | JINST 13 (2018) P10032 | (0.2472p + 0.0414)/(p + 0.2372), 5% TMA only | 0.23 at 1 bar |

The values are applied with λ = 0 and to the noble-gas component only. — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc); [Garfield++ Penning tutorial](https://garfieldpp.docs.cern.ch/tutorials/penning)
- **[fact]** The tutorial states that the "Penning transfer probabilities are not taken into account directly when generating a gas file using Magboltz". It also says r_P "can be parameterised as function of the concentration c of the admixture and the pressure p" when sufficient data exist. — [Garfield++ Penning tutorial](https://garfieldpp.docs.cern.ch/tutorials/penning)

**Gas tables (macroscopic tracking)**
- **[fact]** `MediumGas::AdjustTownsendCoefficient()` (comment "GASSPT", inherited from Garfield-9) computes the corrected coefficient at every (E, B, θ) table point as α' = α₀·(ν_ion + Σ_i r_i·ν_exc,i)/ν_ion. The unadjusted α₀ is kept separately (`m_eAlp0`).

  The global r is applied only to excitation levels whose energy exceeds the *lowest* ionisation potential in the mixture. For Ar/CO₂ that means Ar 3d and up, *and* any CO₂ excitations above 13.77 eV. This is the same Townsend relation Şahin et al. use to fit gain curves (their eq. 3.1). — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc); [Şahin 2010, §3.1](https://doi.org/10.1088/1748-0221/5/05/P05002)

**Simple Penning in microscopic tracking**
- **[fact]** In `MediumMagboltz::ElectronCollision`, if the level energy exceeds `m_minIonPot` and a uniform random number is below `m_rPenning[level]`, a secondary electron is created. Its properties are:
  - energy = E_exc − I_min, which is very small; for example about 0.1–0.5 eV for Ar 3d to CO₂;
  - position uniform in a sphere of radius λ (`λ·cbrt(u)`), or at the collision point if λ = 0;
  - **no time delay**: it is emitted at the excitation time.

  The source comment calls this a "Simplified treatment of Penning ionisation". — [Garfield++ MediumMagboltz.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc)

**Detailed de-excitation model**
- **[fact]** `EnableDeexcitation()` switches simple Penning off, and vice versa ("Penning transfer will be switched off" / "Deexcitation handling will be switched off"). It is implemented for **Ar only** ("currently only Ar") and reads data files from `$GARFIELD_INSTALL/share/Garfield/Data/Deexcitation/`:
  - `OscillatorStrengths_Ar.txt`;
  - `TransitionRates_Ar.txt` (NIST ASD plus Zatsarinny & Bartschat 2006; e.g. 1s₄→ground 0.119 ns⁻¹, 1s₂→ground 0.51 ns⁻¹);
  - `RateConstants_Ar_Ar.txt` (Kolts & Setser excimer and mixing rates; 4p→4s transfer from Sadeghi 2001 and Chang & Setser 1978);
  - `RateConstants_Ar_<Q>.txt` for Q = CO₂, CH₄, C₂H₆, iC₄H₁₀, C₂H₂, CF₄.

  The radiative rates above correspond to 8.4 ns and 2.0 ns, consistent with section 1. — [Garfield++ MediumMagboltz.cc `ComputeDeexcitationTable`](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc); [Data/Deexcitation](https://gitlab.cern.ch/garfield/garfieldpp/-/tree/master/Data/Deexcitation)
- **[fact]** Total quenching rate constants in those files, in cm³ ns⁻¹:

  | Quencher | Source | 1s₅ | 1s₄ | 1s₃ | 1s₂ | 4p |
  |---|---|---|---|---|---|---|
  | CO₂ | Velazco et al. 1978; Sadeghi 2001 | 5.3e-19 | 5.0e-19 | 5.9e-19 | 7.4e-19 | 6.1–6.6e-19 |
  | iC₄H₁₀ | Piper et al. 1973; Velazco 1978 (1s₃); ethane-scaled 4p | 7.1e-19 | 6.1e-19 | 8.5e-19 | 11.0e-19 | ≈5.7e-19 average |
  | CH₄ | Chen & Setser 1991; Velazco | 4.55e-19 | 4.5e-19 | 5.30e-19 | 5.7e-19 | 3.4–9.3e-19 |

  1×10⁻¹⁹ cm³ ns⁻¹ = 1×10⁻¹⁰ cm³ s⁻¹, so these are the classic metastable-quenching values of a few ×10⁻¹⁰ cm³ s⁻¹. — [Garfield++ RateConstants_Ar_*.txt](https://gitlab.cern.ch/garfield/garfieldpp/-/tree/master/Data/Deexcitation)
- **[fact]** The quoted uncertainty of these rate constants is about 10–20%. For 4p, only 2p₈, 2p₆, 2p₅ and 2p₁ are measured (Sadeghi et al.); the other six 4p levels use the average. — [Schindler 2012, §5.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[model]** Where measurements are missing, the code estimates parameters:
  - **Ionising branching ratio**: k_{M+}/k_M ≈ η(ε)^0.4, where η is the photo-ionisation yield of the quencher at the Ar level energy (`pIon = pow(eta, 0.4)`), unless a third column in the data file overrides it.
  - **Higher resonant levels**: the rate constant comes from the Watanabe–Katsuura formula, k = 2.591e-19·(u_A·u_Q)^0.4·(T/μ)^0.3.
  - **Non-resonant 3d and 5s levels**: a hard-sphere estimate with r(Ar 3d) = 436 pm and r(Ar 5s) = 635 pm, and quencher radii of 165 pm (CO₂) to 250 pm (iC₄H₁₀).
  - **Ar 3d/5s → 4p collisional transfer**: hard-coded k = 1e-20 cm³ ns⁻¹, commented "Parameter to be tuned (order of magnitude guess)".
  - **Hornbeck–Molnar**: kHM = 2e-18 cm³ ns⁻¹, commented "This value seems high, to be checked!".
  - **The "Ar_HIGH" lumped level**: split evenly among 5 levels.

  — [Garfield++ MediumMagboltz.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc); [Schindler 2012, eq. 5.3–5.4](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[fact]** Radiation trapping (`EnableRadiationTrapping`) is done by explicit photon transport. Each resonance photon gets a Voigt line shape: Doppler width plus resonance broadening with prefactor k = 1.92π/√3 (Ali & Griem), and the wings are cut off. The photon is re-absorbed by Ar (restarting the cascade) or by the quencher, which ionises with probability η. Natural broadening (about 10⁻⁷ eV) is neglected. — [Schindler 2012, §5.1–5.3](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf); [Garfield++ `LineCollisionRate`](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc)

**Magboltz versions**
- **[fact]** Garfield++ documents MediumMagboltz as an interface to Magboltz 11. The 11.12 update, relayed by H. Schindler from S. Biagi's notes, removed an artificial ionisation cross-section used to speed up attaching mixtures, reverted the Townsend coefficient calculation to the 10.n method, and fixed a magnetic-field bug. — [Garfield++ MR "Magboltz 11.12"](https://gitlab.cern.ch:8443/garfield/garfieldpp/-/merge_requests/222); [MediumMagboltz doxygen](https://garfieldpp.docs.cern.ch/doxygen/classGarfield_1_1MediumMagboltz.html)

### Inferences
- **[inference]** Degrad: I found no source documenting a Penning option inside Degrad, and the Garfield++ TrackDegrade docs expose only clusters (electrons and excitations) plus a switch for the detailed de-excitation cascade. In practice, Penning in primary ionisation is handled downstream, through Garfield's W and Fano adjustment, or by the user converting excitations to electrons with r_P. Verify in the Degrad Fortran or TrackDegrade.cc before relying on this.
- **[inference]** Because of the λ = 0, no-delay simplification, the simple model gets *gain* right by construction (if r_P is correct) but puts Penning electrons in the wrong place and at the wrong time. τ_P is 0.1 ns to tens of ns, and the spatial spread from photon-mediated transfer reaches tens of μm. That matters for Micromegas and GEM timing (e.g. PICOSEC) and for avalanche transverse spread.
- **[inference]** Several defaults are just stand-ins: the default r_P for Ar/iC₄H₁₀ is a single 1 atm, 10% point, and nothing exists for ternary mixtures. Anyone simulating, for example, Ar/CF₄/iC₄H₁₀ or Ar/iC₄H₁₀ at 95:5 with the default gets a fallback warning, or an error for ternaries.

### Gaps
- I could not access the current Garfield++ User Guide PDF (404 at the old URL); the source code was used instead.
- No Biagi Magboltz 12/13 changelog was found. Whether newer Magboltz versions changed the Ar level splitting or the "HIGH" level is unknown.
- I did not confirm whether Degrad itself offers a Penning flag.

---

## 3. Why r_P is a single fitted parameter, and what blocks an ab-initio treatment

### Takeaway
The obstacle is mainly **missing and uncertain data plus model degeneracy**, not unsolvable quantum mechanics in principle. Schindler (2012) built the full kinetic model (now in Garfield++). With all available literature data, he still could not predict r_P: "Literature data on collision-induced deexcitation transfer are less abundant. A 'first-principle' simulation is, therefore, not possible." He tuned about 15 unknown parameters by simulated annealing to reproduce the fitted r_P values, and notes the solution "is clearly not a unique solution".

A genuinely hard quantum-chemistry problem does sit underneath: the Penning-ionisation branching ratio for each Ar* level with each polyatomic molecule depends on autoionising complex potential surfaces and competing dissociative channels. In practice nobody computes it; it is measured, and measurements exist only for the metastable 4s levels and a few 4p levels.

### Cited Findings

**What a microscopic model needs (Schindler 2012, §5.5)**
- **[fact]** The required inputs are:
  1. optical transition rates for the noble gas and line shapes of the resonance lines;
  2. the photo-absorption cross-section and ionisation yield of the quencher;
  3. rate constants *and ionisation yields* for collisional de-activation by the quencher;
  4. Ar–Ar collisional transition rates.

  Items 1 and 2 are available, except no isobutane photo-absorption data were found. Item 3 is sparse. — [Schindler 2012](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Which inputs are measured and which are not**
- **[fact]** Measured: total quenching rate constants for 4s levels (Velazco 1978; Chen & Setser 1991; Piper 1973) and for four of the ten 4p levels (Sadeghi 2001). For 3d, 5s and higher levels, "no measured rate constants could be found".
- **[fact]** Ionising branching fraction: only Ar metastables with C₂H₂ were found, at 0.61 ± 0.14 and 0.74 ± 0.07. For all other levels and quenchers the ratio is *estimated* as η^(2/5). — [Schindler 2012, §5.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**The tuned parameter set (Schindler 2012, Table 5.3)**
- **[fact]** The free parameters were:
  - k(3d,5s→4p) = 10×10⁻²⁰ cm³ ns⁻¹ and k(high→4p) = 0.8×10⁻²⁰ cm³ ns⁻¹;
  - 3d and 5s quenching by CO₂ at 12 and 22×10⁻¹⁹ cm³ ns⁻¹, ionisation yield 0.8;
  - by CH₄ at 1.8 and 3.1×10⁻¹⁹ cm³ ns⁻¹, yield 0.2;
  - by C₂H₆, yield 0.8;
  - 4p ionisation yields of 0.2 (CH₄) and 0.3 (C₂H₆).

  The targets were Şahin's r values for Ar with 5–20% CO₂, 2–10% CH₄ and 10% C₂H₆ at 1 atm, simulated in a 50 µm parallel-plate gap at 70 kV/cm. For Ar/5% C₂H₂ the model gave r = 0.68, against 0.79 (cylindrical) and 0.73 (parallel-plate) from gain fits. — [Schindler 2012, §5.2.4](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[fact]** The 3d→4p transfer constant and the Hornbeck–Molnar rate remain flagged in the live source as "order of magnitude guess" and "seems high, to be checked!". — [Garfield++ MediumMagboltz.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc)

**Why gain fits cannot cleanly constrain the kinetics**
- **[fact]** Collision diameters are ambiguous by factors of 2–3. This biases the extracted f_B⁺ but not the net transfer or the time scale. Fits with an excimer term were unstable, so f_A2 was set to 0. The Ar/Xe four-parameter fit gave "unphysical" values (negative a2) because of "too few data points and too large uncertainties for the number of parameters". — [Şahin 2010, §2.1, §4.1, §4.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Separating 3d from 4p is impossible: "As for the xenon mixtures, 3p⁵4p and 3p⁵3d are produced in almost equal amounts and with a similar electric field dependence. Since we can not separate them in the fits, we assume these states have the same transfer probability." — [Şahin 2010, §4.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Schindler's motivation for the microscopic model: deriving r(p,c) from gain curves "is possible but requires a large set of gain measurements"; in a microscopic simulation "the pressure and concentration dependence of Penning ionization should emerge naturally". — [Schindler 2012, §5.1.1](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

### Inferences
- **[inference]** Ranking the obstacles:
  1. Missing level-resolved quenching rates and ionising branching ratios for the 4p (6 of 10), 3d, 5s and higher levels with each quencher. This matters most for Ar/CO₂, where only 3d and up are eligible and *none* of the relevant rates are measured.
  2. Missing photo-absorption data for iC₄H₁₀ and a few other quenchers.
  3. Uncertain level-resolved Ar excitation cross-sections. The Penning gain is ∝ Σ r_i·ν_exc,i, so a 20% error in the excitation cross-section maps directly into a 20% error in the required r.
  4. Radiation trapping and escape, which depend on geometry.
  5. Ab-initio Penning-ionisation dynamics: complex (autoionising) potential surfaces for Ar* plus a polyatomic molecule, as in the Hotop–Niehaus exchange picture. This is research-grade quantum chemistry and is effectively never attempted for detector gases.

  Items 1–4 are data problems; item 5 is the "hard QM", and it is why item 1 cannot just be computed.
- **[inference]** A single r_P works operationally because τ_P (sub-ns to a few ns for 5–20% quencher at 1 atm) is short compared with avalanche development in wires but not in Micromegas. It also works because gain depends only on Σ r_i·ν_exc,i, which is one number per field. A gain curve therefore fixes essentially one effective combination per (p, c) and can never constrain level-resolved physics.

### Gaps
- No quantum-chemistry calculation of Ar*–CO₂ or Ar*–iC₄H₁₀ Penning branching ratios was found (searched indirectly). Absence of evidence: I did not search the chemistry literature (Yencha, Siska, Hotop reviews) in depth.
- No post-2012 update of the de-excitation data compilation was found.

---

## 4. Şahin / Kowalski / Veenhof fits: methods, values, dependences and their own caveats

### Takeaway
The method is always the same. Take precise gain curves G(V) from a cylindrical proportional counter or a parallel plate, compute α₀, ν_ion and ν_exc,i with Magboltz, and fit with:
- the Penning-corrected Townsend coefficient (eq. 3.1);
- a free gain-calibration factor g;
- a photon-feedback term G/(1 − βG).

Then repeat over (p, c) and fit r(p,c) with a rational model (§1). Asymptotic collisional branching ratios a1 come out around 0.4–0.7. At atmospheric pressure, r_P is about 0.2 for Ar/CH₄ and 0.3–0.6 for Ar/CO₂ and Ar/heavier alkanes. Ar/Xe, Ar/C₂H₂ and Ne mixtures are higher (0.6–1).

### Cited Findings

**Method (Şahin 2010, §3)**
- **[fact]** The fit model is dα'/dr in the form α'(E) = α(E)·[Σν_ion + Σ r_i·ν_exc,i]/Σν_ion, with gain scaling G := g·G and photon feedback G := G/(1 − βG).
  - The gain scale g is *correlated* with r_i; the two are disentangled by their different field dependence. A common g is used across voltages and, where possible, across mixtures.
  - In some literature data sets the absolute gain is uncertain by up to a factor of 10, and applied g values range from 0.23 to 1.25.
  - β is left free and is "not strongly correlated" with r_i.
  - A temperature change of 5 K can change the gain by up to 23%.

  — [Şahin 2010, §3.1–3.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** The authors note that Magboltz ionisation cross-sections for some gases (Xe in particular) were conjectured to be too high, and a uniform scaling of α was rejected as unsatisfactory. — [Şahin 2010, §3.4](https://doi.org/10.1088/1748-0221/5/05/P05002)

**Results, 2010 (literature gain curves)**
- **[fact]** At 1 atm with 10% admixture: CH₄ r = 0.212 ± 0.002, C₂H₆ 0.31 ± 0.01, C₃H₈ 0.43 ± 0.01, iC₄H₁₀ 0.40 ± 0.01. The increase with molecular weight is explained by larger UV absorption (higher f_rad) and possibly a larger f_B⁺. iC₄H₁₀ falls off the trend; the authors suggest because it is compact. — [Şahin 2010, Table 2, §4.5](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Ar/CO₂ (Kowalski data, 5–20% CO₂ at 1070 hPa): asymptote a1 = 0.620 ± 0.057, with a lower limit a2 > 0.042. The fit parameters are a1 = 0.620 ± 0.057, a2 = 0.008 ± 0.034, a3 = −0.011 ± 0.009, and the parametrisation diverges below about 3–5% CO₂. There is no evidence of excimer loss or photon feedback despite G > 10⁵. — [Şahin 2010, §4.2](https://doi.org/10.1088/1748-0221/5/05/P05002); [Schindler 2012, §4.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[fact]** Ar/CH₄: a1 = 0.43 ± 0.05 for 10%, 0.26 ± 0.04 for 2% and 0.27 ± 0.04 for 5%. The transfer probability increases with pressure up to about 20 atm, attributed partly to Hornbeck–Molnar ionisation. The 10% CH₄, 1 atm compilation shows excess scatter (S = 2.6). — [Şahin 2010, §4.4](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Ar/C₂H₂: a1 = 0.72 ± 0.02 (cylindrical) versus 0.67 ± 0.02 (parallel plate), with ±0.1 systematic. This is an explicit geometry difference. — [Şahin 2010, §4.3](https://doi.org/10.1088/1748-0221/5/05/P05002)
- **[fact]** Ar/Xe: a1 = 1.248 ± 0.086, so r exceeds 1 above 15% Xe. The authors discuss possible causes: photo-ionisation farther from the anode, ionisation cross-sections being too small (Xe cross-sections were being updated at the time), and transfer from Ar 4s via Xe associative ionisation. — [Şahin 2010, §4.1](https://doi.org/10.1088/1748-0221/5/05/P05002)

**Later papers**
- **[fact]** 2014, NIM A 768, 104: high-precision gains in Ar/CO₂, 1–50% CO₂, at 400/800/1200/1800 hPa, gains from 1 up to 5×10⁵, extracting both r_P and photon-feedback β. It yields the Garfield++ 1 atm formula r = (0.6643c + 0.0028)/(c + 0.0518). — [Uludağ record](https://acikerisim.uludag.edu.tr/items/270e1f80-f61d-4d35-9e12-d4439abb747b/full); [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)
- **[fact]** 2017, JINST 12 C01035 (Şahin & Kowalski): a "comprehensive model" of Ar/CO₂ at 0.4, 0.8, 1.2 and 1.8 atm, identifying the contributions of several transfer processes from the p and c dependence. It yields the six-parameter (p,c) formula above. In that formula r_P rises with pressure at low p and saturates by about 1 atm (10% CO₂: 0.36 at 0.4 atm, 0.43 at 1–1.8 atm). — [IOP abstract](https://iopscience.iop.org/article/10.1088/1748-0221/12/01/C01035); [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)
- **[fact]** Ne mixtures:
  - 2016 (JINST 11 P01003): Ne/CO₂, 6–60% CO₂, several pressures, with Penning rates and secondary processes.
  - 2020 (JINST 15, issue 12): Ne/N₂, 2.5–20% N₂, 0.05–1.8 atm, fitting r_P and β.
  - 2021 (JINST 16 P03026): a Ne/CO₂ and Ne/N₂ model in which r *decreases* at high quencher fraction.

  — [Uludağ Ne-CO₂ record](https://acikerisim.uludag.edu.tr/items/9682467d-bbe9-4677-a755-f7acb1b73edd/full); [Ne-N₂ 2020](https://avesis.uludag.edu.tr/yayin/55691752-4c45-4fce-b3f7-ed4faa0aa721/measurements-and-calculations-of-gas-gains-in-ne-n-2-mixtures-pressure-and-concentration-scaling); [JINST 16 P03026](https://iopscience.iop.org/article/10.1088/1748-0221/16/03/P03026)
- **[fact]** Xe/TMA (JINST 13 (2018) P10032): r_P depends only on pressure, at 5% TMA; the formula gives about 0.23 at 1 bar. An earlier microbulk model assumed r_p = 0.24 ± 0.04, and a direct measurement found up to about 35%. — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc); [arXiv:1506.05077](https://arxiv.org/pdf/1506.05077)

### Inferences
- **[inference]** For Ar/iC₄H₁₀, the most-used Micromegas gas, the *only* systematic r_P is a single 10%, 1 atm literature-curve point (0.40 ± 0.01, from Agrawal & Ramsey 1988–89 data). There is no concentration or pressure scan comparable to the Ar/CO₂ programme. Values used for 95:5 or 98:2 mixtures are extrapolations.
- **[inference]** Model dependence is real. With Magboltz 2 (lumped levels), Andronic et al. fitted a 37% "Penning fraction" for Ar/CO₂ 85:15. The Şahin 2014 formula gives about 0.52 for that composition with Magboltz 8.x and later level-resolved cross-sections. r_P is defined relative to a specific cross-section set and level eligibility, and is not portable between Magboltz versions. — [Andronic 2004](https://arxiv.org/abs/physics/0402044)

### Gaps
- I could not read the exact 2014 tabulated r(c, p) with errors, the β values, or the 2017 process breakdown (paywalled).
- No Şahin-group Ar/iC₄H₁₀ concentration or pressure scan was found.

---

## 5. Known problems: degeneracies, geometry dependence, and gain versus W/Fano consistency

### Takeaway
A fitted r_P is an *effective* parameter. It absorbs:
- errors in Magboltz excitation and ionisation cross-sections;
- the gain normalisation (g);
- photon feedback (if β is not fitted);
- field-map errors, space charge and charging-up in MPGDs.

It also depends on geometry, through radiation trapping and photon escape and through the different fields and time scales of wires versus parallel plates versus MPGDs. The same r_i in principle sets the Jesse-effect W and Fano reduction, but there is little direct cross-check.

### Cited Findings

**THGEM study (Azevedo et al. 2016)**
- **[fact]** Garfield++ THGEM gain calculations were "known to produce differences of 2 orders of magnitude" versus data. Including Penning resolves this, using r = 0.18 (Ar/5% CH₄, Şahin 2010), r = 0.57 (Ar/30% CO₂, Şahin 2014), and an empirical r = 0.4 for Ne/5% CH₄ ("just an estimative").
- **[fact]** For Ar/30% CO₂, the simulation matches the data up to G ≈ 500 but is a factor of 2 low at about 2×10⁴. This is attributed to photon feedback, which was not simulated. A summary snippet says r ≈ 0.4 restored agreement at high gain; I did not verify this in the full text.
- **[fact]** Charging-up was not included.

— [Azevedo et al., arXiv:1606.04852 (2016)](https://ar5iv.arxiv.org/html/1606.04852)

**Spatial extent of transfers (ALICE TPC example)**
- **[fact]** Microscopic tracking with on-the-spot Penning transfers (r = 0.57) gives an azimuthal avalanche spread σ_φ ≈ 33–37°. The ALICE TPC wire chambers in Ne/CO₂/N₂ show 50–120°. Explaining that needs photon-mediated, spatially extended transfers. The gain-curve fits also need a feedback term, β = 2×10⁻⁵. "Because of a lack of input data" a detailed model reproducing this could not be implemented. — [Schindler 2012, §3.3](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Micromegas single-electron response**
- **[fact]** For Micromegas in Ne/iC₄H₁₀ 95:5 (Zerguerras et al.), the Polya relative width is f ≈ 0.31 ± 0.02 below U_mesh = 500 V and rises above that, attributed to photon feedback. "The dominant source of uncertainty in the calculation of gain spectra for this gas mixture is given by the Penning transfer efficiency: as yet, no estimates for r in Ne/iC₄H₁₀ from gain curve [fits exist]". — [Schindler 2012, §3.4.1](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf); [Zerguerras et al., arXiv:0901.2882](https://arxiv.org/pdf/0901.2882)
- **[fact]** Penning transfer reduces the relative avalanche variance, because it adds more ionising "trials". This is part of why r_P matters for Polya θ as well as for mean gain. — [Schindler 2012, §3.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Geometry dependence**
- **[fact]** Ar/C₂H₂ asymptotic r: 0.72 ± 0.02 in cylindrical chambers versus 0.67 ± 0.02 in a parallel plate, with ±0.1 systematic. Garfield++ hard-codes the cylindrical value. — [Şahin 2010, §4.3](https://doi.org/10.1088/1748-0221/5/05/P05002); [Garfield++ MediumGas.cc comment](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)

**Jesse effect, W and Fano**
- **[fact]** The Jesse effect is W' = W₀/(1 + (1/N_ion)·Σ_j r_j·N_exc,j). With Şahin's Ar/CO₂ r(c), Schindler computes for 1 keV electrons:
  - a W minimum of about 25.3 ± 0.1 eV at c ≈ 8%;
  - a Fano-factor minimum at c ≈ 5%.

  Bortner et al. measured, with alphas, W_min = 26.05 eV at 4.5% CO₂. Schindler attributes the difference partly to alpha versus electron W, and to the r(c) parametrisation being poorly constrained below 5%. The calculation uncertainty from cross-sections is about ±5% for W and ±10% for F. — [Schindler 2012, §4.2](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)
- **[fact]** In Penning mixtures, linear interpolation of W and F between the pure-gas values does not apply. — [Schindler 2012](https://repository.cern/records/a5nhh-0sg66/files/CERN-THESIS-2012-208.pdf)

**Practitioner reports**
- **[fact]** A neon-based PICOSEC/two-stage Micromegas optimisation study (2024) fixed r_P = 0.7 for all Ne mixtures instead of fitting. A Fermilab user talk reported that Garfield's built-in Penning did not reproduce expected results and resorted to a modified transfer coefficient. — [arXiv:2411.04429](https://arxiv.org/pdf/2411.04429); [FNAL indico talk](https://indico.fnal.gov/event/44260/contributions/190497/attachments/130861/159799/GarfieldGain.pdf)

### Inferences
- **[inference]** In primary ionisation (Degrad, Heed), r_P acts on excitations made by low-energy δ-electrons at zero field. In the avalanche it acts on excitations at 10–100 kV/cm. The level mix differs: relatively more 4s at low energy. Using one global r_P for both can therefore be inconsistent, especially for Ar/iC₄H₁₀, where 4s is eligible.

  This is directly relevant to this repo's Stage A primary-ionisation and W checks; see memory note "Primary ionization Heed check". Fitting r_P to gain and then using it for W can double-count or miscount if the gain fit absorbed cross-section errors.
- **[inference]** For Micromegas (gap of about 128 µm, field about 30–40 kV/cm, avalanche lasting about 1 ns), the following are all comparable to the avalanche itself:
  - the electron drift time over the gap;
  - the Penning delay τ_P;
  - photon ranges of 10–100 µm.

  The instantaneous, on-the-spot model is therefore least justified there. Expected consequences, not verified: signal timing (PICOSEC), avalanche lateral size, and possibly a gain bias if delayed electrons are created in the lower-field region near the mesh.
- **[inference]** Any r_P fitted to *your own* Micromegas gain curve with a Garfield++ field map will also absorb mesh-geometry and gap-tolerance errors, charging and space charge. It is a detector-specific tuning constant, not a gas property.

### Gaps
- No published direct comparison was found of r from gain fits versus r from measured W or Fano in the same Ar/iC₄H₁₀ or Ar/CO₂ mixture.
- No Micromegas-specific Ar/iC₄H₁₀ r_P fit was located.

---

## 6. Recent work (2018–2026)

### Takeaway
Since 2018 the field has mostly *extended the fitted-r_P tables*: Xe/TMA (2018), Ne/N₂ (2020), Ne/CO₂ and Ne/N₂ model (2021). Garfield++ has refined but not fundamentally changed its de-excitation model (Ar only). I found no new ab-initio or excited-state-kinetics Penning model for detector gases, no Penning study specific to wet gases, and no ternary-mixture r_P parametrisation in Garfield++.

### Cited Findings
- **[fact]** Garfield++ master (2026) contains default r_P only for the binary mixtures in the §2 table. The newest entries are Ne/CO₂ and Ne/N₂ (2021) and Xe/TMA (2018). — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)
- **[fact]** The current de-excitation code still supports only Ar as the excited species, with quencher data for CO₂, CH₄, C₂H₆, iC₄H₁₀, C₂H₂ and CF₄. It has received code-level fixes (thread-safe counters, explicit resolution of HIGH-level targets, per-file streams) but no new physics data that I could see. — [Garfield++ MediumMagboltz.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumMagboltz.cc)
- **[fact]** A GPU microscopic Monte Carlo (OuroborosBEM, 2021) cites a measured Penning probability of 0.574 for Ar/CO₂ 70:30, i.e. it reuses the Şahin r_P. — [arXiv:2110.09214](https://arxiv.org/pdf/2110.09214)
- **[fact]** 2026 conference abstracts describe kinetic models for Ar–CF₄ and He–CF₄ scintillation using Garfield++ cross-sections. They are about scintillation, not Penning gain. — [indico.nikhef.nl 7584](https://indico.nikhef.nl/event/7584/contributions/31322/); [indico.cern.ch 1473492](https://indico.cern.ch/event/1473492/contributions/7106078/)
- **[fact]** The 2024 neon two-stage Micromegas optimisation used a fixed r_P = 0.7 instead of a fit. — [arXiv:2411.04429](https://arxiv.org/pdf/2411.04429)

### Inferences
- **[inference]** The practical state of the art for an Ar/iC₄H₁₀ or Ar/CO₂ Micromegas (relevant to MX17) is:
  - (a) Magboltz gas table plus `EnablePenningTransfer(r, 0, "ar")`, with r tuned to your own measured gain curve, while fitting a gain scale and a feedback β at the same time;
  - (b) quoting the r_P spread, for example 0.3–0.5 for Ar/iC₄H₁₀ at 1 atm, as a systematic.

  The detailed `EnableDeexcitation` model is available for Ar with iC₄H₁₀ or CO₂ but is not validated beyond the 2012 tuning set. It is useful for exploring timing and spatial effects, not for predicting absolute gain.

### Gaps
- My searches did not cover the RD51/DRD1 indico talks (2022–2026) thoroughly. There may be unpublished Penning work, for example by Veenhof or Schindler, or gas-gain campaigns in DRD1 WG4/WG2.
- No sources found on wet gases (H₂O admixture) and Penning, or on photon-feedback modelling in Garfield++ beyond `EnableRadiationTrapping` and the β fit term.
