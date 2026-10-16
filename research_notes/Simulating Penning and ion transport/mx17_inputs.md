# Quantitative inputs for an MX17 Micromegas simulation in Ar/iC4H10 95/5 (dry and wet): Penning transfer and ion transport

> **Correction 2026-10-16:** MX17's amplification gap is **150 µm** (ESL to mesh underside), not the generic 128 µm assumed in these notes, and the field the ions see is the T6 map's 31.0 kV/cm = 124 Td (effective gap 158 µm, because the woven mesh's equivalent plane sits ~8 µm inside the weave; `response/meshcell/mesh_plane_check.py`). The 38.3 kV/cm / 153 Td figures and the transit tables below that use d = 128 µm do not apply to MX17; the speed ratios relative to Ar⁺ still do. Corrected numbers: `reports/Simulating Penning and ion transport.md`, field reconciliation.

Scope: resistive bulk Micromegas, 128 µm gap (WRONG for MX17, see above), ~490 V on the mesh, Ar/iC4H10 95/5 at about 1 atm and room temperature, with and without a few hundred to a few thousand ppm H2O. Researched 2026-10-15.

Conventions used below:
- N(293.15 K, 1 atm) = p/(k_B T) = 2.503e19 cm^-3.
- K is the mobility at operating density. K0 is the reduced mobility at 273.15 K and 1 atm, so K = K0 · (T/273.15 K) · (1 atm/p).

**A unit check the report writer must carry forward (my calculation, not from a source):**
- E = 490 V / 0.0128 cm = 38.28 kV/cm. At 293 K and 1 atm, E/N = 38281 / 2.503e19 = 1.529e-15 V cm² = **153 Td**, not the 120–125 Td in the brief.
- The project's own close-out note says "the amp gap runs at 123.8 Td". That matches N ≈ 3.09e19 cm^-3, roughly 1.23 atm at 293 K. The gas-table pressure/temperature used in the sim is the likely place to look.
- In a uniform gap the ion transit is t = d²/(K·V), independent of N for a given K. So:
  - Garfield's own Ar+ table at 153 Td and 293 K (K0 = 1.154, K = 1.238 cm²/Vs) gives t = **270 ns**, not 340 ns.
  - K ≈ 0.98 cm²/Vs, the value that gives 340 ns, follows only if the sim's gas density is about 25% higher than 1 atm at 293 K.
- **The 140 ns transit the fits prefer requires K = d²/(V t) = (0.0128)²/(490 × 140e-9) = 2.39 cm²/Vs. That is K0 = 2.23 cm²/Vs at 1 atm/293 K, or proportionally more if the true density is higher.** At the 1.23 atm implied by "123.8 Td" it would be K0 ≈ 2.75, which exceeds the Langevin limit for any hydrocarbon ion in Ar (see Q3).
- The required ratio to Garfield's Ar+ at the same E/N is 2.23/1.154 = **1.9×**. The ×2.4 quoted in HANDOFF.md comes from the K ≈ 1.0 baseline.

---

## Q1. Penning transfer probability r_P for Ar/iC4H10 at 2–10% iso (especially 5%): published values, Garfield++ default, Micromegas fits, and field/gap dependence

### Takeaway
- The only 1-atm literature value is **r_P = 0.40 ± 0.01 for Ar + 10% iC4H10** (Şahin et al. 2010). Garfield++ applies it as a flat constant for any iso fraction.
- Şahin's later fits at 30 kPa, with photon feedback fitted at the same time, give r_P ≈ 0.19–0.30 for 3.3–10% iso, with a non-monotonic 5% point of 0.19.
- Micromegas-specific Garfield++ fits for Ar/iso 90/10 range from **0.40–0.50** (64–220 µm gaps) to **0.80** (128/192 µm). These are effective parameters that absorb other mismodelling.
- For 95/5 at 1 atm, **r_P ≈ 0.35–0.45 is the defensible range, with 0.40 as the central value** (inference). Treat r_P and photon feedback as degenerate unless fitted jointly.

### Cited Findings
- **Garfield++ default for Ar/iC4H10** (the table is in `MediumGas::EnablePenningTransfer()` in MediumGas.cc, not MediumMagboltz.cc; MediumMagboltz just calls it). The code hard-codes `rP = 0.40` with the comment "There is only one value for this mixture: c = 0.1, p = 1 atm". If the iso fraction differs from 10% by more than 1%, or p differs from 1 atm, it prints "Using transfer probability for 10% iC4H10" / "at atmospheric pressure" and still uses 0.40. The default only works for 2-component mixtures (`m_nComponents != 2` returns false), so it cannot be used for the wet ternary. — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc) (local copy, commit 927e5c21, 2026-08-06)
- Other defaults from the same function, for comparison:
  - Ar/CO2 at 1 atm: rP = (0.6643c + 0.0028)/(c + 0.0518), giving **0.405 at 7% CO2** and 0.456 at 10% (my evaluation).
  - Ar/C3H8: rP = 0.4536c/(c + 0.0035), giving **0.424 at 5%** and 0.438 at 10% (my evaluation).
  - Ar/C2H6 10%: 0.31.
  - Ar/CH4: a pressure- and concentration-dependent formula from Şahin 2010.
  - — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)
- Şahin, Tapan, Özmutlu, Veenhof, JINST 5 (2010) P05002, Table 2, "Transfer probabilities in argon with 10% alkane at 1 atm":

  | Alkane | Ionisation potential | r_P |
  |---|---|---|
  | CH4 | 12.65 eV | 0.212 ± 0.002 |
  | C2H6 | 11.52 eV | 0.31 ± 0.01 |
  | C3H8 | 10.95 eV | 0.43 ± 0.01 |
  | iC4H10 | 10.67 eV | **0.40 ± 0.01** |

  The data are Agrawal et al. proportional-counter gain curves, and iso was measured only at 10%. The authors say isobutane "probably falls out of this trend because it is compact and hence has a smaller collision time". — [Şahin et al. 2010, repository copy](https://acikerisim.uludag.edu.tr/items/0c97c615-5554-4cf7-b343-a1557f6ca0ae); DOI 10.1088/1748-0221/5/05/P05002
- Şahin 2010 method:
  - A single r is fitted per gain curve, with a common gain-scale factor g for all anode voltages. r is not field-dependent in the model.
  - Its concentration and pressure dependence comes from collision time versus excited-state lifetime. Table 1 gives the mean time until an Ar* meets an iC4H10 molecule as **τ = 100 ps at 1 atm for c = 1**, so about 2 ns at 5% (my scaling).
  - Ar 3p⁵4p levels (~13 eV) have lifetimes of 21.7–40.5 ns. 3p⁵3d levels (threshold 13.85 eV) have ~50 ns, and their J=1 levels ~3.5 ns.
  - Argon excimer (second continuum, 9.8 eV) photons cannot ionise iC4H10.
  - — [Şahin et al. 2010](https://acikerisim.uludag.edu.tr/items/0c97c615-5554-4cf7-b343-a1557f6ca0ae)
- **Ar/iC4H10 at 3.3–30% iso** (Şahin 2018, Uludağ Univ. J. Fac. Eng. 23(3) 49–60):
  - Data: Bronić & Grosswendt (1998) single-wire counter gain data at **30 kPa, 300 K**.
  - Fit: a common r_Pen for all Ar excited levels, plus a photon-feedback parameter β.
  - Table 2 results:

    | iC4H10 fraction | r_Pen |
    |---|---|
    | 3.3% | 0.295 ± 0.005 |
    | **5%** | **0.189 ± 0.007** |
    | 10% | 0.258 ± 0.004 |
    | 20% | 0.188 ± 0.006 |
    | 30% | 0.147 ± 0.009 |

  - Feedback β falls systematically from 3.3% to 20%, and no β is needed at 30%. Feedback visibly breaks proportionality above gain ~1000.
  - The author reads the falling r with rising iso as excitation energy going into non-ionising dissociation of iso.
  - — [Şahin 2018](https://acikerisim.uludag.edu.tr/items/0be7e437-5772-4e34-a74a-dedfbd782b6f)
- **Bulk Micromegas, Ar/iC4H10 90:10, 128 and 192 µm gaps** (Bhattacharya et al. 2013):
  - The authors "chose to carry out the calculations with 80% transfer rate for Argon-Isobutane mixture (90:10) and 25% for P10 since these values agreed well with the experimental data".
  - The 80% "was much higher than that predicted in [Şahin 2010]", and even with 80% the simulated gains are "still lower than the measured data". They suggest the real gap may differ from nominal.
  - — [arXiv:1304.6228](https://arxiv.org/pdf/1304.6228)
- **Bulk Micromegas, several gaps (64–220 µm), Ar/iC4H10 90:10** (Bhattacharya et al. 2016, with Colas and Attié of Saclay): "we chose 40 − 50% transfer rate ... slightly higher than that estimated (40%) in reference [28] [Şahin 2010] and needs further investigation". Their Fig. 4 uses 40% in two panels and 50% in two others. — [arXiv:1605.02921](https://arxiv.org/pdf/1605.02921)
- A THGEM study reports that Garfield++ without Penning transfer can miss measured gain by up to ~2 orders of magnitude. — [arXiv:1606.04852](https://arxiv.org/pdf/1606.04852) (search snippet only, not read in full)
- The lowest Ar excited level (3p⁵4s) is at 11.55 eV, above the iC4H10 ionisation threshold of 10.67 eV, so all Ar excited levels can in principle Penning-ionise iso. — [Şahin 2018](https://acikerisim.uludag.edu.tr/items/0be7e437-5772-4e34-a74a-dedfbd782b6f); [Şahin 2010](https://acikerisim.uludag.edu.tr/items/0c97c615-5554-4cf7-b343-a1557f6ca0ae)

### Inferences
- **r_P versus iso fraction:**
  - The Ar/C3H8 Garfield fit is the closest molecular analogue with a measured concentration dependence. It changes r_P by only −3% from 10% to 5% (0.438 to 0.424), because collisions (~2 ns at 5%) remain much faster than the radiative loss of the metastable and radiation-trapped 4s states.
  - **So r_P(5% iso, 1 atm) ≈ 0.38–0.40** from scaling Şahin's 10% value. This is an inference, not a measurement.
  - The 30 kPa Şahin 2018 values (0.19–0.30) are probably lower partly because pressure is lower: collision times are 3.4× longer, so more excited states decay radiatively before transfer. The non-monotonic 5% point suggests r_P and β are not cleanly separated at low iso.
- **Micromegas "fitted" values (0.4–0.8) are not clean measurements of r_P.** They absorb gap tolerance, mesh transparency, electronics calibration and, most importantly for MX17, photon feedback. MX17's own campaign (EOS_README, FRESH_LOOK §9) found that changing r_P from 0.30 to 0.80 does not reproduce the measured gain-vs-voltage slope, while avalanche photon feedback G/(1−bG) does. Şahin 2018 likewise needed β at 3.3–10% iso. At 5% iso a quencher-poor mixture shows feedback, so **r_P should be fitted together with a feedback term or fixed at ~0.40.**
- **Field/gap dependence:** none of the sources finds a genuine field dependence of r_P. Şahin's framework treats r as field-independent. The Bhattacharya 40% vs 50% split across gaps is within what one fit's systematics can produce. An explicit gap-dependent r_P in Micromegas is therefore not supported.

### Gaps
- No 1-atm measurement of r_P for Ar/iC4H10 below 10% iso was found.
- No Micromegas Penning fits specific to T2K, CAST, ILD, MIMAC or COMPASS Ar/iso were located. T2K uses Ar/CF4/iC4H10 95/3/2 (ternary), and ATLAS NSW uses Ar/CO2 93/7, where Garfield's default gives 0.405. Zerguerras et al. (NIM A 772 (2015) 76) measured single-electron gain in Ar/iC4H10 95/5 in a Micromegas, but the full text was blocked by the CDS bot-check, so any r_P they used is unverified. — [CDS record](https://cds.cern.ch/record/2153708/files/nima772-76.pdf)
- The Şahin 2010 per-mixture concentration fits for Ar/iC4H10 do not exist, since there is only one data point.

---

## Q2. Effect of water vapour: Penning-energetics check, Ar* quenching by H2O, gain vs humidity data, attachment, and H3O+(H2O)n ions

### Takeaway
- Ar 4s (metastable and resonant) states, at 11.55–11.83 eV, cannot Penning-ionise H2O (IP 12.62 eV). Water quenches them dissociatively, with **k ≈ 4.8e-10 cm³/s**, comparable to iso's collision rate.
- At 5% iso, 300–3000 ppm H2O takes only **~0.5–9% of Ar* quenching**, so it should lower the effective r_P by only a few percent relative. The model choice of 0.40 dry and 0.38–0.40 wet is therefore close (inference).
- Measured Micromegas/GEM gains in Ar/CO2 90/10 change by **≤10–20% over 0–5000 ppm**.
- **H3O+ does not form in significant amounts**, because hydrocarbon ions do not transfer protons to water. Water can still slow ions at low drift fields through clustering, observed at **−5% K0 per +100 ppm** in Ne/CO2. That effect is expected to be weaker at the ~150 Td amplification field.

### Cited Findings
- Ar 3p⁵4s levels: lowest 11.55 eV. 3p⁵4p levels, about 13 eV, decay to 4s in 21.7–40.5 ns. 3p⁵3d starts at 13.85 eV. — [Şahin et al. 2010](https://acikerisim.uludag.edu.tr/items/0c97c615-5554-4cf7-b343-a1557f6ca0ae)
- The H2O ionisation potential, 12.62 eV, is the value given in the brief (NIST WebBook lists 12.621 eV, not re-fetched this session; [NIST WebBook H2O](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=20)). The four 4s levels, 11.55–11.83 eV, therefore lie 0.8–1.07 eV below it, so **ground-state-water Penning ionisation by Ar(4s) is endothermic**.
- **Ar(1s5) [= 3P2] + H2O → Ar + H + OH, k = 4.8e-10 cm³/s**, attributed to de Jong and confirmed by Leblond et al.; used in a chemical-kinetic model of Ar–H2O dielectric-barrier discharges. — [Ar–H2O DBD kinetic model (ResearchGate)](https://www.researchgate.net/publication/340031752_A_chemical-kinetic_model_of_DBDs_in_Ar-H_2_O_mixtures). The primary de Jong / Velazco–Kolts–Setser papers were not retrieved, because the DTIC copies returned HTTP 403.
- Ar(3P2) quenching rate constants by small molecules span 1.5–50 × 10⁻¹¹ cm³/s (Piper, Richardson, Taylor, Setser 1972). — [Faraday Discussions 53 (1972)](https://pubs.rsc.org/en/content/articlelanding/1972/dc/dc9725300100) (abstract only)
- Novicki & Krenos measured an absolute quenching cross-section for Ar(3P2,0) + H2O (J. Chem. Phys. 89, 7031 (1988)). The value was not retrieved. — [Rutgers Krenos page](https://chem.rutgers.edu/people/emeritus-faculty/emeritus-faculty-member/207-krenos-john)
- **Gain vs humidity in MPGDs** (Friberg et al., arXiv:2503.03012, NIM A 2026): GEM, THGEM and Micromegas (including a 128 µm zMMG), Ar/CO2 90/10, 0–5000 ppmV H2O.
  - GEM gain rises by up to 10% up to ~3000 ppmV (less clear on the decreasing scan).
  - GEM current at 360–380 V varies by less than 3%, attributed to T/p.
  - Micromegas gain "remains relatively stable, or no clear dependency can be observed".
  - The authors bound the gain change across humidity at ≤10–20%.
  - Magboltz α−η for Ar/CO2 90/10 at 10–120 kV/cm is "constant" with water content.
  - Penning effects are not discussed. The recommended range is 100–1000 ppmV, and discharge stability improves with humidity.
  - — [arXiv:2503.03012](https://arxiv.org/html/2503.03012v2)
- The same paper notes that prior water-vapour studies are "inconclusive, and no consensus has been reached". — [arXiv:2503.03012](https://arxiv.org/html/2503.03012v2)
- **Ion–molecule rates** (Anicich 1993):
  - **Ar⁺ + H2O**: total 1.62e-9 cm³/s, giving H2O⁺ (0.80) and ArH⁺ (0.20).
  - **H2O⁺ + H2O → H3O⁺ + OH**: 1.85e-9 cm³/s.
  - **H2O⁺ + C2H6 → H3O⁺ + C2H5** (0.83 of 1.6e-9 cm³/s).
  - Ar⁺ + C2H6: 1.15e-9 cm³/s.
  - The 1993 evaluation has **no Ar⁺ + iC4H10 entry**.
  - — [Anicich, J. Phys. Chem. Ref. Data 22, 1469 (1993)](https://srd.nist.gov/jpcrdreprint/1.555940.pdf)
- **Proton affinity of isobutene = 802 ± 4 kJ/mol** (theory), against an experimental 801.7 ± 1.4 kJ/mol. The tert-butyl cation t-C4H9⁺ is protonated isobutene. — [Smith & Radom, via ACS Figshare](https://acs.figshare.com/collections/Heat_of_Formation_of_the_i_tert_i_-Butyl_Radical/3416628)
- **Water and ion mobility** (Deisting et al. 2016, ALICE TPC R&D, 1 atm, low drift fields): "A decrease of K was measured for increasing water content."
  - In Ne/CO2 90/10: "For an increase of 100 ppm in water a decrease by 5 % in K0 was found".
  - Their CO2-ion K0 is up to 10% (Ar) or 20% (Ne) lower than Coimbra's 8 Torr, <20 ppm H2O values, which they attribute to water and pressure.
  - — [arXiv:1603.07638](https://arxiv.org/pdf/1603.07638)
- **Cluster ions** (Kalkan, Arslandok, Cortez, Kaya, Tapan, Veenhof, JINST 10 (2015) P07004): in Ar/CO2 and Ne/CO2 at 1 atm, noble-gas ions transfer charge to CO2 "within a few nanoseconds", and a few ns later CO2⁺ has picked up CO2 to form CO2⁺·(CO2)n. These "drift more slowly", changing the ion-motion signals, which they expect to affect Micromegas. — [Kalkan et al. 2015](https://avesis.uludag.edu.tr/yayin/6a020f55-18d2-48e6-bfed-9e6fd91abe0e/cluster-ions-in-gas-based-detectors)
- **Old warning in a different device**: in a macro-argon ionisation detector, sensitivity fell 10-fold when water went from ~30 to ~1000 ppm. This is a different, non-Penning-optimised device. — [Nature 205, 593 (1965)](https://www.nature.com/articles/205593a0.pdf)

### Inferences
- **Share of Ar* quenching taken by water** (my calculation, using k_w = 4.8e-10 cm³/s):
  - The iso quenching rate coefficient is not directly sourced. Şahin's collision time of 100 ps at c = 1 corresponds to k_iso ≈ 1/(τN) = 4.0e-10 cm³/s; I bracket it at 3–9e-10.
  - Water share = k_w·n_w / (k_w·n_w + k_iso·n_iso), with n_iso = 0.05 N:

    | k_iso | 300 ppm | 1000 ppm | 3000 ppm |
    |---|---|---|---|
    | 3e-10 cm³/s | 0.95% | 3.1% | 8.8% |
    | 6e-10 cm³/s | 0.48% | 1.6% | 4.6% |
    | 9e-10 cm³/s | 0.32% | 1.1% | 3.1% |

  - The effective r_P scales roughly as (1 − share), so **r_P(wet) ≈ 0.40 × (0.91–0.995) ≈ 0.36–0.40**. A 5–10% change in r_P moves the gain by of order 10–20% at G ~ 10³–10⁴, consistent with the ≤10–20% bound measured in Ar/CO2 (inference).
- **Ar 4p/3d levels (12.9–14 eV) are energetically able to Penning-ionise H2O**, but with water ~50× rarer than iso, essentially all of their energy goes to iso first. Garfield's per-gas Penning (`EnablePenningTransfer(r, λ, "H2O")`) applies only to levels above that gas's IP. A global r applied to all levels above the lowest IP, 10.67 eV for iso, already covers this (inference).
- **Ar⁺ diverted to H2O⁺** = 1.62e-9·n_w / (1.62e-9·n_w + 1.37e-9·n_iso), using the Langevin estimate for Ar⁺ + iC4H10 from Q5: **0.7% at 300 ppm, 2.3% at 1000 ppm, 6.6% at 3000 ppm**.
- **H2O⁺ will not reliably become H3O⁺.** With iso 17–170× more abundant than water, H2O⁺ (12.62 eV) is overwhelmingly likely to charge-transfer to iso (10.67 eV) before meeting a water molecule (inference from IPs). So **H3O⁺ is at most a ~1–7% channel, and probably much less.**
- Hydrocarbon ions cannot hand their proton to water. The PA of H2O (~691 kJ/mol) is below that of isobutene (802) and propene (~752). [Values for H2O and propene are standard Hunter & Lias 1998 NIST numbers, quoted from memory and not re-fetched; flag.] So **proton-hydrate H3O⁺(H2O)n chemistry does not take over the ion population at ppm water in Ar/iso.**
- What water can do is **associate** with t-C4H9⁺ or C3H7⁺ in three-body collisions, forming e.g. protonated tert-butanol. Whether this happens and survives at the ~1000–1400 K effective ion temperature in the gap (see Q3) depends on the ion–H2O binding energy, which I did not find.
- The Deisting −5%/100 ppm effect is a low-field drift-region result. **If anything, water lowers the amplification-gap ion mobility**, pushing the real ion further from the 140 ns that the fits prefer.
- **Electron attachment**: Friberg et al.'s Magboltz calculation shows no change in α−η with up to 5000 ppm H2O in Ar/CO2. H2O's dissociative attachment (H⁻, threshold ~4.4 eV) is weak (physics from memory, flag). Attachment is therefore not a significant gain lever at ppm-level water.

### Gaps
- No primary measurement of the Ar(3P2)/Ar(3P0) + iC4H10 quenching rate constant was retrieved. Velazco–Kolts–Setser 1978 (J. Chem. Phys. 69, 4357) would have it, but the DTIC access returned 403.
- No measured gain-vs-humidity data specifically in Ar/iC4H10 (any device) was found, and no Penning model of humidity effects was found.
- Ion-cluster binding energies (t-C4H9⁺·H2O, C3H7⁺·H2O, H3O⁺·(H2O)n) and their survival at E/N ~150 Td were not sourced.
- The cross-section for electron attachment to water clusters at ppm levels was not researched.

---

## Q3. Mobility of hydrocarbon ions (C3H7⁺, C3H5⁺, C4H9⁺, iC4H10⁺) in Ar at ~120–150 Td; Langevin limit; comparison with Ar⁺ in Ar

### Takeaway
- **Ar⁺ in Ar is slow because of resonant charge exchange.** Garfield's table gives K0 = 1.22 at 120 Td and 1.16 at 150 Td.
- Hydrocarbon ions in Ar have no resonant charge exchange and sit close to the Langevin polarisation limit: **K0_L = 2.37 (C3H7⁺, m = 44), 2.23 (C4H9⁺/iC4H10⁺, m = 57–58), 1.99 (C8H19⁺ dimer)**.
- Measured analogues in Ar-rich mixtures are 80–100% of Langevin: CF3⁺ in 90% Ar has K0 ≈ 2.11, and the CO2-derived ion in pure Ar extrapolates to K0 = 1.9.
- Blanc's law with 5% iso (K0 ≈ 0.65 in pure iso) lowers this to **K0 ≈ 1.7–2.0** at low field.
- **K0 ≈ 2.2 (140 ns) is at the very top of what is physically possible. K0 ≈ 1.7–2.0 (155–185 ns) is the realistic expectation**, i.e. 1.5–1.7× Garfield's Ar⁺ at 153 Td, not 1.9× (inference).

### Cited Findings
- **Garfield `IonMobility_Ar+_Ar.txt`** (T = 300 K; Hornbeck 1951, Beaty 1968, Basurto PRE 61, 3053): K0 = 1.53 (0–8 Td), 1.41 (50), 1.27 (100), **1.22 (120)**, **1.16 (150)**, 1.06 (200), 0.95 (300), 0.78 (500) cm²/Vs. Accuracy is quoted as 1% below 250 Td and 3% above. — [Garfield++ Data](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Data/IonMobility_Ar+_Ar.txt)
- Garfield's `LoadIonMobility` treats file values as reduced mobilities at 273.15 K and 1 atm, scaled by (T/273.15)/(p/atm). The E/N axis is in Td, scaled with the Loschmidt number. — [Garfield++ MediumGas.cc](https://gitlab.cern.ch/garfield/garfieldpp/-/blob/master/Source/MediumGas.cc)
- **Garfield `IonMobility_C8Hn+_iC4H10.txt`**: K0 = **0.65 cm²/Vs** (flat, 0–45 Td) for C8Hn⁺ in pure isobutane at 293 K, from Cortez et al. JINST 14 (2019) P04010. — [Garfield++ Data](https://gitlab.cern.ch/garfield/garfieldpp/-/tree/master/Data)
- **Cortez et al. 2019** (Coimbra): ion mobilities in pure iC4H10 and Ar/iC4H10 with 5–95% Ar, at 8 Torr (10.6 mbar), room temperature, 10–45 Td, taken from time-of-arrival peak centroids. The numerical Ar-rich values could not be retrieved, because the full text is paywalled. — [IOP abstract](https://iopscience.iop.org/article/10.1088/1748-0221/14/04/P04010); [avesis record](https://avesis.uludag.edu.tr/yayin/f2ec6582-4ec7-4f20-8cc5-92abdb2bbb59/experimental-ion-mobility-measurements-in-pure-ic4h10-and-ar-ic4h10-mixtures)
- Coimbra's method uses the **Langevin limit K0 = 13.88/√(α μ)** (α in Å³, μ in amu) and **Blanc's law 1/K0_mix = f1/K0_g1 + f2/K0_g2** to identify ions. Measurements at 8 Torr:
  - In pure Ar: K0 ≈ 1.57 (assigned Ar⁺) and **≈ 1.92 (assigned Ar2⁺)** at 30 Td.
  - **In 90% Ar/10% CF4: K0 ≈ 2.11 cm²/Vs, assigned CF3⁺ (m = 69)**, at 15 Td. Ar⁺ + CF4 → CF3⁺ "prevents the formation of Ar2⁺".
  - In pure CF4, CF3⁺ has K0 ≈ 1.10 against a Langevin value of 1.13.
  - — [Cortez, LCTPC meeting, Nov 2017 (slides)](https://agenda.linearcollider.org/event/7795/contributions/40342/attachments/32468/49376/ACortez_-_LCTPC_Collaboration_Meeting_-_November_2017.pdf)
- In argon–ethane mixtures (Coimbra), two arrival peaks were attributed to 3-carbon and 4-carbon ions. — search snippet for [Cortez et al., Xe–C2H6 / Ar–C2H6 papers](https://www.osti.gov/etdeweb/biblio/22335751) (not read in full)
- **Deisting et al. 2016** (1 atm, Ar/CO2): Blanc fit gives the CO2-derived ion **K0 = 1.9(1) cm²/Vs in pure Ar** and 0.9(1) in CO2. In Ne the values are 4.0(2) and 1.0(1). — [arXiv:1603.07638](https://arxiv.org/pdf/1603.07638)
- Blanc's law holds when charge transfer is negligible compared with polarisation attraction and short-range repulsion. — search snippet, Coimbra group (see [arXiv:1603.07638](https://arxiv.org/pdf/1603.07638) for the Blanc usage)
- MX17's own close-out note used Blanc to put the iso ion at K0 ≈ 1.16 against Ar⁺'s 1.21. That calculation evidently put **Ar⁺'s own mobility in Ar** in as K_Ar for the iso ion. — local file `design/report/S3_ION_CLOSEOUT_2026-08-09.md`

**Langevin calculation (explicit):**
- K0_L = 13.876 / √(α[Å³] · μ[amu]) cm²/Vs, with α(Ar) = 11.08 a0³ = 1.6411 Å³ and m(Ar) = 39.948.

| Ion | m (amu) | μ = mM/(m+M) | √(α μ) | K0_L (cm²/Vs) |
|---|---|---|---|---|
| Ar⁺ (no charge exchange; hypothetical) | 39.95 | 19.97 | 5.725 | 2.42 |
| C3H5⁺ | 41.07 | 20.25 | 5.765 | 2.41 |
| C3H7⁺ | 44.10 | 20.96 | 5.865 | **2.37** |
| t-C4H9⁺ | 57.11 | 23.51 | 6.211 | **2.23** |
| iC4H10⁺ | 58.12 | 23.68 | 6.234 | **2.23** |
| C8H19⁺ (iso dimer) | 115.2 | 29.66 | 6.977 | 1.99 |
| CF3⁺ (check vs Coimbra 2.11) | 69.0 | 25.30 | 6.444 | 2.15 |

- Langevin in pure iso, using α(iC4H10) = 8.14 Å³ (CRC-type value from memory, flag): C4H9⁺ 0.91, C8H19⁺ 0.78 cm²/Vs. Coimbra measured 0.65 for C8Hn⁺, about 83% of Langevin.

**Blanc for 95% Ar / 5% iso (explicit):**
- 1/K0 = 0.95/K_Ar + 0.05/K_iso.

| K_Ar | K_iso = 0.65 | K_iso = 0.75 |
|---|---|---|
| 2.2 | 1.97 | 2.01 |
| 2.0 | 1.81 | 1.85 |
| 1.9 | 1.73 | 1.77 |

**Transit at 490 V, 128 µm, 293 K, 1 atm (t = d²/(K V)):**

| K0 | Transit |
|---|---|
| 1.154 (Ar⁺, 153 Td) | 270 ns |
| 1.53 | 204 ns |
| 1.8 | 173 ns |
| 1.9 | 164 ns |
| 2.0 | 156 ns |
| 2.22 | 140 ns |

**Ion heating (Wannier, explicit):**
- (3/2)kT_eff = (3/2)kT + (1/2)M v_d², with M = m(Ar).
- For K0 = 2.0: v_d = K·E = 2.146 × 38281 = 8.2e4 cm/s, so (1/2)M v² = 0.14 eV and **T_eff ≈ 1370 K**.
- For K0 = 1.2: T_eff ≈ 680 K.

### Inferences
- A non-resonant ion of mass 43–58 in Ar should drift **1.5–1.9× faster than Ar⁺ in Ar at 150 Td**. The physical reason is that Ar⁺'s mobility is suppressed by resonant charge transfer. This supports the direction of the waveform fits.
- **The 140 ns fit (K0 = 2.23 at 1 atm/293 K) demands that the drifting ion be at or above the pure-Ar Langevin ceiling for m ≥ 57, with no Blanc dilution from the 5% iso and no high-field reduction.** That is not plausible for t-C4H9⁺ (the expected terminal ion; see Q5).
  - It is only marginally possible for C3H7⁺ (K0_L = 2.37), and Blanc already brings that to about 2.1.
  - The measured analogues in Ar-rich mixtures sit at 80–100% of Langevin (CF3⁺ 2.11 against 2.15; CO2-type ion 1.9 against 2.37 for m = 44, or 2.07 for a dimer).
  - **Best estimate: K0 ≈ 1.7–2.0, giving a transit of ~155–185 ns** (inference with ~±10% uncertainty).
- At T_eff ~1000–1400 K, the polarisation (1/r⁴) interaction gives a field-independent mobility in the pure Maxwell-model limit. Short-range repulsion usually lowers K by roughly 5–20% at a mean collision energy of ~0.1–0.2 eV, so the high-field value is unlikely to exceed the low-field one (inference, not sourced).
- The same ion heating argues against weakly bound clusters surviving in the amplification gap. The C8Hn⁺ dimer seen in pure iso at 8 Torr and low field may well not form or survive in 5% iso at 150 Td. If it did, Langevin falls to 1.99 and Blanc gives ~1.65 (inference).
- **If the sim's gas density really is ~1.23 atm** (the 123.8 Td reading), the required K0 becomes ~2.75, which **exceeds Langevin for every candidate ion**. The density and E/N assumption must be resolved before using the mobility argument quantitatively.

### Gaps
- No measured mobility for C3H7⁺, C3H5⁺, C4H9⁺ or iC4H10⁺ in Ar at any E/N above ~45 Td was found. The Cortez 2019 Ar-rich numbers (10–45 Td) are behind a paywall.
- No Viehland/Ellis-type (LXCat) high-E/N data for alkyl ions in Ar was located this session.
- α(iC4H10) = 8.14 Å³ and the PA values for H2O and propene are quoted from memory (CRC / Hunter & Lias) and should be checked.

---

## Q4. Micromegas signal and ion-tail studies that determined the effective ion mobility or transit from waveforms

### Takeaway
- The literature gives only **qualitative or order-of-magnitude ion-tail durations**, with no fitted effective mobility for Ar/iso:
  - Giomataris: "typically within 100 ns", and an ion tail "developed within … 80 ns" in Ar + 10% iso with a fast current preamplifier, gap not stated for that figure.
  - Generic ATLAS NSW / SPD statements: ~100–150 ns for 128 µm Ar/CO2.
- None of these determines K, though they are compatible with faster-than-Ar⁺ ions.

### Cited Findings
- Giomataris (ICFA Instrum. Bull., Fall 1999):
  - "The charge signal is mainly due to the positive ion drifting to the micromesh electrode, which takes place typically within 100ns, depending on the amplification gap and the gas mixture."
  - For Ar + 10% iso, going from a 100 to a 50 µm gap "reduces the signal rise time by a factor 3", and "Using a gap of 30 µm the rise time is only 17ns".
  - Fig. 7, current preamplifier in Ar + 10% iso: "the electron signal and the ion tail are developed within 5 ns and 80 ns respectively". The gap and voltage are not given in the text. The electron peak is ~10× the ion tail amplitude.
  - — [Giomataris, "MICROMEGAS: results and prospects"](https://www.slac.stanford.edu/pubs/icfa/fall99/paper1/paper1.pdf)
- ATLAS NSW Micromegas (128 µm, Ar/CO2 93/7): a CERN presentation says ions take "about 100 ns to reach the mesh". — [Indico slides](https://indico.cern.ch/event/235764/material/slides/1?contribId=0) (search snippet only)
- The NICA SPD TDR gives an "ion collection time for a single-cluster avalanche" of about 150 ns. — [arXiv:2404.08317](https://arxiv.org/pdf/2404.08317) (search snippet only; gap and gas not confirmed)
- Kalkan et al. 2015 explicitly expect cluster-ion formation (slower ions) to alter Micromegas ion signals in Ar/CO2. — [Kalkan et al. 2015](https://avesis.uludag.edu.tr/yayin/6a020f55-18d2-48e6-bfed-9e6fd91abe0e/cluster-ions-in-gas-based-detectors)
- Zerguerras et al. (NIM A 772, 2015) studied Micromegas single-electron response with a laser in **Ar/iC4H10 95/5** (and Ne/He with 5% iso), the MX17 gas. The full text was blocked, so it is unknown whether they report an ion-tail duration. — [CDS](https://cds.cern.ch/record/2153708/files/nima772-76.pdf)

### Inferences
- **Checking Giomataris's numbers against K** (t = d²/(K V)):
  - Fig. 7, if it is a 100 µm gap at ~400–450 V: K = d²/(V t) = 1e-4/(420 × 80e-9) ≈ 3.0 cm²/Vs. That is above Langevin, so the gap was probably smaller or the "80 ns" is not a full transit.
  - At 50 µm and ~350 V it would be K ≈ 0.9.
  - **The figure cannot be used quantitatively without the gap and voltage** (inference).
- **Rise-time scaling**: going from 50 to 30 µm gives ~3× in d² at modestly lower voltage, consistent with the d²/V scaling, so the data are self-consistent but not mobility-calibrating.
- The "~100 ns" NSW statement for 128 µm Ar/CO2 93/7 at ~550 V would imply K ≈ (0.0128)²/(550 × 100e-9) ≈ 3.0 cm²/Vs, again above any Langevin value. These round-number statements are best read as order-of-magnitude, or as the time for the bulk of the induced charge (inference).
- **No published Micromegas waveform fit of the ion mobility in Ar/iso was found.** MX17's own fit (T_ion ≈ 140 ns at the ~166 ns shaper) is therefore new information, and it is degenerate with the shaper model and the ion fraction f_eff, as HANDOFF §10.2 notes.

### Gaps
- No T2K ND280, ILD/LCTPC, CAST or MIMAC paper reporting an ion-tail-derived mobility was located. Their electronics (AFTER/AGET, shaping ≥100 ns) are mostly too slow to resolve the ion tail.
- Bortfeldt (floating-strip Micromegas) and NSW signal-shape papers were not retrieved.

---

## Q5. At 5% iso, does any Ar⁺ survive? Ar⁺ + iC4H10 rate and time constant at 1 atm; terminal ion

### Takeaway
- No Ar⁺ survives to matter. The Langevin capture estimate is k ≈ 1.4e-9 cm³/s, consistent with the measured Ar⁺ + C2H6 rate of 1.15e-9 in Anicich. That gives **τ(Ar⁺) ≈ 0.6 ns at 5% iso and 1 atm**, during which Ar⁺ drifts only ~0.3 µm.
- Ar2⁺ formation is slower (τ ≈ 7 ns), and Ar2⁺ is also destroyed by iso.
- The drifting ion is a hydrocarbon ion from the start, most plausibly ending as t-C4H9⁺ (protonated isobutene), the most stable carbocation in the system.
- **Using the Ar⁺/Ar table for the entire transit is physically wrong for this gas** (inference).

### Cited Findings
- The Anicich 1993 evaluation lists Ar⁺ + C2H6 → products at 1.15e-9 cm³/s and Ar⁺ + H2O at 1.62e-9 cm³/s. **No Ar⁺ + iC4H10 entry** exists in that compilation. — [Anicich 1993](https://srd.nist.gov/jpcrdreprint/1.555940.pdf)
- In argon–CO2, noble-gas ions transfer their charge to the molecular admixture "within a few nanoseconds" at 1 atm. — [Kalkan et al. 2015](https://avesis.uludag.edu.tr/yayin/6a020f55-18d2-48e6-bfed-9e6fd91abe0e/cluster-ions-in-gas-based-detectors)
- Coimbra's 90% Ar/10% CF4 data show that a charge-transfer product (CF3⁺) replaces Ar⁺ and Ar2⁺ entirely. — [Cortez LCTPC slides](https://agenda.linearcollider.org/event/7795/contributions/40342/attachments/32468/49376/ACortez_-_LCTPC_Collaboration_Meeting_-_November_2017.pdf)
- In isobutane ion chemistry, C3H7⁺ is the most abundant fragment at most ionising energies. Fragment ions such as C4H7⁺ and smaller react readily with the parent, mainly producing C5H11⁺ and C4H9⁺. — search snippet ([ResearchGate, "Electron impact ionization and ion reactions in n-butane"](https://www.researchgate.net/publication/230927357_Electron_impact_ionization_and_ion_reactions_in_n-butane); not read in full, and the chemistry may refer to n-butane)
- Isobutene PA is 802 kJ/mol (t-C4H9⁺ is protonated isobutene). — [ACS Figshare](https://acs.figshare.com/collections/Heat_of_Formation_of_the_i_tert_i_-Butyl_Radical/3416628)

**Calculation (explicit):**
- **Langevin capture**: k_L = 2.342e-9 · √(α[Å³]/μ[amu]) cm³/s. For Ar⁺ + iC4H10, α = 8.14 Å³ (flag: from memory) and μ = 39.948 × 58.12/98.07 = 23.67 amu, so **k_L = 2.342e-9 × √0.3439 = 1.37e-9 cm³/s**.
- n_iso = 0.05 × 2.503e19 = 1.25e18 cm⁻³, so ν = k n = 1.72e9 s⁻¹ and **τ = 0.58 ns**.
- **Distance drifted by Ar⁺ in that time**: K ≈ 1.24 cm²/Vs at 38.3 kV/cm gives v ≈ 4.7e4 cm/s, so 0.58 ns × v ≈ 0.27 µm.
- **Survival fraction after the full transit**: exp(−270 ns / 0.58 ns) ≈ e⁻⁴⁶⁰, effectively zero.
- **Competing Ar2⁺ formation**: k3 ≈ 2.5e-31 cm⁶/s (typical literature value, from memory; flag) × (0.95 N)² gives τ ≈ 7 ns, so iso wins by ~12×.
- **Fraction of Ar⁺ that goes to H2O instead**: 0.7% at 300 ppm, 2.3% at 1000 ppm, 6.6% at 3000 ppm (see Q2).

### Inferences
- **Hydrocarbon ion sequence** (inference from energetics; specific rate constants for iso were not found):
  1. Ar⁺ (15.76 eV) dissociatively charge-transfers to iso within <1 ns, giving C3H7⁺, C3H6⁺, C3H5⁺, C2H5⁺ and similar fragments.
  2. The fragments react further with iso at near-capture rates (~1e-9 cm³/s, so again ~ns at 5%), funnelling toward t-C4H9⁺ (m = 57).
- So for ~99% of the ~155–270 ns transit the drifting ion has mass 43–57, and possibly a cluster.
- **Simulation recommendation (inference):**
  - Replace `IonMobility_Ar+_Ar.txt` with a single-species table for a mass-57 ion in 95/5. No measured high-E/N table exists.
  - A defensible placeholder is a flat **K0 ≈ 1.8–2.0 cm²/Vs**: Langevin 2.23, reduced by Blanc with K_iso ≈ 0.65–0.9 to about 1.9–2.0, then −0–10% for high-field and short-range effects.
  - Bracket it with 1.6 (dimer or partly clustered, or wet) and 2.2 (Langevin ceiling).
  - Expected transit at 490 V/128 µm/1 atm/293 K: **~165–185 ns (bracket 140–200 ns)**.
- The fitted 140 ns therefore sits at the edge of the physically allowed range. If the fit holds, part of the waveform discrepancy should be sought elsewhere:
  - shaper response;
  - the ion signal fraction f_eff;
  - field non-uniformity near the mesh, where ions start ~0–20 µm above the anode and the field near the mesh holes is higher;
  - the gas density assumption.

### Gaps
- No measured rate constant or product distribution for Ar⁺ + iC4H10, Ar2⁺ + iC4H10, or C3H7⁺ + iC4H10 → t-C4H9⁺ + C3H8 was retrieved. Anicich 2003 (JPL Publication 03-19) or NIST kinetics would be the next places to look.
- The three-body association rates of t-C4H9⁺ with iC4H10 or H2O at 1 atm, and their reversibility at T_eff ~1000+ K, are unknown.
