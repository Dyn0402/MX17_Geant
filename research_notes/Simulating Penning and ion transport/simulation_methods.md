# Simulation methods for coupled electron transport, excited-state kinetics and ion–molecule chemistry in weakly ionized gases (context for Magboltz/Garfield++)

Notes compiled 2026-10-15. Searches were web-only (about 18 tool calls). Most findings rest on abstracts, arXiv landing pages and conference abstracts, not full texts; this is flagged where it matters. Items marked "[prior knowledge, not re-verified]" come from the researcher's background knowledge, with a DOI given but not fetched this session. Treat them as leads.

## 1. Electron transport solvers (Boltzmann and Monte Carlo) and how they handle excited states, superelastic collisions and Penning transfer

### Takeaway
Electron-swarm solvers are mature and benchmarked against each other. Two-term Boltzmann solvers (BOLSIG+, LoKI-B) are fast but bias diffusion coefficients. Monte Carlo codes (Magboltz/PyBoltz, METHES, Betaboltz, ThunderBoltz) are exact for the cross sections they are given. The plasma-physics codes differ from Magboltz in one way: they accept user-specified excited-state populations, which gives superelastic collisions, and they couple iteratively to a chemistry solver (LoKI-B+C). Magboltz/Garfield++ treat Penning transfer as one ad hoc transfer probability, with no excited-state population or kinetics.

### Cited findings
- **LoKI-B+C (IST Lisbon)** couples a two-term Boltzmann solver (LoKI-B) with a chemistry solver (LoKI-C). The coupling runs nested iterations, one of them over "the densities of the most relevant excited states affecting the electron Boltzmann equation", until the EEDF converges globally. LoKI-B is open source under GPL v3.0. — [GEC2025 LoKI-B+C poster](https://nprime.tecnico.ulisboa.pt/psi.com/img/Publications/poster_loki-B+C_GEC2025.pdf); [IPFN news: LoKI-B released open-source](https://www.ipfn.tecnico.ulisboa.pt/news-and-events/news/281663955271717)
- In 2025, LoKI-B+C was applied self-consistently to a low-pressure Ar ICP (0.2–1.5 Pa, ne ~ 1e16 m^-3) with a combined stochastic plus Ohmic Fokker–Planck heating operator. It outputs the EEDF, the densities of the main Ar species, E/N, the sheath potential and the coupled power, and was compared with probe EEDFs. The group lists time-dependent B+C coupling as work still in progress. — [APS GEC 2025 abstract DT4.71](https://archive.aps.org/gec/2025/dt4/71); [GEC2025 poster](https://nprime.tecnico.ulisboa.pt/psi.com/img/Publications/poster_loki-B+C_GEC2025.pdf)
- **LoKI-MC (Lisbon, 2021):** a Monte Carlo code built on the LoKI-B approach. It handles "electron-neutral collisions with any target state (electronic, vibrational and rotational), characterized by any user-prescribed population, allowing to include the effects of superelastic collisions in a general manner." — [GEC2021 LoKI-MC abstract](https://nprime.tecnico.ulisboa.pt/loki/img/publications/GEC2021_abstractLoKI-MC.pdf)
- **ThunderBoltz (LANL; Park, Scheiner, Zammit; arXiv Oct 2023; published in PSST, DOI 10.1088/1361-6595/ad6fcc):** an open-source 0D DSMC code for "a generalized combination of species and arbitrary cross sections". It evolves electrons, ions and background neutrals self-consistently, supports particle-based chemical kinetics and rate-constant post-processing, and has a Python interface to LXCat data. It was benchmarked against BOLSIG+ and swarm data. The abstract does not mention superelastic collisions. — [arXiv:2310.07913](https://arxiv.org/abs/2310.07913)
- **PyBoltz** is a Python/Cython refactoring of Magboltz. It matches the original at comparable speed and adds a Modified Effective Range Theory interface. Magboltz is called "the industry standard", with about 20 years of development. — [arXiv:1910.06983](https://ar5iv.arxiv.org/html/1910.06983)
- **Betaboltz** is a C++17 MC code for electron and ion transport in arbitrary mixtures with static E and B fields. Its drift velocities matched literature well only for atomic and non-polar, spherically symmetric molecular gases. — [Betaboltz summary, via arXiv:1901.08140 / AIP Advances 2019](https://doaj.org/article/a27e8e59f8ae40c1b09ed02caa32ec73)
- **METHES** is a MATLAB MC code for steady-state electron transport in uniform fields. It reads LXCat files, outputs transport coefficients, rates and the EEDF, and was benchmarked on Ar, N2, O2, CF4, SF6 and N2/O2. — [METHES on Mendeley Data](https://data.mendeley.com/datasets/rcp4ssbc23); [pymethes docs](https://pymethes.readthedocs.io/en/stable/search.html)
- **Two-term vs MC accuracy:** for Ar, two-term Boltzmann solvers overestimate diffusion coefficients by about 30% in the 5–100 Td range compared with Monte Carlo. — [GEC 2011 QRP1.69 (Pitchford et al.)](https://meetings-archive.aps.org/gec/2011/qrp1/69/)
- **Penning in Garfield++:** Şahin et al. (JINST 5 P05002, 2010) fitted Penning transfer probabilities to published gain curves for Ar-based mixtures. They found the probability depends on composition and pressure. — [Penning transfer in argon-based gas mixtures](https://acikerisim.uludag.edu.tr/items/0c97c615-5554-4cf7-b343-a1557f6ca0ae)
- Without Penning transfer, Garfield++ underestimated thick-GEM gains by about 2 orders of magnitude. Including it (for example r = 0.18 for Ar/5%CH4) reproduced data for Ne/5%CH4, Ar/5%CH4 and Ar/30%CO2. The authors note that the Penning option "is not often used … since the data in Penning fractions is sparse and not reliable." — [Azevedo et al., arXiv:1606.04852](https://arxiv.org/pdf/1606.04852)
- A Fermilab talk reports that Garfield's Penning implementation did not give the expected results, so the author modified Magboltz to change the transfer coefficient. It quotes r ≈ 0.29 for Ar/10%CH4 at 10 atm from Şahin. — [Korb, FNAL Indico](https://indico.fnal.gov/event/44260/contributions/190497/attachments/130861/159799/GarfieldGain.pdf)

### Inferences
- In plasma codes, excited-state densities are dynamic variables that feed back into the electron kinetics. In Garfield++, an excitation is a terminal event with a fixed probability of converting to ionization. The tools to do it properly exist (LoKI-MC populations, LoKI-B+C loops, ThunderBoltz species chemistry). What is missing is input data for Ar* + quencher channels, not algorithms.
- In a detector avalanche, the excited-state density is tiny compared with the neutral density, so superelastic feedback is negligible. The real gap is the branching of each Ar* (Penning vs. dissociation vs. radiation vs. collisional quenching) and its time and space delay. That is a chemistry and data problem, not a transport one.

### Gaps
- No public documentation of MultiBolt's features was retrieved. It appears only in lists of open swarm codes.
- I could not confirm from the abstract whether ThunderBoltz implements superelastic collisions.
- BOLSIG+ (Hagelaar & Pitchford, PSST 14, 722, 2005, DOI 10.1088/0963-0252/14/4/011) can take excited-state populations and superelastic collisions as input [prior knowledge, not re-verified].

## 2. Plasma chemistry / global models (ZDPlasKin, LoKI-C, QDB, PLASIMO, CHEMKIN-type) and their treatment of Ar metastables, Penning ionization, radiation trapping and three-body processes

### Takeaway
0D chemistry codes solve species-density rate equations. Electron-impact rates come from an embedded Boltzmann solver, and heavy-particle rates come from Arrhenius or tabulated data. Ar metastable kinetics are routine in them: stepwise ionization, Penning ionization, and resonance-state radiation trapping via escape factors. The physics is simple once the rate coefficients are known. Total quenching rates of Ar* by hydrocarbons are generally available; product branching ratios are not.

### Cited findings
- **ZDPlasKin** (Pancheshnyi, Eismann, Hagelaar, Pitchford; GEC 2008) is freeware. A preprocessor turns a text list of species, reactions and rates into a Fortran 90 module with a DVODE ODE solver. It links automatically to BOLSIG+: the EEDF is computed at a given E/N and integrated over cross sections to get electron-impact rates. Heavy-particle rates are Arrhenius or tabulated. — [GEC 2008 abstract](https://meetings-archive.aps.org/gec/2008/mwp1/60)
- **Successors:** a 2024 Julia-compiled ZDPlasKin coupled to BOLSIG+ and Cantera with operator splitting ([GEC 2024 HW6.10](https://meetings-archive.aps.org/gec/2024/hw6/10)); ChemPlasKin, which couples CppBOLOS with Cantera ([arXiv:2405.04224](https://arxiv.org/html/2405.04224v1)); and RePlaChem, a mechanism-reduction library ([arXiv:2507.04401](https://arxiv.org/pdf/2507.04401)).
- **QDB (Quantemol):** the 2021 release holds 28,917 processes among 2,485 species, plus surface processes. It exports to HPEM, COMSOL, ChemKIN, CFD-ACE+ and VisGlow, and includes a Boltzmann solver and a 0D model. Unlike LXCat, it aims to *recommend* datasets for an application. — [QDB 2017, arXiv:1704.04088](https://arxiv.org/pdf/1704.04088); [2021 release, PSST 31 095020](https://pure.york.ac.uk/portal/en/publications/the-2021-release-of-the-quantemol-database-qdb-of-plasma-chemistr/)
- **Ar metastables in global models:** metastables contribute to ionization through stepwise and Penning routes, and the contribution grows with plasma density ([Plasma 5(1):3, MDPI 2022](https://mdpi.com/2571-6182/5/1/3)). A global model for the Ar metastable density during the E–H transition needed a non-Maxwellian EEDF to reproduce a measured trend ([KAIST](https://koasas.kaist.ac.kr/handle/10203/99582)). An afterglow global model for Ar* matched data once electron energy loss to walls was included ([HAL-01005941](https://hal.archives-ouvertes.fr/hal-01005941)).
- **Radiation trapping:** in an Ar ICP, the escape factor of resonance radiation is set by pressure broadening of the line wings, and models often lump the 4s states into one effective lifetime. — [DTIC ADA136633](https://apps.dtic.mil/sti/pdfs/ADA136633.pdf)

### Inferences
- All the parts a detector "Penning + chemistry" model would need already exist in these codes as standard features: Ar(4s) resonant and metastable states with escape factors, Ar* + M Penning and quenching with rate constants, and three-body excimer formation (Ar2*). A detector-gas chemistry set for Ar/iC4H10 could be written in ZDPlasKin syntax in a day. The rate data for the products, and their dependence on the quencher's internal energy, could not.
- Global models assume spatial uniformity and timescales of μs to ms. An avalanche lasts about 1 ns over about 100 μm with steep gradients, so 0D chemistry fits only as a per-cell or per-track post-processor, not as a direct replacement.

### Gaps
- No source was retrieved on PLASIMO or on LoKI-C's specific Ar metastable kinetic scheme.
- No published Ar/hydrocarbon (isobutane, CH4) plasma-chemistry reaction set with Penning branching ratios was found in these searches.
- Holstein–Biberman theory proper was not retrieved. Only an escape-factor application was found.

## 3. PIC/MCC, fluid/hybrid and streamer codes: could they be applied to avalanches or streamers in detectors?

### Takeaway
Yes, in principle, and partly in practice for RPCs. Streamer codes (Afivo-streamer fluid; PIC-MCC; Itô/kinetic MC) already model space charge plus photoionization and agree well with each other. The detector community has its own GPU and space-charge microscopic MC codes and RPC fluid models. The limit is not the field solver or the particle push. It is the lack of photoionization and photon-feedback data for detector gas mixtures.

### Cited findings
- **Afivo-streamer** (Teunissen & Ebert, J. Phys. D 2017) is an open-source plasma fluid code with adaptive mesh refinement and geometric multigrid Poisson solving. It is used for 3D streamer branching. — [GEC 2018 TF1.8](https://meetings-archive.aps.org/gec/2018/tf1/8)
- **Photoionization in streamer codes:** Zhelezniak's model, solved either in the Helmholtz approximation (three-term expansion) or as a Monte Carlo method with discrete photons. Stochastic photoionization changes branching, and gives more radio emission above 100 MHz. — [arXiv:2401.12353](https://arxiv.org/pdf/2401.12353); [Malla et al. 2024, arXiv:2404.09772](https://arxiv.org/abs/2404.09772v2)
- **PIC-MCC vs fluid** (Wang, Sun, Teunissen; PSST 31(1)): for positive streamers in air, mean velocities agree to within 4%. Flux transport coefficients match the particle results better than bulk ones. Photoionization is the main source of stochastic fluctuations. — [CWI](https://ir.cwi.nl/pub/31317/); [CWI](https://ir.cwi.nl/pub/31318/)
- **Kinetic Monte Carlo / Itô model** (Marskar): PIC-style electrons with Itô diffusion, plus particle photoionization, using fluid-model inputs. — [arXiv:2307.01797](https://ar5iv.labs.arxiv.org/html/2307.01797)
- **RPC streamer modelling:** a 1.5D fluid model with photoionization (Dujko, Bošnjaković, Petrović) found that electron absorption at the anode strongly affects space charge and positive-streamer formation ([GEC 2020 KT3.1](https://archive.aps.org/gec/2020/kt3/1)). A COMSOL RPC streamer model used Magboltz transport data ([arXiv:2005.13911](https://arxiv.org/pdf/2005.13911)). A 2024 paper states that "accurate photoionization models for complex RPC gas mixtures are not yet available" and so assumes background ionization ([MOTIONER 2024](https://aob.rs/images/doi/konferencije/2024/PDF/2024-motioner/35-37.pdf)).
- **Detector-side particle codes:** OuroborosBEM is a multi-GPU microscopic MC with space charge and recombination, which Garfield++ lacks ([arXiv:2110.09214](https://arxiv.org/abs/2110.09214)). MATOQ is an MC for RPC eco-gases that includes the space-charge field ([arXiv:2304.10307](https://arxiv.org/pdf/2304.10307)). Garfield++ was parallelized with neBEM for RPC space charge ([arXiv:2211.06361](https://ar5iv.arxiv.org/html/2211.06361)).

### Inferences
- Large PIC frameworks (XOOPIC, EDIPIC, PICLas, BIT1, VSim, WarpX-MCC) are built for bounded plasmas with self-consistent fields. A sub-ns avalanche with 1e4–1e7 electrons is small enough for microscopic MC (Garfield++) or a streamer fluid code. Gaseous detectors do not need PIC for its own sake; space charge only matters near the streamer and saturation regime.
- The streamer codes' photoionization models are parametrized for air (N2 to O2). Mixtures such as Ar/iC4H10, Ar/CO2 or RPC freon gases lack absorption and ionization-yield data. This is the same data gap that blocks feedback modelling in Garfield++.

### Gaps
- No PIC study of an MPGD avalanche was found, and none using the named PIC codes (XOOPIC, EDIPIC, PICLas, BIT1, VSim, WarpX) for a detector. I am not aware of any such study.

## 4. Cross-section and rate-data infrastructure (LXCat, QDB, NIST, IAEA) and validation

### Takeaway
LXCat is the community hub. It hosts many sets (Biagi = Magboltz transcription, Phelps, Hayashi, IST-Lisbon, Morgan, BSR, Viehland ions). It does not rank sets; users compare them. The noble gases were systematically cross-validated against swarm data around 2013, and Ar sets agree with each other to within modelling accuracy. LXCat 3 (versioned JSON, chemistry mechanisms) was planned for 2024 release, and I could not confirm its launch. Molecular detector quenchers (isobutane, etc.) are covered almost only by Biagi's swarm-adjusted sets.

### Cited findings
- **Pitchford et al.**, "Comparisons of sets of electron–neutral scattering cross sections and swarm parameters in noble gases: I. Argon", J. Phys. D 46 (2013) 334001. It compared nine independent Ar sets on LXCat via Boltzmann solvers against measured swarm data. — [summary of search result; APS GEC 2011 precursor](https://meetings-archive.aps.org/gec/2011/qrp1/69/)
- The swarm parameters calculated from the various Ar sets agreed with experiment and with each other to within the accuracy needed for plasma modelling. How cross sections are interpolated can shift swarm parameters by about 6% at intermediate E/N. — [GEC 2011 QRP1.69](https://meetings-archive.aps.org/gec/2011/qrp1/69/); [GEC 2011 QRP1.64](https://meetings-archive.aps.org/gec/2011/qrp1/64/)
- LXCat "make[s] no judgment about the relative merit" of its datasets. More than 60 contributors had contributed by 2024. — [LXCat status 2024, GEC HT4.3](https://meetings-archive.aps.org/gec/2024/ht4/3); [Pitchford IAEA DCN23 presentation](https://amdis.iaea.org/media/presentations/DCN23-Pitchford-presentation.pdf)
- **LXCat 3:** a JSON-Schema data model with ArangoDB, versioned data and permalinks, and support for plasma-chemistry mechanisms. Public release was "scheduled" for Sept/Oct 2024. A 2025 abstract describes it as supporting both state-to-state data and "complete cross section sets". — [GEC 2024 EM2.14](https://meetings-archive.aps.org/gec/2024/em2/14); [ESCAMPIG 2024 abstract](https://nprime.tecnico.ulisboa.pt/psi.com/img/Publications/abstract_lxcat_ESCAMPIG2024.pdf); [GEC 2025 GF2.1](https://archive.aps.org/gec/2025/gf2/1)
- **QDB** recommends datasets, unlike LXCat (see section 2). — [arXiv:1704.04088](https://arxiv.org/pdf/1704.04088)
- A 2021 comparison of the Hayashi, Triniti, Morgan and BSR Ar sets focused on the electron mean energy. — [Acta Sci. Biol. (asb.pmf.unibl.org)](https://asb.pmf.unibl.org/index.php/asb/article/view/29)

### Inferences
- Detector-gas modelling shares the LXCat ecosystem: Magboltz's data appear there as the "Biagi" database. Its blind spot is the same as plasma physics': validation covers electron swarms in the ground-state gas. There is no equivalent validated compilation for heavy-particle Ar* + M branching or for ion mobilities in detector mixtures (Viehland covers mainly atomic ions).

### Gaps
- LXCat 3 was not confirmed live as of 2026. Check lxcat.net.
- Part II and later of the Pitchford noble-gas comparisons (He, Ne, Kr, Xe), and any comparable exercise for detector quenchers, were not retrieved.
- NIST and IAEA (e.g. ALADDIN, CollisionDB) holdings were not searched in detail.

## 5. First-principles cross-section theory: what is computable (Ar excitation) vs hard (polyatomics, Penning dynamics)

### Takeaway
Electron–Ar excitation is a solved problem for theory. B-spline R-matrix (BSR) cross sections agree with differential measurements near threshold and are in LXCat. Electron–polyatomic scattering (isobutane) is computed only at modest levels: static-exchange Schwinger multichannel, elastic only, ≥10 eV. Penning ionization dynamics is computed ab initio, via complex optical potentials V(R) − iΓ(R)/2, only for atom–atom or atom–diatom systems (He*/Ne* + Ar, H2, HD). Merged-beam experiments benchmark those calculations at mK energies. Nothing comparable exists for Ar* + hydrocarbons.

### Cited findings
- BSR (B-spline R-matrix with non-orthogonal, term-dependent orbitals) predictions for Ar near-threshold differential excitation at 135° agree well with measurement (Allan, Zatsarinny, Bartschat, PRA 74, 030701, 2006). The BSR sets were used to predict swarm parameters with Boltzmann solvers, and are on LXCat. — [NIST record](https://www.nist.gov/node/1589816); [GEC 2014 MW1.4](https://meetings-archive.aps.org/gec/2014/mw1/4); [Bartschat NIST slides](https://bartschat.wp.drake.edu/wp-content/uploads/sites/22/2019/12/Bartschat-NIST-slides.pdf)
- The BSR method has also been applied to Be (with CCC), N, Kr and Yb. For Be, coupling to the target continuum matters most for ground-state transitions. — [arXiv:1606.01326](https://arxiv.org/pdf/1606.01326)
- **Isobutane:** Lopes, Bettega, Lima, Ferreira (J. Phys. B, 2004) used the Schwinger multichannel method in static-exchange for *elastic* scattering from isobutane and n-butane at 10–50 eV. — [UChile repository record](https://repositorioslatinoamericanos.uchile.cl/handle/2250/1242182)
- **Penning theory (optical potential):** the standard model is W(R) = V(R) + (i/2)Γ(R), with the real part from ab initio RCCSD(T), e.g. for He(3P)–Ar. — [Skomorowski et al., arXiv:1601.00433 / JPCA 120, 3309 (2016)](https://arxiv.org/pdf/1601.00433)
- Ne*(3P2) + Ar Penning and associative ionization were modelled with combined empirical and ab initio ionization widths in a "Rotating Atom Approximation" ([arXiv:physics/0610131](https://arxiv.org/pdf/physics/0610131)). Such calculations reproduce experimental product ratios down to about 0.02 K ([arXiv:2101.03486](https://arxiv.org/pdf/2101.03486)).
- **Merged-beam experiments (Narevicius group):** Henson et al., Science 2012, reached collision energies down to 10 mK and saw orbiting resonances in He* + Ar and He* + H2 Penning ionization. A He* + ortho-H2 shape resonance near 270 mK benchmarks quantum-dynamics and electron-correlation treatments ([arXiv:1606.04384](https://arxiv.org/pdf/1606.04384)). A Ne* + HD p-wave resonance near 22 mK is sensitive to the van der Waals C6 term ([arXiv:2201.07716](https://arxiv.org/abs/2201.07716v1)).
- Even for these benchmark systems, "the highly excited entrance channel and the limited precision of ab initio potential surfaces remain major obstacles." — [Skomorowski et al. 2016](https://physik.fu-berlin.de/en/einrichtungen/ag/ag-koch/Publications/paper-pdf/2016_Photoassociation-Spectroscopy_JPhysChemA_120_3309_2016.pdf)
- For Ne*–Ne*, ab initio well depths disagree with experiment because of spin–orbit coupling, so relativistic effects matter. — [arXiv:1904.05672](https://arxiv.org/pdf/1904.05672)

### Inferences
- The hardness ladder for theory runs: e + Ar excitation (solved; BSR/CCC) < e + small molecule elastic/vibrational (R-matrix/SMC fine) < e + polyatomic electronic excitation and dissociation (SMC/R-matrix with many open channels; isobutane has 14 atoms and dense, dissociative Rydberg/valence states) < heavy-particle Ar* + polyatomic Penning branching. The last needs Γ(R) over a many-dimensional surface, many spin–orbit-coupled entrance channels (Ar 4s/4p with J resolution), non-adiabatic coupling to neutral dissociation (the competing channel that sets the Penning fraction), and a final state in the ionization continuum.
- The thermal-energy Penning *fraction* relevant to detectors (about 0.1–0.5) is the ratio of Penning to total quenching. Accurate theory has only been demonstrated for atom/diatom targets at sub-K energies. I found no indication that ab initio work extends to Ar* + C4H10.

### Gaps
- The main Zatsarinny–Bartschat Ar paper (PRA/J. Phys. B, about 2006–2013) was not retrieved directly.
- No ab initio Penning calculation for Ar* + any hydrocarbon was found. Hickman, Miller and Gordon's classic optical-potential work, and Osterwalder's Penning electron spectroscopy studies, were not retrieved this session.
- No UK R-matrix or Quantemol-N calculation for isobutane electronic excitation was found.

## 6. Machine learning and data-driven approaches (inverse swarm problem)

### Takeaway
Neural-network inversion of swarm data into cross sections has existed since about 2019–2020 (JCU/ANU/Flinders: Stokes, Cocks, Brunger, White). It works well for He, poorly for Ar (well understood), and has been applied to biomolecular gases (THF). Later variants exist: iterative and multi-branch networks (2024), PINN Boltzmann solvers, and DeepCSNet. I found no application to detector gases (isobutane, CF4 mixtures) or to Penning-transfer coefficients.

### Cited findings
- Stokes, Cocks, Brunger, White, PSST 29 055009 (2020). The network was trained on LXCat cross sections paired with Boltzmann-computed transport coefficients. He momentum-transfer and ionization were recovered to about 4%. Simultaneous He elastic, n=2 excitation and ionization were recovered to about 10%, 20% and 25%, and singlet and triplet n=2 could not be separated. The Ar extension "was less successful" for well-understood reasons. — [arXiv:1912.05842](https://arxiv.org/abs/1912.05842v2); [JCU record](https://researchonline.jcu.edu.au/63563/)
- The approach was applied to electron–THF cross sections from pulsed-Townsend data, refining the dissociation and attachment channels. — [arXiv:2007.02762](https://arxiv.org/abs/2007.02762v1); [TIB](https://oa.tib.eu/renate/items/9093842b-ea43-4fd6-b4c5-87f03088cce2/full)
- Iterative deep-learning procedure (JCU). — [JCU 86023](https://researchonline.jcu.edu.au/86023/)
- Muccignat et al. (2024, Machine Learning: Science and Technology) use a multi-branch network, one branch per cross section. They stress that the inverse problem is non-unique. — [DOAJ](https://doaj.org/article/efa870389a434aacb0ca0e4bc17e2336)
- DeepCSNet (Wang & Zhong 2024, PSST) predicts electron-impact doubly differential ionization cross sections. — [arXiv:2407.09778](https://arxiv.org/html/2407.09778v2)
- PINNs have been used for the electron Boltzmann equation (Kawaguchi et al., SPIG 2024). — [SPIG 2024 proc.](https://aob.rs/images/doi/konferencije/2024/PDF/2024-32nd-SPIG/193.pdf)

### Inferences
- Inverse swarm methods constrain only what swarm data are sensitive to (momentum transfer, total inelastic, ionization, attachment). The Penning fraction is a property of the heavy-particle chemistry, not of the electron cross sections. Gain curves constrain only the product (excitation × r_Penning). This is the same fitting Şahin et al. did by hand, and ML adds degeneracy-breaking only if extra observables (pressure scaling, scintillation yields, time structure) are included.

### Gaps
- No 2025–2026 papers applying ML to detector gases or Penning fractions were found.

## 7. Work bridging plasma-chemistry / excited-state physics and gaseous detectors

### Takeaway
The closest bridges come from the electroluminescence (EL) TPC community: Coimbra/Aveiro (Oliveira, Santos, Veloso) plus Schindler/Veenhof/Biagi, and NEXT. Garfield++ records each excitation's level, position and time, and post-processing converts them to VUV photons (assuming one photon per excitation, or excimers via three-body formation). Penning-fluorescent Xe-TMA chemistry was explicitly declared out of scope of the microscopic simulations (2017). Experiments measured about 35% Penning transfer and about 3% fluorescence transfer, but no time-dependent kinetic model was found.

### Cited findings
- **Oliveira et al. (arXiv:1103.6237)** built a toolkit on Magboltz and Garfield that computes light emission from electron impact in Ne, Ar, Kr and Xe, giving excitation and EL efficiencies vs field. Co-authors include Schindler, Veenhof, Biagi, dos Santos, Ferreira and Veloso. — [arXiv:1103.6237](https://ar5iv.arxiv.org/html/1103.6237)
- For VUV EL in GEM/MHSP with Garfield++ (arXiv:1206.1646), Xe excited levels are parametrized in 50 groups. Each excitation's position, time, level and energy are retrievable, and each excited atom is assumed to emit one isotropic VUV photon. — [arXiv:1206.1646](https://ar5iv.arxiv.org/html/1206.1646)
- **NEXT energy resolution** (arXiv:1105.2954) used microscopic Garfield++ with Magboltz 8.9.3 cross sections for EL in a 5 mm Xe gap at 10 bar. The EL mechanism is excimer formation via three-body collisions. — [arXiv:1105.2954](https://arxiv.org/pdf/1105.2954)
- **Molecular additives in Xe optical TPCs** (Azevedo et al. 2017, arXiv:1705.09481): modelling "akin to the Penning-Fluorescent Xe-TMA TPC falls outside the considerations of the present work, but will be performed in the future." — [arXiv:1705.09481](https://arxiv.org/pdf/1705.09481)
- **Xe-TMA experiments:** Penning transfer up to about 35%, about 3% fluorescence transfer to TMA (about 300 nm), and strong suppression of primary scintillation. — [arXiv:1511.02257](https://www.arxiv.org/pdf/1511.02257v1)
- NEXT-MM measured Xe-TMA gain, Penning effects and diffusion. — [ICHEP 2014 slides](https://indico.ific.uv.es/event/2025/contributions/1186/attachments/1142/1283/tdafni_mm_ichep2014.pdf)
- Additives (CO2, CH4, CF4) reduce diffusion but lower the EL yield. — [JHEP 2020 (DOAJ)](https://doaj.org/article/3ae2750321094dc3a67e7a1caa41d58c)

### Inferences
- The detector community's "excited-state" modelling is bookkeeping on Garfield++/Degrad excitation records plus constant transfer or branching factors. It is not a coupled kinetic model with rate equations, pressure-dependent three-body excimer formation, quenching and radiation transport in one framework. The exceptions are analytic excimer-kinetics treatments in EL studies, not covered here.
- Why Garfield++ stays simple (synthesis): (i) the timescales separate, since avalanche electrons see a negligible excited-state density and the chemistry is a downstream branching problem; (ii) the missing inputs are Ar*/Xe* + polyatomic branching and state-resolved quenching rates, which theory cannot yet compute (section 5) and which swarm inversion cannot isolate (section 6); (iii) a single fitted r_Penning per mixture and pressure already reproduces gain curves (section 1). The extra fidelity of LoKI-B+C/ZDPlasKin-style models would be limited by the same missing rates.

### Gaps
- No paper was found that couples Garfield++/Degrad to a chemistry solver (ZDPlasKin/LoKI-C/Cantera) for Penning or excimer kinetics in a detector avalanche.
- No Ar/CF4 time-dependent scintillation-kinetics simulation was retrieved. ArDM-specific EL modelling was not searched.
- Escada et al. (2011) Monte Carlo of Xe EL with CH4/CF4 dopants was cited but not retrieved.
