#!/usr/bin/env python3
"""Slide note: the two candidate mechanisms behind the T14 sim/data gap (photon feedback,
hydrocarbon ions), the archive evidence for each, and what the literature says.

Reads this directory's outputs (ladders_all.csv, fits_all.csv, width_vs_v.log,
det0_toy_scan.log, det0_curves.npz) plus Garfield's ion-mobility file, so a rerun of the
analysis updates the deck. Run with the nTof_x17 venv:

    ~/PycharmProjects/nTof_x17/.venv/bin/python make_lit_deck.py
    python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py out/mx17-feedback-and-ions.html \
        --slug mx17-feedback-and-ions --force --deploy
"""
import os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser(os.environ.get('SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd

HERE = os.path.dirname(os.path.abspath(__file__))
L = pd.read_csv(f'{HERE}/ladders_all.csv')
F = pd.read_csv(f'{HERE}/fits_all.csv')
C = np.load(f'{HERE}/det0_curves.npz')

# ------------------------------------------------------------------ numbers from the analysis
def frange(s, col):
    v = F[F.set == s][col]
    return v.min(), v.max()
RATIO = {s: frange(s, 'data_over_garf') for s in F.set.unique()}
VSTAR = {s: frange(s, 'Vstar_gfb') for s in F.set.unique()}
ARISO = ['det3 (June, T14) Ar/iso', 'det3o Ar/iso', 'det0 Ar/iso']
ar_lo = min(RATIO[s][0] for s in ARISO); ar_hi = max(RATIO[s][1] for s in ARISO)

# det0 toy scan table: rms[(tp, f)] = list over T
scan = open(f'{HERE}/det0_toy_scan.log').read()
TS = [100, 140, 180, 220, 340]
RMS = {}
for m in re.finditer(r'tp (\d+) f ([\d.]+): ([\d. ]+)', scan):
    RMS[(int(m.group(1)), float(m.group(2)))] = [float(x) for x in m.group(3).split()]
HE = {}
for m in re.finditer(r'tp (\d+): ([\d.]+)\s+([\d.]+)\s+([\d.]+)', scan):
    HE[int(m.group(1))] = [float(m.group(i)) for i in (2, 3, 4)]

# ion transit predictions: T = g^2 / (K V), K = K0 N0/N (745.83 Torr, 293.15 K)
N0N = (760 / 745.83) * (293.15 / 273.15)
NTd = 745.83 / 760 * 2.6868e19 * 273.15 / 293.15      # gas density [cm^-3]
gar = np.loadtxt(os.path.expanduser('~/Software/garfield/Data/IonMobility_Ar+_Ar.txt'), comments=['#', '*'])
def T_ion(g_um, K0, V=490.):
    g = g_um * 1e-4
    return g * g / (K0 * N0N * V) * 1e9
def K0_ar(g_um, V=490.):
    td = V / (g_um * 1e-4) / NTd * 1e17
    return float(np.interp(td, gar[:, 0], gar[:, 1])), td
PRED = []
for g in (150, 128):
    k, td = K0_ar(g)
    PRED.append((f'Ar⁺ (Garfield), {g} µm', T_ion(g, k), sd.RED,
                 f'K0 = {k:.2f} cm²/V·s at {td:.0f} Td (Garfield IonMobility_Ar+_Ar.txt)\ngap {g} µm, 490 V'))
for name, k in (('C₈Hₙ⁺', 1.742), ('C₄Hₙ⁺', 2.084)):
    for g in (150, 128):
        PRED.append((f'{name}, {g} µm', T_ion(g, k), sd.GREEN,
                     f'K0 = {k} cm²/V·s: Cortez et al. 2019, Table 5, 95 % Ar, 15 Td, 8 Torr (low field)\n'
                     f'gap {g} µm, 490 V. High-field value not measured.'))
chain_T = T_ion(150, K0_ar(150)[0])

D = sd.Deck('MX17: feedback and ions',
            'Two mechanisms behind the T14 sim/data gap, the archive evidence for each, and the literature.')

# ------------------------------------------------------------------ 1 cover
tip_ar = '\n'.join(f'{s}: ×{RATIO[s][0]:.2f}–{RATIO[s][1]:.2f}' for s in ARISO)
cover = (sd.kicker('MX17 response simulation · T14 · 13 Oct 2026') +
         sd.p('The sim/data gap points to two known gas-physics effects the chain leaves out: '
              'avalanche photon feedback, and a drifting ion that is a hydrocarbon fragment, not Ar⁺',
              54, sd.DINK, 600) +
         sd.p('Evidence: 3 chambers × 3 gases from the April–June cosmic bench, plus SPS. '
              'Literature: both effects are described before, and the gain form used here is the field’s standard one.',
              30, sd.DMUT) +
         sd.row(sd.bignum(f'×{ar_lo:.1f}–{ar_hi:.1f}', 'data / Garfield gain slope, Ar/iso',
                          sd.DBLUE, 'same on all three Ar/iso chambers; ×3–7 in Ar/CF₄ and He/C₂H₆', tip=tip_ar),
                sd.bignum('G/(1−βG)', 'fits every ladder', sd.DGREEN,
                          'the photon-feedback form of Şahin, Kowalski & Veenhof (2016)'),
                sd.bignum(f'{chain_T:.0f} → 140–220 ns', 'ion transit, chain → data', sd.DRED,
                          'Cortez et al. 2019 hydrocarbon-ion mobilities give 150–240 ns'), gap=64))
D.slide('cover', cover, dark=True, short='Summary', notes='''
<p>Two separate effects, each now constrained by more than one chamber and more than one gas, account for T14's oldest rows: the gain slope and amplitude, and the rise time / ion signal. Neither is a new idea. Both have a literature, summarised on slides 5 and 8.</p>
<p>Source of every number: <code>MX17_Geant/design/report/FRESH_LOOK_2026-10-12.md</code> §12–13; scripts and logs in <code>MX17_Geant/scratchpad/2026-10-13c/</code>.</p>''')

# ------------------------------------------------------------------ 2 setup
steps1 = [dict(label='Primary avalanche', sub='Townsend, G ≈ 2·10⁴ (Garfield models this)', color=sd.BLUE,
               tip='Garfield AvalancheMicroscopic + Magboltz cross sections, Penning as configured'),
          dict(label='UV photons', sub='Ar excimer ~10 eV; He₂* 16–20 eV', color=sd.GOLD,
               tip='Not followed to the electrodes by Garfield'),
          dict(label='Photoelectron at the mesh / in gas', sub='fraction β per avalanche electron', color=sd.GOLD,
               tip='Micromegas: the mesh (cathode) sees the avalanche directly, unshielded'),
          dict(label='Secondary full avalanche', sub='prompt (~ns): G_eff = G/(1−βG)', color=sd.GREEN,
               tip='Same pulse shape, more charge; diverges as βG → 1')]
steps2 = [dict(label='Ar⁺ made in avalanche', sub='Garfield/chain drifts it as Ar⁺', color=sd.RED,
               tip='chain: IonMobility_Ar+_Ar.txt, T ≈ 340 ns over 150 µm at 490 V'),
          dict(label='Charge transfer to iC₄H₁₀', sub='k ≈ 1.4·10⁻⁹ cm³/s → ~1 ns at 1 atm, 5 %', color=sd.GOLD,
               tip='Cortez et al. 2019 (rate assumed equal to n-butane, Jiao et al. 2007)'),
          dict(label='Hydrocarbon ion drifts', sub='C₄Hₙ⁺ / C₈Hₙ⁺, no resonant charge exchange', color=sd.GREEN,
               tip='K0 1.7–2.1 at 95 % Ar vs Ar⁺ 1.2 at high field'),
          dict(label='Faster ion tail', sub='T 150–240 ns: rise, FWHM, peak time, +8–31 % peak', color=sd.GREEN)]
setup = (sd.title('What the chain leaves out: two links, one on each side of the avalanche',
                  'Top: the gain side (slope, curvature, amplitude). Bottom: the signal side (rise, FWHM, peak time).') +
         sd.p('1 · Gain: photon feedback', 30, sd.INK, 600) + sd.flow(steps1, size=24) +
         sd.p('2 · Signal: which ion drifts', 30, sd.INK, 600) + sd.flow(steps2, size=24) +
         sd.p('Hover the dotted terms, boxes and data points for definitions, counts and sources.', 22, sd.MUT))
D.slide('setup', setup, short='The two links', notes='''
<p>T14 compares the full chain (Geant4 → drift → Garfield avalanche → induction kernels → ion template → DREAM shaper) with det3 cosmic data. Its open rows were: amplitude ×1.8 low, gain slope ×1.5 shallow and curved, rise and FWHM too slow (needing either a ×2.4 faster ion or a ×5 weaker ion signal), and the neighbour-strip delay. The first three are addressed here; the neighbour delay is not.</p>
<p>The two mechanisms act on different observables, so they can be tested separately: feedback is prompt, so it changes the charge and not the shape; the ion species changes the shape and only indirectly the peak.</p>''')

# ------------------------------------------------------------------ 3 ladders
COLS = {'det3 (June, T14) Ar/iso': (sd.BLUE, 'p10', 'circle'), 'det3o Ar/iso': (sd.BLUE, 'p25', 'diamond'),
        'det0 Ar/iso': (sd.BLUE, 'p75', 'open'), 'det3o Ar/CF4': (sd.ORANGE, 'p50', 'diamond'),
        'det0 He/C2H6': (sd.GREEN, 'p75', 'open')}
GAS_S = {'Ar_iC4H10_95_5': 0.318, 'Ar_CF4_90_10': 0.208, 'He_C2H6_96p5_3p5': 0.190}
GCOL = {'Ar_iC4H10_95_5': sd.BLUE, 'Ar_CF4_90_10': sd.ORANGE, 'He_C2H6_96p5_3p5': sd.GREEN}
P = sd.Plot(1060, 640, x=(415, 640, 'lin'), y=(0, 2.0, 'lin'), xlabel='mesh (resist) voltage [V]',
            ylabel='local slope d ln A / dV  [per 10 V]')
P.xticks([(v, str(v)) for v in range(425, 641, 25)]).yticks([(v / 2, f'{v / 2:.1f}') for v in range(0, 5)])
for gas, s in GAS_S.items():
    P.hline(s, GCOL[gas], dash='10 8', w=2.5,
            tip=f'Garfield slope {s:.2f} per 10 V for {gas} (uniform 150 µm gap; no curvature in any gas)')
for s, (c, q, mk) in COLS.items():
    g = L[(L.set == s) & (L.q == q) & L.ok & (L.view == 'x')]
    if s.startswith('det3 (June'):
        g = g[g.scan == 'scan1']
    g = g.groupby('volt').A.mean()
    V = g.index.values.astype(float); y = np.log(g.values)
    xm = 0.5 * (V[1:] + V[:-1]); sl = np.diff(y) / np.diff(V) * 10
    tips = [f'{s}, X {q}\n{V[i]:.0f} → {V[i + 1]:.0f} V: slope {sl[i]:.2f} per 10 V' for i in range(len(sl))]
    P.line(xm, sl, c, w=3, marker=mk, r=7, tips=tips)
P.text(432, 0.93, 'Ar/iso: det3 June (●), old det3 (◆), det0 (○)', 21, sd.BLUE)
P.text(590, 1.02, 'Ar/CF₄ old det3', 21, sd.ORANGE)
P.text(510, 1.2, 'He/C₂H₆ det0', 21, sd.GREEN)
P.text(520, 0.38, 'dashed: Garfield 0.32 / 0.21 / 0.19', 21, sd.MUT)
side = sd.col(
    sd.callout(f'Ar/iso: data slope <b>×{ar_lo:.2f}–{ar_hi:.2f}</b> Garfield on three chambers bulked in April, May and June. '
               'It is the gas, not T14’s chamber.', sd.BLUE, 25),
    sd.callout('Same chamber, gas changed: old det3 ×1.4 → <b>×3.1 in Ar/CF₄</b>; det0 ×1.5 → <b>×3–7 in He/C₂H₆</b>.',
               sd.ORANGE, 25),
    sd.callout('Every data ladder curves up; Garfield’s never does. No α(E) can give this curvature at real gain.',
               sd.GREEN, 25), w=540, gap=22)
D.slide('ladders', sd.title('Every gain ladder is steeper than Garfield, and the excess is set by the gas',
                            'Local slope between neighbouring voltages, X view, clean points only; dashed: Garfield’s slope per gas') +
        sd.row(P.svg('local slope vs voltage'), side, gap=56),
        foot='Threshold-free peak-strip max on M3-selected muons (same observable as T14). FRESH_LOOK §13.2, ladders_all.py.',
        short='Ladders', notes=f'''
<p><b>Data.</b> June det3 (T14 chamber) Ar/iso 95/5; old det3 (bulked 5 May) Ar/iso then Ar/CF₄ 90/10 the next day; det0 (April; 6 mm drift gap) He/C₂H₆ 96.5/3.5 then Ar/iso. Threshold-free peak-strip maximum on M3-good rays, so the event population is the same at every voltage; fixed fiducial per chamber; quantile kept only if ≥ 5× the noise-only maximum and ≤ 0.7 × rail. Points excluded for time drift: det0 Ar/iso 490 V, old det3 Ar/iso 510 V, det0 He 505 V (out of sequence), Ar/CF₄ above 630 V (plateau).</p>
<p><b>Bias direction.</b> All scans ran top-down in voltage and det0's gain rose with time (He: ×2–5 between the first and last 505 V points), which flattens measured slopes: det0's ratios are lower bounds.</p>
<p><b>Fits</b> (<code>fits_all.csv</code>): exponential, G/(1−bG) with Garfield’s slope (2 parameters), and with a free slope. Runaway V* from the Garfield-slope fit: June det3 {VSTAR['det3 (June, T14) Ar/iso'][0]:.0f}–{VSTAR['det3 (June, T14) Ar/iso'][1]:.0f} V, old det3 Ar/iso {VSTAR['det3o Ar/iso'][0]:.0f}–{VSTAR['det3o Ar/iso'][1]:.0f} V, Ar/CF₄ {VSTAR['det3o Ar/CF4'][0]:.0f}–{VSTAR['det3o Ar/CF4'][1]:.0f} V, He {VSTAR['det0 He/C2H6'][0]:.0f}–{VSTAR['det0 He/C2H6'][1]:.0f} V.</p>
<p><b>Why α(E) cannot do it.</b> For any Townsend form the curvature is bounded by (s/V)(sV/ln G − 2); at real gains (ln G ≈ 10) that is 3–7× below what is measured. Penning, gas tables, water and the iso fraction are all α(E) changes.</p>''')

# ------------------------------------------------------------------ 4 width
wlog = open(f'{HERE}/width_vs_v.log').read()
def wseries(block_title, view='x'):
    i = wlog.find(block_title); seg = wlog[i:].split('\n')
    for ln in seg[1:4]:
        if ln.strip().startswith(view + ' '):
            return [(int(a), float(b)) for a, b, _ in re.findall(r'(\d+):([\d.]+)\|([\d.na]+)', ln)]
    return []
W = [('Ar/iso, June det3 (p50/p02)', None, sd.BLUE), ]
P2 = sd.Plot(1000, 600, x=(-100, 5, 'lin'), y=(-0.3, 0.6, 'lin'), xlabel='V − V*  (fitted runaway voltage)  [V]',
             ylabel='Δ width  ln(p_hi/p_lo) − first point')
P2.xticks([(v, str(v)) for v in range(-100, 1, 20)]).yticks([(v / 10, f'{v / 10:+.1f}') for v in range(-2, 7, 2)])
P2.hline(0, sd.MUT, dash='4 6')
Jl = L[(L.set == 'det3 (June, T14) Ar/iso') & (L.view == 'x') & L.ok]
jw = Jl.pivot_table(index='volt', columns='q', values='A').dropna()
jw = np.log(jw.p50 / jw.p02)
SER = [('June det3 Ar/iso', list(zip(jw.index, jw.values)), 520, sd.BLUE, 'circle', 'ln(p50/p02), no floor subtraction'),
       ('old det3 Ar/iso', wseries('det3o Ar/iso: W = ln((p50-f)/(p10-f))'), 525, sd.BLUE, 'diamond', 'ln((p50−f)/(p10−f))'),
       ('det0 Ar/iso', wseries('det0 Ar/iso: W = ln((p90-f)/(p50-f))'), 495, sd.BLUE, 'open', 'ln((p90−f)/(p50−f))'),
       ('old det3 Ar/CF₄', wseries('det3o Ar/CF4: W = ln((p50-f)/(p10-f))'), 638, sd.ORANGE, 'diamond', 'ln((p50−f)/(p10−f))'),
       ('det0 He/C₂H₆', wseries('det0 He/C2H6: W = ln((p75-f)/(p50-f))'), 515, sd.GREEN, 'open', 'ln((p75−f)/(p50−f))')]
for name, ser, vs, c, mk, defn in SER:
    if len(ser) < 2: continue
    v = np.array([a for a, _ in ser], float); w = np.array([b for _, b in ser])
    P2.line(v - vs, w - w[0], c, w=3, marker=mk, r=7,
            tips=[f'{name}, {int(a)} V: width {b:.2f} ({defn})\nV* ≈ {vs} V' for a, b in ser])
side2 = sd.col(sd.p('<b>Pure gain scaling keeps the width fixed.</b> A wrong α(E), a sagging gap or a field-map error multiply every avalanche alike.', 25),
               sd.p('<b>Feedback broadens it near runaway.</b> A random number of secondaries per avalanche (a branching process, mean βG → 1) spreads the amplitudes.', 25),
               sd.p('Seen in Ar/CF₄ and He/C₂H₆ in the last 10–20 V; flat in Ar/iso, which runs at βG ≈ 0.3–0.5. Also true inside 50 mm cells, so not spatial gain spread.', 25),
               w=560, gap=20)
D.slide('width', sd.title('Near runaway the amplitude distribution broadens, as a branching process must',
                          'Change of relative width of the peak-amplitude distribution vs distance from the fitted runaway voltage V*, X view') +
        sd.row(P2.svg('width vs V'), side2, gap=56),
        foot='Noise-only floor subtracted where noted (det0 22 ADC, old det3 38 ADC). width_vs_v.py → width_vs_v.log.',
        short='Width signature', notes='''
<p>Space charge would narrow the distribution at high gain, the opposite sign. The higher quantiles also run away first (He: V* = 510 V on p90, 515 V on p75, 525 V on p50; June det3: βG higher on p50 than on p02), which is the same branching signature seen another way.</p>
<p>Caveat: near the noise floor the lowest quantile is lifted, which narrows the low-voltage widths and could fake broadening; hence the floor subtraction and the ≥ 3× floor cut. The Ar/CF₄ and He broadening survives both.</p>''')

# ------------------------------------------------------------------ 5 feedback literature
refs1 = [
    ('Iguaz et al., microbulk Micromegas (arXiv:1110.2641, 1201.3012)',
     '“the gain curves for low quantities of isobutane show a deviation from the Rose and Korff gain model at high amplification fields. This over-exponential behaviour is due to the low quencher concentrations, which cannot avoid the production of secondary avalanches by the UV photons generated in the primary avalanche.” Ar/iC₄H₁₀ 1–35 %, 50 µm gap; also: resolution worsens below 5 % iso “as not all the UV photons … are absorbed”.',
     'Same gas family, same detector type: the effect we fit is reported, with the same attribution.'),
    ('Şahin, Kowalski, Veenhof, JINST 11 P01003 (2016); NIM A 768 (2014); NIM A 718 (2013)',
     '“Production of the additional electron can be described with a single parameter β … G_T = G/(1 − βG)”, fitted together with Magboltz + Penning gains; β falls with quencher fraction and pressure, tied to the quencher’s photo-absorption cross section.',
     'Our gain model is exactly theirs: Garfield for G, plus one feedback parameter. β ≈ 2·10⁻⁵ here is in their range.'),
    ('Peskov & Fonte, review (arXiv:0911.0463)',
     '“Detectors the cathode of which are not geometrically shielded from the avalanche emission, for example MICROMEGAS … can operate at high gains only in quenched gases”; feedback shows as “gain curve deviation from linear behaviour in log scale” because the feedback pulses are prompt (τ ≪ shaping time); photoemission rises steeply with photon energy.',
     'Explains why the pulse shape is untouched and why He (16–20 eV photons) feeds back more than Ar (~10 eV).'),
    ('Garfield++ MPGD studies (arXiv:1606.04852; MPGD WG4 twiki; 2025 GEM pitch study arXiv:2505.03454)',
     'Simulated gains below measurement are common; Penning fixes much of it for THGEM; photon feedback is named as a missing process (Garfield does not emit electrons from electrodes).',
     'Our ×1.3–1.5 slope gap in Ar/iso is the same kind of known shortfall, with the curvature pointing at feedback rather than Penning.')]
rows = [[f'<b>{a}</b>', b, c] for a, b, c in refs1]
D.slide('lit1', sd.title('Literature: over-exponential Micromegas gain from photon feedback is known, and G/(1−βG) is the standard fit',
                         None) + sd.table(['source', 'what it says', 'what it means here'], rows, size=21,
                                          widths=[330, 820, 500], align=['left', 'left', 'left']),
        short='Lit: feedback', notes='''
<p><b>What the literature does not give us.</b> (1) A measured β for Ar/iC₄H₁₀ 95/5 in a bulk Micromegas: Iguaz et al. show the over-exponential rise at “low quantities” of isobutane but the excerpt does not say where it stops; Şahin et al.’s β tables are for Ar–CO₂ and Ne–CO₂ in wire chambers. (2) A first-principles β: Garfield’s de-excitation model covers Ar levels, its optical data lacks isobutane, and it does not do photoemission from the mesh.</p>
<p><b>Not new, then, is the mechanism.</b> What is new for us is its size at 5 % isobutane (βG ≈ 0.3–0.5 at G ≈ 2·10⁴), and the gas ordering Ar/iso &lt; Ar/CF₄ ≈ He/C₂H₆ seen on the same chambers. My earlier prediction that He/C₂H₆ would feed back less was wrong; the photon-energy dependence of photoemission (Peskov & Fonte) is the likely reason.</p>
<p>Links: <a href="https://arxiv.org/pdf/1110.2641">arXiv:1110.2641</a>, <a href="https://arxiv.org/pdf/1201.3012">arXiv:1201.3012</a>, <a href="https://www.researchgate.net/publication/289584320">Şahin 2016 (Ne–CO₂)</a>, <a href="https://www.sciencedirect.com/science/article/pii/S0168900214010924">Şahin 2014 (Ar–CO₂)</a>, <a href="https://arxiv.org/pdf/0911.0463">Peskov &amp; Fonte</a>, <a href="https://arxiv.org/abs/1606.04852">THGEM Garfield++</a>, <a href="https://arxiv.org/html/2505.03454v1">GEM pitch study</a>, <a href="https://twiki.cern.ch/twiki/bin/view/MPGD/WG4-Simulation">MPGD WG4</a>.</p>''')

# ------------------------------------------------------------------ 6 det0 pulses
t = C['t']; m = (t >= -200) & (t <= 800)
P3 = sd.Plot(800, 560, x=(-200, 800, 'lin'), y=(-0.1, 1.05, 'lin'), xlabel='t from leading 50 % [ns]',
             ylabel='mean peak-strip pulse (norm.)', title='Vertical tracks, det0, 6 mm gap')
P3.xticks([(v, str(v)) for v in range(-200, 801, 200)]).yticks([(v / 4, f'{v / 4:.2f}') for v in range(0, 5)])
P3.line(t[m][::3], C['He'][m][::3], sd.GREEN, w=3, markers=False, tip=f'det0 He/C2H6 data, n = {int(C["nHe"])}')
P3.line(t[m][::3], C['Ar'][m][::3], sd.BLUE, w=4, markers=False, tip=f'det0 Ar/iso data, n = {int(C["nAr"])}')
P3.line(t[m][::3], C['toy_Ar_tp140_f0.9_T180'][m][::3], sd.INK, w=2.5, dash='8 6', markers=False,
        tip='toy Ar/iso, shaper tp 140 (He-calibrated), f 0.9, T 180 ns')
P3.line(t[m][::3], C['toy_Ar_tp140_f0.9_T340'][m][::3], sd.RED, w=2.5, dash='3 5', markers=False,
        tip='toy Ar/iso, tp 140, f 0.9, T 340 ns (the chain’s Ar+)')
P4 = sd.Plot(800, 560, x=(80, 360, 'lin'), y=(0, 0.2, 'lin'), xlabel='ion transit T in the toy [ns]',
             ylabel='rms (toy − data), Ar/iso', title='Ar/iso fit, shaper calibrated on He (tp 130–150)')
P4.xticks([(v, str(v)) for v in TS]).yticks([(v / 20, f'{v / 20:.2f}') for v in range(0, 5)])
FC = {0.2: sd.GREY, 0.5: sd.GOLD, 0.9: sd.INK}
for tp, dash in ((130, '8 6'), (140, None), (150, '2 6')):
    for f in (0.2, 0.5, 0.9):
        r = RMS.get((tp, f))
        if r:
            P4.line(TS, r, FC[f], w=3 if dash is None else 2, dash=dash, r=6 if dash is None else 4,
                    tips=[f'tp {tp} ns, f_ion {f}, T {T} ns: rms {v:.3f}' for T, v in zip(TS, r)])
P4.vline(chain_T, sd.RED, label='chain (Ar⁺)', tip=f'chain ion transit {chain_T:.0f} ns')
he_best = min((min(v), tp) for tp, v in HE.items())
D.slide('det0', sd.title('det0: the Ar/iso pulse rules out the chain’s slow Ar⁺ ions',
                         'Left: mean pulses. Right: Ar/iso rms vs ion transit for f_ion 0.2 / 0.5 / 0.9 (grey / gold / ink), shaper tp 130 (dashed), 140 (solid), 150 (dotted)') +
        sd.row(P3.svg('det0 pulses'), P4.svg('rms vs T'), gap=44) +
        sd.legend([('He/C₂H₆ data', sd.GREEN), ('Ar/iso data', sd.BLUE), ('toy f 0.9, T 180', sd.INK, 'dash'),
                   ('toy f 0.9, T 340 (chain)', sd.RED, 'dash')], 22),
        foot='det0 6 mm drift (Feb logbook run 69+), 800 V ≈ 1.33 kV/cm. det0_toy.py / det0_toy_scan.py; FRESH_LOOK §13.3 + correction.',
        short='det0 ions', notes=f'''
<p><b>Why det0.</b> Its 6 mm gap gives short columns (Ar/iso ~218 ns at ~27.5 µm/ns, He/C₂H₆ ~400 ns at ~15 µm/ns, Magboltz extrapolated to 1.33 kV/cm), and the same chamber and electronics ran He/C₂H₆, whose ions are ~5–10× faster than any ion in Ar. The toy models both gases (cluster statistics, drift, diffusion, gain spread, strips, sampling, selection, 50 %-edge alignment). The shaper time scale is fitted on He first (best rms {he_best[0]:.3f} at tp {he_best[1]} ns), then Ar’s (f, T) is scanned.</p>
<p><b>Result.</b> At every shaper He allows, Ar⁺ with full signal (f 0.9, T 340 ns) is excluded (rms 0.17–0.19 vs a best of ~0.03). The best region is T ≈ 140–220 ns with f ≥ 0.5. Suppressed slow ions (f 0.2, T 340) are disfavoured (0.049–0.059) but not excluded, and at tp 150 a weak fast ion fits about as well: <b>det0 rules out the chain’s ion, it does not by itself fix f.</b></p>
<p><b>A new tension.</b> He calibrates the shaper at tp 130–150 ns, faster than the nominal 180 ns and the 166 ns from the noise autocorrelation. Either the shaper is faster, or the He column or He ion is mis-modelled. Open.</p>
<p>Earlier today I quoted a template fit (Ar = He pulse ⊗ ion) giving f 0.8–0.96, T 120–150 ns; that treated He as a point response, which the 6 mm gap and He’s slow drift make wrong. It is superseded by this toy.</p>''')

# ------------------------------------------------------------------ 7 ion literature + prediction
bars = [(lab, v, c, tp) for lab, v, c, tp in PRED]
P5 = sd.Plot(1000, 600, x=(0, 7, 'lin'), y=(0, 400, 'lin'), ylabel='ion transit T at 490 V [ns]')
P5.yticks([(v, str(v)) for v in range(0, 401, 100)])
P5.band([0, 7], [140, 140], [220, 220], sd.GREEN, 0.12, tip='data: det0 toy fit 140–220 ns; June det3 rise needs ~140 ns at 166 ns shaping')
P5.text(6.95, 232, 'data band (det0 + June det3)', 21, sd.GREEN, 'end')
for i, (lab, v, c, tp) in enumerate(bars):
    P5.vbar(i + 0.6, v, 70, c, tip=f'{lab}: T = {v:.0f} ns\n{tp}', label=f'{v:.0f}')
    P5.text(i + 0.6, -14, lab.split(',')[0], 19, sd.INK, 'middle')
    P5.text(i + 0.6, -38, lab.split(',')[1].strip(), 19, sd.MUT, 'middle')
refs2 = sd.col(
    sd.callout('<b>Cortez et al., JINST 14 P04010 (2019)</b>: in Ar–iC₄H₁₀, Ar⁺ charge-transfers to isobutane; the drifting ions are C₄Hₙ⁺ / C₈Hₙ⁺. K₀ at 95 % Ar: <b>2.08 (C₄)</b>, <b>1.74 (C₈)</b> cm²/V·s (15 Td).', sd.GREEN, 22),
    sd.callout('<b>Giomataris (ICFA 1999)</b>: in Ar + 10 % iC₄H₁₀ “the electron signal and the ion tail are developed within 5 ns and <b>80 ns</b>”, ion drift “typically within 100 ns” (100 µm gap).', sd.GREEN, 22),
    sd.callout(f'<b>Garfield</b> Ar⁺ in Ar: K₀ 1.2 at our E/N (resonant charge exchange slows it). With the chain’s 150 µm gap this is its {chain_T:.0f} ns.', sd.RED, 22),
    w=600, gap=18)
D.slide('ions', sd.title('Literature: the drifting ion is a hydrocarbon, and its mobility predicts the measured transit',
                         'T = g²/(K·V) at 490 V, for the chain’s Ar⁺ and for the measured hydrocarbon-ion mobilities, two gap values') +
        sd.row(P5.svg('ion transit predictions'), refs2, gap=48),
        foot='K = K₀·N₀/N at 745.8 Torr, 293 K. Hydrocarbon K₀ are low-field (15 Td, 8 Torr); their high-field values are not measured.',
        short='Lit: ions', notes='''
<p><b>Ion identity.</b> Cortez et al. assume the Ar⁺ → iC₄H₁₀ charge-transfer rate equals Ar⁺ → n-butane (1.4×10⁻⁹ cm³/s, Jiao et al. 2007), giving C₂–C₃ fragments that react on to C₄Hₙ⁺ and then C₈Hₙ⁺. At 1 atm with 5 % isobutane (n ≈ 1.2×10¹⁸ cm⁻³) the transfer takes ~1 ns, so essentially the whole transit is by hydrocarbon ions. At 8 Torr they see C₄Hₙ⁺ survive only above 80 % Ar; at 1 atm the C₄ → C₈ step is likely complete, so C₈Hₙ⁺ (K₀ 1.74) is the more probable carrier.</p>
<p><b>Gap.</b> The chain uses 150 µm (the effective-gap choice that closes part of the gain); a standard bulk Micromegas has 128 µm. T ∝ g², so the gap is worth ×0.73.</p>
<p><b>Not settled.</b> No measured high-field (~130–160 Td) mobility for these ions; non-resonant ions keep their mobility better than Ar⁺ at high E/N, but by how much is a guess. Giomataris’ 80–100 ns is a 100 µm gap at an unspecified voltage; scaled to 128 µm and 490 V it is ~120–150 ns.</p>
<p>Links: <a href="https://iopscience.iop.org/article/10.1088/1748-0221/14/04/P04010">Cortez 2019</a>, <a href="https://www.slac.stanford.edu/pubs/icfa/fall99/paper1/paper1.pdf">Giomataris 1999</a>, <a href="https://indico.cern.ch/event/709670/contributions/3019642/attachments/1670094/2678865/Ion_Mobility_studies_in_Ar-CF4-iC4H10_mixtures_for_the_LCTPC_Collaboration.pdf">Coimbra LCTPC slides</a>, <a href="https://arxiv.org/pdf/1603.07638">ion mobility in Ar/Ne mixtures</a>, <a href="https://www.nist.gov/publications/resonant-charge-exchange-and-transport-ions-high-electric-field-gas-density-ratios-en">Rao et al. (Ar⁺ at high E/N)</a>.</p>''')

# ------------------------------------------------------------------ 8 ledger
led = sd.hbars([('ion speed', 1.2, sd.RED, 'toy, f 0.91: T 340 → 140 ns raises the peak per charge ×1.08 (vertical), ×1.23–1.26 (10.7°), ×1.27–1.31 (20°); mixed angles ~×1.15–1.25', '×1.15–1.25, mixed angles'),
                ('photon feedback', 1.5, sd.GREEN, '1/(1−βG) at 490 V from the June ladder shape alone: ×1.5–2.1', 'shape fit allows ×1.5–2.1'),
                ('product', 1.8, sd.INK, '×1.2 × ×1.5 ≈ ×1.8', ''),
                ('T14 measured', 1.8, sd.BLUE, 'X, 490 V, rail-bias corrected ×1.82; published ×1.77', 'data / sim amplitude'),
                ('if bench gas ~1 % wet', 2.3, sd.GREY, 'Garfield gain ×0.77 at 1 % H₂O (wet bracket, condor 16705137)', 'need rises; still inside the feedback band × ion')],
               2.6, width=700, label_w=330, fmt=lambda v: f'×{v:.2f}', size=26)
D.slide('ledger', sd.title('Together the two effects close the amplitude row within its errors',
                           'Factors on the det3 June amplitude at 490 V (data over chain)') + led +
        sd.callout('The ADC/fC scale is still a datasheet chain (no pulser shot is possible now), so this closes within calibration, not to a few %.', sd.GOLD, 24),
        short='Amplitude ledger', notes='''
<p>Ion speed does not change the gain slope; feedback does not change the pulse shape. So the slope and curvature rows are feedback alone, the rise/FWHM rows are the ion alone, and only the amplitude row needs both.</p>
<p>The ion factor assumes f ≈ 0.9 (full ion signal). If the ion signal is weaker (det0 does not exclude f 0.5), the ion factor changes and feedback carries more of the amplitude.</p>''')

# ------------------------------------------------------------------ 9 open
end = (sd.kicker('What this does not rule out, and what is next') +
       sd.row(sd.col(sd.p('Not ruled out', 34, sd.DINK, 600),
                     sd.p('• A weak but fast ion signal (f ≈ 0.5, or 0.2 at a 150 ns shaper) instead of a full fast one.', 26, sd.DMUT),
                     sd.p('• A shaper faster than nominal: He says tp 130–150 ns, the noise ACF 166 ns, the cfg 180 ns.', 26, sd.DMUT),
                     sd.p('• Feedback size from first principles: not computed; β is fitted, as in the literature.', 26, sd.DMUT),
                     sd.p('• Hydrocarbon-ion mobility at 130–160 Td: not measured anywhere I found.', 26, sd.DMUT),
                     sd.p('• Neighbour-strip delay and Y lobe: untouched by either effect.', 26, sd.DMUT), gap=14),
              sd.col(sd.p('Next (simulation only)', 34, sd.DINK, 600),
                     sd.p('1. Chain: ion transit 140–220 ns (an exact time rescale) and rerun T14 shapes.', 26, sd.DMUT),
                     sd.p('2. Chain: adopt G/(1−βG), β from the June ladder; re-quote amplitude.', 26, sd.DMUT),
                     sd.p('3. Old det3 Ar/CF₄ shapes with the new Magboltz table: a second ion species (CF₃⁺).', 26, sd.DMUT),
                     sd.p('4. Settle the shaper: He calibration vs noise ACF vs cfg.', 26, sd.DMUT), gap=14), gap=80))
D.slide('open', end, dark=True, short='Open + next')

out = D.write(f'{HERE}/out/mx17-feedback-and-ions.html', note_meta=dict(
    title='MX17 response sim: photon feedback and hydrocarbon ions',
    summary='Archive gain ladders (3 chambers × 3 gases) and det0 pulse shapes against the literature: '
            'G/(1−βG) feedback and a hydrocarbon (not Ar⁺) drifting ion explain the T14 slope, curvature and rise gaps.',
    tags='x17, mx17, simulation, micromegas, literature', date='2026-10-13'))
print(out)
