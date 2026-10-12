#!/usr/bin/env python3
"""Slide note: state of the MX17 response simulation vs det3 data (T14), 2026-10-12.

Figure-first summary of design/report/FRESH_LOOK_2026-10-12.md and
design/DREAM_SETTINGS.md, for the site (dylan-neff.web.cern.ch/notes/).

The waveform figures are computed here from the real peak-strip waveforms
(no quality cut) and from the independent toy model:

    waveforms_angled_pm3.npz, waveforms_pm3.npz
        laptop ~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/
        (made by scratchpad/2026-10-12/laptop_dump_*.py)
    scratchpad/2026-10-12/toy_pulse.py   (numpy only, no response/ imports)

The ledger numbers (amplitude, slope, scans) are the outputs of the
scratchpad scripts as recorded in FRESH_LOOK §1-7.6; each LEDGER entry names
its section. Usage:

    python3 make_sim_status_deck.py --wf-dir <dir with both npz> --out note.html
    python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py note.html \
        --slug mx17-response-sim-status --force --deploy
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / 'scratchpad' / '2026-10-12'))
sys.path.insert(0, os.path.expanduser(os.environ.get(
    'SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd  # noqa: E402
from slidedoc import (BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, INK, MUT,  # noqa: E402
                      DBLUE, DRED, DGREEN, DMUT, DINK)
import toy_pulse as TOY  # noqa: E402

DT = 60.0
GRID = np.arange(-4, 16, 0.5)          # samples from the leading 50 % crossing
ANIM_URL = 'mx17-response-chain.html'

# --------------------------------------------------------------------------- #
# ledger: numbers from FRESH_LOOK_2026-10-12.md (section in each key's comment)
# --------------------------------------------------------------------------- #
LEDGER = dict(
    # §6: q_sum median ratio sim/data within ±3° of the sim's median angle
    amp=[  # selection, data kept X, Y, ratio X, ratio Y
        ('quality_ok (frozen T14)', 0.72, 0.51, 0.7344, 0.6929),
        ('no quality cut', 1.0, 1.0, 0.6111, 0.4836),
        ('gain-normalised χ² < data p95', 0.95, 0.95, 0.5967, 0.4692),
        ('single track candidate', 0.95, 0.92, 0.6180, 0.4858),
    ],
    sat_frac=dict(data_x=0.44, data_y=0.41, sim_x=0.11, sim_y=0.07),
    # §1, §3: gain slope d lnA / dV per 10 V
    slope=[
        ('sim (Garfield S3, meshfield)', 0.311, RED, 'avalanche calib, mesh field map; slope 0.3106 (mesh) vs 0.3107 (uniform)'),
        ('data: T14 HV test', 0.449, INK, 'mx17_sim_wft/hv_slope, det3 27 June scan'),
        ('data: peak sample', 0.419, INK, '425–505 V, 27 June det3 scan (2026-08-29)'),
        ('data: deconvolved charge', 0.4185, INK, '0.4185 ± 0.0023, same scan'),
        ('data: model-free window sum', 0.428, INK, '0.4280 ± 0.0024, same scan'),
    ],
    # §7.3: rise 10–90 % medians (ns) at 0°, 10°, 20°
    rise_toy=dict(full=(370, 279, 265), noion=(303, 176, 144)),
    # §7.2: X edge 10→50 % (ns) at matched peak delay bins
    matched=dict(bins=((350, 450), (450, 550), (550, 700)),
                 data=(159, 184, 210), chain=(154, 189, 219),
                 toy_full=(164, 184, 208), toy_noion=(126, 141, 141)),
    early=dict(data_x=0.32, data_y=0.49, chain_x=0.04, chain_y=0.15, toy=0.0),
    # §7.5: shape rms (units of peak) of the averaged 20° / 10° peak strip, f_ion 0.908
    tscan_T=(340, 260, 227, 200, 170, 110, 80),
    tscan={
        'manual 180 ns': ((0.183, 0.109, 0.085, 0.067, 0.051, 0.029, 0.024),
                          (0.172, 0.107, 0.086, 0.070, 0.056, 0.035, 0.028)),
        'CR-RC² 180 ns': ((0.233, 0.169, 0.148, 0.133, 0.120, 0.102, 0.097),
                          (0.225, 0.169, 0.150, 0.137, 0.125, 0.108, 0.103)),
        'manual 123 ns': ((0.133, 0.037, 0.036, 0.056, 0.080, 0.116, 0.129),
                          (0.115, 0.032, 0.030, 0.046, 0.066, 0.099, 0.110)),
        'CR-RC² 123 ns': ((0.164, 0.069, 0.037, 0.020, 0.027, 0.059, 0.071),
                          (0.151, 0.069, 0.042, 0.024, 0.019, 0.041, 0.050)),
    },
    noion_rms=(0.021, 0.023),
)


# --------------------------------------------------------------------------- #
# waveform analysis
# --------------------------------------------------------------------------- #
def cross(w, frac, side):
    ipk = int(np.argmax(w))
    lvl = frac * w[ipk]
    if side < 0:
        b = np.flatnonzero(w[:ipk + 1] < lvl)
        if not len(b):
            return np.nan
        i = b[-1]
        j = i + 1
    else:
        b = np.flatnonzero(w[ipk:] < lvl)
        if not len(b):
            return np.nan
        j = ipk + b[0]
        i = j - 1
    return i + (lvl - w[i]) / (w[j] - w[i]) if w[j] != w[i] else float(i)


def avg_shape(W):
    """Peak-strip pulses (unsaturated), normalised to peak, aligned at the leading 50 %."""
    W = W[W.max(1) < 3500]
    acc, n = np.zeros_like(GRID), 0
    for w in W:
        l50 = cross(w, .5, -1)
        if np.isnan(l50):
            continue
        acc += np.interp(GRID + l50, np.arange(len(w)), w / w.max())
        n += 1
    return acc / n, n


def rise_median(W):
    W = W[W.max(1) < 3500]
    r = [(cross(w, .9, -1) - cross(w, .1, -1)) * DT for w in W]
    return float(np.nanmedian(r))


TOYS = {  # name -> overrides on toy_pulse.BASE, all through the manual 180 ns shaper
    'noion': dict(shaper='manual_1', f_ion=0.0),
    'full': dict(shaper='manual_1', f_ion=0.908, t_ion=340.0),
    'fast': dict(shaper='manual_1', f_ion=0.908, t_ion=110.0),
}


def _toy_job(arg):
    name, th = arg
    B = TOY.sampled_peak_strips({**TOY.BASE, **TOYS[name], 'theta': th}, 800, seed=4, target_peak=2000)[:, 3]
    return name, th, avg_shape(B)[0].tolist(), rise_median(B)


def analyse(wf_dir):
    za = np.load(Path(wf_dir) / 'waveforms_angled_pm3.npz')
    out = dict(shape={}, rise={}, n={})
    for tag, th in (('th00', 0.0), ('thx10', 10.7), ('thx20', 20.0)):
        for leg in ('data', 'sim'):
            W = za[f'{tag}_{leg}_wf'][:, 3]
            a, n = avg_shape(W)
            out['shape'][(leg, th)] = a
            out['n'][(leg, th)] = n
            out['rise'][(leg, th)] = rise_median(W)
    with ProcessPoolExecutor(min(6, os.cpu_count() or 1)) as ex:
        for name, th, a, r in ex.map(_toy_job, [(n, t) for n in TOYS for t in (10.7, 20.0)]):
            out['shape'][(name, th)] = np.array(a)
            out['rise'][(name, th)] = r
    m = GRID >= -2
    out['rms'] = {(k, th): float(np.sqrt(np.mean((out['shape'][(k, th)][m] - out['shape'][('data', th)][m]) ** 2)))
                  for k in ('sim', *TOYS) for th in (10.7, 20.0)}
    return out


# --------------------------------------------------------------------------- #
# slides
# --------------------------------------------------------------------------- #
C = dict(data=INK, sim=RED, noion=BLUE, full=ORANGE, fast=GREEN)
LBL = dict(data='data (det3 cosmics)', sim='simulation chain', noion='toy: no ion signal',
           full='toy: full ions, Ar⁺ (340 ns)', fast='toy: full ions, 110 ns transit')


def chain_diagram():
    """Two-row block diagram of the chain with a status colour per stage."""
    ok, warn, bad = GREEN, GOLD, RED
    st = [
        ('Geant4', 'muon track, primary clusters', ok,
         'Geant4 + PAI: ClusterTree of primary ionisation in Ar/iC4H10 95/5.\nPrimary ionisation checked: 91 e⁻/cm, matches data (primary_ionisation_2026-08-10).'),
        ('Drift', 'Magboltz v, D_T, D_L', warn,
         'Electrons drift 30 mm at 333 V/cm. Dry 95/5 table: 39.14 µm/ns.\nBundle (data-fitted) v = 36.60 µm/ns: 6.9 % apart, humidity never measured.'),
        ('Mesh', 'transparency (Garfield, T6)', ok,
         'Fraction of drift electrons reaching the amplification gap, from the Garfield mesh field map.\nFixed 2026-08-11 (×1.1455 had been counted twice).'),
        ('Avalanche', 'Garfield S3, 150 µm gap', bad,
         'Garfield++ microscopic avalanche in the meshfield map, 490 V.\nSim gain 24 094; data needs ≈ 4×10⁴ (×1.64 low, X).\nSlope d lnG/dV 0.311 vs data 0.42–0.45 per 10 V.'),
        ('Ions', 'drift back to mesh (S3)', bad,
         'Ions drift 150 µm back to the mesh with the Ar⁺-in-Ar mobility table: ~340 ns rectangle,\nf_ion = 0.908 of the induced charge.\nThe inclined data want this signal largely absent (f_eff ≈ 0.2) at the true 180 ns shaper.'),
        ('Induction', 'resistive ESL → strips (S1/W2)', warn,
         'Semi-spectral weighting-potential solver for the resistive strip stack (ρ_s = 2 MΩ/sq, kapton + glue).\nDelivers the per-strip induced current incl. resistive spreading.\nY/X peak ratio 18 % low; ±1 neighbour delay ~×2 the data.'),
        ('DREAM', 'CSA + 180 ns shaper', warn,
         'DREAM ASIC: CSA (5 µs) + PZC + Sallen-Key shaper, peaking 180 ns (code 2).\nSetting read back from the run cfg 2026-10-12. Actual analog response never measured.'),
        ('FEU', '60 ns samples, ADC, noise', ok,
         'FEU sampling 32 × 60 ns, 200 fC range, pedestal noise from det3.\nADC scale end-to-end 1.2 % from the datasheet chain; noise round-trip 0.2 %.'),
        ('wft', 'same reconstruction as data', ok,
         'Output written as sim_decoded_*.root and read by wft FeuReader unmodified:\nthe identical reconstruction runs on both legs.'),
    ]
    W, H = 1664, 500
    bw, bh, gx = 160, 178, 22
    o = []
    x = 0
    for i, (name, sub, c, tip) in enumerate(st):
        bx = 6 + i * (bw + gx)
        by = 120
        o.append(f'<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" rx="14" fill="#fffdf9" '
                 f'stroke="{c}" stroke-width="4"{sd.tipattr(tip)}/>')
        o.append(f'<rect x="{bx}" y="{by}" width="{bw}" height="14" rx="7" fill="{c}"/>')
        o.append(sd.T(bx + bw / 2, by + 62, name, 27, INK, weight=600))
        words = sub.split(' ')
        lines, cur = [], ''
        for w_ in words:
            if len(cur + ' ' + w_) > 14 and cur:
                lines.append(cur)
                cur = w_
            else:
                cur = (cur + ' ' + w_).strip()
        lines.append(cur)
        for k, ln in enumerate(lines[:3]):
            o.append(sd.T(bx + bw / 2, by + 96 + 24 * k, ln, 19, MUT))
        if i < len(st) - 1:
            o.append(sd.arrow(bx + bw + 2, by + bh / 2, bx + bw + gx - 2, by + bh / 2, MUT, 2.5, 9))
    # group labels
    groups = [(0, 2, 'Geant4 + Magboltz', GREY), (2, 5, 'Garfield++', GREY), (5, 6, 'own solver', GREY),
              (6, 8, 'electronics model', GREY), (8, 9, 'analysis', GREY)]
    for a, b, lab, c in groups:
        x0, x1 = 6 + a * (bw + gx), 6 + b * (bw + gx) - gx
        o.append(sd.line(x0, 92, x1, 92, c, 3))
        o.append(sd.T((x0 + x1) / 2, 76, lab, 21, MUT))
    # observed mismatch annotations below
    notes = [(3, 'gain ×1.6 low,\nslope ×1.5 low', RED), (4, 'pulse too slow /\ntoo wide', RED),
             (5, 'Y/X, ±1 delay', GOLD), (6, 'nominal only', GOLD)]
    for i, txt, c in notes:
        cx = 6 + i * (bw + gx) + bw / 2
        o.append(sd.line(cx, 120 + bh + 6, cx, 120 + bh + 40, c, 2.5))
        for k, ln in enumerate(txt.split('\n')):
            o.append(sd.T(cx, 120 + bh + 70 + 26 * k, ln, 21, c, weight=600))
    return sd.svg(W, H, ''.join(o), 'simulation chain block diagram')


def build(res, out_path):
    L = LEDGER
    D = sd.Deck('MX17 response simulation: state of play',
                'Where the end-to-end MX17 Micromegas response simulation (Geant4 → Garfield → induction → DREAM) '
                'stands against det3 cosmic data on 2026-10-12, and the discrepancies nobody can explain yet.')

    amp_x_nocut = 1 / L['amp'][1][3]
    amp_y_nocut = 1 / L['amp'][1][4]
    s_sim = L['slope'][0][1]
    s_dat = [r[1] for r in L['slope'][1:]]
    r20d, r20s = res['rise'][('data', 20.0)], res['rise'][('sim', 20.0)]

    # 1 cover --------------------------------------------------------------
    D.slide('cover', sd.kicker('MX17 response simulation · det3 cosmics · 2026-10-12')
            + '<h1 style="font-size:76px;font-weight:600;line-height:1.08;letter-spacing:-1px;max-width:1600px">'
              'The full simulation runs end to end and is in the right ballpark, but it fails on three counts nobody can explain yet</h1>'
            + sd.p('Charge too low, gain rising too slowly with HV, and pulses too slow. A clean independent check says the '
                   'simulation code is not the culprit; the physics inputs or a missing mechanism are.', 30, DMUT)
            + sd.row(
                sd.bignum(f'×{amp_x_nocut:.2f}', 'charge deficit (X)',
                          DRED, f'sim/data q_sum {L["amp"][1][3]:.2f} with no quality cut; Y is ×{amp_y_nocut:.2f}',
                          tip='q_sum median, ±3° of the sim median angle, no quality cut (FRESH_LOOK §6).\n'
                              'Lower bound: 44 % of data peak strips hit the ADC ceiling.'),
                sd.bignum(f'×{np.mean(s_dat) / s_sim:.1f}', 'gain slope too shallow',
                          DRED, f'sim {s_sim:.2f} vs data {min(s_dat):.2f}–{max(s_dat):.2f} per 10 V; ≈ ×1.6 at 490 V (curved)',
                          tip='d lnA / dV over 425–505 V, four data estimators vs the Garfield S3 calibration.'),
                sd.bignum(f'+{r20s - r20d:.0f} ns', 'rise time too slow at 20°',
                          DRED, f'sim {r20s:.0f} vs data {r20d:.0f} ns (10–90 %)',
                          tip='Median 10–90 % rise, peak strip, unsaturated, no quality cut.\n'
                              'At 20° the waveform is close to the single-electron response.'),
                gap=64),
            dark=True, short='Summary')

    # 2 the chain ------------------------------------------------------------
    D.slide('chain', sd.title('One muon, nine stages: each piece of the chain is a separate code',
                              'Colour = status. Green: checked against something outside the code. Amber: nominal / unmeasured. '
                              'Red: where the mismatch with data shows up.')
            + chain_diagram()
            + sd.row(sd.callout('Hover a box for what it does and what was checked. '
                                f'An animated walk-through of the chain is <a href="{ANIM_URL}">here</a>.', BLUE, 26),
                     sd.legend([('anchored', GREEN, 'box'), ('nominal / unmeasured', GOLD, 'box'),
                                ('mismatch appears here', RED, 'box')], 24)),
            notes='<p>Target run: det3 cosmic bench, Ar/iC4H10 95/5, 490 V mesh / 1000 V drift over 30 mm '
                  '(<code>cosmic_bench/june_tests/Run/mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V</code>). '
                  'The defaults were frozen before the comparison: ρ_s = 2 MΩ/sq, dry 95/5 gas table, det3 bundle as analysed.</p>'
                  '<p>The chain writes <code>sim_decoded_*.root</code>, which wft reads unmodified, so both legs go through '
                  'the identical reconstruction. Code: <code>MX17_Geant/response/</code> (digitizer, solver, avalanche, dream), '
                  'comparison harness <code>nTof_x17/mx17_sim_wft/t14_compare.py</code>.</p>',
            short='The chain')

    # 3 scoreboard -----------------------------------------------------------
    rows = [
        ['<b>Amplitude</b> q_sum, X / Y', f'×{amp_x_nocut:.2f} / ×{amp_y_nocut:.2f} low', 'open',
         'needs gain ≈ 4×10⁴, above the stable Ar/iso band'],
        ['<b>Gain slope</b> d lnA/dV', f'{s_sim:.2f} vs {min(s_dat):.2f}–{max(s_dat):.2f} /10 V', 'open', 'sharpest number; data slope is curved'],
        ['<b>Rise time</b> 10–90 %, 20°', f'{r20s:.0f} vs {r20d:.0f} ns', 'open', 'needs ion signal ≈ absent at 180 ns'],
        ['FWHM', '×1.18 wide (X)', 'open', 'same family as rise'],
        ['±1 neighbour peak delay', '+93–142 vs +60 ns', 'open', 'points at induction / sharing'],
        ['Y/X peak ratio', 'sim 18 % low on Y', 'open', 'measured under the biased cut, redo'],
        ['Early foot before the pulse', 'data only', 'unexplained', '+0.04–0.16 of peak at −360…−240 ns'],
        ['Negative lobe on Y', 'sim only (−0.2)', 'unexplained', 'data −0.11, sum ≈ 0'],
    ]
    D.slide('score', sd.title('Six discrepancies, two oddities: none is closed',
                              'det3, 490 V, X = FEU 7, Y = FEU 8. Sim ÷ data unless stated.')
            + sd.table(['observable', 'sim vs data', 'status', 'note'], rows, 29,
                       widths=[520, 440, 200, 500], align=['left', 'left', 'left', 'left'],
                       tips=['FRESH_LOOK §6, no quality cut', 'FRESH_LOOK §1, §3', 'computed here from the angled waveform dump',
                             'T14 frozen comparison', 'invariant across 6 SPS arms and 2 gases in data',
                             'OVERNIGHT §5; the quality cut removed 49 % of data Y vs 28 % of X',
                             'FRESH_LOOK §7.2, wf_sumshape.py', 'FRESH_LOOK §7.2']),
            notes='<p>These are the frozen-default numbers, after both 2026-08-11 corrections (transparency double count, '
                  'rise-estimator censoring) and, for amplitude, after removing the biased quality cut (slide 4).</p>'
                  '<p>Amplitude and rise are treated as separate threads: a faster shaper with a long ion signal increases '
                  'the ballistic deficit, which would make the amplitude worse, so no single knob found so far moves both the right way.</p>',
            short='Scoreboard')

    # 4 amplitude --------------------------------------------------------------
    bars = []
    for sel, kx, ky, rx, ry in L['amp']:
        bars.append((sel, 1 / rx, BLUE, f'{sel}\ndata kept: X {kx:.0%}, Y {ky:.0%}\nsim/data X {rx:.4f} → deficit ×{1 / rx:.2f}', 'X'))
        bars.append(('', 1 / ry, PURPLE, f'{sel}\nsim/data Y {ry:.4f} → deficit ×{1 / ry:.2f}', 'Y'))
    rows_ = [(b[0], b[1], b[2], b[3]) for b in bars]
    D.slide('amp', sd.title('The quality cut was hiding a third of the charge deficit',
                            'Data ÷ sim of the median event charge q_sum, by event selection (±3° of the sim\'s median angle).')
            + sd.row(sd.col(sd.hbars(rows_, 2.3, width=600, h=50, label_w=440, fmt=lambda v: f'×{v:.2f}', size=27),
                            sd.legend([('X view', BLUE, 'box'), ('Y view', PURPLE, 'box')]), gap=28, w=1150),
                     sd.col(sd.callout(f'quality_ok is χ²/dof &lt; 300 with pedestal-only weights, so χ² ∝ amplitude². It removed '
                                       '<b>28 % (X) / 49 % (Y)</b> of data tracks but 1 % of sim tracks: the brightest data.', RED, 25),
                            sd.callout('Three fair selections agree: <b>×1.6 (X), ×2.1 (Y)</b>. The rise time does not move (≤ 7 ns).', BLUE, 25),
                            sd.callout(f'Still a lower bound: {L["sat_frac"]["data_x"]:.0%} of data peak strips hit the ADC ceiling '
                                       f'vs {L["sat_frac"]["sim_x"]:.0%} in sim.', GREY, 25), gap=26)),
            notes='<p>Method: <code>t14_compare.py</code> rerun unmodified with <code>{x,y}_quality_ok</code> forced True '
                  '(laptop <code>fresh_look_2026-10-12/t14_compare_transpfix_noq/</code>); ledger re-derived with '
                  '<code>scratchpad/2026-10-12/quality_cut_amplitude.py</code>, which reproduces the frozen numbers exactly with the cut in.</p>'
                  '<p>Consequences: a pure gain fix would need ≈ 4×10⁴, above the stable bulk-Micromegas band of 3–4×10⁴. '
                  'The Y/X asymmetry grows from 1.06× to 1.26×. The q_sum ratio is a fitted quantity and the data leg is on the '
                  'pre-r06 wft bundle, so it still needs a re-quote on r06.</p>',
            foot='FRESH_LOOK_2026-10-12 §6 · TRANSPFIX sim leg vs data_leg_satdet3', short='Amplitude')

    # 5 slope ------------------------------------------------------------------
    srows = [(lab, v, c, f'{lab}: {v:.3f} per 10 V\n{t}') for lab, v, c, t in L['slope']]
    D.slide('slope', sd.title('The gain rises 1.4–1.6× faster with HV in data than in Garfield',
                              'd ln(amplitude) / dV per 10 V around 490 V; det3, 27 June HV scan.')
            + sd.row(sd.col(sd.hbars(srows, 0.5, width=560, h=60, label_w=440, fmt=lambda v: f'{v:.3f}', size=28), w=1150),
                     sd.col(sd.callout('Three data estimators that share nothing but the events agree to 2 %. '
                                       'It is not a data artifact.', INK, 25),
                            sd.callout('The data slope is curved: <b>0.37</b> over 425–450 V, <b>0.50</b> near 490–500 V, '
                                       'so at the operating point the ratio is ≈ 1.6.', BLUE, 25),
                            sd.callout('Penning transfer can give the gain or the slope, never both. '
                                       'An effective gap closes the gain but only 9 % of the slope.', RED, 25), gap=26)),
            notes='<p>Data: <code>mpgd26/slides/HANDOFF_hv_window.md</code> (nTof_x17 commit d8d1101, 2026-08-29). '
                  'Sim: S3 avalanche calibration on the meshfield map (<code>response/avalanche/aval_calib_meshfield_hvscan.json</code>); '
                  'the T6 field-map shape is excluded (mesh 0.3106 vs uniform-field 0.3107).</p>'
                  '<p>Not yet done: a Garfield literature benchmark (a plain 128 µm bulk Micromegas in Ar/iso 95/5 against published '
                  'gain curves). If it under-slopes too, both amplitude and slope become a known Magboltz cross-section limitation '
                  'to calibrate rather than derive.</p>',
            foot='FRESH_LOOK §1, §3', short='Gain slope')

    # 6 what's been excluded -----------------------------------------------------
    def lst(items):
        return '<ul style="font-size:27px;line-height:1.4;padding-left:26px;display:flex;flex-direction:column;gap:10px">' + \
               ''.join(f'<li>{i}</li>' for i in items) + '</ul>'
    D.slide('excluded', sd.title('A long list has been excluded, so the cause is outside the obvious knobs',
                                 'Each was a pre-registered test in August; "do not re-open without new evidence".')
            + sd.row(sd.card('<p style="font-size:30px;font-weight:600;color:#c63d4f">Amplitude &amp; slope</p>' + lst([
                sd.term('Penning transfer', 'Outcome C of the slope hunt: it can close the gain or the slope, never both.'),
                sd.term('Effective gap', 'Closes the gain (×1.69) but only 9 % of the slope.'),
                'Field-map shape (mesh vs uniform)', 'Water, iso fraction (both worse)',
                sd.term('Primary ionisation', 'Measured: 91 e⁻/cm, as simulated.'),
                sd.term('Mesh transparency', 'Fixed 2026-08-11, now counted once (×1.1455).')])),
                sd.card('<p style="font-size:30px;font-weight:600;color:#c63d4f">Rise time</p>' + lst([
                    'Track inclination', 't0 and frame alignment', 'The rise estimator (fixed)',
                    sd.term('Resistive-sheet screening', 'Retired twice by independent derivations: the sheet time response factorises out; ≤ 8 %.'),
                    sd.term('Ion mobility ±30 %', 'Superseded 2026-10-12: a wrong ion species (×1.3–1.6) was not covered; at 180 ns it is not enough anyway.'),
                    sd.term('Shaper 123 ns', 'Closed 2026-10-12 by the register read-back: 180 ns (slide 10).'),
                    sd.term('A chain bug', 'Closed 2026-10-12 by the independent toy (slide 7).')])),
                sd.card('<p style="font-size:30px;font-weight:600;color:#c99318">Seven self-inflicted defects</p>' + lst([
                    'LUT window at 1 µs', 'Mislabelled voltage scan', 'Duplicated transparency',
                    'Railed-waveform rise', 'Bundle overwritten by a variant run', 'Cert harness 101 ns off',
                    'Collect grouping merging arms']) + sd.p('Each moved a headline number. Every exclusion is only as good as the code behind it.', 22, MUT)),
                gap=32),
            notes='<p>Source: FRESH_LOOK §2, <code>DEEP_DIVE_2026-08-11.md</code>, <code>OVERNIGHT_2026-08-10.md</code>. '
                  'The two last rise items were closed today by the fresh look.</p>',
            short='Excluded')

    # 7 toy vs angle --------------------------------------------------------------
    ang = (0, 10.7, 20)
    P = sd.Plot(980, 640, x=(-2, 22), y=(100, 400), xlabel='track angle (deg)', ylabel='median rise 10–90 % (ns)')
    P.xticks([(0, '0°'), (10, '10°'), (20, '20°')]).yticks([(v, str(v)) for v in range(100, 401, 50)])
    dr = [res['rise'][('data', a if a else 0.0)] for a in ang]
    sr = [res['rise'][('sim', a if a else 0.0)] for a in ang]
    for ys, c, lab, d in ((L['rise_toy']['full'], ORANGE, 'toy: full ions', '10 6'),
                          (L['rise_toy']['noion'], BLUE, 'toy: no ions', '10 6')):
        P.line(ang, ys, c, dash=d, tips=[f'{lab}, {a}°: {v} ns' for a, v in zip(ang, ys)])
    P.line(ang, sr, RED, tips=[f'chain, {a}°: {v:.0f} ns' for a, v in zip(ang, sr)])
    P.line(ang, dr, INK, r=8, tips=[f'data, {a}°: {v:.0f} ns' for a, v in zip(ang, dr)])
    D.slide('toy', sd.title('An independent toy agrees with the chain, so the slowness is not a code bug',
                            'Median peak-strip rise vs track angle. Toy: numpy only, shares no code with the chain.')
            + sd.row(P.svg('rise vs angle'),
                     sd.col(sd.legend([('data', INK), ('simulation chain', RED), ('toy, full ions (chain inputs)', ORANGE, 'dash'),
                                       ('toy, no ion signal', BLUE, 'dash')], 24),
                            sd.callout('With the chain\'s own inputs the toy is even <b>slower</b> than the chain '
                                       '(it lacks the PZC undershoot and the spread of ion birth heights).', ORANGE, 25),
                            sd.callout('Data sit on the <b>no-ion</b> toy at every angle. Neither the toy nor the chain gets near the data '
                                       'with the ion signal in.', INK, 25), gap=28, w=600), gap=60),
            notes='<p>Toy (<code>scratchpad/2026-10-12/toy_pulse.py</code>): straight track, Magboltz v and diffusion, '
                  '91.2 e/cm in 3.83-e clusters, Gaussian prompt induction (σ 222 µm), parallel-plate ion current '
                  'Q(t) = (1−f) + f·min(t/T,1), DREAM shaper scaled to 180 ns (5→100 %), 60 ns samples, 11 ADC noise, and the '
                  'identical rise/FWHM measurement of <code>t14_compare.py</code>. The chain\'s settings: f_ion 0.908, T 340 ns.</p>'
                  '<p>Data and chain points are computed here from the angled dumps (unsaturated, no quality cut); the toy points '
                  'are FRESH_LOOK §7.3.</p>',
            foot='Data/chain: waveforms_angled_pm3.npz (no quality cut) · toy: FRESH_LOOK §7.3', short='Toy check')

    # 8 vertical tracks: peak position --------------------------------------------
    mb = L['matched']
    xc = [(a + b) / 2 for a, b in mb['bins']]
    P2 = sd.Plot(900, 600, x=(330, 720), y=(100, 240), xlabel='peak delay after the 10 % onset (ns)',
                 ylabel='leading edge 10→50 % (ns)')
    P2.xticks([(400, '350–450'), (500, '450–550'), (625, '550–700')]).yticks([(v, str(v)) for v in range(100, 241, 20)])
    for k, c, lab, d in (('toy_noion', BLUE, 'toy, no ions', '10 6'), ('toy_full', ORANGE, 'toy, full ions', '10 6'),
                         ('chain', RED, 'chain', None), ('data', INK, 'data', None)):
        P2.line(xc, mb[k], c, dash=d, tips=[f'{lab}, delay {a}–{b} ns: {v} ns' for (a, b), v in zip(mb['bins'], mb[k])])
    e = L['early']
    ebars = [('data X', e['data_x'], INK, 'data X: 32 % peak within 350 ns of onset'),
             ('data Y', e['data_y'], INK, 'data Y: 49 %'),
             ('chain X', e['chain_x'], RED, 'chain X: 4 %'), ('chain Y', e['chain_y'], RED, 'chain Y: 15 %'),
             ('toy', 0.002, ORANGE, 'toy (full ions): 0 %')]
    D.slide('vertical', sd.title('Vertical tracks measure where the peak lands, not the ion signal',
                                 'Left: at matched peak delay the data edge follows the full-ion model. '
                                 'Right: the data simply peak early more often.')
            + sd.row(P2.svg('edge at matched delay'),
                     sd.col(sd.p('Share of pulses peaking &lt; 350 ns after onset', 26, INK, 600),
                            sd.hbars(ebars, 0.55, width=360, h=38, label_w=140, fmt=lambda v: f'{v:.0%}', size=25),
                            sd.callout('So the August "no ions" reading from vertical tracks was confounded. '
                                       'Why data peak early (fast hand-off to neighbours in Y, mild loss with depth) is itself open.',
                                       RED, 24), gap=24, w=700), gap=48),
            notes='<p>For a vertical track the peak strip sees a ~770 ns flat-topped input, so "10–90 % of max" mostly '
                  'measures where on the flat top the maximum falls. Binned by that delay, data (X) and the full-ion chain '
                  'and toy agree within ~10 ns; the no-ion toy is 30–70 ns too fast.</p>'
                  '<p>What makes data pulses peak early: the data peak strip passes late charge to its neighbours faster '
                  '(Y: 0.37 of peak at +480 ns vs the chain\'s 0.66, ±3 sums agree better); summed charge is mildly front-loaded '
                  'in X (0.57 vs 0.68 at +600 ns). Drift velocity × gap is consistent: summed signals end together. '
                  'Scripts: <code>wf_matched.py</code>, <code>wf_sumshape.py</code>.</p>',
            foot='FRESH_LOOK §7.2 · X view · no quality cut', short='Vertical tracks')

    # 9 inclined shapes -----------------------------------------------------------
    t_ns = GRID * DT
    panels = []
    for th, lab in ((20.0, '20° tracks'), (10.7, '10° tracks')):
        P3 = sd.Plot(810, 560, x=(-240, 900), y=(-0.15, 1.05), title=lab,
                     xlabel='time from the leading 50 % point (ns)', ylabel='amplitude / peak')
        P3.xticks([(v, str(v)) for v in range(-200, 901, 200)]).yticks([(v, f'{v:g}') for v in (0, 0.25, 0.5, 0.75, 1)])
        P3.hline(0, MUT, None, 1.5)
        for k in ('full', 'fast', 'noion', 'sim', 'data'):
            a = res['shape'][(k, th)]
            dash = '10 6' if k in TOYS else None
            rms = '' if k == 'data' else f"\nrms vs data {res['rms'][(k, th)]:.3f} of peak"
            nn = f"\nn = {res['n'][(k, th)]}" if k in ('data', 'sim') else ''
            P3.line(t_ns, a, C[k], w=4 if k == 'data' else 3.5, dash=dash, markers=False,
                    tip=f'{LBL[k]}, {lab}\nrise median {res["rise"][(k, th)]:.0f} ns{rms}{nn}')
        panels.append(P3.svg(lab))
    D.slide('inclined', sd.title('At 10° and 20° the data pulse is as narrow as the shaper alone',
                                 'Average peak-strip pulse, normalised and aligned at the leading 50 %. Toys use the 180 ns DREAM shaper.')
            + sd.row(*panels, gap=44)
            + sd.legend([(LBL['data'], INK), (LBL['sim'], RED), (LBL['full'], ORANGE, 'dash'),
                         (LBL['noion'], BLUE, 'dash'), (LBL['fast'], GREEN, 'dash')], 23),
            notes='<p>At 10° / 20° each strip sees only ~110 / ~55 ns of drift, so the waveform is close to the single-electron '
                  'response and the early-peaking escape of slide 8 does not apply (peak-delay distributions are narrow).</p>'
                  f'<p>Shape rms against data, 20° / 10° (units of peak, t ≥ −120 ns): chain {res["rms"][("sim", 20.0)]:.3f} / '
                  f'{res["rms"][("sim", 10.7)]:.3f}; toy full ions {res["rms"][("full", 20.0)]:.3f} / {res["rms"][("full", 10.7)]:.3f}; '
                  f'toy no ions {res["rms"][("noion", 20.0)]:.3f} / {res["rms"][("noion", 10.7)]:.3f}; '
                  f'toy 110 ns ion {res["rms"][("fast", 20.0)]:.3f} / {res["rms"][("fast", 10.7)]:.3f}.</p>'
                  '<p>A front-loaded ion current (Ψ = (1 − z/g)^p, p ≈ 5–8) fixes the rise (20°: 168 vs 165 ns) but leaves the tail too long. '
                  'The "10°" sample is the 10.7° ladder point. Scripts: <code>wf_tail20.py</code>, <code>wf_peak_scan.py</code>.</p>',
            short='Inclined shapes')

    # 10 transit scan -------------------------------------------------------------
    Ts = L['tscan_T']
    P4 = sd.Plot(1060, 640, x=(60, 360), y=(0, 0.25), xlabel='ion transit time T (ns), f_ion 0.908',
                 ylabel='shape rms vs data at 20° (× peak)')
    P4.xticks([(v, str(v)) for v in (80, 110, 170, 227, 260, 340)]).yticks([(v, f'{v:.2f}') for v in (0, .05, .1, .15, .2, .25)])
    P4.raw(f'<rect x="{P4.X(212):.1f}" y="{P4.y0}" width="{P4.X(340) - P4.X(212):.1f}" height="{P4.ph}" fill="{GREEN}" '
           f'fill-opacity="0.08"{sd.tipattr("Plausible: Ar⁺ (340 ns) to an isobutane-like ion ×1.6 faster (~210 ns)")}/>', back=True)
    P4.text(220, 0.238, 'plausible ions', 21, GREEN)
    P4.hline(L['noion_rms'][0], BLUE, '10 6', 2.5, label='no ion signal, 180 ns: 0.021',
             tip='Manual 180 ns shaper with f_ion = 0: rms 0.021 (20°) / 0.023 (10°)')
    sty = {'manual 180 ns': (RED, None, 5), 'CR-RC² 180 ns': (RED, '6 5', 3), 'manual 123 ns': (GREY, None, 3),
           'CR-RC² 123 ns': (GREY, '6 5', 3)}
    for k, (r20, r10) in L['tscan'].items():
        c, d, w = sty[k]
        P4.line(Ts, r20, c, w=w, dash=d, r=5,
                tips=[f'{k}, T {t} ns (×{340 / t:.1f} Ar⁺ mobility): rms {a:.3f} (20°) / {b:.3f} (10°)'
                      for t, a, b in zip(Ts, r20, r10)])
    D.slide('tscan', sd.title('At the real 180 ns shaper, keeping the ions needs an impossibly fast ion',
                              'Fit quality of the averaged 10°/20° shapes vs ion transit time (= mobility scaling), four shaper settings.')
            + sd.row(P4.svg('transit scan'),
                     sd.col(sd.legend([('manual 180 ns (the real setting)', RED), ('CR-RC² 180 ns', RED, 'dash'),
                                       ('123 ns (ruled out)', GREY)], 23),
                            sd.callout('At 180 ns the fit only reaches the no-ion level for <b>T ≲ 110 ns</b>: '
                                       'an ion 3–4× faster than Ar⁺. No plausible ion does that.', RED, 25),
                            sd.callout('A 123 ns shaper with an isobutane-like ion (×1.3–1.6) would have fit with no suppression, '
                                       'but the register says 180 ns.', GREY, 25), gap=26, w=540), gap=56),
            notes='<p>S3 drifts ions with AvalancheMC and no ion diffusion (v = K(E)·E along fixed field lines), so scaling '
                  'K by s rescales the ion current exactly, I(t) → s·I(st). No Garfield rerun is needed for the bracket.</p>'
                  '<p>The August "species: 4 % slower" check is invalid: <code>response/meshcell/ion_template_check.py:193</code> '
                  'applies Blanc\'s law to two different ions. In 95/5 the Ar⁺ hands its charge to isobutane within ~ns, '
                  'and the drifting ion\'s mobility in argon is ~1.6–1.9 at low field vs Ar⁺\'s 1.21 (estimate; Cortez et al., '
                  'JINST 14 (2019) P04010 for low field, no ~124 Td measurement found).</p>',
            foot='FRESH_LOOK §7.5 · wf_ion_T_scan.py', short='Ion transit scan')

    # 11 shaper settled ---------------------------------------------------------------
    rows = [
        ['peaking time', '<b>180 ns</b> (code 2)', '123 ns'],
        ['input polarity / VICM', 'negative / 1.35 V ✓', 'positive / 1.25 V (mismatched) ✗'],
        ['CSA', 'in the path ✓', 'bypassed ✗'],
        ['discriminator range', '17.5 %, ThDAC −127', '5 %, HIT off'],
        ['FEU 1 override 0x881F 0xD043', '180 → 283 ns, 5 → 50 µs', 'peaking unchanged ✗'],
    ]
    D.slide('shaper', sd.title('The shaper is 180 ns: read back from the run\'s own archived config',
                               'Register 1 of FEUs 7/8: Feu * Dream * 1 0x081F 0xD023, identical in all 50 cfgs of the 27 June scan.')
            + sd.row(sd.col(sd.table(['field', 'high word first → 0x081FD023', 'low word first → 0xD023081F'], rows, 25,
                                     widths=[420, 460, 480], align=['left', 'left', 'left']), w=1160),
                     sd.col(sd.callout('Only high-word-first gives a working Micromegas front end, and only it makes '
                                       'the deliberate overrides decode as deliberate changes.', GREEN, 25),
                            sd.callout('Every data FEU in the cosmic, n_TOF and SPS runs (3574 cfgs): <b>180 ns</b>.', INK, 25),
                            sd.callout('Still unmeasured: the chip\'s <i>actual</i> response (spread, PZC trim). '
                                       'One FEU pulser shot measures it.', GOLD, 25), gap=24)),
            notes='<p>August\'s audit said the target run "archived no .cfg". It did: '
                  '<code>/eos/experiment/ntof/data/x17/cosmic_bench/june_tests/Run/mx17_det3_saturday_scan_6-27-26/&lt;subrun&gt;/raw_daq_data/CosmicTb_MX17.cfg</code>. '
                  'Decoder: <code>nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py</code>. Full record: <code>MX17_Geant/design/DREAM_SETTINGS.md</code>.</p>'
                  '<p>A third witness: the July and SPS self-trigger threshold scans encoded ThDAC (which straddles the two words) '
                  'high word first and produced working scans.</p>',
            foot='design/DREAM_SETTINGS.md · per-run table design/dream_cfg/dream_settings_by_run.csv', short='Shaper = 180 ns')

    # 12 close ------------------------------------------------------------------------
    def dl(items, c):
        return '<ul style="font-size:26px;line-height:1.4;padding-left:28px;display:flex;flex-direction:column;gap:10px;color:#eef0f3">' + \
               ''.join(f'<li>{i}</li>' for i in items) + '</ul>'
    D.slide('open', sd.kicker('Open problems and next steps')
            + '<h2 style="font-size:58px;font-weight:600;line-height:1.1">No mechanism yet for any of the three. '
              'Two measurements and one benchmark would narrow it most.</h2>'
            + sd.row(sd.col(sd.p('Unexplained', 30, DRED, 600), dl([
                'Ion signal ≈ 5× suppressed at the true shaper (f_eff ≈ 0.2 vs 0.91)',
                'Charge ×1.6 (X) / ×2.1 (Y) low', 'Gain slope ×1.4–1.6 shallow (data slope curved)',
                'Early foot before data pulses; −0.2 lobe on sim Y', '±1 neighbour delay ×2, Y/X 18 %'], DRED)),
                sd.col(sd.p('Next', 30, DGREEN, 600), dl([
                    '<b>FEU pulser shot</b> through the test capacitor: the shaper\'s real response and the absolute ADC/fC scale',
                    'Look for an ion-suppression mechanism specific to the <b>resistive</b> bulk Micromegas (only sheet screening, ≤ 8 %, sized so far)',
                    '<b>Garfield literature benchmark</b>: plain bulk MM in Ar/iso 95/5 vs published gain curves (does Magboltz under-slope?)',
                    'Template-to-template: wft calibrate on sim output vs data', 'Re-quote q_sum on the r06 bundle; redo Y/X without the cut'], DGREEN)),
                sd.col(sd.p('Not ruled out', 30, DMUT, 600), dl([
                    'ρ_s off 2 MΩ/sq (band 1.4–2.6)', 'Gas contamination (no hygrometer on any run)',
                    'Bench conditions differing from the cfg (HV, temperature)', 'A shared assumption all chain pieces make'], DMUT)),
                gap=48),
            dark=True, short='Open problems')

    D.write(out_path, note_meta=dict(
        title='MX17 response simulation: state of play (Oct 2026)',
        summary='The Geant4 → Garfield → induction → DREAM chain vs det3 cosmics: charge ×1.6 low, gain slope ×1.5 shallow, '
                'pulses too slow. An independent toy clears the code; the shaper is confirmed 180 ns, so the ion signal must be '
                '≈ 5× suppressed, with no mechanism found.',
        tags='MX17,simulation,T14', date='2026-10-12'),
        footer='Built by MX17_Geant/design/report/deck/make_sim_status_deck.py from FRESH_LOOK_2026-10-12.md and the waveform dumps.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--wf-dir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = analyse(a.wf_dir)
    for k, v in sorted(res['rise'].items(), key=str):
        print('rise', k, round(v))
    for k, v in sorted(res['rms'].items(), key=str):
        print('rms', k, round(v, 3))
    build(res, a.out)
    print('wrote', a.out)
