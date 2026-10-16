#!/usr/bin/env python3
"""Slide note: the MX17 response simulation against det3 data, state on 2026-10-17.

Follow-up to make_sim_status_deck.py (2026-10-12, "six discrepancies, none closed"). It summarises
design/report/FRESH_LOOK_2026-10-12.md §8-20: what now matches the det3 target run, what still
disagrees, and ideas to fix each remaining row.

Inputs:
    scratchpad/2026-10-16/results/      ledgers, shape and Y-profile logs of the 2026-10-16 chain points
                                        (ledger_m3.py / shape_all.py / profile_tails.py on yoga), and the
                                        pooled nc0 / nc1 avalanche calibrations (from EOS response_sim/avalanche/)
    LEDGER below                        earlier numbers, each with its FRESH_LOOK section

Usage:
    python3 make_sim_state_deck.py --out note.html
    python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py note.html \\
        --slug mx17-response-sim-state --force --deploy
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = REPO / 'scratchpad' / '2026-10-16' / 'results'
sys.path.insert(0, os.path.expanduser(os.environ.get(
    'SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd  # noqa: E402
from slidedoc import (BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, INK, MUT,  # noqa: E402
                      DBLUE, DRED, DGREEN, DMUT, DINK)

STATUS_URL = 'mx17-response-sim-status.html'
ANIM_URL = 'mx17-response-chain.html'
ANGLES = ['th00', 'thx10', 'thx20', 'thy10', 'thy20']
ANG_LABEL = {'th00': 'vertical', 'thx10': '10° in x', 'thx20': '20° in x', 'thy10': '10° in y', 'thy20': '20° in y'}

# one colour per configuration, used on every slide
C_DATA, C_OLD, C_V1, C_E17, C_R5 = INK, GREY, BLUE, ORANGE, GREEN
C_X, C_Y = BLUE, PURPLE

# --------------------------------------------------------------------------- #
# earlier numbers (FRESH_LOOK section in each comment)
# --------------------------------------------------------------------------- #
LEDGER = dict(
    # §18.2: unbiased (M3-angle) ledger, vertical tracks, peak p50 data/sim, X and Y
    journey=[
        ('August', 'Ar⁺ ion, no feedback, 180 ns shaper (cluster 4450898, s=1)', 2.07, 2.43),
        ('+ fast ion', 'ion ×1.75 Ar⁺ mobility, no feedback (4450898, s=1.75)', 1.66, 1.90),
        ('+ feedback', 'photon feedback G/(1−bG), bG fitted to the June ladder shape (4450899)', 0.90, 1.01),
    ],
    # §16 table, X view, T14 data legs: rise 10–90 % median, sim − data [ns], Ar⁺ / 180 ns default
    rise_aug=dict(x10=251 - 176, x20=231 - 162),
    # §10.4: ±1 neighbour peak delay after the central strip, vertical, data vs chain [ns] (pre-verdict-1)
    nb_delay=dict(x=(90, 135), y=(75, 120)),
    # §1 / old deck: gain slope d lnA/dV per 10 V
    slope=dict(sim=0.311, data=(0.419, 0.449)),
    # §18.4-18.5: strips >= 5σ, sim X/Y vs data X/Y, before the common-mode fix
    nover_before=((12, 13), (6, 8)),
    sigma_cm=(271, 5.7),
    # §18.2 Y vertical: the T14 data-leg selection vs true angle
    selection_y=1.37,
    # §20.2 / §18.10: rho_s band from T2b
    t2b=(1.42, 2.56),
    # §19: Irakli TPOT bench gain vs nominal (mV/fC)
    adc_fc=(18.05, 20.0),
    # e17c72c Stage A fix: nPrimDrift medians new/old
    prim=(0.898, 0.915),
)


# --------------------------------------------------------------------------- #
# readers
# --------------------------------------------------------------------------- #
def read_ledger(name):
    """ledger_m3.py output -> {(point, angle, view): dict(p50, p25, sat_d, sat_s)}."""
    out = {}
    for ln in open(RES / name):
        m = re.match(r'(\S+)/(\S+)\s+([xy])\s+\d+\s+\S+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)/([\d.]+)', ln)
        if m:
            pt = m[1].replace('w2_rho5M_', '')
            out[(pt, m[2], m[3])] = dict(p50=float(m[4]), p25=float(m[5]), sat_d=float(m[6]), sat_s=float(m[7]))
    return out


def read_shape():
    """shape_all.py output -> {(point, view, angle): (rise_sim, rise_diff, fwhm_sim, fwhm_diff)}."""
    out = {}
    for ln in open(RES / 'shape_1016.log'):
        m = re.match(r'(\S+)\s+([xy])\s+(\d+)\s+rise\s+(\d+)\s+\(\s*([+-]?\d+)\)\s+fwhm\s+(\d+)\s+\(\s*([+-]?\d+)\)', ln)
        if m:
            pt = m[1].replace('w2_rho5M_', '')
            out[(pt, m[2], int(m[3]))] = tuple(int(m[i]) for i in (4, 5, 6, 7))
    return out


def read_profile(name):
    """profile_tails.py output -> {'x sim'|'x data'|'y sim'|'y data': {offset: amp/peak}}."""
    out, cur, offs = {}, None, None
    for ln in open(RES / name):
        m = re.match(r'([xy]) (sim|data)\s+n=', ln)
        if m:
            cur = f'{m[1]} {m[2]}'
            continue
        f = ln.split()
        if f and f[0] == 'off':
            offs = [int(v) for v in f[1:]]
        elif cur and f[:1] == ['amp/pk']:
            out[cur] = dict(zip(offs, [float(v) for v in f[1:]]))
    return out


def read_calib(name):
    p = json.load(open(RES / name))['point']
    q = p['polya']
    return dict(surv=q['survival'], surv_err=q['survival_err'], gain=q['gain_mean'],
                charge=q['gain_mean_uncond'], theta=q['theta'], nev=p['nev_total'], n=p['n_slices'])


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs)


def load():
    L = {}
    for f in ('v1', 'ionK', 'r5', 'e17'):
        L.update(read_ledger(f'ledger_1016_{f}.log'))
    return dict(L=L, S=read_shape(),
                prof={(r, a): read_profile(f'profile_1016_{r}_{a}.log')
                      for r in ('rho2M_verdict1', 'rho5M_r5verdict1') for a in ('th00', 'thy10')},
                nc0=read_calib('aval_calib_meshfield_nc0_pooled.json'),
                nc1=read_calib('aval_calib_meshfield_nc1_pooled.json'))


# --------------------------------------------------------------------------- #
# derived numbers
# --------------------------------------------------------------------------- #
def view_mean(L, pt, v, key='p50'):
    return mean(L[(pt, a, v)][key] for a in ANGLES if (pt, a, v) in L)


def interp_zero(xs, ys):
    """x where ln(y) = 0, linear in ln(y) between the first two points that bracket it."""
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        a, b = math.log(y0), math.log(y1)
        if a == b:
            continue
        t = -a / (b - a)
        if 0 <= t <= 1 or (x1 == xs[-1]):
            return x0 + t * (x1 - x0)
    return float('nan')


def derived(R):
    L = R['L']
    d = {}
    d['bg'] = [0.45, 0.50, 0.55]
    d['bx'] = [view_mean(L, f'e17b{int(b * 100)}', 'x') for b in d['bg']]
    d['by'] = [view_mean(L, f'e17b{int(b * 100)}', 'y') for b in d['bg']]
    d['bj'] = [math.sqrt(x * y) for x, y in zip(d['bx'], d['by'])]
    d['bg_x'], d['bg_y'], d['bg_j'] = (interp_zero(d['bg'], d[k]) for k in ('bx', 'by', 'bj'))
    d['v1x'], d['v1y'] = view_mean(L, 'verdict1', 'x'), view_mean(L, 'verdict1', 'y')
    # amplitude ranges of the current best configurations
    d['v1_rng'] = (min(L[k]['p50'] for k in L if k[0] == 'verdict1'), max(L[k]['p50'] for k in L if k[0] == 'verdict1'))
    # Y/X: (Y data/sim)/(X data/sim) per angle, verdict 1 and rho 5M
    d['yx'] = {}
    for a in ANGLES:
        row = {}
        for pt in ('verdict1', 'r5verdict1'):
            if (pt, a, 'x') in L and (pt, a, 'y') in L:
                row[pt] = L[(pt, a, 'y')]['p50'] / L[(pt, a, 'x')]['p50']
        d['yx'][a] = row
    # Y wings at vertical: mean of both sides, data / sim 2M / sim 5M
    P2, P5 = R['prof'][('rho2M_verdict1', 'th00')], R['prof'][('rho5M_r5verdict1', 'th00')]
    w = {}
    for o in (1, 2, 3):
        w[o] = dict(data=(P2['y data'][o] + P2['y data'][-o]) / 2, s2=(P2['y sim'][o] + P2['y sim'][-o]) / 2,
                    s5=(P5['y sim'][o] + P5['y sim'][-o]) / 2)
        # log-log interpolation in rho of sim/data to 1
        r2, r5 = w[o]['s2'] / w[o]['data'], w[o]['s5'] / w[o]['data']
        w[o]['rho'] = 2 * 2.5 ** (math.log(r2) / (math.log(r2) - math.log(r5)))
    d['wing'] = w
    yx0 = d['yx']['th00']
    d['rho_yx0'] = 2 * 2.5 ** (math.log(yx0['verdict1']) / (math.log(yx0['verdict1']) - math.log(yx0['r5verdict1'])))
    d['yx_move'] = [d['yx'][a]['verdict1'] / d['yx'][a]['r5verdict1'] - 1 for a in ANGLES if 'r5verdict1' in d['yx'][a]]
    return d


# --------------------------------------------------------------------------- #
# slides
# --------------------------------------------------------------------------- #
def chain_flow():
    st = dict(ok=GREEN, warn=GOLD, bad=RED)
    a = [
        dict(label='Geant4', sub='primary ionisation<br><b>fixed −11 %</b> (e17c72c)', color=st['ok'],
             tip='mm_sim, PAI in the gas. Fixed 2026-10-15 (e17c72c): gas densities, W value and PAI were wrong;\n'
                 'primaries were ~11 % high. Heed vs Geant4 agreement confirmed by Alexandra.'),
        dict(label='Drift', sub='Magboltz table<br>dry gas assumed', color=st['warn'],
             tip='Dry 95/5 table. Data-fitted drift speed is 6.9 % lower (wet gas?). Humidity never measured.'),
        dict(label='Mesh', sub='transparency ε<br><b>0.955 → 0.993</b>, net +0.7 %', color=st['ok'],
             tip='0.955 was a Garfield stepping artefact. Null-collision steps give 0.993 ± 0.001,\n'
                 'but the gain per surviving electron is 3 % lower, so ε × gain moves +0.7 % (condor 4474926).'),
        dict(label='Avalanche', sub='Garfield S3<br><b>+ feedback bG 0.45–0.47</b>', color=st['warn'],
             tip='Garfield Townsend gain plus photon feedback G/(1 − bG). bG 0.45 fitted to the June ladder shape;\n'
                 'with the corrected primaries the joint X+Y value is 0.47. Not tested out of sample.'),
        dict(label='Ions', sub='Ar⁺ template ×1.75<br><b>= hydrocarbon K0 ≈ 2.2</b>', color=st['ok'],
             tip='Ion signal from the calib\'s birth heights rebuilt in the T6 field map (--ion-model reconstructed).\n'
                 'K0 2.05 ≈ verdict 1\'s 1.75 × Ar⁺ template; 1.8 / 2.25 bracket it.'),
    ]
    b = [
        dict(label='Induction', sub='resistive anode, ρ_s 2 MΩ/sq<br><b>Y/X and Y wings off</b>', color=st['bad'],
             tip='COMSOL-validated kernels (2D, axisymmetric, real-strip 3D). ρ_s = 2 MΩ/sq default.\n'
                 'The Y wings want 4–5.6 MΩ/sq; the Y/X ratio does not follow ρ_s at all.'),
        dict(label='DREAM', sub='<b>shaper 166 ns</b> (noise ACF)<br>ADC/fC from datasheet', color=st['warn'],
             tip='Shaping time from the target run\'s own noise autocorrelation (FEU 7/8: 166 ns; FEU 3: 158 ns).\n'
                 'ADC/fC is the datasheet chain; Irakli\'s TPOT bench gives −10 %.'),
        dict(label='FEU noise', sub='<b>common mode fixed</b><br>σ_cm 271 → 5.7 ADC', color=st['ok'],
             tip='σ_cm was 271 ADC (std incl. railed events); quiet-event value 5.7 / 4.85 ADC.\n'
                 'Fake hits gone; sim noise 7.41 ADC as in data.'),
        dict(label='wft reco', sub='same code both legs<br><b>data by true M3 angle</b>', color=st['ok'],
             tip='The chain writes sim_decoded_*.root, which wft reads unmodified.\n'
                 'Comparison: ledger_m3.py (data by true M3 angle), shape_all.py, profile_tails.py.'),
    ]
    return sd.flow(a, size=29) + sd.flow(b, size=29)


def amp_journey(L, d):
    steps = list(LEDGER['journey'])
    v1 = ('verdict 1', 'hydrocarbon ion (s 1.75), bG 0.45, shaper 166 ns, noise fixed (4451436)',
          L[('verdict1', 'th00', 'x')]['p50'], L[('verdict1', 'th00', 'y')]['p50'])
    e45 = ('+ primaries', 'Stage A fix (−9 % primaries), bG 0.45 (4472245)',
           L[('e17b45', 'th00', 'x')]['p50'], L[('e17b45', 'th00', 'y')]['p50'])
    e50 = ('bG 0.50', 'e17 primaries, bG 0.50 (4472245)',
           L[('e17b50', 'th00', 'x')]['p50'], L[('e17b50', 'th00', 'y')]['p50'])
    steps += [v1, e45, e50]
    xs = list(range(len(steps)))
    P = sd.Plot(1060, 620, x=(-0.5, len(steps) - 0.5), y=(0.7, 2.6), xlabel='',
                ylabel='peak-strip median, data ÷ sim (vertical)', margin=(24, 30, 70, 104))
    P.band([-0.5, len(steps) - 0.5], [0.9, 0.9], [1.1, 1.1], GREEN, 0.12, tip='±10 % of the data')
    P.hline(1.0, INK, '6 6', 2, label='sim = data')
    P.xticks([(i, s[0]) for i, s in enumerate(steps)])
    P.yticks([(v, f'{v:.1f}') for v in (0.8, 1.0, 1.4, 1.8, 2.2, 2.6)])
    tx = [f'{s[0]}: {s[1]}\nX data/sim {s[2]:.2f}' for s in steps]
    ty = [f'{s[0]}: {s[1]}\nY data/sim {s[3]:.2f}' for s in steps]
    P.line(xs, [s[2] for s in steps], C_X, tips=tx, r=8)
    P.line(xs, [s[3] for s in steps], C_Y, tips=ty, r=8, marker='diamond')
    return P.svg('peak amplitude data/sim through the fixes')


def bg_plot(d, L):
    P = sd.Plot(1000, 620, x=(0.425, 0.565), y=(0.74, 1.18), xlabel='feedback bG (corrected primaries)',
                ylabel='peak median data ÷ sim (mean of 5 angles)')
    P.hline(1.0, INK, '6 6', 2, label='sim = data')
    P.xticks([(b, f'{b:.2f}') for b in (0.43, 0.45, 0.47, 0.50, 0.53, 0.55)])
    P.yticks([(v, f'{v:.1f}') for v in (0.8, 0.9, 1.0, 1.1)])
    per = {v: ['\n'.join(f'  {ANG_LABEL[a]}: {L[(f"e17b{int(b * 100)}", a, v)]["p50"]:.3f}' for a in ANGLES)
               for b in d['bg']] for v in 'xy'}
    P.line(d['bg'], d['bx'], C_X, r=8, tips=[f'X, bG {b:.2f}: mean {y:.3f}\n{t}' for b, y, t in zip(d['bg'], d['bx'], per['x'])])
    P.line(d['bg'], d['by'], C_Y, r=8, marker='diamond',
           tips=[f'Y, bG {b:.2f}: mean {y:.3f}\n{t}' for b, y, t in zip(d['bg'], d['by'], per['y'])])
    P.line(d['bg'], d['bj'], INK, w=2.5, dash='3 6', r=5,
           tips=[f'joint (geometric mean of X and Y), bG {b:.2f}: {y:.3f}' for b, y in zip(d['bg'], d['bj'])])
    P.points([0.45], [d['v1x']], C_X, r=9, marker='open',
             tips=[f'verdict 1 (old primaries, bG 0.45): X mean {d["v1x"]:.3f}'])
    P.points([0.45], [d['v1y']], C_Y, r=9, marker='open',
             tips=[f'verdict 1 (old primaries, bG 0.45): Y mean {d["v1y"]:.3f}'])
    for k, c, lab in (('bg_x', C_X, 'X'), ('bg_y', C_Y, 'Y'), ('bg_j', INK, 'joint')):
        P.vline(d[k], c, '2 5', 2, tip=f'{lab} alone matches the data at bG ≈ {d[k]:.3f} (linear in ln between points)')
    P.text(d['bg_x'] - 0.002, 0.765, f'X {d["bg_x"]:.2f}', 22, C_X, 'end')
    P.text(d['bg_y'] + 0.002, 0.765, f'Y {d["bg_y"]:.2f}', 22, C_Y)
    P.text(d['bg_j'] + 0.002, 1.15, f'joint {d["bg_j"]:.2f}', 22, INK)
    return P.svg('amplitude vs feedback bG')


def rise_plot(S):
    cats = [('x', 10), ('x', 20), ('y', 10), ('y', 20)]
    series = [('Ar⁺, 180 ns (Aug default)', C_OLD, None),
              ('verdict 1 (template ×1.75)', C_V1, 'verdict1'),
              ('K0 1.80', '#9bb7d9', 'ionK180'),
              ('K0 2.05', '#4f86c6', 'ionK205'),
              ('K0 2.25', '#163f73', 'ionK225')]
    P = sd.Plot(1080, 600, x=(-0.5, 3.5), y=(-10, 80), xlabel='',
                ylabel='rise 10–90 %, sim − data [ns]', margin=(24, 30, 70, 104))
    P.hline(0, INK, None, 2)
    P.band([-0.5, 3.5], [-5, -5], [5, 5], GREEN, 0.12, tip='±5 ns')
    P.xticks([(i, f'{v.upper()} {a}°') for i, (v, a) in enumerate(cats)])
    P.yticks([(v, f'{v:+d}' if v else '0') for v in (0, 20, 40, 60, 80)])
    bw = 34
    for j, (name, c, pt) in enumerate(series):
        for i, (v, a) in enumerate(cats):
            px = P.X(i) + (j - 2) * (bw + 4)
            if pt is None:
                val = LEDGER['rise_aug'].get(f'{v}{a}')
                if val is None:
                    continue
                tip = f'{name}, {v.upper()} {a}°: +{val} ns (T14 data leg, FRESH_LOOK §16)'
            else:
                rs, val, fs, fd = S[(pt, v, a)]
                tip = (f'{name}, {v.upper()} {a}°\nrise sim {rs} ns, data {rs - val} ns → {val:+d} ns\n'
                       f'FWHM sim {fs} ns → {fd:+d} ns vs data')
            P.vbar((px,), val, bw, c, tip=tip, base=0)
    return P.svg('rise residual by configuration'), series


def profile_plot(R, view):
    P2, P5 = R['prof'][('rho2M_verdict1', 'th00')], R['prof'][('rho5M_r5verdict1', 'th00')]
    P = sd.Plot(810, 560, x=(-4.5, 4.5), y=(0.005, 1.3, 'log'), xlabel='strip offset from the peak strip',
                ylabel='amplitude ÷ peak (median)' if view == 'y' else '', title=f'{view.upper()} view, vertical tracks')
    P.xticks([(o, f'{o:+d}' if o else '0') for o in range(-4, 5)])
    P.yticks([(0.01, '0.01'), (0.03, '0.03'), (0.1, '0.1'), (0.3, '0.3'), (1, '1')])
    offs = list(range(-4, 5))
    for lab, prof, key, c, dash in (('data', P2, f'{view} data', C_DATA, None),
                                    ('sim ρ_s 2 MΩ/sq (verdict 1)', P2, f'{view} sim', C_V1, None),
                                    ('sim ρ_s 5 MΩ/sq', P5, f'{view} sim', C_R5, '10 6')):
        ys = [prof[key][o] for o in offs]
        P.line(offs, ys, c, dash=dash, r=6, tips=[f'{lab}, offset {o:+d}: {y:.3f} of peak' for o, y in zip(offs, ys)])
    return P.svg(f'{view} strip profile')


def rho_plots(d):
    rmin, rmax = 1.0, 30.0
    t2b = LEDGER['t2b']
    # left: Y/X excess
    A = sd.Plot(810, 560, x=(rmin, rmax, 'log'), y=(0.95, 1.30), xlabel='ρ_s [MΩ/sq]',
                ylabel='(Y data/sim) ÷ (X data/sim)', title='Y/X: barely moves')
    B = sd.Plot(810, 560, x=(rmin, rmax, 'log'), y=(0.5, 2.2, 'log'), xlabel='ρ_s [MΩ/sq]',
                ylabel='Y wing, sim ÷ data', title='Y wings: close near 4–6 MΩ/sq')
    for P in (A, B):
        P.band(list(t2b), [P.ys[0]] * 2, [P.ys[1]] * 2, GOLD, 0.14,
               tip=f'T2b ρ_s band {t2b[0]}–{t2b[1]} MΩ/sq (measurement on the board type)')
        P.xticks([(v, f'{v:g}') for v in (1, 2, 5, 10, 20)])
        P.hline(1.0, INK, '6 6', 2, label='matches data', anchor='start')
        P.vline(2, MUT, '2 5', 1.5)
        P.vline(5, MUT, '2 5', 1.5)
    A.yticks([(v, f'{v:.2f}') for v in (1.0, 1.1, 1.2, 1.3)])
    B.yticks([(v, f'{v:g}') for v in (0.5, 0.7, 1, 1.5, 2)])
    shades = {'th00': INK, 'thx10': BLUE, 'thy10': PURPLE, 'thy20': ORANGE}
    for a, c in shades.items():
        r = d['yx'][a]
        ys = [r['verdict1'], r['r5verdict1']]
        A.line([2, 5], ys, c, r=7, tips=[f'{ANG_LABEL[a]}, ρ_s {x} MΩ/sq: {y:.3f}' for x, y in zip((2, 5), ys)])
        A.text(5.4, ys[1] - 0.005, ANG_LABEL[a], 21, c)
    # extrapolate the vertical one (log-linear in rho) to where it reaches 1
    r0 = d['yx']['th00']
    k = math.log(r0['r5verdict1'] / r0['verdict1']) / math.log(2.5)
    xs = [5, d['rho_yx0']]
    A.line(xs, [r0['verdict1'] * (x / 2) ** k for x in xs], INK, w=2, dash='4 6', markers=False)
    A.points([d['rho_yx0']], [1.0], INK, r=7, marker='open',
             tips=[f'vertical: log-extrapolated Y/X match at ρ_s ≈ {d["rho_yx0"]:.1f} MΩ/sq;\n'
                   'the inclined angles need ≫ 10 MΩ/sq'])
    w = d['wing']
    for o, c in ((2, C_V1), (3, C_R5)):
        ys = [w[o]['s2'] / w[o]['data'], w[o]['s5'] / w[o]['data']]
        k = math.log(ys[1] / ys[0]) / math.log(2.5)
        lo_x = max(1.0, 2 * (B.ys[1] / ys[0]) ** (1 / k)) if k < 0 else 1.0
        hi_x = min(30.0, 2 * (B.ys[0] / ys[0]) ** (1 / k)) if k < 0 else 30.0
        xs = [lo_x, 2, 5, hi_x]
        B.line(xs, [ys[0] * (x / 2) ** k for x in xs], c, w=2, dash='4 6', markers=False)
        B.line([2, 5], ys, c, r=7, tips=[f'±{o} strips, ρ_s {x} MΩ/sq: sim {s:.3f} vs data {w[o]["data"]:.3f} → {y:.2f}'
                                         for x, s, y in zip((2, 5), (w[o]['s2'], w[o]['s5']), ys)])
        B.points([w[o]['rho']], [1.0], c, r=7, marker='open',
                 tips=[f'±{o} wing matches the data at ρ_s ≈ {w[o]["rho"]:.1f} MΩ/sq (log-log interpolation)'])
        B.text(1.08, ys[0] * (0.55 ** k) * 1.04, f'±{o} strips', 21, c)
    return A.svg('Y/X vs rho_s'), B.svg('Y wings vs rho_s')


def mesh_plot(nc0, nc1):
    rows = [('mesh transparency ε', nc1['surv'] / nc0['surv'],
             f'survival nc0 {nc0["surv"]:.4f} ± {nc0["surv_err"]:.4f} → nc1 {nc1["surv"]:.4f} ± {nc1["surv_err"]:.4f}'),
            ('gain | survived', nc1['gain'] / nc0['gain'],
             f'Polya mean of surviving seeds: {nc0["gain"]:.0f} → {nc1["gain"]:.0f} (θ {nc0["theta"]:.2f} → {nc1["theta"]:.2f});\n'
             'statistical error on the ratio ≈ ±1.2 %'),
            ('ε × gain (what the chain sees)', nc1['charge'] / nc0['charge'],
             f'mean charge per drifted electron {nc0["charge"]:.0f} → {nc1["charge"]:.0f}')]
    P = sd.Plot(900, 600, x=(-0.5, 2.5), y=(0.95, 1.05), xlabel='', ylabel='null-collision ÷ default stepping',
                margin=(24, 30, 70, 110))
    P.hline(1.0, INK, None, 2)
    P.xticks([(i, r[0]) for i, r in enumerate(rows)])
    P.yticks([(v, f'{v:.2f}') for v in (0.95, 0.97, 0.99, 1.0, 1.01, 1.03, 1.05)])
    for i, (lab, v, tip) in enumerate(rows):
        c = GREEN if i == 2 else (C_V1 if v > 1 else ORANGE)
        P.vbar(i, v, 150, c, tip=tip, base=1.0, label=f'×{v:.3f}')
    return P.svg('mesh nc1 vs nc0')


def build(R, out):
    L, S, d = R['L'], R['S'], derived(R)
    nc0, nc1 = R['nc0'], R['nc1']
    lo, hi = d['v1_rng']
    yx_v = [d['yx'][a]['verdict1'] for a in ANGLES]
    yx_lo, yx_hi = (min(yx_v) - 1) * 100, (max(yx_v) - 1) * 100
    v1_rise = [S[('verdict1', v, a)][1] for v in 'xy' for a in (10, 20)]
    wing_lo, wing_hi = sorted((d['wing'][2]['rho'], d['wing'][3]['rho']))
    yxm_lo, yxm_hi = min(d['yx_move']) * 100, max(d['yx_move']) * 100
    D = sd.Deck('MX17 response simulation: where it stands',
                'The MX17 Micromegas response simulation against the det3 target run on 2026-10-17: what now matches, '
                'what still disagrees, and ideas to close the rest.')

    # 1 cover ----------------------------------------------------------------
    D.slide('cover', sd.kicker('MX17 response simulation · det3 cosmics, 490 V · 2026-10-17')
            + '<h1 style="font-size:72px;font-weight:600;line-height:1.08;letter-spacing:-1px;max-width:1640px">'
              'The simulation now matches det3 in amplitude to within 15 % and in pulse shape to ~10 ns; '
              'what is left is how the charge splits between the X and Y strips</h1>'
            + sd.p('Five days ago it was ×1.6–2.4 too dim and 70 ns too slow. Two physics additions (fast hydrocarbon ions, '
                   'avalanche feedback) and six bug fixes closed that. The Y view and the neighbour strips are the open rows, '
                   'and none of the tuning has yet been tested on a second gas or chamber.', 30, DMUT)
            + sd.row(
                sd.bignum(f'{lo:.2f}–{hi:.2f}', 'peak amplitude, data ÷ sim', DGREEN,
                          'all angles, X and Y, verdict 1 (was 2.1–3.0)',
                          tip='Peak-strip median, data selected by true M3 angle, sim all events (ledger_m3.py).\n'
                              'Verdict 1: hydrocarbon ion, bG 0.45, shaper 166 ns, fixed noise.'),
                sd.bignum(f'+{min(v1_rise)}…+{max(v1_rise)} ns', 'inclined rise time, sim − data', DGREEN,
                          f'10° and 20°, X and Y (was +{LEDGER["rise_aug"]["x20"]}…+{LEDGER["rise_aug"]["x10"]} ns in X)',
                          tip='Median 10–90 % rise, peak strip, inclined tracks (shape_all.py).'),
                sd.bignum(f'+{yx_lo:.0f}…+{yx_hi:.0f} %', 'Y/X amplitude ratio still off', DRED,
                          'sim puts too little in Y; ρ_s does not fix it',
                          tip='(Y data/sim) ÷ (X data/sim) at the five angles, verdict 1.\n'
                              'A ρ_s of 5 MΩ/sq moves it by only ' + f'{yxm_lo:.0f}–{yxm_hi:.0f} %.'),
                gap=64),
            dark=True, short='Summary',
            notes='<p>Target run: det3 cosmic bench, Ar/iC4H10 95/5, 490 V mesh, 1000 V drift over 30 mm, the June 27 long run '
                  '(<code>long_run_resist_490V_drift_1000V</code>). Data are selected by the true angle of the M3 reference track, '
                  'with no det3 reco or quality cut; the sim takes all events at its fixed angle (FRESH_LOOK §18.2).</p>'
                  f'<p>The previous state, with the six open discrepancies, is the <a href="{STATUS_URL}">2026-10-12 note</a>. '
                  'The full record is <code>MX17_Geant/design/report/FRESH_LOOK_2026-10-12.md</code> §8–20.</p>')

    # 2 the chain ------------------------------------------------------------
    D.slide('chain', sd.title('Nine stages, one muon: five are anchored, three are nominal, one is where it fails',
                              'Colour = status today. Hover a box for what changed since 2026-10-12.')
            + chain_flow()
            + sd.row(sd.callout('Hover the dotted terms, the boxes and every point or bar for numbers and sources. '
                                f'An animated walk through the chain is <a href="{ANIM_URL}">here</a>.', BLUE, 26),
                     sd.legend([('checked against something outside the code', GREEN, 'box'),
                                ('nominal or fitted on this run', GOLD, 'box'),
                                ('the disagreement shows here', RED, 'box')], 24), gap=64),
            notes='<p>Each stage is a separate code: Geant4 (<code>mm_sim</code>) → Stage A clusters → Stage B digitiser '
                  '(<code>response/digitizer</code>: drift, mesh, avalanche sampling, ion signal, induction kernels) → '
                  'Stage C DREAM/FEU (<code>response/dream</code>) → <code>sim_decoded_*.root</code> → the same wft reconstruction as the data.</p>'
                  '<p>The offline inputs (gas table, transparency, avalanche calibration, kernels) are computed once by Magboltz, '
                  'Garfield++ and COMSOL. Per event, only Geant4 and the sampling run.</p>',
            short='The chain')

    # 3 scorecard ------------------------------------------------------------
    ok, part, bad = (f'<b style="color:{c}">{t}</b>' for c, t in ((GREEN, 'closed'), (GOLD, 'explained'), (RED, 'open')))
    rise20 = S[('verdict1', 'x', 20)][1]
    rows = [
        ['<b>Amplitude</b>', '×1.6–2.4 too dim', f'{lo:.2f}–{hi:.2f}', ok, 'hydrocarbon ion + feedback'],
        ['<b>Gain slope</b> d lnA/dV', f'{LEDGER["slope"]["sim"]:.2f} vs {LEDGER["slope"]["data"][0]:.2f}–{LEDGER["slope"]["data"][1]:.2f} /10 V',
         'same bG', part, 'fitted on the ladder shape'],
        ['<b>Rise</b> 10–90 %, 20°', f'+{LEDGER["rise_aug"]["x20"]} ns', f'X {rise20:+d}, Y {S[("verdict1", "y", 20)][1]:+d} ns', ok, 'ion speed + 166 ns shaper'],
        ['<b>FWHM</b>', '×1.18 wide', f'+{min(S[("verdict1", v, a)][3] for v in "xy" for a in (10, 20))}…'
         f'+{max(S[("verdict1", v, a)][3] for v in "xy" for a in (10, 20))} ns', ok, 'same'],
        ['<b>Multiplicity</b> strips ≥ 5σ', '12/13 vs 6/8', 'equal', ok, 'common-mode noise bug'],
        ['<b>Y/X</b> amplitude ratio', 'biased by the cut', f'+{yx_lo:.0f}…+{yx_hi:.0f} %', bad, 'not ρ_s (slide 8)'],
        ['<b>Y wings</b> ±2–3 strips', 'too wide', f'match at {wing_lo:.0f}–{wing_hi:.1f} MΩ/sq', part, 'outside the T2b ρ_s band'],
        ['<b>±1 neighbour delay</b>', '+135 vs +90 ns (X)', 'not re-measured', bad, 'not ion speed, not ρ_s'],
        ['Early foot, Y lobe', 'data only / sim deeper', 'unchanged', bad, 'lobe is real board physics'],
    ]
    D.slide('score', sd.title('Five of the nine rows are closed; the Y view and the neighbours are not',
                              'det3, 490 V. X = FEU 7, Y = FEU 8. "Then" is 2026-10-12, "now" the verdict-1 configuration.')
            + sd.table(['observable', 'then', 'now', 'status', 'how / why'], rows, 27,
                       widths=[440, 360, 330, 170, 360], align=['left'] * 5,
                       tips=['peak p50 data/sim, unbiased ledger (§18.2); then: §18.2 s=1 row and T14',
                             'feedback multiplies gain and slope by the same 1/(1−bG) (§9)',
                             'shape_all.py on the M3-selected data (§18.8)', 'same', '§18.4–18.5, A/B condor 4451433',
                             '§18.3, §20.2', '§18.10, §20.2 (profile_tails.py)', '§10.3–10.4, vertical tracks',
                             '§7.2, §11']),
            notes='<p>"Closed" means within the pre-registered band on the target run. It does not mean tested on a second run: '
                  'the ion speed, the feedback bG and the shaper were all read from this same run (slide 13).</p>'
                  '<p>The gain slope is "explained" rather than "closed": the feedback model fitted to the shape of the June HV '
                  'ladder gives both the slope and the amplitude with one parameter, but the decisive test (a ladder at a second '
                  'isobutane fraction) needs bench time.</p>',
            short='Scorecard')

    # 4 amplitude journey -----------------------------------------------------
    D.slide('amp', sd.title('Amplitude: two physics additions took the sim from ×2 too dim to the data',
                            'Peak-strip median, data ÷ sim, vertical tracks, data by true M3 angle. Each point adds one change.')
            + sd.row(amp_journey(L, d),
                     sd.col(sd.legend([('X (FEU 7)', C_X, 'line'), ('Y (FEU 8)', C_Y, 'line')], 24),
                            sd.callout('A ' + sd.term('hydrocarbon ion', 'The drifting ion in Ar/iC4H10 is a C4Hn⁺ cluster, '
                                                      '~1.5–1.85× faster than Ar⁺ (FRESH_LOOK §14.2, §15).')
                                       + ' makes the pulse faster, so more of the charge reaches the peak: ×1.2.', BLUE, 25),
                            sd.callout(sd.term('Feedback', 'G_eff = G/(1 − bG): each avalanche seeds further ones, '
                                                'most likely through photons freeing electrons at the mesh (§9).')
                                       + ' multiplies the gain by 1/(1 − bG) ≈ ×1.8. bG was fitted to the shape of the HV ladder, '
                                         'not to this amplitude.', BLUE, 25),
                            sd.callout('The corrected primaries (−9 %) then over-shoot at bG 0.45, so bG moves up to 0.47.', GOLD, 25),
                            w=540), gap=56),
            notes='<p>Points 1–3 are from the M3-angle ledger in FRESH_LOOK §18.2 (clusters 4450898/4450899, 180 ns shaper, '
                  'old noise model). Points 4–6 are the 2026-10-16 ledgers (<code>scratchpad/2026-10-16/results/ledger_1016_*.log</code>).</p>'
                  '<p>The T14 harness had selected the data leg with the det3 reco angle and the quality cut, which dropped bright '
                  f'events; on Y vertical it understated the data by ×{LEDGER["selection_y"]} (§18.1–18.2). '
                  'Every number here uses the unbiased selection.</p>',
            foot='det3 June 27 long run, 490 V; data |view angle| < 3°, |other| < 6°, fiducial M3 box; sim 3000 events per angle.',
            short='Amplitude')

    # 5 bG with corrected primaries ------------------------------------------
    D.slide('bg', sd.title(f'With the corrected primaries, X wants bG {d["bg_x"]:.2f} and Y wants {d["bg_y"]:.2f}',
                           'Mean over five angles of the peak median, data ÷ sim, at three feedback values (cluster 4472245).')
            + sd.row(bg_plot(d, L),
                     sd.col(sd.legend([('X', C_X, 'line'), ('Y', C_Y, 'line'), ('joint', INK, 'dash'),
                                       ('verdict 1 (old primaries)', MUT, 'dot')], 24),
                            sd.p(f'The joint optimum is <b>bG ≈ {d["bg_j"]:.2f}</b>, a small step from verdict 1\'s 0.45. '
                                 'The pre-registration said 0.49–0.53: the amplitude falls faster with bG than predicted.', 26),
                            sd.p(f'The gap between the X and Y optima <b>is the Y/X problem</b>: '
                                 f'one gain cannot serve both views until the split is fixed. A joint fit now just averages over it.', 26),
                            sd.callout('Use bG ≈ 0.47 for now; refit once Y/X has a cause.', GREEN, 26),
                            w=580), gap=56),
            notes='<p>Pre-registered in <code>scripts/condor/stageBC_w2_e17.sub</code>: bG 0.45 raises every p50 data/sim by '
                  '7–12 % (falsifier: < 4 % or > 15 %). Observed +5.2…+10.3 %, so it passes; the X saturated fraction falls '
                  '0.025–0.075 as predicted, and rise and FWHM stay within 1 ns of verdict 1.</p>'
                  '<p>Optima: linear interpolation in ln(data/sim) between bG 0.45 and 0.50. '
                  'Stage A: <code>mm_sim -g ArIso -p muon -e 4000 -n 3000 --beam-spread 31.2</code>, built at e17c72c+.</p>',
            short='Feedback bG')

    # 6 rise -----------------------------------------------------------------
    rp, series = rise_plot(S)
    D.slide('rise', sd.title('Pulse shape: the physical ion mobility reproduces the fitted one within 6 ns',
                             'Median 10–90 % rise of the peak strip, sim − data, inclined tracks. Ion rebuilt from Garfield '
                             'birth heights in the mesh field map, at three mobilities.')
            + sd.row(rp, sd.col(sd.legend([(s[0], s[1], 'box') for s in series], 23),
                                sd.p('K0 is the reduced mobility in cm²/(V·s). Verdict 1 compressed the Ar⁺ template ×1.75; '
                                     'that matches <b>K0 ≈ 2.2</b>.', 25),
                                sd.p('The literature value K0 2.05 is 5–6 ns slower than verdict 1. The rise is already '
                                     '+8…+16 ns over the data, so that costs a little agreement for a physical parameter.', 25),
                                sd.callout('Decision: adopt K0 2.05 as default, with 1.80 / 2.25 as the systematic?', GOLD, 25),
                                w=560), gap=40),
            notes='<p>Condor 4472241 (<code>stageBC_w2_ionK.sub</code>), <code>--ion-model reconstructed --ion-mobility K0=…</code>. '
                  'Pre-registered: K0 2.05 within +2…+8 ns of verdict 1 in rise and 0…+4 % in amplitude; falsifier > 10 ns or > 6 %. '
                  'Observed +5…+6 ns and +0.9…+2.2 %. FWHM is 1 ns over its 0…+5 band in three of four views.</p>'
                  '<p>The grey bars are the August default (Ar⁺, 180 ns shaper) on the T14 data legs (FRESH_LOOK §16), '
                  'available for X only.</p>',
            foot='Hover a bar for sim and data medians and the FWHM. Data inclined: X 10° n = 1202, 20° n = 70; Y 10° n = 782, 20° n = 33.',
            short='Rise time')

    # 7 Y profile --------------------------------------------------------------
    w2 = d['wing'][2]
    D.slide('profile', sd.title('The sim spreads Y over too many strips; a higher ρ_s narrows it to the data',
                                'Median amplitude ÷ peak by strip offset, vertical tracks. X for comparison.')
            + sd.row(profile_plot(R, 'y'), profile_plot(R, 'x'), gap=44)
            + sd.legend([('data', C_DATA, 'line'), ('sim ρ_s 2 MΩ/sq (verdict 1)', C_V1, 'line'),
                         ('sim ρ_s 5 MΩ/sq', C_R5, 'dash')], 25),
            notes=f'<p>At ±2 strips the Y amplitude is {w2["s2"]:.3f} of the peak in the sim at 2 MΩ/sq, {w2["s5"]:.3f} at 5 MΩ/sq, '
                  f'and {w2["data"]:.3f} in the data. X barely changes with ρ_s and already matched.</p>'
                  '<p><code>profile_tails.py</code> on the verdict-1 and r5verdict1 points (<code>profile_1016_*.log</code>). '
                  'The kernels were checked against COMSOL 3D on all Y channels (§18.10), so this is the input ρ_s, not a solver error.</p>',
            short='Y profile')

    # 8 rho_s does not fix Y/X ---------------------------------------------------
    pa, pb = rho_plots(d)
    D.slide('rho', sd.title('ρ_s explains the Y wings but not the Y/X ratio: two different effects',
                            'Verdict 1 (2 MΩ/sq) and the same chain at 5 MΩ/sq. Dashed: log-extrapolation. '
                            'Gold: the measured ρ_s band.')
            + sd.row(pa, pb, gap=44),
            notes='<p>Pre-registered in <code>scripts/condor/stageBC_w2_rho5M.sub</code>: Y/X rises 20–50 % and the Y-wing and Y/X '
                  'interpolations agree within a factor 1.5 in ρ_s. Falsifier: they disagree by more, or Y/X moves < 10 %. '
                  f'<b>Triggered</b>: Y/X moves {yxm_lo:.0f}–{yxm_hi:.0f} %, because X brightens with ρ_s too '
                  '(X data/sim falls 5–9 % instead of rising).</p>'
                  f'<p>The Y wings match at {wing_lo:.1f}–{wing_hi:.1f} MΩ/sq (±3 and ±2 strips). The vertical Y/X would need '
                  f'~{d["rho_yx0"]:.1f} MΩ/sq and the inclined angles far more. The 2026-10-14 point-source estimate (2.5 MΩ/sq for both) '
                  'is wrong once the electron cloud is included.</p>',
            short='ρ_s')

    # 9 mesh -------------------------------------------------------------------
    D.slide('mesh', sd.title('The mesh transparency was 4 % low, but correcting it changes the signal by only 0.7 %',
                             f'Garfield single electrons through the woven-mesh map, default vs null-collision stepping, '
                             f'{nc1["n"]} slices, {nc1["nev"]} seeds each (condor 4474926).')
            + sd.row(mesh_plot(nc0, nc1),
                     sd.col(sd.p('With default stepping, electrons jump past field gradients at the wires and are lost: '
                                 f'ε = {nc0["surv"]:.3f}. Null-collision steps give <b>{nc1["surv"]:.3f}</b>, as the G7 check predicted.', 26),
                            sd.p('The rescued electrons pass close to the wires and come with smaller avalanches, so the gain per '
                                 'surviving electron drops 3 %.', 26),
                            sd.callout('The chain thins by ε, then draws the gain: the product is what matters, and it moves '
                                       f'×{nc1["charge"] / nc0["charge"]:.3f}. No bG refit is needed.', GREEN, 26),
                            w=640), gap=56),
            notes='<p>Pre-registered in <code>mx17_aval_nc1.sub</code>: survival 0.990–0.997 (observed 0.9928, passes) and '
                  'gain|survived within ±3 % of 24 094 (observed ×0.970, on the edge, about 2.5σ).</p>'
                  '<p>Pooled file: EOS <code>response_sim/avalanche/aval_calib_meshfield_nc1_pooled.json</code> '
                  '(<code>scratchpad/2026-10-16/pool_nc1.py</code>). The chain still runs the nc0 file.</p>'
                  '<p>On lxplus this campaign first failed completely: the LCG_109 view ships its own Garfield headers, and '
                  'ROOT built our classes from them (ComponentGrid 688 vs 720 bytes), which corrupted memory. Fixed in '
                  '<code>setup_garfield.sh</code>; a same-setup rerun of an August calibration agrees within 1.4σ.</p>',
            short='Mesh')

    # 10 fixes ---------------------------------------------------------------
    def fix(num, head, body, c, tip):
        return sd.card(f'<p style="font-size:52px;font-weight:600;color:{c};line-height:1">{num}</p>'
                       + sd.p(f'<b>{head}</b>', 27) + sd.p(body, 23, MUT), tip=tip, pad=28)
    sc = LEDGER['sigma_cm']
    D.slide('fixes', sd.title('Six bugs and biases were found by checking the chain, not by tuning it',
                              'Each one changed what the comparison said. Hover a card for the source.')
            + sd.row(fix(f'×{LEDGER["selection_y"]}', 'Data selection',
                         'T14 picked the data with an amplitude-dependent cut; Y was understated by a third.', RED,
                         '§18.1–18.2: data now selected by true M3 angle, no det3 reco or quality cut.'),
                     fix(f'{sc[0]}→{sc[1]}', 'Common-mode noise [ADC]',
                         'Railed events inflated σ_cm 50×, making fake hits and wrong multiplicity.', RED,
                         '§18.5, response/dream/noise.py; A/B condor 4451433: amplitude −1 %, shape +3 ns.'),
                     fix('−11 %', 'Primary ionisation',
                         'Wrong gas densities, W value and missing PAI in Geant4 Stage A.', GOLD,
                         'e17c72c; nPrimDrift medians 0.898–0.915 of the old files; Heed agreement confirmed by Alexandra.'),
                     gap=32)
            + sd.row(fix('0.955→0.993', 'Mesh transparency', 'A Garfield stepping artefact (slide 9).', GOLD,
                         'G7 and condor 4474926; net effect on the signal +0.7 %.'),
                     fix('688≠720', 'lxplus Garfield [bytes]',
                         'ROOT read the wrong headers; jobs crashed or hung. August results check out.', GOLD,
                         'nTof_x17 faefdac; audit condor 4474927: 45.0k vs 42.8k / 46.4k at 490 V.'),
                     fix('0 → 1.6 %', 'Mesh-map inactive region',
                         'The lxplus Garfield lacked the region-flag patch: electrons crossed the wires.', GOLD,
                         'Rebuilt as lcg109-927e5c21p; check_grid_active.py guards every job.'),
                     gap=32)
            + sd.callout('None of these is a tuning parameter. Before them the Y amplitude was biased by a third and the '
                         'multiplicity rows were noise artefacts; after them the comparison can be trusted at the 5–10 % level.', GREEN, 27),
            notes='<p>None of these is a tuning parameter. They are why the comparison can now be trusted at the 5–10 % level: '
                  'before them, the amplitude ledger was biased by a third in Y and the multiplicity rows were noise artefacts.</p>',
            short='Fixes')

    # 11 what still disagrees --------------------------------------------------
    p25x = L[('verdict1', 'th00', 'x')]['p25']
    p50x = L[('verdict1', 'th00', 'x')]['p50']
    nb = LEDGER['nb_delay']
    rows = [
        ['<b>Y/X amplitude</b>', f'sim Y {yx_lo:.0f}–{yx_hi:.0f} % low vs X', 'every angle; not ρ_s', 'charge division, strip C, FEU 7 vs 8 gain'],
        ['<b>Y wings vs ρ_s</b>', f'need {wing_lo:.0f}–{wing_hi:.1f} MΩ/sq', f'T2b says {LEDGER["t2b"][0]}–{LEDGER["t2b"][1]}', 'ρ_s of det3 itself, or the stack-up'],
        ['<b>±1 neighbour delay</b>', f'X +{nb["x"][1]} vs +{nb["x"][0]} ns', 'Y the same, ×0.67', 'crosstalk or prompt sharing'],
        ['<b>Inclined rise</b>', f'+{min(v1_rise)}…+{max(v1_rise)} ns', 'Y 20° is the largest', 'shaper 158 vs 166 ns, ion z-lift'],
        ['<b>Low-side width</b>', f'X p25 {p25x:.2f} vs p50 {p50x:.2f}', 'data broader than sim', 'gain map over the area'],
        ['Early foot', '+0.04–0.16 of peak', 'data only, −360…−240 ns', 'not looked at since 10-12'],
        ['Y negative lobe', 'sim −0.2, data −0.11', 'COMSOL reproduces the shape', 'tied to the Y-kernel inputs'],
        ['ADC/fC scale', 'datasheet chain', 'TPOT bench −10 %', 'absorbed in bG today'],
    ]
    D.slide('open', sd.title('What still disagrees: mostly the Y view and the neighbour strips',
                             'Verdict-1 configuration unless noted. The first three share the Y/induction side of the chain.')
            + sd.table(['observable', 'size now', 'where', 'leading suspect'], rows, 27,
                       widths=[380, 420, 420, 440], align=['left'] * 4,
                       tips=['§18.3, §20.2', '§20.2 vs T2b (§18.10)', '§10.4: event-aligned on the central strip, before verdict 1',
                             'shape_all.py, verdict 1', 'ledger_m3, verdict 1 vertical X', '§7.2', '§7.2, §11',
                             '§19: Irakli TPOT 18.05 vs 20 mV/fC']),
            notes='<p>The neighbour delay was last measured on 2026-10-13, before verdict 1. §10.4 showed it does not depend on '
                  'the ion speed, so it is expected to be unchanged, but it should be re-measured on the verdict-1 output first.</p>'
                  '<p>The low-side width: the data p25/p50 ratio is lower than the sim\'s, so the data are broader on the low side. bG '
                  'hardly changes that (§18.7), which points to something the sim does not have, such as gain varying across the '
                  'chamber.</p>',
            short='Open rows')

    # 12 ideas ---------------------------------------------------------------
    rows = [
        ['<b>Compute the X and Y strip capacitances</b> from the board stack-up (2D field solve)',
         'Y/X', 'sim only', 'If C_Y < C_X by ~75–100 pF, the protection network\'s charge division explains the 15 % (§19).'],
        ['<b>Ask Irakli</b> for the protection network values (C_c, BOM)', 'Y/X', 'email', 'Needed for the same estimate.'],
        ['<b>Check the X–Y stack-up</b> (layer thickness, ε_r) against the board', 'Y/X, Y wings, Y lobe', 'drawings, photos',
         'The Y strips sit deeper; their induction depends on it. Peel photos of a board exist (design/figures).'],
        ['<b>Compare FEU 7 and FEU 8 gains</b> in pulser or calibration runs in the archive', 'Y/X', 'archive',
         'X and Y are on different FEUs; a 10 % gain difference looks exactly like Y/X.'],
        ['<b>Re-measure the neighbour delay</b> on verdict 1, then pulser crosstalk from the archive', '±1 delay', 'sim + archive',
         'Separates crosstalk from extra prompt sharing (§10.4).'],
        ['<b>Gain map vs M3 position</b> in data, applied in the sim', 'low-side width', 'data only', 'Tests gain non-uniformity.'],
        ['<b>Ar/CF₄ out-of-sample run</b> on the archived old-det3 data', 'everything tuned', 'condor', 'The first test not fitted on this run.'],
        ['<b>Measure det3\'s ρ_s</b> bus to bus', 'Y wings', 'ohmmeter, minutes', 'Decides whether ρ_s 4–5.6 is real.'],
    ]
    D.slide('ideas', sd.title('Ideas, roughly in order of information per effort',
                              'Most need no bench time. Y/X is first because it blocks the bG fit and three other rows.')
            + sd.table(['idea', 'closes', 'needs', 'why'], rows, 23, widths=[640, 220, 230, 570], align=['left'] * 4),
            notes='<p>The Y/X candidates are not exclusive. Charge division at a protected FEU input scales each view by '
                  'C_c/(C_c + C_strip); with C_c ≈ 510 pF (Irakli\'s NSW fit, rough) a 15 % Y/X needs C_Y roughly 75–100 pF below C_X. '
                  'MX17 strip capacitances have never been measured or computed.</p>'
                  '<p>Applying Irakli\'s −10 % ADC/fC would raise the needed gain, so bG goes up by ~0.03; it changes no shape row.</p>',
            short='Ideas')

    # 13 close -----------------------------------------------------------------
    def dl(items, c):
        return ''.join(f'<div style="border-left:5px solid {c};padding:4px 0 4px 22px">{sd.p(t, 27, DINK)}</div>' for t in items)
    D.slide('close', sd.kicker('What this does not rule out, and what needs deciding')
            + sd.row(sd.col(sd.p('<b>Not ruled out</b>', 32, DINK),
                            dl(['<b>One run.</b> The ion speed, bG and the shaper were all read from the det3 490 V long run. '
                                'Agreement here is a fit, not a prediction, until a second gas or chamber is checked.',
                                '<b>The gain scale.</b> ADC/fC is a datasheet chain; a 10 % error moves bG, not the shape.',
                                '<b>Wet gas.</b> The table is dry; the drift speed says the gas may not be (×0.77 gain if wet, §12.2).',
                                '<b>Feedback is a model.</b> G/(1 − bG) fits; the carrier (photons) is the leading candidate, '
                                'not a measurement.'], DRED), gap=18),
                     sd.col(sd.p('<b>Decisions</b>', 32, DINK),
                            dl(['Ion default: <code>--ion-model reconstructed --ion-mobility K0=2.05</code>, '
                                '1.80 / 2.25 as the systematic?',
                                'Next full chain pass ("verdict 2"): corrected primaries, bG 0.47, the null-collision '
                                'avalanche calibration, ρ_s 2 or 5 MΩ/sq?',
                                'Where to spend effort on Y/X: strip capacitance and protection values first (no hardware needed).'],
                               DGREEN), gap=18), gap=72),
            dark=True, short='Open questions')

    D.write(out, note_meta=dict(
        title='MX17 response simulation: where it stands (2026-10-17)',
        summary='The MX17 Micromegas response simulation against det3 cosmics: amplitude now matches to ~10 % and '
                'inclined pulse shape to ~10 ns; the Y/X charge split, the Y wings and the neighbour strips still disagree. '
                'What was fixed, what is open, and ideas to close it.',
        tags='MX17,simulation,micromegas,Garfield', date='2026-10-17'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    print('wrote', build(load(), a.out))


if __name__ == '__main__':
    main()
