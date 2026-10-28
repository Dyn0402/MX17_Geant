#!/usr/bin/env python3
"""Slide note: mechanical model of the MX17 readout-board bow under gas overpressure (2026-10-28), figure first.

Reads results/*.json (scenarios.py, real_edge.py, high_pressure.py, edge1d.py, window_gauge.py, validate.py) and the
bench measurement from nTof_x17 rerun/model_audit/followups/fb_bow_maps/results, so re-running those and then this
rebuilds every number and chart.

    ~/PycharmProjects/nTof_x17/.venv/bin/python mechanics/pcb_bow/make_mech_deck.py
    python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py mechanics/pcb_bow/deck/mx17-board-mechanics.html \
        --slug mx17-board-mechanics --force --deploy
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser(os.environ.get('SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd  # noqa: E402
from slidedoc import BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, INK, MUT  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')
FB = os.path.expanduser('~/PycharmProjects/nTof_x17/rerun/model_audit/followups/fb_bow_maps/results')
OUT = os.path.join(HERE, 'deck', 'mx17-board-mechanics.html')
SRC = 'MX17_Geant/mechanics/pcb_bow/'
term = sd.term
DETS = ['det2', 'det3', 'det4', 'det6', 'det7']
BAND = (1.2, 2.3)
LOOSE, BOND, BARE, WIN, MEAS = BLUE, GREEN, GREY, ORANGE, GOLD


def J(n, d=RES):
    return json.load(open(os.path.join(d, n)))


SC = {r['name']: r for r in J('scenarios.json')}
SC.update({r['name']: r for r in J('real_edge.json')})
SC.update({r['name']: r for r in J('high_pressure.json')})
WG, VAL, E1 = J('window_gauge.json'), J('validate.json'), J('edge1d.json')
SUMM = J('summary.json', FB)
PROF = {d: J(f'profiles_{d}.json', FB) for d in DETS}
SKEPTIC = dict(det2=1.86, det3=1.41, det4=1.53, det6=1.15, det7=1.54)
NOM_L, NOM_B = 'G_t1.6_loose_Keff', 'G_t1.6_bonded_sandwich_Keff'
NU_MY, T_MY, AH = 0.38, 60e-6, 0.2001


def cur(n, key='fit_depth_mm'):
    r = SC[n]['rows']
    return [x['p_mbar'] for x in r], [x[key] for x in r]


def at(n, p, key='fit_depth_mm'):
    ps, v = cur(n, key)
    return float(np.interp(p, ps, v))


def p_for(n, d):
    ps, v = cur(n)
    return float(np.interp(d, v, ps)) if v[-1] >= d else float('inf')


def win_sag(p_mbar, E=4.5e9):
    """membrane centre sag (mm) of the 60 um window; matches the solver to 0.1 mm (window_gauge.py)."""
    return (100 * p_mbar * AH ** 4 * (1 - NU_MY) * (0.8 + 0.062 * NU_MY) ** 3 / (E * T_MY)) ** (1 / 3) * 1e3


def ticks(vals, fmt='{:g}'):
    return [(v, fmt.format(v)) for v in vals]


def band(P, x0, x1, label=True):
    P.raw(f'<rect x="{P.X(x0):.1f}" y="{P.Y(BAND[1]):.1f}" width="{P.X(x1) - P.X(x0):.1f}" '
          f'height="{P.Y(BAND[0]) - P.Y(BAND[1]):.1f}" fill="{MEAS}" fill-opacity="0.16"'
          f'{sd.tipattr(MEAS_TIP)}/>', back=True)
    if label:
        P.text(x0 + 0.02 * (x1 - x0), BAND[1] - 0.12, 'measured dish, 5 chambers (incl. systematics)', 21, '#8a6a10')


MEAS_TIP = ('Bench dish, depth at 150 mm from the quadratic surface fit (nTof_x17 fb_bow_maps summary.json):\n'
            + '\n'.join(f'{d}: {SUMM[d]["depth_mm"]:.2f} ± {SUMM[d]["depth_err"]:.2f} mm (pivot estimator {SKEPTIC[d]:.2f})'
                        for d in DETS) + '\nBand 1.2–2.3 mm includes the ±0.4 mm per-chamber systematic.')


def tipsof(n, extra=''):
    out = []
    for x in SC[n]['rows']:
        out.append(f'{x["p_mbar"]:g} mbar: bench depth {x["fit_depth_mm"]:.2f} mm\ntrue dip at 150 mm '
                   f'{x["dip150_mm"]:.2f} mm, centre {x["centre_mm"]:.2f} mm{extra}')
    return out


# --------------------------------------------------------------------------- #
def s_cover(D):
    lo, hi = at(NOM_B, 3), at(NOM_L, 5)
    frac_lo = lo / 1.7; frac_hi = hi / 1.7
    n_sc = len(SC)
    pL = p_for('HG_t1.6_loose_Keff', 1.7)
    body = (sd.kicker('MX17 · readout board · mechanics · 28 Oct 2026') +
            '<h1 style="font-size:66px;font-weight:600;line-height:1.1;letter-spacing:-1.5px;color:#eef0f3">'
            f'A few mbar of gas bows the readout board by only {lo:.1f}–{hi:.1f} mm: '
            f'pressure is ~{frac_lo * 100:.0f}–{frac_hi * 100:.0f} % of the 1.2–2.3 mm dish we measure</h1>' +
            sd.row(sd.bignum(f'{n_sc}', 'mechanical scenarios', sd.DBLUE,
                             'PCB 1.2–2.4 mm × Rohacell loose / taped / bonded × edge × mesh tension × 0–100 mbar'),
                   sd.bignum('402 mm', 'free span', sd.DGREEN,
                             'the support-plate aperture, from the STEP file; the O-ring and frame sit outside it'),
                   sd.bignum(f'~{pL:.0f} mbar', 'needed for the measured 1.7 mm', sd.DRED,
                             f'1.6 mm board, Rohacell loose; the window would then bulge ~{win_sag(pL):.0f} mm'),
                   sd.bignum('5', 'chambers measured', sd.DGREY, 'det2, det3, det4, det6, det7, June cosmic bench'),
                   gap=44))
    notes = f"""
<p><b>Question.</b> On the cosmic bench every MX17 chamber shows its mesh centre 1.2–2.3 mm below a point 150 mm out,
with the cathode roughly flat (nTof_x17 <code>rerun/model_audit/followups/fb_bow_maps</code>). The chamber runs a few
mbar over atmosphere, which visibly bulges the 60 µm window. Does the same pressure bow the readout PCB by the measured
amount, and how much does the Rohacell taped behind the board help?</p>
<p><b>Answer.</b> No. A large-deflection two-layer plate model of the board over the 402 mm aperture gives, with the
mounting's own edge stiffness, a bench-equivalent dish of {at(NOM_B, 3):.2f}–{at(NOM_L, 3):.2f} mm at 3 mbar and
{at(NOM_B, 5):.2f}–{at(NOM_L, 5):.2f} mm at 5 mbar (bonded – loose Rohacell, 1.6 mm PCB). No combination of
thickness (1.2–2.4 mm), bond, edge or mesh tension reaches 1.2 mm below ~15 mbar. The measured profile shape is also not
a pressure shape. The rest of the dish needs another origin (board warp, mounting) or the bench depth is over-read.</p>
<p><b>What would settle it:</b> measure the window sag (it is a pressure gauge, slide 6), change the overpressure and
re-measure the tilt route, and put a straightedge on an unpressurised board.</p>
<p>Source: <code>{SRC}REPORT.md</code>; generator <code>{SRC}make_mech_deck.py</code>.</p>"""
    D.slide('cover', body, notes, dark=True, short='cover')


# --------------------------------------------------------------------------- #
def edge_schematic(w=1664, h=600):
    o = []
    x0mm, x1mm = 160.0, 240.0
    X = lambda v: 120 + (v - x0mm) / (x1mm - x0mm) * (w - 220)
    base, zs = 540, 11.5
    Z = lambda z: base - z * zs
    def rect(a, b, z0, z1, fill, tip, stroke='none', op=1.0):
        return (f'<rect x="{X(a):.1f}" y="{Z(z1):.1f}" width="{X(b) - X(a):.1f}" height="{Z(z0) - Z(z1):.1f}" '
                f'fill="{fill}" fill-opacity="{op}" stroke="{stroke}"{sd.tipattr(tip)}/>')
    # plate, foam, PCB, frame
    o.append(rect(201, 240, 0, 8, '#b9bfc8', '8 mm Al support plate. 402 × 402 mm through-aperture: its edge at 201 mm '
                                             'from the active-area axis is where the board stops having anything behind it.'))
    o.append(rect(x0mm, 240, 8, 13, '#efe3b0', '5 mm Rohacell (51 assumed), 470 × 470 mm, taped to the back of the board. '
                                              'Over the aperture it is a backing layer; outside it, a foundation on the plate '
                                              '(k = E/t ≈ 1.4×10¹⁰ N/m³, compression only).'))
    o.append(rect(x0mm, 240, 13, 14.6, '#3f8f5f', 'Readout PCB, 1.6 mm structural laminate (FR4 + 5 Cu layers + coverlay). '
                                                 'The CAD body (1.70 mm) also counts the mesh and the 150 µm gap.'))
    o.append(rect(x0mm, 240, 14.6, 14.9, '#7d8796', 'Micromesh on bulk pillars, 150 µm above the board.'))
    fr_tip = ('30 mm Al gas frame: 410 mm aperture (inner face at 205 mm), 15 mm wall to 220 mm. '
              'Its lower face pinches the board onto the foam over 205–206.6 and 210.3–220 mm.')
    o.append(rect(205, 206.55, 14.9, 40, '#cfd3d9', fr_tip, stroke=MUT))
    o.append(rect(206.55, 210.25, 18.9, 40, '#cfd3d9', fr_tip, stroke=MUT))
    o.append(rect(210.25, 220, 14.9, 40, '#cfd3d9', fr_tip, stroke=MUT))
    # break marks on the frame top
    o.append(sd.line(X(204), Z(40), X(221), Z(40), '#fff', 6))
    o.append(sd.T(X(221.5), Z(30), 'Al gas frame (30 mm, cut)', 21, INK, 'start'))
    # O-ring
    ortip = ("O-ring ('joint' in the STEP): 3.7 mm wide, 206.55–210.25 mm, seated in a groove in the frame's lower face; "
             'it seals on the board at ≈208.4 mm, so the gas pushes on the board out to there.')
    o.append(f'<ellipse cx="{X(208.4):.1f}" cy="{Z(16.9):.1f}" rx="{(X(210.25) - X(206.55)) / 2:.1f}" ry="{2.0 * zs:.1f}" '
             f'fill="#1c2230"{sd.tipattr(ortip)}/>')
    # screw
    sctip = 'Screws (M4-size in the STEP) at 214 mm on a 55 mm pitch, through window flange, frame and board.'
    o.append(rect(213.4, 214.6, 8, 40, '#5a6170', sctip))
    # pressure arrows
    for xa in np.linspace(163, 207.5, 12):
        o.append(sd.arrow(X(xa), Z(23.5), X(xa), Z(15.4), ORANGE, 3, 12))
    o.append(sd.T(X(176), Z(25.6), 'gas overpressure Δp, out to the O-ring seal', 22, ORANGE, 'middle', 600))
    # labels
    o.append(sd.T(X(162), Z(4), '← no aluminium behind the board for 402 mm (free span) →', 21, INK, 'start'))
    o.append(sd.T(X(225), Z(4), 'Al support plate', 21, INK, 'middle'))
    o.append(sd.T(X(235), Z(10.5), 'Rohacell', 21, INK, 'middle'))
    o.append(sd.T(X(229), Z(16.2), 'PCB + mesh', 21, INK, 'middle'))
    # dimension ticks
    for v, lab in ((201, '201'), (205, '205'), (208.4, '208.4'), (214, '214'), (220, '220')):
        o.append(sd.line(X(v), base + 6, X(v), base + 20, MUT, 2))
        o.append(sd.T(X(v), base + 46, lab, 21, MUT))
    o.append(sd.T(X(240), base + 46, 'mm from the active-area axis', 21, MUT, 'end'))
    return sd.svg(w, h, ''.join(o), 'cross-section of the MX17 board edge')


def s_setup(D):
    body = (sd.title('What holds the board: it spans the 402 mm plate aperture and is pinched just outside the O-ring',
                     'Edge cross-section from the CAD (MX17 2026 v1.3 STEP, global placements); heights ×4. '
                     'Hover the parts and dotted terms.') + edge_schematic())
    notes = f"""
<p><b>Geometry</b> (distances from the active-area axis, read from the STEP file with gmsh's OCC importer):</p>
<table><tr><th>mm</th><th>what</th></tr>
<tr><td>0 – 201</td><td>free span: the 8 mm Al plate has a 402 × 402 mm through-aperture; only the Rohacell is behind the board</td></tr>
<tr><td>201 – 205</td><td>board on the Rohacell over the plate, gas above (inside the 410 mm frame aperture)</td></tr>
<tr><td>205 – 206.55</td><td>inner lip of the frame face</td></tr>
<tr><td>206.55 – 210.25</td><td>O-ring groove; the O-ring ("joint", 3.7 mm wide) seals on the board, seal ≈ 208.4</td></tr>
<tr><td>210.25 – 220</td><td>frame face on the board; screws at 214 mm, 55 mm pitch</td></tr>
<tr><td>220 – 250</td><td>board overhang (the 470 mm board, Rohacell and plate are offset by +15, +15 mm)</td></tr></table>
<p>The active area (399.4 mm) is a few mm smaller than the aperture; the board (470 mm) is much larger, but nothing
outside 201 mm can bend much because the board is pressed onto the foam and the plate there, and pinched by the frame
beyond the O-ring. So the bending span is the plate aperture, not the board, the frame or the O-ring line. Slide 4
quantifies how firmly the edge is held.</p>
<p>The screws stop 1.5 mm into the Rohacell in the CAD (their tips at z = 267.7 mm, plate top at 266.2): the model does
not depend on what they thread into, only on the frame bearing on the board.</p>"""
    D.slide('setup', body, notes, short='geometry')


# --------------------------------------------------------------------------- #
def model_cartoon(w=860, h=560):
    o = []
    xa, xb = 40, 820
    o.append(sd.T(xa, 36, 'two layers, one deflection w(x, y), each with its own in-plane stretch', 22, INK, 'start', 600))
    # PCB
    o.append(f'<rect x="{xa}" y="90" width="{xb - xa}" height="46" fill="#3f8f5f" rx="6"/>')
    o.append(sd.T(xa + 20, 121, 'PCB  (D₁, A₁)', 24, '#fff', 'start', 600))
    # springs
    for xs in np.linspace(xa + 30, xb - 30, 14):
        pts_x = [xs, xs + 8, xs - 8, xs + 8, xs - 8, xs]
        pts_y = [136, 142, 150, 158, 166, 172]
        o.append(sd.poly(pts_x, pts_y, GOLD, 2.5))
    o.append(sd.T(xb + 4, 160, '', 20))
    o.append(f'<rect x="{xa}" y="172" width="{xb - xa}" height="110" fill="#efe3b0" rx="6"/>')
    o.append(sd.T(xa + 20, 236, 'Rohacell  (D₂, A₂)', 24, INK, 'start', 600))
    o.append(f'<rect x="{xa}" y="282" width="{xb - xa}" height="8" fill="#9aa3b0"/>')
    o.append(sd.T(xa + 20, 318, 'mylar skin 25 µm (bonded to the foam, or absent)', 20, MUT, 'start'))
    o.append(sd.T(xa + 20, 346, 'gold springs: interface shear k(x, y); 0 = loose, ∞ = glued, a map = a tape pattern', 21, '#8a6a10', 'start'))
    # legend of physics
    lines_ = ['large deflection (von Kármán): stretching stiffens the board once w ≈ t',
              'edge: rotational spring K from the real mounting (slide 4); in-plane held',
              'loads: Δp out to the seal, board weight on the bench, mesh-tension moment',
              'output: the bench observable, a quadratic surface fit, depth at 150 mm']
    for i, L in enumerate(lines_):
        o.append(f'<circle cx="{xa + 10}" cy="{400 + 42 * i}" r="6" fill="{INK}"/>')
        o.append(sd.T(xa + 30, 408 + 42 * i, L, 22, INK, 'start'))
    return sd.svg(w, h, ''.join(o), 'two-layer plate model')


def s_model(D):
    D_k0, D_kinf = VAL['interface_k0'], VAL['interface_k1e+12']
    rows = [
        ['linear plate, simply supported', f'{VAL["linear_N61"]["ss"] * 1e3:.3f}', '4.06', 'Timoshenko'],
        ['linear plate, clamped', f'{VAL["linear_N61"]["clamped"] * 1e3:.3f}', '1.26', 'Timoshenko'],
        ['interface k = 0 (D, N·m)', f'{D_k0:.2f}', f'{VAL["D_sum"]:.2f}', 'D₁ + D₂'],
        ['interface k → ∞ (D, N·m)', f'{D_kinf:.2f}', f'{VAL["D_composite"]:.2f}', 'composite section'],
        ['edge moment, SS (mm)', f'{abs(VAL["mesh_moment_ss_mm"]):.3f}', '1.323', 'exact (Poisson)'],
        ['membrane limit, 10 mbar (mm)', f'{VAL["nl_t0.2_q1000"]["w0"] * 1e3:.2f}', f'{VAL["nl_t0.2_q1000"]["membrane"] * 1e3:.2f}',
         'Vlassak–Nix'],
    ]
    tab = sd.table(['check', 'model', 'reference', 'source'], rows, size=23, widths=[330, 110, 120, 200],
                   align=['left', 'right', 'right', 'left'])
    body = (sd.title('The model: a stretching two-layer plate, checked against six closed-form results',
                     'plate.py, finite differences + Q4 membrane, Newton iteration with pressure continuation') +
            sd.row(model_cartoon(), sd.col(sd.p('Validation (validate.py)', 26, INK, 600), tab,
                                           sd.p('α×10³ = w₀D/(q a⁴); grid used for the scans: 41×41 (≈ 1 % error).', 21, MUT),
                                           w=760), gap=44))
    notes = f"""
<p><b>Energy.</b> Bending ½ΣD<sub>i</sub>[w<sub>xx</sub>² + w<sub>yy</sub>² + 2ν w<sub>xx</sub>w<sub>yy</sub> + 2(1−ν)w<sub>xy</sub>²]
(each layer about its own neutral surface); membrane ½Σ ε<sub>i</sub>ᵀA<sub>i</sub>ε<sub>i</sub> with von Kármán strains
ε = sym∇u<sub>i</sub> + ½∇w∇wᵀ; interface ½k(s<sub>x</sub>² + s<sub>y</sub>²), slip s = u₁ − u₂ + d∇w with d the distance
between the layers' neutral surfaces; loads −qw − m∇²w. Slip uses one-point quadrature (the 2×2 rule locks as k → ∞).</p>
<p><b>Materials.</b> FR4 20 GPa (17–24 as a variant), Cu 117 GPa smeared by each layer's in-window coverage (0, 0.76,
0.458, 0.457, 0.53), 50 µm coverlay; Rohacell 51: E 70 MPa (31: 36 MPa as a variant), 52 kg/m³; mylar 4.5 GPa.
D of the 1.6 mm laminate: {SC['A_t1.6_none']['D_pcb']:.1f} N·m; with Rohacell + skin bonded: {SC['A_t1.6_bonded_sandwich']['D_composite']:.1f} N·m.</p>
<p><b>Bench observable.</b> Each model surface is fitted with the bench basis [1, U, W, U²/2, W²/2, UW] (U, W in 100 mm)
over the efficient area 399 × 360 mm, and the depth is 0.5·k̄·1.5², exactly nTof_x17 <code>fb_bow_maps/summary.py</code>.
For a board held at its edges this reads about half of the true dip w(150) − w(0), because the shape is far from a
parabola (slide 5 shows both).</p>
<p>Files: <code>{SRC}plate.py</code>, <code>validate.py</code>, results in <code>results/validate.json</code>.</p>"""
    D.slide('model', body, notes, short='model')


# --------------------------------------------------------------------------- #
def s_edge(D):
    P = sd.Plot(1080, 680, x=(150, 235), y=(-3.2, 0.7), xlabel='distance from the active-area axis (mm)',
                ylabel='board height (mm, + toward the cathode)', margin=(24, 30, 92, 110))
    P.xticks(ticks([150, 170, 190, 201, 210, 220, 235])).yticks(ticks([-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5]))
    zones = [(201, 205, '#efe3b0', 'foam below, gas above'), (205, 206.55, '#cfd3d9', 'frame lip'),
             (206.55, 210.25, '#9aa3b0', 'O-ring groove'), (210.25, 220, '#cfd3d9', 'frame face + screws (214)')]
    for a, b, c, lab in zones:
        P.raw(f'<rect x="{P.X(a):.1f}" y="{P.y0}" width="{P.X(b) - P.X(a):.1f}" height="{P.ph}" fill="{c}" '
              f'fill-opacity="0.45"{sd.tipattr(f"{a}–{b} mm: {lab}")}/>', back=True)
    P.text(152, 0.5, 'free span (aperture)', 21, MUT)
    P.text(203, 0.55, 'clamp zone', 21, MUT)
    cols = {'bottomed_5mbar': (BLUE, 'frame bottomed, 5 mbar'), 'bottomed_10mbar': (GREEN, 'frame bottomed, 10 mbar'),
            'oring_10mbar': (RED, 'frame riding on the O-ring, 10 mbar')}
    for k, (c, lab) in cols.items():
        pr = E1['profiles'][k]
        x, w = np.array(pr['x_mm']), np.array(pr['w_um']) / 1e3
        P.line(list(x), list(w), c, 4, markers=False,
               tip=f'{lab}: K_eff = {pr["K_eff"]:.0f} N·m/m at the aperture edge')
    P.text(153, -3.0, '1-D cross-section, 1.6 mm PCB, 0.25 mm cells', 21, MUT)
    rows = []
    for t in (1.3, 1.6, 2.0):
        r1 = [r for r in E1['rows'] if r['t_mm'] == t and r['mode'] == 'bottomed' and r['q_Pa'] == 500 and r.get('F_or') == 2000.0][0]
        r2 = [r for r in E1['rows'] if r['t_mm'] == t and r['mode'] == 'oring' and r['q_Pa'] == 1000 and r.get('k_or') == 1e6][0]
        rows.append([f'{t} mm', f'{r1["K_eff"]:.0f}', f'{100 * r1["frac_to_pinned"]:.0f} %', f'{r2["K_eff"]:.0f}',
                     f'{100 * r2["frac_to_pinned"]:.0f} %'])
    tab = sd.table(['PCB', 'K bottomed', '→ pinned', 'K on O-ring', '→ pinned'], rows, size=22)
    fb = [float(r[2].rstrip(' %')) for r in rows]; fo = [float(r[4].rstrip(' %')) for r in rows]
    fb_lo, fb_hi, fo_lo, fo_hi = min(fb), max(fb), min(fo), max(fo)
    side = sd.col(sd.legend([(v[1], v[0]) for v in cols.values()], 22),
                  sd.p('Edge stiffness K (N·m/m) and how far the centre deflection moves from a perfect clamp '
                       '(0 %) toward a pinned edge (100 %):', 22, INK), tab,
                  sd.callout(f'Frame bottomed: the centre moves only {fb_lo:.0f}–{fb_hi:.0f} % of the way to a pinned edge '
                             f'({fo_lo:.0f}–{fo_hi:.0f} % if the frame rides on the O-ring at 10 mbar).', BLUE, 24), w=540, gap=20)
    body = (sd.title('The O-ring sits outside the span: the board is pinched there, so its edge is nearly clamped',
                     'Board height near the edge under pressure, fine 1-D model of the real mounting') +
            sd.row(P.svg('edge deflection'), side, gap=40))
    notes = f"""
<p><b>Model</b> (<code>{SRC}edge1d.py</code>): the board in cylindrical bending, 0.25 mm cells, from the centre to the
overhang end at 235 mm. Supports: Rohacell foundation outside 201 mm (k = E/t = 1.4×10¹⁰ N/m³, compression only, so the
board may lift); the frame face as a rigid one-sided contact over 205–206.55 and 210.25–220 mm ("bottomed"); the O-ring's
sealing reaction (2000 N/m assumed) pushing down at 208.4 mm; gas inside the seal. Alternative: the frame rides on the
O-ring (no metal contact), which is then a preloaded line spring (10⁶–10⁷ N/m²).</p>
<p><b>Result.</b> K_eff = M/θ at the aperture edge is 700 / 1110 / 1830 N·m/m for a 1.3 / 1.6 / 2.0 mm board with the
frame bottomed (K·a/D ≈ 50–90), almost independent of pressure, the sealing force (0 or 2000 N/m), a 20 µm foam
precompression (+40 %) or a softer foam grade (−15 %). The board lifts by < 40 µm anywhere outside the aperture. With the
frame riding on the O-ring the edge softens with pressure (K ≈ 460 N·m/m at 10 mbar) and the rim pivots about the seal,
lifting ~0.5 mm at the overhang.</p>
<p><b>In 2-D</b> the K_eff values are used as the spring edge of the plate model ("real edge" from here on). Earlier
estimate K = E b³/12t ≈ 3900 N·m/m (15 mm foam pinch) was ~3.5× too stiff; the dish changes by +5–7 %. Ideal clamped
vs simply supported, for scale: 0.55 vs 0.76 mm at 5 mbar (1.6 mm, loose Rohacell).</p>"""
    D.slide('edge', body, notes, foot='STEP geometry: frame 410/440 mm, O-ring 206.55–210.25 mm, screws at 214 mm, '
                                      'plate aperture 402 mm (all from the active-area axis).', short='edge')


# --------------------------------------------------------------------------- #
def s_pcb_window(D):
    pw = WG['p_mbar']
    A = sd.Plot(810, 600, x=(0, 10), y=(0, 14), xlabel='gas overpressure Δp (mbar)', ylabel='centre deflection (mm)',
                title='60 µm mylar window, 400 mm')
    A.xticks(ticks([0, 2, 4, 6, 8, 10])).yticks(ticks([0, 2, 4, 6, 8, 10, 12, 14]))
    ps = [0] + pw
    s38 = [0] + WG['rows']['3.8GPa']['solver_mm']; s50 = [0] + WG['rows']['5.0GPa']['solver_mm']
    s45 = [0] + WG['rows']['4.5GPa']['solver_mm']
    pd = list(np.linspace(0, 10, 101) ** 1) ; pd = [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75] + pd[9:]
    A.band(pd, [win_sag(p, 5.0e9) for p in pd], [win_sag(p, 3.8e9) for p in pd], WIN, 0.18, tip='mylar E = 3.8–5.0 GPa')
    A.line(pd, [win_sag(p) for p in pd], WIN, 4, markers=False)
    A.points(pw, s45[1:], WIN, 7, tips=[f'{p:g} mbar: window sag {s:.1f} mm (solver, E = 4.5 GPa)\nband: E 3.8–5.0 GPa'
                                        for p, s in zip(pw, s45[1:])])
    xs = [x['p_mbar'] for x in SC[NOM_L]['rows']]
    A.line(xs, [x['centre_mm'] for x in SC[NOM_L]['rows']], LOOSE, 3, markers=False)
    A.text(5.3, 2.3, 'board centre, same scale', 21, LOOSE)
    A.text(0.4, 12.6, 'sag ∝ Δp^⅓: 1 → 10 mbar only doubles it', 23, INK, weight=600)
    B = sd.Plot(810, 600, x=(0, 10), y=(0, 2.6), xlabel='gas overpressure Δp (mbar)', ylabel='board deflection (mm)',
                title='readout board, 1.6 mm, real edge')
    B.xticks(ticks([0, 2, 4, 6, 8, 10])).yticks(ticks([0, 0.5, 1, 1.5, 2, 2.5]))
    band(B, 0, 10)
    for n, c, lab in ((NOM_L, LOOSE, 'Rohacell loose'), (NOM_B, BOND, 'Rohacell bonded')):
        B.line(xs, [x['centre_mm'] for x in SC[n]['rows']], c, 3, dash='10 7', markers=False,
               tip=f'{lab}: true centre deflection')
        B.line(xs, [x['fit_depth_mm'] for x in SC[n]['rows']], c, 4, tips=tipsof(n, f'\n{lab}'))
    B.text(4.6, 0.12, 'solid: bench depth · dashed: true centre', 21, MUT)
    body = (sd.title(f'At 3 mbar the window bulges {win_sag(3):.1f} mm, the board only {at(NOM_B, 3):.1f}–{at(NOM_L, 3):.1f} mm',
                     'Deflection vs overpressure: window (solver, no pre-tension) and board (bench-equivalent depth and '
                     'true centre deflection)') +
            sd.row(A.svg('window sag'), B.svg('board bow'), gap=44) +
            sd.legend([('window, E = 4.5 GPa (band 3.8–5.0)', WIN), ('board, Rohacell loose', LOOSE),
                       ('board, Rohacell + skin bonded', BOND), ('measured board dish', MEAS, 'box')], 22))
    notes = f"""
<p><b>Window</b> (<code>{SRC}window_gauge.py</code>): the same solver in its membrane limit (60 µm PET, ν 0.38, clamped
and held on the 400.2 mm flange aperture, no pre-tension). It agrees with the Vlassak–Nix square-membrane law
Δp = E t w³ / ((1−ν)(0.8+0.062ν)³ a⁴) to 0.1 mm. Sag at 1 / 3 / 5 / 10 mbar: {', '.join(f'{v:.1f}' for v in WG['rows']['4.5GPa']['solver_mm'][:1] + WG['rows']['4.5GPa']['solver_mm'][2:4] + WG['rows']['4.5GPa']['solver_mm'][-1:])} mm.
A pre-tensioned (taut-glued) window sags less at low pressure, so this is an upper bound on sag for a given Δp, i.e. a
lower bound on Δp for a given sag.</p>
<p><b>Board</b>: 1.6 mm PCB with the mounting's edge stiffness, bench orientation (its weight adds ~0.35 mbar). Solid:
the bench observable (quadratic fit, depth at 150 mm); dashed: the true centre deflection. The board's response is also
sub-linear (membrane stiffening once w ≈ t), roughly ∝ Δp^0.5–0.6 at 5–10 mbar.</p>"""
    D.slide('pcbwin', body, notes, short='board & window')


# --------------------------------------------------------------------------- #
def s_gauge(D):
    P = sd.Plot(1100, 680, x=(5, 28), y=(0, 2.6), xlabel='window centre sag (mm); top ticks: the Δp (mbar) that produces it',
                ylabel='expected board dish, bench depth at 150 mm (mm)')
    P.xticks(ticks([6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28])).yticks(ticks([0, 0.5, 1, 1.5, 2, 2.5]))
    band(P, 5, 28, label=False)
    P.text(10.6, 2.2, 'measured dish, 5 chambers (incl. systematics)', 21, '#8a6a10')
    cases = [('HG_t1.3_loose_Keff', RED, '1.3 mm PCB, loose (softest)'), ('HG_t1.6_loose_Keff', LOOSE, '1.6 mm, loose'),
             ('HG_t1.6_bonded_sandwich_Keff', BOND, '1.6 mm, bonded'), ('HG_t2.0_bonded_sandwich_Keff', PURPLE, '2.0 mm, bonded (stiffest)')]
    for n, c, lab in cases:
        if n not in SC:
            continue
        g = n.replace('HG_', 'G_')
        rows = [x for x in SC[g]['rows'] if 0 < x['p_mbar'] <= 10] + [x for x in SC[n]['rows'] if x['p_mbar'] > 10]
        xs = [win_sag(x['p_mbar']) for x in rows]; ys = [x['fit_depth_mm'] for x in rows]
        P.line(xs, ys, c, 4, tips=[f'{lab}\nΔp {x["p_mbar"]:g} mbar → window {win_sag(x["p_mbar"]):.1f} mm, board '
                                   f'{x["fit_depth_mm"]:.2f} mm (true dip {x["dip150_mm"]:.2f})' for x in rows])
    for pm in (1, 3, 5, 10, 20, 40, 100):
        P.raw(sd.line(P.X(win_sag(pm)), P.y0, P.X(win_sag(pm)), P.y0 + 14, INK, 2))
        P.text(win_sag(pm), 2.40, f'{pm:g}', 20, INK, 'middle')
    P.raw(f'<rect x="{P.X(5):.1f}" y="{P.y0}" width="{P.X(win_sag(5)) - P.X(5):.1f}" height="{P.ph}" fill="{WIN}" '
          f'fill-opacity="0.08"/>', back=True)
    P.text(5.25, 1.0, '“a few mbar”', 21, WIN)
    d5max = max(at(n.replace('HG_', 'G_'), 5) for n, c, lab in cases)
    P.text(14.5, 0.42, f'a ~10 mm window (5 mbar) means ≤ {d5max:.1f} mm of board dish', 23, INK, weight=600)
    pl = p_for('HG_t1.6_loose_Keff', 1.7)
    side = sd.col(sd.legend([(lab, c) for n, c, lab in cases], 22), sd.legend([('measured dish', MEAS, 'box')], 22),
                  sd.callout(f'For the measured 1.7 mm from pressure alone: Δp ≈ {pl:.0f} mbar (1.6 mm, loose), '
                             f'a window bulge of ≈ {win_sag(pl):.0f} mm.', RED, 24),
                  sd.callout('Measure the window sag at the centre with a straightedge while flowing; '
                             'read the expected board dish off this curve.', ORANGE, 24), w=500, gap=22)
    body = (sd.title('Read the window, predict the board: the window sag is the pressure gauge',
                     'Expected board dish against the window sag the same overpressure produces, 1–100 mbar, real edge') +
            sd.row(P.svg('gauge'), side, gap=40))
    notes = f"""
<p>Each curve is one board scenario swept to 100 mbar (<code>{SRC}high_pressure.py</code>, real edge stiffness); its x
coordinate is the window sag at the same Δp (membrane law with E = 4.5 GPa, which the solver reproduces to 0.1 mm;
ticks along the bottom mark Δp). Because window sag ∝ Δp^⅓, a ±1 mm reading of a ~10 mm sag pins Δp to ±30 %, which
moves the predicted board dish by only ±15 %.</p>
<p>Overpressure needed for the measured dish from pressure alone (bench depth):</p>
<table><tr><th>scenario</th><th>1.2 mm</th><th>1.7 mm</th><th>2.3 mm</th></tr>
{''.join(f'<tr><td>{lab}</td>' + ''.join(f'<td>{p_for(n, d):.0f} mbar</td>' if math.isfinite(p_for(n, d)) else '<td>&gt;100</td>' for d in (1.2, 1.7, 2.3)) + '</tr>' for n, c, lab in cases if n in SC)}
</table>
<p>Caveats: a pre-tensioned window sags less, so a measured sag converts to a slightly higher Δp than shown; the bench
overpressure may differ chamber to chamber and with the flow.</p>"""
    D.slide('gauge', body, notes, short='window gauge')


# --------------------------------------------------------------------------- #
def s_thick(D):
    P = sd.Plot(1100, 680, x=(0, 10), y=(0, 2.6), xlabel='gas overpressure Δp (mbar)',
                ylabel='board dish, bench depth at 150 mm (mm)')
    P.xticks(ticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])).yticks(ticks([0, 0.5, 1, 1.5, 2, 2.5]))
    band(P, 0, 10)
    for t, dash in ((1.3, None), (1.6, '12 7'), (2.0, '3 7')):
        for bk, c in (('loose', LOOSE), ('bonded_sandwich', BOND)):
            n = f'G_t{t}_{bk}_Keff'
            xs, ys = cur(n)
            P.line(xs, ys, c, 4, dash=dash, markers=dash is None, tips=tipsof(n, f'\nPCB {t} mm, Rohacell {bk}'))
    P.text(0.3, 1.0, 'nothing reaches the band below 10 mbar', 24, INK, weight=600)
    side = sd.col(sd.legend([('Rohacell loose', LOOSE), ('Rohacell + mylar skin bonded', BOND)], 22),
                  sd.legend([('PCB 1.3 mm (solid)', INK, 'line'), ('1.6 mm (dashed)', INK, 'dash'), ('2.0 mm (dotted)', INK, 'dash')], 22),
                  sd.table(['PCB', 'loose', 'bonded'],
                           [[f'{t} mm', f'{at(f"G_t{t}_loose_Keff", 5):.2f}', f'{at(f"G_t{t}_bonded_sandwich_Keff", 5):.2f}']
                            for t in (1.3, 1.6, 2.0)], size=24),
                  sd.p('depth at 5 mbar (mm)', 21, MUT), w=460, gap=22)
    body = (sd.title('Board thickness and Rohacell bond move the bow by ×2 at most; none reaches the band',
                     'Bench-equivalent dish vs overpressure, real edge stiffness, bench orientation (board weight included)') +
            sd.row(P.svg('thickness x backing'), side, gap=40))
    rowsA = ''.join(f'<tr><td>{r["t"]}</td><td>{r["backing"]}</td><td>{r["D"]:.1f}</td><td>{r["d2"]:.2f}</td>'
                    f'<td>{r["d5"]:.2f}</td><td>{r["d10"]:.2f}</td></tr>'
                    for r in J('compare.json')['table_thickness_backing'])
    notes = f"""
<p>Thickness is the structural laminate (CAD 1.70 mm body minus mesh and 150 µm gap ≈ 1.5–1.6 mm; a 1.6 mm
standard board fits as well). The board's own weight gives 0.03–0.08 mm. At n_TOF the boards stand vertical: drop that
term (−0.03 to −0.05 mm at a few mbar).</p>
<p>Full grid with the earlier foam-pinch edge (K = 3900 N·m/m, ~5 % stiffer than the real edge), from
<code>results/compare.md</code>:</p>
<table><tr><th>PCB mm</th><th>backing</th><th>D N·m</th><th>2 mbar</th><th>5 mbar</th><th>10 mbar</th></tr>{rowsA}</table>
<p>FR4 modulus 17 / 24 GPa: +7 / −9 % (loose); Rohacell 31 instead of 51: +3 % loose, +13 % bonded.</p>"""
    D.slide('thick', body, notes, short='thickness × backing')


# --------------------------------------------------------------------------- #
def tile(kind, k, size=150):
    """mini map of where the Rohacell is bonded (green) inside the aperture."""
    o = [f'<rect x="0" y="0" width="{size}" height="{size}" fill="#f7f2df" stroke="{MUT}" stroke-width="1.5"/>']
    s = size / 0.402
    c = size / 2
    G_ = BOND
    if kind == 'none':
        pass
    elif kind == 'uniform':
        op = {1e8: 0.35, 1e9: 0.6, 1e10: 0.85, 1e12: 1.0}.get(k, 1.0)
        o.append(f'<rect x="0" y="0" width="{size}" height="{size}" fill="{G_}" fill-opacity="{op}"/>')
    elif kind == 'band':
        b = k * s
        o.append(f'<path d="M0 0H{size}V{size}H0Z M{b:.1f} {b:.1f}V{size - b:.1f}H{size - b:.1f}V{b:.1f}Z" '
                 f'fill="{G_}" fill-rule="evenodd"/>')
    elif kind == 'debond':
        o.append(f'<rect x="0" y="0" width="{size}" height="{size}" fill="{G_}"/>')
        o.append(f'<circle cx="{c}" cy="{c}" r="{k * s:.1f}" fill="#f7f2df"/>')
    elif kind == 'patch':
        o.append(f'<circle cx="{c}" cy="{c}" r="{k * s:.1f}" fill="{G_}"/>')
    elif kind == 'strips':
        wd, pitch = k
        xs = np.arange(-0.201, 0.2011, 0.0005)
        on = np.abs(((xs + pitch / 2) % pitch) - pitch / 2) < wd / 2
        i = 0
        while i < len(xs):
            if on[i]:
                j = i
                while j < len(xs) and on[j]:
                    j += 1
                o.append(f'<rect x="{(xs[i] + 0.201) * s:.1f}" y="0" width="{(xs[j - 1] - xs[i]) * s:.1f}" '
                         f'height="{size}" fill="{G_}"/>')
                i = j
            else:
                i += 1
    return ''.join(o)


def s_tape(D):
    C = J('compare.json')['table_bond']
    kinds = {'A_t1.6_loose': ('none', 0), 'B_psa_1e8': ('uniform', 1e8), 'B_psa_1e9': ('uniform', 1e9),
             'B_psa_1e10': ('uniform', 1e10), 'B_band25': ('band', 0.025), 'B_band60': ('band', 0.060),
             'B_debond_R100': ('debond', 0.100), 'B_debond_R150': ('debond', 0.150), 'B_patch_R100': ('patch', 0.100),
             'B_strips50_p100': ('strips', (0.050, 0.100)), 'B_strips25_p130': ('strips', (0.025, 0.134)),
             'A_t1.6_bonded_sandwich': ('uniform', 1e12)}
    short = {'A_t1.6_loose': 'loose', 'B_psa_1e8': 'soft foam tape, full face', 'B_psa_1e9': 'thin PSA, full face',
             'B_psa_1e10': 'stiff PSA, full face', 'B_band25': '25 mm edge band', 'B_band60': '60 mm edge band',
             'B_debond_R100': 'centre unstuck, r < 100', 'B_debond_R150': 'centre unstuck, r < 150',
             'B_patch_R100': 'centre patch only, r < 100', 'B_strips50_p100': 'strips 50 mm / 100 mm',
             'B_strips25_p130': 'strips 25 mm / 134 mm', 'A_t1.6_bonded_sandwich': 'glued, full face'}
    C = sorted(C, key=lambda r: r['frac'])
    size, gx, gy = 150, 268, 300
    o = []
    for i, r in enumerate(C):
        kx, ky = i % 6, i // 6
        x0, y0 = kx * gx + 40, ky * gy + 10
        kind, k = kinds[r['name']]
        tp = (f'{r["label"]}\ncaptures {100 * r["frac"]:.0f} % of the glued stiffening\n'
              f'bench depth at 2 / 5 / 10 mbar: {r["d2"]:.2f} / {r["d5"]:.2f} / {r["d10"]:.2f} mm')
        o.append(f'<g transform="translate({x0},{y0})"{sd.tipattr(tp)}>{tile(kind, k, size)}</g>')
        col = BOND if r['frac'] > 0.8 else (GOLD if r['frac'] > 0.5 else RED)
        o.append(sd.T(x0 + size / 2, y0 + size + 44, f'{100 * r["frac"]:.0f} %', 36, col, 'middle', 600))
        o.append(sd.T(x0 + size / 2, y0 + size + 76, short[r['name']], 20, INK, 'middle'))
        o.append(sd.T(x0 + size / 2, y0 + size + 102, f'{r["d5"]:.2f} mm at 5 mbar', 19, MUT, 'middle'))
    svg = sd.svg(1664, 620, ''.join(o), 'where the Rohacell is stuck')
    body = (sd.title('Where the Rohacell is stuck matters more than how strongly: bond the outer rim, the centre can let go',
                     'Bonded area (green) inside the 402 mm aperture → share of the full-glue stiffening it captures '
                     '(1.6 mm PCB, 5 mbar). Sorted.') + svg +
            sd.legend([('bonded to the PCB', BOND, 'box'), ('loose (contact only)', '#e6dcb8', 'box')], 22))
    notes = f"""
<p><b>Why.</b> Rohacell helps by turning the board into a sandwich: the PCB is one skin, the 25 µm mylar on the far face
of the foam the other, and the foam has to carry the shear between them. That shear flow is largest near the supports
and zero at the centre (as in any uniformly loaded beam), so a bond near the edge does most of the work: a 60 mm rim band
captures {100 * [r for r in C if r['name'] == 'B_band60'][0]['frac']:.0f} %, a centre that has come unstuck out to
r = 100 mm costs nothing ({100 * [r for r in C if r['name'] == 'B_debond_R100'][0]['frac']:.0f} %), while a
central patch alone captures only {100 * [r for r in C if r['name'] == 'B_patch_R100'][0]['frac']:.0f} %.</p>
<p><b>Size of the effect.</b> Loose → glued: {at('A_t1.6_loose', 5):.2f} → {at('A_t1.6_bonded_sandwich', 5):.2f} mm at
5 mbar (−{100 * (1 - at('A_t1.6_bonded_sandwich', 5) / at('A_t1.6_loose', 5)):.0f} %). Without the mylar skin (foam
glued, nothing on its far face) the gain is about half of that. Tape that only joins the foam's edge to the board
outside the aperture is the "loose" case. A bad bond therefore makes the bow at most ~1.5× the well-bonded value.</p>
<p><b>Practical:</b> if re-taping, cover the outer 60 mm all round with a continuous tape (or full-face transfer
adhesive), and make sure the back mylar is bonded to the foam: without that second skin the sandwich has nothing to work
with. Patterns were run with the earlier foam-pinch edge (K = 3900 N·m/m); the fractions are insensitive to the edge.</p>
<p>Interface model: shear spring k(x, y) = 10¹² N/m³ (glued) inside the green area, 0 elsewhere; full-face tapes
k = G/t: 10⁸ (soft foam tape, ~0.1 MPa over 1 mm), 10⁹ (acrylic PSA 0.1 MPa / 100 µm), 10¹⁰.</p>"""
    D.slide('tape', body, notes, short='where to tape')


def s_adhesive(D):
    ks = [0, 1e7, 1e8, 1e9, 1e10, 1e12]
    Dk = [VAL[f'interface_k{k:g}'] for k in ks]
    d0, d1 = Dk[0], Dk[-1]
    P = sd.Plot(1100, 660, x=(1e6, 1e13, 'log'), y=(0, 105), xlabel='interface shear stiffness k = G/t of the bond (N/m³)',
                ylabel='share of the glued stiffening captured (%)')
    P.xticks([(10 ** e, f'10{sd.sup(e)}') for e in range(6, 14)]).yticks(ticks([0, 20, 40, 60, 80, 100]))
    for (lo, hi, lab, y) in ((3e7, 3e8, 'foam tapes', 12), (5e8, 2e10, 'acrylic PSA / transfer tape', 30),
                             (3e11, 1e13, 'epoxy', 48)):
        P.raw(f'<rect x="{P.X(lo):.1f}" y="{P.y0}" width="{P.X(hi) - P.X(lo):.1f}" height="{P.ph}" fill="{GREY}" '
              f'fill-opacity="0.12"/>', back=True)
        P.text(math.sqrt(lo * hi), y, lab, 20, MUT, 'middle')
    xs = [1e6] + ks[1:]
    ys = [100 * (d - d0) / (d1 - d0) for d in Dk]
    P.line(xs, ys, BOND, 4, tips=[f'k = {sd.sci(k) if k else 0} N/m³: D_eff = {d:.2f} N·m (linear plate, 1.7 mm PCB)'
                                  for k, d in zip(ks, Dk)])
    C = {r['name']: r for r in J('compare.json')['table_bond']}
    for n, k in (('B_psa_1e8', 1e8), ('B_psa_1e9', 1e9), ('B_psa_1e10', 1e10)):
        P.points([k], [100 * C[n]['frac']], ORANGE, 10,
                 tips=[f'{C[n]["label"]}: {100 * C[n]["frac"]:.0f} % at 5 mbar (large deflection, 1.6 mm)'])
    P.text(1.3e6, 92, 'any full-face adhesive tape gets ≥ 85 %', 24, INK, weight=600)
    side = sd.col(sd.legend([('linear plate, D_eff', BOND), ('full model at 5 mbar', ORANGE, 'dot')], 22),
                  sd.callout('The foam itself is soft in-plane (70 MPa), so even a weak bond can drag it along. '
                             'Coverage, not grip, is what can go wrong.', BOND, 24), w=500, gap=22)
    body = (sd.title('How strongly it sticks hardly matters: a soft tape is already 85 % of glue',
                     'Board stiffness gained from a full-face bond vs the bond\'s shear stiffness') +
            sd.row(P.svg('adhesive'), side, gap=40))
    notes = f"""
<p>Green: effective flexural rigidity of PCB (1.7 mm) + Rohacell + mylar skin, clamped, linear, from
<code>validate.py</code> check 3: k = 0 gives D₁ + D₂ = {d0:.2f} N·m, k = 10¹² gives {d1:.2f} N·m (composite section
{VAL['D_composite']:.2f}). Orange: the full large-deflection model at 5 mbar (1.6 mm PCB, slide 8).</p>
<p>A rough partial-composite estimate (Newmark/γ-method) gives γ = 1/(1 + π² E₂A₂/(k L²)): with the foam's in-plane
stiffness E₂A₂ ≈ 0.5 MN/m and L ≈ 0.2 m, γ > 0.9 for k > 3×10⁸ N/m³, as the curve shows.</p>"""
    D.slide('adhesive', body, notes, short='bond strength')


# --------------------------------------------------------------------------- #
def s_mesh(D):
    P = sd.Plot(1100, 660, x=(0, 10), y=(0, 1.2), xlabel='gas overpressure Δp (mbar)', ylabel='board dish, bench depth (mm)')
    P.xticks(ticks([0, 2, 4, 6, 8, 10])).yticks(ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]))
    series = [('G_t1.6_loose_Keff', LOOSE, None, 'loose, no mesh tension'),
              ('M_t1.6_loose_T10', LOOSE, '3 6', 'loose, T = 10 N/m'),
              ('M_t1.6_loose_T1000', RED, None, 'loose, T = 10 N/cm'),
              ('C_loose_ss', GREY, '12 7', 'loose, SS edge, no tension'),
              ('M_t1.6_loose_T1000_ss', PURPLE, '12 7', 'loose, SS edge, T = 10 N/cm'),
              ('G_t1.6_bonded_sandwich_Keff', BOND, None, 'bonded, no tension'),
              ('M_t1.6_bonded_sandwich_T1000', GOLD, None, 'bonded, T = 10 N/cm')]
    for n, c, dash, lab in series:
        if n not in SC:
            continue
        xs, ys = cur(n)
        P.line(xs, ys, c, 4 if dash is None else 3, dash=dash, markers=False, tip=lab)
    d10 = at('M_t1.6_loose_T10', 3) - at('G_t1.6_loose_Keff', 3)
    d1k = at('M_t1.6_loose_T1000', 3) - at('G_t1.6_loose_Keff', 3)
    dss = at('M_t1.6_loose_T1000_ss', 3) - at('C_loose_ss', 3)
    P.text(0.3, 1.1, f'10 N/cm: +{d1k * 1e3:.0f} µm at 3 mbar (real edge)', 24, INK, weight=600)
    side = sd.col(sd.legend([(lab, c, 'dash' if dash else 'line') for n, c, dash, lab in series if n in SC], 21),
                  sd.table(['at 3 mbar', 'added dish'],
                           [['T = 10 N/m, real edge', f'{d10 * 1e3:+.1f} µm'], ['T = 10 N/cm, real edge', f'{d1k * 1e3:+.0f} µm'],
                            ['T = 10 N/cm, SS edge', f'{dss * 1e3:+.0f} µm']], size=22), w=520, gap=20)
    body = (sd.title(f'Mesh tension (~10 N/cm) adds only {d1k * 1e3:.0f}–{(at("M_t1.6_loose_T1000", 0) - at("G_t1.6_loose_Keff", 0)) * 1e3:.0f} µm with the real edge',
                     'Board dish vs overpressure with the mesh-tension moment, 1.6 mm PCB') + sd.row(P.svg('mesh'), side, gap=40))
    notes = f"""
<p><b>Mechanism.</b> The mesh is laminated on the board (bulk pillars). Its tension T is balanced inside the board by an
equal compression, so the pair carries no net in-plane force: no stress stiffening or softening of the bend. What is
left is the moment of the mesh about the board's neutral surface, m = T·e with e ≈ 0.8 mm (PCB neutral surface 0.66 mm
below the top + 0.15 mm gap). It cups the board toward the mesh side (centre away from the cathode, the measured sign).
A uniform moment does nothing to a perfectly clamped plate; with the real edge stiffness a little of it gets through.</p>
<p><b>Numbers.</b> T = 10 N/m gives m ≈ 0.008 N·m/m: +{d10 * 1e3:.1f} µm. T = 10 N/cm (1000 N/m, the range quoted for
stretched meshes) gives +{d1k * 1e3:.0f} µm on the real edge and +{dss * 1e3:.0f} µm on a simply supported one, close
to additive with the pressure (the two couple only through the small membrane stiffening).</p>
<p><b>Floating-mesh case</b> (mesh held by the frame and resting on pillars, not laminated): its tension then acts on
the mesh alone. A dish of curvature κ needs a pressure T·κ to keep the mesh on the pillars, against the electrostatic
pull ε₀E²/2 ≈ 48 Pa at 500 V over 150 µm. At 10 N/m, T·κ ≈ 0.6 Pa: the mesh follows the board. At 10 N/cm it would
not follow beyond ~0.5 mm of dish. Not our case if the pillars are bulk, but it is the one way mesh tension could matter.</p>
<p>Files: <code>{SRC}real_edge.py</code> (M_ scenarios), <code>results/real_edge.json</code>.</p>"""
    D.slide('mesh', body, notes, short='mesh tension')


# --------------------------------------------------------------------------- #
POS = (np.linspace(30, 370, 10)[:-1] + np.linspace(30, 370, 10)[1:]) / 2 - 200.0


def shape(n, ax):
    if n == 'parabola':
        return POS / 150.0
    rows = SC[n]['rows']
    i = int(np.argmin([abs(r['fit_depth_mm'] - 1.7) for r in rows])) if not n.startswith('F_') else 0
    return np.array(rows[i]['slope'][ax])


def wfit(y, e, s):
    W = 1 / e ** 2
    X = np.c_[np.ones(len(y)), s]
    b = np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (W * y))
    return b, float(np.sum(W * (y - X @ b) ** 2))


def s_shape(D):
    SH = [('parabola', GREY, 'parabolic dish'), ('G_t1.6_loose_Keff', BLUE, 'pressure, real edge'),
          ('C_loose_ss', GREEN, 'pressure, simply supported'), ('F_edge_rot_m5mrad', PURPLE, 'imposed edge rotation')]
    chi = {n: 0.0 for n, _, _ in SH}
    panels = []
    for i, ax in enumerate('xy'):
        rowp = []
        for d in DETS:
            pr = PROF[d][ax]
            y, e = np.array(pr['w0']), np.array(pr['w0_err'])
            lo, hi = float(np.min(y - e)), float(np.max(y + e))
            pad = 0.1 * (hi - lo)
            P = sd.Plot(330, 330, x=(-180, 180), y=(lo - pad, hi + pad), margin=(44, 10, 50 if i else 20, 60),
                        xlabel='u − centre (mm)' if i else '')
            P.xticks([(-150, '−150'), (0, '0'), (150, '150')] if i else [])
            P.yticks([(round(v, 1), f'{v:.1f}'.replace('-0.0', '0.0')) for v in np.linspace(lo, hi, 3)])
            for n, c, lab in SH:
                s = shape(n, ax)
                b, x2 = wfit(y, e, s)
                if (d, ax) != ('det7', 'x'):
                    chi[n] += x2
                P.line(list(POS), list(b[0] + b[1] * s), c, 3, markers=False, tip=f'{lab}: χ² = {x2:.0f} (7 dof)')
            P.points(list(POS), list(y), INK, 6,
                     tips=[f'{d} {ax.upper()}, {p + 200:.0f} mm: w0 = {a:.3f} ± {b_:.3f} µm/ns' for p, a, b_ in zip(POS, y, e)])
            for p_, a, b_ in zip(POS, y, e):
                P.raw(sd.line(P.X(p_), P.Y(a - b_), P.X(p_), P.Y(a + b_), INK, 2))
            P.raw(sd.T(P.x0 + 6, P.y0 - 14, f'{d} {ax.upper()}' + (' (faulty half)' if (d, ax) == ('det7', 'x') else ''),
                       22, INK, 'start', 600))
            rowp.append(P.svg(f'{d} {ax}'))
        panels.append(sd.row(*rowp, gap=4))
    best = min(chi, key=chi.get)
    leg = sd.legend([(f'{lab} (Σχ² {chi[n]:.0f})', c) for n, c, lab in SH] + [('measured w0', INK, 'dot')], 21)
    body = (sd.title('The measured tilt steepens toward the rim; a pressure-loaded, edge-held board does the opposite',
                     'Drift-tilt offset w0 (µm/ns) along each view, fitted with each model\'s slope shape (+ a rigid tilt)') +
            sd.col(*panels, gap=0) + leg)
    notes = f"""
<p>w0 = v·r_f·dh/du is the solid route of the bench study (fit errors ±0.02–0.1 µm/ns). For each view: weighted fit
w0 = c + α·s(u), with s the model's mean slope in the same 9 bins (other coordinate averaged over the efficient area),
taken at the pressure where that model's dish is 1.7 mm. Σχ² over 9 views (det7 X, with its faulty half, excluded),
63 dof: {', '.join(f'{lab} {chi[n]:.0f}' for n, c, lab in SH)}.</p>
<p>An edge-held plate's slope turns over before the rim (curvature reverses near a clamp); the data steepen toward it,
most clearly on det2 Y, det3 X/Y, det4 X/Y and det6 Y. None of the shapes fits (best χ²/dof ≈ {chi[best] / 63:.0f}).
Part of that is reconstruction near the strip ends (the skeptic review flagged own-axis edge effects). Read as
geometry, curvature concentrated near the rim points to a load at the edge (mounting, frame) or a built-in board shape,
not to a uniform pressure. The amplitude α, read through each shape with a flat cathode (r_f = ½), implies a dish of
2.8–3.0 mm.</p>"""
    D.slide('shape', body, notes, short='shape test')


# --------------------------------------------------------------------------- #
def s_close(D):
    lo5, hi5 = at(NOM_B, 5), at(NOM_L, 5)
    body = (sd.kicker('summary · what this does not rule out · what to check') +
            '<h2 style="font-size:54px;font-weight:600;line-height:1.12;color:#eef0f3">'
            f'Pressure gives {lo5:.1f}–{hi5:.1f} mm at 5 mbar. The rest of the 1.2–2.3 mm needs another cause, '
            'or ~20–60 mbar of overpressure.</h2>' +
            sd.row(sd.col(sd.p('Not ruled out', 30, '#eef0f3', 600),
                          sd.p('• the bench overpressure is really ≥ 20 mbar (the window would bulge ≥ 16 mm)', 24, sd.DGREY),
                          sd.p('• built-in board warp: normal fab tolerance 0.75 % ≈ 3 mm over 400 mm; it would survive '
                               'vertical mounting at n_TOF', 24, sd.DGREY),
                          sd.p('• a mounting load at the edge (frame flatness, screw torque, cables, M1 cards); the rim-'
                               'steepening shape points there', 24, sd.DGREY),
                          sd.p('• the bench depth over-read: ±0.4 mm reconstruction systematic and ≤ 0.6 mm M3 bow put the '
                               'low end at 0.8–1.2 mm', 24, sd.DGREY), w=800, gap=14),
                   sd.col(sd.p('Cheap checks', 30, '#eef0f3', 600),
                          sd.p('1. Window sag with a straightedge while flowing: ~8–10 mm → 3–5 mbar; ~20 mm → ~40 mbar', 24, sd.DGREY),
                          sd.p('2. Double the overpressure, re-measure the w0 route: a pressure dish grows ×1.5, a warp '
                               'does not move', 24, sd.DGREY),
                          sd.p('3. Straightedge or dial gauge on an open, unpressurised board', 24, sd.DGREY),
                          sd.p('4. The fab stack-up (t enters as t³) and how the Rohacell and its mylar are taped', 24, sd.DGREY),
                          w=800, gap=14), gap=64))
    notes = f"""
<p><b>Model limits.</b> Kirchhoff kinematics per layer (no core shear in the Rohacell; estimated &lt; 1 % at this span);
41×41 grid (~1 %); the board stays in contact with the foam over the aperture (true under pressure); isotropic FR4 and an
ideal square aperture; no initial warp, so the model gives the change from the unloaded shape; the window has no
pre-tension.</p>
<p><b>Open questions for the hardware.</b> Actual overpressure on the bench (bubbler height or outlet restriction),
any observed window sag, the PCB stack-up from the fab, how the Rohacell is taped (full face, strips or edge) and
whether the back mylar is bonded to it, and whether the frame bottoms on the board or rides on the O-ring.</p>
<p>Everything: <code>{SRC}</code> (REPORT.md, results/compare.md, figures/).</p>"""
    D.slide('close', body, notes, dark=True, short='summary')


def main():
    D = sd.Deck('MX17 board mechanics', 'Mechanical model of the MX17 readout-board bow under gas overpressure, '
                                        'Rohacell bonding, edge mounting and mesh tension, vs the bench dish')
    for f in (s_cover, s_setup, s_model, s_edge, s_pcb_window, s_gauge, s_thick, s_tape, s_adhesive, s_mesh, s_shape, s_close):
        f(D)
    p = D.write(OUT, note_meta=dict(
        title='MX17 board mechanics: does gas pressure bow the readout board?',
        summary='Large-deflection two-layer plate model of the readout PCB + Rohacell over the 402 mm aperture: a few '
                'mbar gives 0.3–0.6 mm of dish, a third of the 1.2–2.3 mm measured. Where to tape the Rohacell, the '
                'window as a pressure gauge, the real edge mounting, mesh tension.',
        tags='x17, mx17, model audit, mechanics, Rohacell', date='2026-10-28'))
    print('wrote', p)


if __name__ == '__main__':
    main()
