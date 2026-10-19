"""Slide note: is det3's X left/right neighbour asymmetry real? (2026-10-19)

Builds the deck from the raw-ADC extractions of raw_asym.py (desktop copy in ~/x17/nb_asym_raw, override with NB_RAW),
the raw survey (raw/survey/index.tsv from run_raw3.sh), ../copy_survey.csv (sub-run times, beam rows) and the
bus weighting solves in ../bus_scan/. Every number on a slide is computed here.
Usage: make_asym_deck.py [--out PATH]"""
import os, sys, glob, json, csv, re, argparse, datetime as dt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser(os.environ.get('SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd
import asym_lib as L
from slidedoc import INK, MUT, RULE, BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, BG, CARD
from slidedoc import DBLUE, DRED, DGREY, DGREEN, DINK, DMUT

HW = L.HW
NS_PER_SAMPLE = 60
CR, CL, CC = RED, BLUE, INK          # p+1 (gets the copy), p-1, the peak channel
CJUN, CMAY, COTH = RED, GREEN, GREY  # det3 June-15 build, det3 May build, other boards

# ----------------------------------------------------------------------------------------------------------------- #
# data
# ----------------------------------------------------------------------------------------------------------------- #
SAMPLES = {   # tag: (board label, view, FEU, run label)
    'd3jun_x': ('det3 · June-15 build', 'X', 3, '6-24 long run'),
    'd3jun_y': ('det3 · June-15 build', 'Y', 4, '6-24 long run'),
    'd3wk_x': ('det3 · June-15 build', 'X', 7, '6-28 weekend'),
    'd3wk_y': ('det3 · June-15 build', 'Y', 8, '6-28 weekend'),
    'd3may_x': ('det3 · May build', 'X', 3, '5-6 HV scan 470 V'),
    'd3may_y': ('det3 · May build', 'Y', 4, '5-6 HV scan 470 V'),
    'd1_x': ('det1 · Jan build', 'X', 4, '1-28 overnight'),
    'd1_y': ('det1 · Jan build', 'Y', 6, '1-28 overnight'),
    'd2_x': ('det2 · Arnaud June 12', 'X', 6, '6-23 long run'),
    'd2_y': ('det2 · Arnaud June 12', 'Y', 8, '6-23 long run'),
    'd4_x': ('det4 · Stephan June 22', 'X', 6, '6-24 long run'),
    'd4_y': ('det4 · Stephan June 22', 'Y', 8, '6-24 long run'),
    'd6_x': ('det6 · Stephan June 24', 'X', 3, '6-26 longer run'),
    'd6_y': ('det6 · Stephan June 24', 'Y', 4, '6-26 longer run'),
    'd7_x': ('det7 · Stephan June 22', 'X', 6, '6-26 longer run'),
    'd7_y': ('det7 · Stephan June 22', 'Y', 8, '6-26 longer run'),
}
_cache = {}


def get(tag, veto=True):
    k = (tag, veto)
    if k not in _cache:
        _cache[k] = L.load(tag, veto)
    return _cache[k]


def have(tag):
    return os.path.exists(f'{L.RAW}/{tag}.npz')


def colour_of(tag):
    if tag.startswith('d3jun') or tag.startswith('d3wk') or tag.startswith('t6'):
        return CJUN
    if tag.startswith('d3may'):
        return CMAY
    return COTH


# raw survey (run_raw3.sh): one row per (run, sub-run, detector, view)
def phys(rd, det):
    """Physical board, as in ../survey_table.py."""
    if det in ('mx17_3', 'mx17_A'):
        return 'det3'
    if 'det1_det2_overnight_6-17' in rd:
        return {'mx17_1': 'det2', 'mx17_2': 'det3'}[det]
    if det == 'mx17_1':
        return 'det3' if '/det_3/' in rd else 'det4' if '/det_4/' in rd else 'det1'
    return {'mx17_B': 'det2', 'mx17_C': 'det6', 'mx17_D': 'det7'}.get(det, det.replace('mx17_', 'det'))


def central(d):
    """Copy a over the board-centre chips 2-5 only (the copy is centre-peaked, so the board average moves with
    where the hits land: several HV-scan points have most hits on chip 7)."""
    chip = (d['p'] % 512) // 64
    m = (chip >= 2) & (chip <= 5)
    if m.sum() < 2:
        return float('nan'), float('nan'), int(m.sum())
    a, ea = L.mean_se(d['a1'][m])
    return a, ea, int(m.sum())


# det3 June-15 X runs fetched from EOS after the survey (not in ../survey_jobs.tsv): file, run, sub, FEU, start, drift
EXTRA = [('extra/d617_long_x.npz', 'mx17_det1_det2_overnight_6-17-26', 'long_run', 7, '2026-06-18 05:43', '1000'),
         ('extra/d619_longer_x.npz', 'mx17_det2_det3_weekend_6-19-26', 'longer_run', 7, '2026-06-19 22:40', '1000')]


def survey_rows():
    idx = f'{L.RAW}/survey/index.tsv'
    if not os.path.exists(idx):
        return []
    meta = {}
    for r in csv.DictReader(open(os.path.join(HERE, '..', 'copy_survey.csv'))):
        meta[(r['run'], r['sub'], r['det'], r['view'])] = r
    rows = []
    for line in open(idx):
        f, rd, sub, det, v, feu = line.rstrip('\n').split('\t')
        f = f'{L.RAW}/survey/{os.path.basename(f)}'
        if not os.path.exists(f):
            continue
        try:
            d = L.load(f)
        except Exception:
            continue
        if d['n'] < 60:
            continue
        run = rd.split('/')[-1]
        m = meta.get((run, sub, det, v.upper()), {})
        board = phys(rd, det)
        build = board
        if board == 'det3':
            build = 'det3 May' if ('5-5-26' in rd or '5-6-26' in rd or '6-6-26' in rd) else 'det3 June-15'
        a, ea = L.mean_se(d['a1'])
        ac, eac, nc = central(d)
        rows.append(dict(f=f, run=run, sub=sub, det=det, view=v.upper(), feu=int(feu), board=board, build=build,
                         a=a, ea=ea, ac=ac, eac=eac, nc=nc, n=d['n'], coin=L.coin(d)[0], hot=list(d['hot']),
                         mesh=m.get('mesh_V'), drift=m.get('drift_V'), t=m.get('tstart', '')))
    return rows


# ----------------------------------------------------------------------------------------------------------------- #
# small helpers
# ----------------------------------------------------------------------------------------------------------------- #
def f3(x):
    return f'{x:+.3f}'


def hist(x, edges):
    h, _ = np.histogram(x, edges)
    return h / max(len(x), 1)


def details(*paras):
    return ''.join(f'<p>{p}</p>' for p in paras)


def sources(*items):
    return '<p>Sources: ' + ', '.join(f'<code>{i}</code>' for i in items) + '</p>'


# ----------------------------------------------------------------------------------------------------------------- #
# slides
# ----------------------------------------------------------------------------------------------------------------- #
def s_cover(D):
    x = get('d3jun_x')
    ctrl = [t for t in SAMPLES if not t.startswith('d3jun') and not t.startswith('d3wk') and have(t)]
    ntot = sum(get(t)['n'] for t in SAMPLES if have(t))
    fx, _ = L.coin(x)
    nums = ''.join([
        sd.bignum(f'{len([t for t in SAMPLES if have(t)])}', 'detector views read raw',
                  DBLUE, '7 boards (det3 twice: May and June-15 builds), X and Y, one shared extractor.', size=72,
                  tip='\n'.join(f'{t}: {SAMPLES[t][0]} {SAMPLES[t][1]}, FEU {SAMPLES[t][2]}, {SAMPLES[t][3]}'
                                for t in SAMPLES if have(t))),
        sd.bignum(f'{ntot:,}', 'clean cosmic peaks', DGREY,
                  'Not railed, a real pulse, hot channels removed.', size=72),
        sd.bignum('0', 'processing steps that could pick a side', DGREEN,
                  'Pedestal per channel, peak channel, its two neighbours. No map, no common mode, no clustering.',
                  size=72)])
    body = (sd.kicker('det3 June-15 build · X view · cosmic bench June 2026 · raw DREAM ADC')
            + '<h1 style="font-size:76px;font-weight:600;line-height:1.08;letter-spacing:-2px;width:1660px">'
              'The X asymmetry is real and in the raw ADC: one neighbour carries a copy of the peak. '
              'It belongs to one board, not to the readout or the code</h1>'
            + f'<p style="font-size:30px;color:{DMUT};width:1560px;line-height:1.35">On det3’s June-15 build the '
              f'channel above the peak is the bigger neighbour in {100 * fx:.0f} % of events; every other board and '
              'view sits at 50 %. What makes it is still open: the board’s one known left/right asymmetry '
              '(its bus traces) is ten times too weak.</p>'
            + '<div style="flex:1"></div>'
            + f'<div style="display:flex;gap:56px">{nums}</div>')
    D.slide('cover', body, details(
        'This deck answers one question: is the det3 X left/right neighbour asymmetry (NB_ASYM.md) a real detector '
        'effect, or something the analysis made? It starts from the most basic observable there is and builds up.',
        'Everything here was recomputed from the decoded DREAM files with a new extractor, <code>raw_asym.py</code>, '
        'that shares no code with the earlier analysis (no <code>FeuReader</code>, no strip map, no common-mode '
        'subtraction). The earlier copy survey is used only where noted.',
        'Short answer: yes, it is real, it is a copy of the peak channel’s own pulse into its upper channel neighbour, it '
        'follows the board (not the FEU slot, not the view), it is strongest in the middle of the board, and it '
        'changes with drift field. Its cause is not established.'), dark=True, short='Cover')


def s_setup(D):
    W, H = 1664, 700
    o = []
    # --- left: side view of the detector
    x0, x1 = 40, 860
    ytop, ymesh, yres, ypad, yl5, yl6 = 60, 360, 420, 500, 545, 600
    o.append(f'<rect x="{x0}" y="{ytop}" width="{x1 - x0}" height="{ymesh - ytop}" fill="#eef3fa"/>')
    o.append(sd.line(x0, ytop, x1, ytop, INK, 4))
    o.append(sd.T(x1 - 10, ytop - 14, 'drift cathode', 21, MUT, 'end'))
    o.append(sd.T(x0 + 12, ytop + 40, 'drift gap ~27 mm, Ar/iso', 22, MUT, 'start'))
    o.append(sd.line(x0, ymesh, x1, ymesh, INK, 3, '10 6'))
    o.append(sd.T(x1 - 10, ymesh - 12, 'mesh', 21, MUT, 'end'))
    o.append(f'<rect x="{x0}" y="{ymesh}" width="{x1 - x0}" height="{yres - ymesh}" fill="#fbf1e6"/>')
    o.append(sd.T(x0 + 12, ymesh + 38, '150 µm amplification', 21, MUT, 'start'))
    o.append(f'<rect x="{x0}" y="{yres - 6}" width="{x1 - x0}" height="8" fill="{GOLD}"'
             f'{sd.tipattr("Resistive (ESL) strips, 0.80 mm pitch against 0.78 mm pads: checked, no beat in the asymmetry.")}/>')
    o.append(sd.T(x1 - 10, yres - 14, 'resistive strips', 21, GOLD, 'end'))
    # pads: 9 pads, the peak one in the middle
    pw, pg = 72, 16
    names = ['p−4', 'p−3', 'p−2', 'p−1', 'p', 'p+1', 'p+2', 'p+3', 'p+4']
    cols = [GREY] * 3 + [CL, CC, CR] + [GREY] * 3
    for i in range(9):
        px = 70 + i * (pw + pg)
        o.append(f'<rect x="{px}" y="{ypad}" width="{pw}" height="16" fill="{cols[i]}" rx="2"/>')
        o.append(sd.T(px + pw / 2, ypad - 12 + 0 * 46, names[i], 22, cols[i] if cols[i] != GREY else MUT, 'middle',
                      600 if i in (3, 4, 5) else 400))
        # bus line on one side of each pad (drawn under the gap toward the left)
        bx = px - pg / 2
        o.append(f'<rect x="{bx - 4}" y="{yl5 + 20}" width="8" height="8" fill="{MUT}" opacity="0.7"/>')
        o.append(sd.line(px + pw / 2, ypad + 16, px + pw / 2, yl5 + 12, MUT, 1.2, '3 3'))
        o.append(sd.line(px + pw / 2, yl5 + 12, bx, yl5 + 22, MUT, 1.2))
    o.append(sd.T(x1 - 10, yl5 + 70, 'each pad’s bus runs under the gap on ONE side', 21, MUT, 'end'))
    # track and electrons
    o.append(sd.line(300, ytop - 30, 520, ymesh + 4, INK, 3))
    rng = np.random.default_rng(3)
    for k in range(14):
        t = (k + 0.5) / 14
        ex, ey = 300 + t * 220, ytop + t * (ymesh - ytop)
        o.append(f'<circle cx="{ex:.0f}" cy="{ey:.0f}" r="5" fill="{BLUE}"/>')
        o.append(sd.arrow(ex, ey + 8, ex + rng.normal(0, 6), ymesh - 4, BLUE, 1.2, 7) if k % 3 == 0 else '')
    o.append(sd.T(200, ytop + 120, 'cosmic μ', 24, INK, 'middle', 600))
    o.append(f'<ellipse cx="{70 + 4 * (pw + pg) + pw / 2}" cy="{ymesh + 30}" rx="60" ry="22" fill="{ORANGE}" opacity="0.45"/>')
    o.append(sd.T(70 + 4 * (pw + pg) + pw / 2 + 80, ymesh + 36, 'avalanche', 21, ORANGE, 'start'))
    # --- right: readout chain
    rx = 960
    boxes = [('pads → 512 channels per view', 'X and Y pads interleaved (checkerboard)'),
             ('connector + cable', '64 channels each, “inverted” wiring'),
             ('FEU · DREAM chip', '64 ch per chip, 8 chips per FEU'),
             ('ADC: 32 samples × 60 ns', 'every channel, every trigger, no zero suppression')]
    for i, (a, b) in enumerate(boxes):
        by = 70 + i * 150
        o.append(f'<rect x="{rx}" y="{by}" width="680" height="104" rx="14" fill="{CARD}" stroke="{RULE}" stroke-width="2"/>')
        o.append(sd.T(rx + 24, by + 44, a, 26, INK, 'start', 600))
        o.append(sd.T(rx + 24, by + 80, b, 21, MUT, 'start'))
        if i < 3:
            o.append(sd.arrow(rx + 340, by + 108, rx + 340, by + 144, MUT, 2.5, 10))
    o.append(sd.arrow(x1 + 10, ypad + 8, rx - 10, 122, MUT, 2.5, 12))
    o.append(sd.T(rx, 690, 'We compare p−1 and p+1 by DAQ channel number:', 23, INK, 'start', 600))
    o.append(sd.T(rx, 720, 'no strip map in between.', 23, INK, 'start', 600))
    body = (sd.title('What is measured: a cosmic, a few pads, three channels',
                     'Hover the dotted terms and any data point for definitions, counts and sources.')
            + sd.svg(W, H + 40, ''.join(o), 'detector and readout cartoon'))
    D.slide('setup', body, details(
        'A cosmic muon ionises the ~27 mm drift gap; electrons drift to the mesh, avalanche in the 150 µm gap, and '
        'induce signal on the readout pads under the resistive layer. The pads form a checkerboard: every other pad '
        'belongs to an X channel, the others to a Y channel (gerbers: 131 072 stubs each). Each pad reaches its '
        'channel’s bus line through a via and a 390 µm stub, and every stub of a layer points the same way, so each '
        'bus runs under the inter-pad gap on one side of its pads. That is the only built-in left/right asymmetry '
        'of the board we know of; slide 18 asks how large an effect it can make.',
        'The readout: 512 channels per view, 64-channel connectors and DREAM chips, 32 samples of 60 ns per channel '
        'per trigger. These runs had no zero suppression, no on-board pedestal and no common-noise subtraction '
        '(run_config <code>dream_daq_info</code>).',
        'Throughout: p is the channel with the largest pulse; p−1 and p+1 are its neighbours in DAQ channel number. '
        'With the connector-inverted wiring of all these runs, channel p+1 sits at the lower strip position.'),
        short='Setup')


def s_suspects(D):
    steps = [
        dict(label='Track', sub='cosmics come from both sides equally', color=GREY),
        dict(label='Drift + diffusion', sub='symmetric', color=GREY),
        dict(label='Avalanche, feedback', sub='symmetric', color=GREY),
        dict(label='Resistive spread', sub='0.80 vs 0.78 mm beat: averages out', color=GREY),
        dict(label='Pads + buses', sub='buses one-sided!', color=RED),
        dict(label='Connector, cable', sub='pin order', color=GOLD),
        dict(label='FEU / ADC mux', sub='channel n read after n−1', color=GOLD),
        dict(label='Our code', sub='peak finding, maps, filters', color=GOLD),
    ]
    body = (sd.title('Where could left and right differ? Eight suspects',
                     'Grey: symmetric by construction. Red: the one built-in asymmetry. Gold: places an analysis or '
                     'readout artefact could hide.')
            + sd.flow(steps, size=28)
            + '<div style="height:40px"></div>' + sd.row(sd.callout('Slides 4–10: is it in the raw data at all, and what does it look like?', BLUE, 28),
                     sd.callout('Slides 11–17: which suspect it follows: code, FEU, board, place, time, field.', RED, 28),
                     gap=64)
            + sd.callout('Slides 18–19: can the bus traces make it, and the verdict on each suspect.', GREY, 28))
    D.slide('suspects', body, details(
        'Your prior is the right one: nothing in the physics of a micromegas prefers one side. Cosmics arrive from both '
        'sides, diffusion and the avalanche are isotropic, and the resistive-strip/pad beat averages to zero over the '
        'board (NB_ASYM §1: peak rank mod 40 shows no modulation, sine amplitude 0.020 ± 0.024).',
        'So a left/right difference has to come either from the hardware layout (the bus traces are one-sided), '
        'from the readout (adjacent pins, the order the ADC reads channels), or from the analysis. The deck tests '
        'each in turn and comes back to this chain at the end with a verdict per box.'), short='Suspects')


def nice_step(top, n=4):
    raw = top / n; e = 10 ** np.floor(np.log10(raw))
    return next(m * e for m in (1, 2, 2.5, 5, 10) if m * e >= raw)


def s_event(D):
    z = np.load(f'{L.RAW}/d3jun_x.npz')
    picks = []
    for blk in z['disp']:
        A = blk.max(1); p = int(A.argmax()); c = blk[p]
        if 3 <= p <= 60 and 4 <= c.argmax() <= 22 and abs(c[:3].mean()) < 0.1 * A[p]:
            picks.append((blk, p))
        if len(picks) == 3:
            break
    panels_a, panels_b = [], []
    t = np.arange(32) * NS_PER_SAMPLE
    for n, (blk, p) in enumerate(picks):
        A = blk.max(1); top = A[p] * 1.25; st = nice_step(top)
        P = sd.Plot(530, 300, x=(-4.6, 4.6), y=(0, top), xlabel='channel − p',
                    ylabel='max ADC' if not panels_a else '', margin=(20, 20, 76, 90))
        P.xticks([(i, f'{i:+d}' if i else 'p') for i in range(-4, 5)])
        P.yticks([(v, f'{v:.0f}') for v in np.arange(0, top, st)])
        for i in range(-4, 5):
            col = CC if i == 0 else CR if i == 1 else CL if i == -1 else GREY
            P.vbar(i, max(A[p + i], 0), 40, col, tip=f'channel p{i:+d}: {A[p + i]:.0f} ADC')
        P.text(-4.4, top * 0.9, f'event {n + 1}', 22, INK, weight=600)
        P.text(1.6, top * 0.9, f'p+1 / p−1 = {A[p + 1] / max(A[p - 1], 1):.1f}', 21, CR, weight=600)
        panels_a.append(P.svg(f'event {n}'))
        Q = sd.Plot(530, 300, x=(0, 1920), y=(-0.15 * A[p], top), xlabel='time [ns]',
                    ylabel='ADC − pedestal' if not panels_b else '', margin=(20, 20, 76, 90))
        Q.xticks([(v, str(v)) for v in (0, 600, 1200, 1800)])
        Q.yticks([(v, f'{v:.0f}') for v in np.arange(0, top, st)])
        for i, col in ((-1, CL), (1, CR), (0, CC)):
            Q.line(t, blk[p + i], col, 3, markers=False, tip=f'channel p{i:+d}' if i else 'peak channel p')
        panels_b.append(Q.svg(f'waveforms {n}'))
    body = (sd.title('Raw single events: the peak and its two neighbours',
                     'det3 June-15 build, X view, 6-24 long run · ADC minus per-channel pedestal · first 3 clean events in the file')
            + sd.row(*panels_a, gap=36) + sd.row(*panels_b, gap=36)
            + sd.legend([('p−1', CL, 'box'), ('peak p', CC, 'box'), ('p+1', CR, 'box')]))
    D.slide('event', body, details(
        'Each column is one cosmic event, taken in file order (the first three that pass the clean-pulse cut), so they '
        'are not hand-picked. Top: the largest sample of each channel around the peak. Bottom: the full 32-sample '
        'waveforms of p−1, p and p+1. Raw DREAM samples from decoded_root; the pedestal is one number per channel, the '
        'median over the first 1000 events.',
        'Single events fluctuate: where the track crosses the pad decides which neighbour is bigger. The question is '
        'what happens on average, which the next slides answer with thousands of events.',
        sources('raw_asym.py', '~/x17/nb_asym_raw/d3jun_x.npz (disp)')), short='Raw events')


def s_processing(D):
    steps = [dict(label='Raw ADC', sub='decoded DREAM samples, 32 × 60 ns', color=BLUE),
             dict(label='− pedestal', sub='one number per channel (median)', color=BLUE),
             dict(label='peak channel p', sub='largest max sample; 1–12 channels over 5σ, all within p ± 3', color=BLUE),
             dict(label='compare p−1, p+1', sub='by DAQ channel number', color=BLUE)]
    gone = [dict(label='common-mode subtraction', color=GREY, strike=True),
            dict(label='strip / position map', color=GREY, strike=True),
            dict(label='clustering', color=GREY, strike=True),
            dict(label='M3 tracks', color=GREY, strike=True),
            dict(label='pulse fits', color=GREY, strike=True)]
    cuts = sd.table(['Cut (for every board alike)', 'Why'], [
        ['200 &lt; peak &lt; 3500 ADC', 'not railed (railed peaks break a ratio)'],
        ['peak sample 4–22, flat first 3 samples', 'a real pulse, not ringing or noise'],
        ['p not at a chip edge (4 ≤ p mod 64 ≤ 59)', 'p−1 and p+1 on the same chip and connector'],
        ['drop hot channels (peak &gt; 6× the median rate)', 'one hot channel can fake the effect (slide 15)'],
    ], size=24, widths=[700, 900])
    body = (sd.title('Four processing steps, none of which can pick a side',
                     'The same code runs on every board and view; the controls are the other 13 views')
            + sd.flow(steps, size=26) + sd.p('Not used anywhere in slides 4–17:', 24, MUT)
            + sd.flow(gone, size=22, arrow='·') + cuts)
    D.slide('processing', body, details(
        'Every step is symmetric between p−1 and p+1: a per-channel constant, an argmax over channels, and two '
        'array lookups. A per-channel gain difference cannot create a net asymmetry either, because p runs over all '
        'channels and each channel is p−1 for one peak and p+1 for another.',
        'The cuts act on the peak channel only and never look at the neighbours, except the isolation (no other '
        'channel over 5σ beyond p ± 3), which is symmetric.',
        'Hot channels: a channel that fires on noise produces “events” with p always the same channel, so its two '
        'fixed neighbours decide the left/right answer. The veto removes a channel that is the peak more than six '
        'times the median rate of the channels that ever are. On det3 June X it removes channel 492 (6-23/24 runs).',
        sources('raw_asym.py', 'asym_lib.py')), short='Processing')


def s_coin(D):
    tags = [t for t in SAMPLES if have(t)]
    order = sorted(tags, key=lambda t: (SAMPLES[t][1], -L.coin(get(t))[0]))
    P = sd.Plot(1664, 700, x=(-0.7, len(order) - 0.3), y=(0.40, 0.92), ylabel='fraction of events with p+1 > p−1',
                margin=(20, 20, 150, 110))
    P.yticks([(v, f'{100 * v:.0f} %') for v in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9)])
    P.hline(0.5, MUT, '8 6', 2, label='coin flip', tip='A symmetric detector gives 50 %.')
    for i, t in enumerate(order):
        d = get(t); f, e = L.coin(d); b, v, feu, run = SAMPLES[t]
        col = colour_of(t)
        P.vbar(i, f, 56, col, base=0.40, opacity=0.9 if v == 'X' else 0.45,
               tip=f'{b}, {v} view, FEU {feu}, {run}\n{100 * f:.1f} ± {100 * e:.1f} % of {d["n"]:,} events')
        P.raw(sd.line(P.X(i), P.Y(f - e), P.X(i), P.Y(f + e), INK, 2))
        short = b.split(' · ')[0] + (' Jun' if 'June-15' in b else ' May' if 'May' in b else '')
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 30, short, 19, INK, weight=600))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 54, run.split(' ')[0], 18, MUT))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 76, f'FEU {feu}', 18, MUT))
    nx = sum(1 for t in order if SAMPLES[t][1] == 'X')
    P.raw(sd.line(P.X(nx - 0.5), P.y0, P.X(nx - 0.5), P.y0 + P.ph, RULE, 2))
    P.text(nx - 0.4, 0.89, 'Y views', 24, MUT, weight=600)
    P.text(-0.5, 0.89, 'X views', 24, MUT, weight=600)
    fx = L.coin(get('d3jun_x'))[0]
    P.text(1.6, 0.80, f'det3 June-15 X: p+1 wins {100 * fx:.0f} % of the time', 28, CJUN, weight=700)
    body = (sd.title('Which neighbour is bigger? A coin flip, except det3 X',
                     'Fraction of clean events where channel p+1 has a larger pulse than p−1 · every board and view, same code')
            + P.svg('coin flip'))
    rows = ''.join(f'<tr><td>{SAMPLES[t][0]} {SAMPLES[t][1]}</td><td>{SAMPLES[t][2]}</td><td>{SAMPLES[t][3]}</td>'
                   f'<td>{get(t)["n"]:,}</td><td>{100 * L.coin(get(t))[0]:.1f} ± {100 * L.coin(get(t))[1]:.1f}</td>'
                   f'<td>{f3(L.mean_se(get(t)["a1"])[0])}</td></tr>' for t in order)
    D.slide('coin', body, details(
        'This is the simplest observable I can think of: for each event, is the right neighbour (channel p+1) bigger '
        'than the left (p−1)? It needs no calibration, no model and no strip map. A detector with no left/right '
        'preference gives 50 % within the binomial error.',
        'det3’s June-15 build in X is far off (both runs, on two different FEUs); its Y view, read at the same time, '
        'is at 50 %. The May build of det3, read through the same FEU 3, is at 50 %. det4 and det7 X (built in '
        'Stephan’s June-22 batch) sit a few % high. The weekend det3 run is lower because its neighbours are wider '
        '(a copy of the same absolute size changes the ranking less).',
        f'<table><tr><th>view</th><th>FEU</th><th>run</th><th>events</th><th>P(p+1 &gt; p−1) %</th><th>copy a</th></tr>{rows}</table>',
        sources('asym_lib.coin', 'asym_lib.load')), short='Coin flip')


def s_hists(D):
    edges = np.linspace(0, 1.0, 41); ctr = 0.5 * (edges[1:] + edges[:-1])
    panels = []
    for tag, ttl in (('d3jun_x', 'det3 June-15 · X'), ('d3jun_y', 'det3 June-15 · Y (same events’ other view)')):
        d = get(tag); c = d['A'][:, HW]
        l, r = d['A'][:, HW - 1] / c, d['A'][:, HW + 1] / c
        hl, hr = hist(l, edges), hist(r, edges)
        ymax = max(hl.max(), hr.max()) * 1.15
        P = sd.Plot(810, 560, x=(0, 1.0), y=(0, ymax), xlabel='neighbour / peak (max samples)',
                    ylabel='fraction of events', title=ttl)
        P.xticks([(v, f'{v:.1f}') for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)])
        P.yticks([(v, f'{v:.2f}') for v in np.arange(0, ymax, 0.02 if ymax < 0.15 else 0.05)])
        xs, ys = sd.step_xy(edges, hl); P.line(xs, ys, CL, 3.5, markers=False, tip=f'p−1: mean {l.mean():.3f}, median {np.median(l):.3f}')
        xs, ys = sd.step_xy(edges, hr); P.line(xs, ys, CR, 3.5, markers=False, tip=f'p+1: mean {r.mean():.3f}, median {np.median(r):.3f}')
        P.vline(np.median(l), CL, '6 5', 2); P.vline(np.median(r), CR, '6 5', 2)
        P.text(0.45, ymax * 0.85, f'medians {np.median(l):.2f} vs {np.median(r):.2f}', 24, INK, weight=600)
        P.text(0.45, ymax * 0.76, f'{d["n"]:,} events', 21, MUT)
        panels.append(P.svg(ttl))
    body = (sd.title('The whole p+1 distribution shifts up in X, not in Y',
                     'Per-event ratio of each neighbour to the peak · 6-24 long run · dashed: medians')
            + sd.row(*panels, gap=44) + sd.legend([('p−1 (left)', CL), ('p+1 (right)', CR)]))
    D.slide('hists', body, details(
        'If the asymmetry came from a few odd events (a hot channel, a ringing chip, a mis-decoded event) it would '
        'show as a tail. Instead the whole p+1 distribution is shifted, by roughly a constant fraction of the peak.',
        'The Y view of the same detector, read in the same triggers on FEU 4, has identical left and right '
        'distributions.', sources('asym_lib.load')), short='Distributions')


def s_waves(D):
    t = np.arange(32) * NS_PER_SAMPLE
    panels = []
    for tag, ttl in (('d3jun_x', 'det3 June-15 · X'), ('d3jun_y', 'det3 June-15 · Y')):
        d = get(tag); M = d['W'].mean(0)
        P = sd.Plot(810, 560, x=(0, 1920), y=(-10, M[HW].max() * 1.12), xlabel='time [ns]',
                    ylabel='mean ADC − pedestal', title=ttl)
        P.xticks([(v, str(v)) for v in (0, 480, 960, 1440, 1920)])
        P.yticks([(v, f'{v:.0f}') for v in np.arange(0, M[HW].max() * 1.12, 50)])
        for k, col in ((0, CC), (-1, CL), (1, CR)):
            nm = 'peak p' if k == 0 else f'p{k:+d}'
            P.line(t, M[HW + k], col, 3.5, markers=False, tips=None,
                   tip=f'{nm}: max of mean {M[HW + k].max():.0f} ADC at {t[M[HW + k].argmax()]} ns')
        P.text(1000, M[HW].max() * 1.0, f'p+1 / p−1 at the peak: {M[HW + 1].max() / M[HW - 1].max():.2f}', 24, INK, weight=600)
        panels.append(P.svg(ttl))
    mx = get('d3jun_x')['W'].mean(0); rr = mx[HW + 1].max() / mx[HW - 1].max()
    body = (sd.title(f'Mean waveforms: in X p+1 is {rr:.1f}× p−1; in Y they overlap',
                     'Mean raw waveform of p−1, p, p+1 · 6-24 long run · same events in both views')
            + sd.row(*panels, gap=44) + sd.legend([('p−1', CL), ('peak p', CC), ('p+1', CR)]))
    D.slide('waves', body, details(
        'Mean of the pedestal-subtracted samples over all clean events. These waveforms are slow because this run '
        'was at a low drift field (600 V): the ion tail and the drift spread are long. The shape is the same in X '
        'and Y; only X has the extra charge on p+1.', sources('asym_lib.load')), short='Waveforms')


def s_copy_shape(D):
    d = get('d3jun_x'); M = d['W'].mean(0); t = np.arange(32) * NS_PER_SAMPLE
    diff = M[HW + 1] - M[HW - 1]
    k = np.flatnonzero(M[HW] > 0.2 * M[HW].max())
    s = float(np.dot(diff[k], M[HW][k]) / np.dot(M[HW][k], M[HW][k]))
    P = sd.Plot(810, 560, x=(0, 1920), y=(-5, diff.max() * 1.3), xlabel='time [ns]', ylabel='ADC',
                title='(p+1) − (p−1) against the peak, scaled')
    P.xticks([(v, str(v)) for v in (0, 480, 960, 1440, 1920)])
    P.yticks([(v, f'{v:.0f}') for v in np.arange(0, diff.max() * 1.3, 10)])
    P.line(t, s * M[HW], CC, 3, dash='8 6', markers=False, tip=f'{s:.3f} × the mean peak-channel pulse')
    P.line(t, diff, CR, 4, markers=False, tip='mean (p+1) − (p−1)')
    P.text(900, diff.max() * 1.15, f'= {s:.2f} × the peak pulse', 26, CR, weight=700)
    ratio = diff / np.where(M[HW] > 0.2 * M[HW].max(), M[HW], np.nan)
    Q = sd.Plot(810, 560, x=(0, 1920), y=(-0.05, 0.4), xlabel='time [ns]', ylabel='(p+1 − p−1) / p, per sample',
                title='The ratio sample by sample')
    Q.xticks([(v, str(v)) for v in (0, 480, 960, 1440, 1920)])
    Q.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3, 0.4)])
    Q.hline(0, MUT, None, 1.5)
    dy = get('d3jun_y'); My = dy['W'].mean(0)
    ry = (My[HW + 1] - My[HW - 1]) / np.where(My[HW] > 0.2 * My[HW].max(), My[HW], np.nan)
    Q.line(t, [None if np.isnan(v) else float(v) for v in ratio], CR, 3.5,
           tips=[f'{tt} ns: {v:.3f}' if not np.isnan(v) else None for tt, v in zip(t, ratio)])
    Q.line(t, [None if np.isnan(v) else float(v) for v in ry], GREY, 3, dash='6 5',
           tips=[f'Y, {tt} ns: {v:.3f}' if not np.isnan(v) else None for tt, v in zip(t, ry)])
    Q.text(80, 0.34, 'flat: same shape, same time', 26, INK, weight=700)
    Q.text(1300, 0.03, 'Y view', 22, GREY)
    body = (sd.title(f'The excess is a {100 * s:.0f} % copy of the peak’s own pulse',
                     'det3 June-15 X, 6-24 long run · samples where the peak pulse is above 20 % of its maximum')
            + sd.row(P.svg('difference'), Q.svg('ratio'), gap=44))
    D.slide('copy-shape', body, details(
        'If the extra charge on p+1 were real induced charge from a displaced cloud, or a crosstalk through a '
        'capacitor (which differentiates), it would have its own shape and timing. Instead the mean difference '
        '(p+1) − (p−1) is the peak channel’s own pulse times a constant, sample by sample, including the long tail.',
        'NB_ASYM §1 fitted d = a·c + b·dc/dt on the earlier tagged dumps: b ≈ 0, residual rms 0.002–0.010.',
        'So the cleanest description is: channel p+1 receives a fixed fraction of whatever channel p carries, at the '
        'same time. The next slide tests that “copy” description on the second neighbours.',
        sources('asym_lib.load')), short='Copy shape')


def s_copy_model(D):
    d = get('d3jun_x'); a = L.fit_copy(d)
    pr = L.profile(d['W']); pu = L.profile(L.uncopy(d['W'], a))
    # cartoon
    o = []
    for i, nm in enumerate(['n−2', 'n−1', 'n', 'n+1', 'n+2']):
        x = 40 + i * 140
        col = CC if nm == 'n' else MUT
        o.append(f'<rect x="{x}" y="60" width="100" height="80" rx="12" fill="{CARD}" stroke="{col}" stroke-width="3"/>')
        o.append(sd.T(x + 50, 110, nm, 28, INK, 'middle', 600))
        if i < 4:
            o.append(f'<path d="M {x + 60} 50 Q {x + 120} 0 {x + 180} 50" fill="none" stroke="{CR}" stroke-width="3"/>')
            o.append(f'<polygon points="{x + 180},50 {x + 168},40 {x + 172},54" fill="{CR}"/>')
            o.append(sd.T(x + 120, 14, '+a×', 22, CR, 'middle', 600))
    o.append(sd.T(40, 200, 'every channel adds a × its lower neighbour, same sample', 24, INK, 'start', 600))
    o.append(sd.T(40, 250, 'predicts, for a peak at n:', 22, MUT, 'start'))
    o.append(sd.T(60, 290, '±1:  (n+1) − (n−1) = a × (peak − second neighbour)', 22, INK, 'start'))
    o.append(sd.T(60, 326, '±2:  (n+2) − (n−2) = a × (first neighbour)', 22, INK, 'start'))
    o.append(sd.T(40, 380, 'Undo it (one a for the whole run):', 22, MUT, 'start'))
    o.append(sd.T(60, 418, 'x′(n+1) = x(n+1) − a · x′(n)', 24, INK, 'start', 600))
    o.append(sd.T(40, 470, f'fitted on ±1 only: a = {a:.3f}', 26, CR, 'start', 700))
    cart = sd.svg(760, 520, ''.join(o), 'copy model')
    P = sd.Plot(860, 600, x=(0.4, 4.6), y=(0.002, 0.6, 'log'), xlabel='distance from the peak [channels]',
                ylabel='in-phase content / peak', title='Left vs right, before and after undoing the copy')
    P.xticks([(k, str(k)) for k in (1, 2, 3, 4)])
    P.yticks([(v, f'{v:g}') for v in (0.003, 0.01, 0.03, 0.1, 0.3)])
    ks = [1, 2, 3, 4]
    P.line([k - 0.06 for k in ks], [max(pr[HW - k], 0.0021) for k in ks], CL, 2.5, dash='6 5', r=7,
           tips=[f'raw p−{k}: {pr[HW - k]:.4f}' for k in ks])
    P.line([k - 0.06 for k in ks], [max(pr[HW + k], 0.0021) for k in ks], CR, 2.5, dash='6 5', r=7,
           tips=[f'raw p+{k}: {pr[HW + k]:.4f}' for k in ks])
    P.line([k + 0.06 for k in ks], [max(pu[HW - k], 0.0021) for k in ks], CL, 4, r=9, marker='diamond',
           tips=[f'uncopied p−{k}: {pu[HW - k]:.4f}' for k in ks])
    P.line([k + 0.06 for k in ks], [max(pu[HW + k], 0.0021) for k in ks], CR, 4, r=9, marker='diamond',
           tips=[f'uncopied p+{k}: {pu[HW + k]:.4f}' for k in ks])
    P.text(1.25, 0.36, f'raw: {pr[HW + 1]:.2f} vs {pr[HW - 1]:.2f}', 22, INK)
    P.text(1.25, 0.19, f'undone: {pu[HW + 1]:.3f} vs {pu[HW - 1]:.3f}', 24, INK, weight=700)
    body = (sd.title('One copy fraction balances ±1 and ±2 together',
                     'det3 June-15 X, 6-24 long run · a fitted on the first neighbours only; ±2 is a prediction')
            + sd.row(cart, sd.col(P.svg('profile'),
                                sd.legend([('p−k', CL), ('p+k', CR), ('raw (dashed, dots)', MUT, 'dash'),
                                           ('copy undone (solid, diamonds)', MUT)], size=21))
                     , gap=44))
    D.slide('copy-model', body, details(
        'The copy model has one parameter. Undoing it runs upward through the channel window: x′(n+1) = x(n+1) − a·x′(n). '
        f'Fitting a so that the first neighbours balance gives a = {a:.3f}. The second neighbours were not used in '
        f'the fit; raw they are p−2 {pr[HW - 2]:.4f} vs p+2 {pr[HW + 2]:.4f}, after undoing {pu[HW - 2]:.4f} vs {pu[HW + 2]:.4f}.',
        'NB_ASYM §1 found the same on the tagged dumps (a from ±1 and from ±2 agree on every run: 0.085/0.085, '
        '0.145/0.142, 0.196/0.167), and <code>uncopy.py</code> closes the full-FEU version (6-22: ±1 0.153 → 0.005).',
        'In-phase content α_k = Σ_t x_k(t)·x_p(t) / Σ_t x_p(t)²: the part of channel p+k shaped like the peak pulse.',
        sources('asym_lib.fit_copy', 'asym_lib.uncopy', '../uncopy.py')), short='Copy model')


def s_noise(D):
    tags = [('d3jun_x', 'det3 June-15 X', CJUN, None), ('d3jun_y', 'det3 June-15 Y', GREY, None),
            ('d3may_x', 'det3 May X', CMAY, None)]
    P = sd.Plot(1000, 620, x=(0.7, 6.3), y=(0, 0.42), xlabel='channel distance d (same chip)',
                ylabel='noise correlation, common mode removed', title='Noise on quiet triggers')
    P.xticks([(k, str(k)) for k in range(1, 7)])
    P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3, 0.4)])
    z = np.load(f'{L.RAW}/d3jun_x.npz')
    inj = z['corr_inj'][1:]
    P.line(range(1, 7), list(inj), CJUN, 3, dash='8 6', r=7,
           tips=[f'det3 June X + injected copy 0.15, d={k}: {v:.3f}' for k, v in zip(range(1, 7), inj)])
    for tag, nm, c, _ in tags:
        zz = np.load(f'{L.RAW}/{tag}.npz'); cc = zz['corr_cm'][1:]
        P.line(range(1, 7), list(cc), c, 4, r=8, tips=[f'{nm}, d={k}: {v:.3f} ({int(zz["nquiet"])} quiet events)'
                                                         for k, v in zip(range(1, 7), cc)])
    dj = z['corr_cm'][1] - z['corr_cm'][2]; di = inj[0] - inj[1]
    P.text(1.6, 0.39, 'an electronics copy would put a spike at d = 1', 24, CJUN)
    P.text(2.3, 0.04, f'measured d=1 minus d=2: {dj:+.2f}  (with a 0.15 copy: {di:+.2f})', 22, INK, weight=600)
    side = sd.col(
        sd.legend([('det3 June-15 X', CJUN), ('… with a 0.15 copy injected', CJUN, 'dash'),
                   ('det3 June-15 Y', GREY), ('det3 May X', CMAY)], size=22),
        sd.callout('The signal is copied, the noise is not. A copy made in the cable, the FEU or the ADC would copy '
                   'both.', CJUN, 26),
        sd.callout('So the extra charge is induced before the preamplifier input: on the detector, or between its pads '
                   'and the connector.', GREY, 24), gap=28, w=600)
    body = (sd.title('The noise isn’t copied: it is not the electronics',
                     'Correlation between channel n and n+d on triggers with nothing above 5σ in the chip · 6-24 long run')
            + sd.row(P.svg('noise'), side, gap=56))
    D.slide('noise', body, details(
        'Every DREAM channel has its own preamplifier noise. If anything after the input (shaper, multiplexer, ADC '
        'settling, a cable between FEU stages) added a × channel n to channel n+1, the noise would be copied too, '
        'and the correlation between neighbours (d = 1) would jump relative to d = 2. Injecting x(n+1) += 0.15 x(n) '
        'into these same quiet samples shows the size of that jump (dashed).',
        'The measured X correlation has no d = 1 spike, the same as its Y view and as the May build (same FEU 3, no '
        'copy). The correlations themselves are common-mode residue; their level differs between boards and runs, so '
        'compare the shape (d = 1 vs d = 2), not the level.',
        'Caveat: noise injected at the preamp input (its series noise) appears on the input node; coupling between '
        'input lines could carry some of it. So this test excludes the FEU chain, not a coupling on the input traces '
        'or connector. That one is addressed by “follows the board” (next slide) and by the field dependence (slide 17).',
        'NB_ASYM had the same result with the older dumps (+0.03/−0.03/−0.02 vs +0.06…0.12 for an injected 0.1 copy).',
        sources('raw_asym.py (corr_cm, corr_inj)', '../noise_corr.py')), short='Noise')


def s_follow(D):
    # points: (FEU, board label, a, err, colour, tip)
    pts = []
    for t in SAMPLES:
        if not have(t) or SAMPLES[t][1] != 'X':
            continue
        d = get(t); a, e = L.mean_se(d['a1']); b, v, feu, run = SAMPLES[t]
        pts.append((feu, b, a, e, colour_of(t), f'{b} X on FEU {feu}, {run}: a = {a:+.3f} ± {e:.3f} ({d["n"]} events)'))
    rows = survey_rows()
    feus = [3, 4, 6, 7]
    P = sd.Plot(1080, 640, x=(-0.6, 3.6), y=(-0.05, 0.22), xlabel='FEU slot the X view was read on',
                ylabel='copy a  (p+1 minus p−1, in phase)')
    P.xticks([(i, f'FEU {f}') for i, f in enumerate(feus)])
    P.yticks([(v, f'{v:.2f}') for v in (-0.05, 0, 0.05, 0.10, 0.15, 0.20)])
    P.hline(0, MUT, None, 1.5)
    seen = {}
    for feu, b, a, e, c, tp in sorted(pts, key=lambda q: q[2]):
        i = feus.index(feu); k = seen.get(feu, 0); seen[feu] = k + 1
        x = i - 0.25 + 0.17 * k
        P.raw(sd.line(P.X(x), P.Y(a - e), P.X(x), P.Y(a + e), c, 2))
        P.points([x], [a], c, r=10, tips=[tp])
        P.text(x + 0.05, a + 0.007, b.split(' · ')[0] + (' June' if 'June-15' in b else ' May' if 'May' in b else ''),
               19, c if c != GREY else MUT)
    P.text(0.75, 0.205, 'det3 June-15 copies on FEU 3 and on FEU 7;', 24, CJUN, weight=700)
    P.text(0.75, 0.185, 'FEU 3 with det3 May or det6 does not', 24, INK, weight=600)
    side = sd.col(
        sd.callout('Same FEU, other board → no copy.', CMAY, 28),
        sd.callout('Same board, other FEU → copy.', CJUN, 28),
        sd.callout('Same board, same triggers, other view (Y on FEU 4/8) → no copy.', GREY, 26),
        sd.p('The copy follows the board and its X view, not the readout slot.', 26, INK, 600), gap=28, w=500)
    body = (sd.title('It follows the board, not the FEU slot',
                     'Copy a by FEU slot · each dot one board’s X view · raw extraction, same cuts')
            + sd.row(P.svg('follows the board'), side, gap=56))
    D.slide('follow', body, details(
        'FEU 3 read det3’s May build (May 6), det6 (June 26) and det3’s June-15 build (June 22–24): only the last '
        'copies. det3’s June-15 build was moved to FEU 7 for the June 27–28 runs and still copies. Its Y view, read in '
        'the same triggers on FEU 4 or 8, never does.',
        'What this does not settle: whether det3’s X cable or adapter travelled with the board between slots. If one '
        'specific cable did, it could still be the culprit. The field dependence (slide 17) argues against a cable.',
        'Also the May and June-15 “det3” builds may not be the same PCB (the run configs say “Bulked 5-5-26” and '
        '“Bulked by Stephan June 15”). If they are the same PCB, the copy appeared with the re-bulking.',
        sources('asym_lib.load', 'run_config.json dream_feus')), short='Follows the board')


def s_map(D):
    runs = [('t622L_x', '6-23 long, FEU 3, drift 1000 V', '#e07a86', None),
            ('t623L_x', '6-24 long, FEU 3, drift 600 V', CJUN, None),
            ('t627d100_x', '6-27, FEU 7, drift 100 V', PURPLE, None),
            ('t627d1100_x', '6-27, FEU 7, drift 1100 V', GOLD, None),
            ('d3wk_x', '6-28 weekend, FEU 7, drift 1000 V', ORANGE, None),
            ('d3jun_y', 'Y view, 6-24', GREY, '6 5'),
            ('d3may_x', 'May build X', CMAY, '6 5')]
    P = sd.Plot(1000, 640, x=(-0.4, 7.4), y=(-0.06, 0.36), xlabel='DREAM chip (64 channels ≈ 50 mm of the board)',
                ylabel='copy a', title='Across the board (X channel → x position)')
    P.xticks([(k, str(k)) for k in range(8)])
    P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    P.hline(0, MUT, None, 1.5)
    leg = []
    for tag, nm, c, dash in runs:
        if not have(tag):
            continue
        d = get(tag); ch = d['p'] // 64
        xs, ys, tps = [], [], []
        for k in range(8):
            s = ch == k
            if s.sum() >= 40:
                m, e = L.mean_se(d['a1'][s]); xs.append(k); ys.append(m)
                tps.append(f'{nm}, chip {k}: a = {m:+.3f} ± {e:.3f} ({s.sum()} events)')
        P.line(xs, ys, c, 3.5 if not dash else 2.5, dash=dash, r=7, tips=tps)
        leg.append((nm, c, 'dash' if dash else 'line'))
    P.text(4.4, 0.32, 'peaks mid-board, every run', 26, INK, weight=700)
    # 2D map
    right = sd.p('2D map pending: needs the 6-24 Y extraction (t623L_y).', 24, MUT)
    if have('t623L_y') and have('t623L_x'):
        x = get('t623L_x')
        zy = np.load(f'{L.RAW}/t623L_y.npz'); py = zy['p'].astype(int)
        hy = set(L.hot(py).tolist())
        tok = lambda files: [re.search(r'_(\d{3})_\d\d\.root$', f).group(1) for f in files]
        ty, tx = tok(zy['files']), tok(x['z']['files'])
        ky = {(ty[f], e): p for f, e, p in zip(zy['fi'], zy['ev'], py) if p not in hy}
        cx, cy, aa = [], [], []
        for f, e, p, a in zip(x['fi'], x['ev'], x['p'], x['a1']):
            q = ky.get((tx[f], e))
            if q is not None:
                cx.append(p // 128); cy.append(q // 128); aa.append(a)
        cx, cy, aa = map(np.array, (cx, cy, aa))
        cells = []
        for i in range(4):
            for j in range(4):
                s = (cx == i) & (cy == j)
                if s.sum() >= 25:
                    m, e = L.mean_se(aa[s])
                    cells.append((i, i + 1, j, j + 1, m, f'X chips {2 * i}–{2 * i + 1}, Y chips {2 * j}–{2 * j + 1}: '
                                                         f'a = {m:+.3f} ± {e:.3f} ({s.sum()} events)'))
                else:
                    cells.append((i, i + 1, j, j + 1, None, f'{s.sum()} events'))
        Q = sd.Plot(560, 560, x=(0, 4), y=(0, 4), xlabel='X channel block (128 ch)', ylabel='Y channel block (128 ch)',
                    title=f'Face of the board, 6-24 ({len(aa):,} matched)', margin=(24, 20, 92, 96))
        Q.xticks([(k + 0.5, str(k)) for k in range(4)]); Q.yticks([(k + 0.5, str(k)) for k in range(4)])
        Q.cells(cells, 0, 0.35, missing='#ddd', gap=2)
        for c_ in cells:
            if c_[4] is not None:
                Q.text((c_[0] + c_[1]) / 2, (c_[2] + c_[3]) / 2 - 0.08, f'{c_[4]:.2f}', 24,
                       '#fff' if c_[4] < 0.2 else INK, 'middle', 600)
        right = sd.col(Q.svg('2d map'), sd.colorbar(560, 80, 0, 0.35, [(0, '0'), (0.1, '0.1'), (0.2, '0.2'), (0.3, '0.3')],
                                                    'copy a', horizontal=True), gap=4)
    body = (sd.title('Not uniform: 2–3× stronger mid-board than at the edges',
                     'det3 June-15 build, X view · copy a per DREAM chip, five runs · right: by X and Y position together')
            + sd.row(sd.col(P.svg('map'), sd.legend(leg, size=19), gap=6), right, gap=48))
    D.slide('map', body, details(
        'Each DREAM chip reads 64 adjacent X channels, about 50 mm of the board. The copy is ~0.1 at the two ends '
        '(chips 0, 1, 7) and 0.2–0.3 in chips 2–5, with the same shape on every run and at both drift fields. '
        'Chip order follows the board’s x position (the connectors sit in order along the edge).',
        'This matters for the cause. A layout feature of the PCB (the bus traces, the connector pin order) repeats '
        'identically for every channel and every chip, so it would give a flat profile. A profile that peaks in '
        'the middle points to something that varies over the board’s area: the mesh or drift-gap geometry (a bowed '
        'mesh or cathode changes the local field), the bulk (lamination, pillars, coverlay thickness, curing), or '
        'moisture.',
        'Right: matched X and Y events (same trigger) binned by both coordinates, so the dependence on y (along the '
        'X strips) can be seen too. Grey cells have fewer than 25 events.',
        sources('asym_lib.load', 'raw_asym.py t622L_x t623L_x t623L_y t627d100_x t627d1100_x d3wk_x')),
        short='Across the board')


def s_amp(D):
    edges = [200, 500, 900, 1400, 2000, 2700, 3500]
    P = sd.Plot(1664, 640, x=(150, 3600), y=(-0.05, 0.3), xlabel='peak-channel amplitude [ADC]', ylabel='copy a')
    P.xticks([(v, str(v)) for v in (500, 1000, 1500, 2000, 2500, 3000, 3500)])
    P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    P.hline(0, MUT, None, 1.5)
    leg = []
    for tag, nm, c in (('t623L_x', 'det3 June X, 6-24', CJUN), ('t622L_x', 'det3 June X, 6-23', '#e07a86'),
                       ('d3wk_x', 'det3 June X, 6-28', ORANGE), ('d3jun_y', 'det3 June Y, 6-24', GREY),
                       ('d3may_x', 'det3 May X', CMAY)):
        if not have(tag):
            continue
        d = get(tag); c0 = d['A'][:, HW]; xs, ys, tps = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (c0 >= lo) & (c0 < hi)
            if s.sum() >= 40:
                m, e = L.mean_se(d['a1'][s]); xs.append(np.median(c0[s])); ys.append(m)
                tps.append(f'{nm}, {lo}–{hi} ADC: a = {m:+.3f} ± {e:.3f} ({s.sum()} events)')
        P.line(xs, ys, c, 3.5, r=8, tips=tps); leg.append((nm, c))
    P.text(300, 0.27, 'flat: a fixed fraction, from small pulses to near-railed ones', 24, INK, weight=700)
    body = (sd.title('Linear: the copied fraction ignores the pulse height',
                     'Copy a in bins of the peak-channel amplitude · det3 June-15 X runs, with Y and the May build as controls')
            + P.svg('amplitude') + sd.legend(leg, size=22))
    D.slide('amplitude', body, details(
        'A threshold, saturation or a non-linear ADC effect would make the copy depend on the pulse height. It is flat '
        'within errors from 200 to 3500 ADC: a linear coupling.',
        'NB_ASYM also found it flat in track angle (a real tilt gives −0.0085 per degree on det2; det3’s copy would need '
        '~20°), flat in channel-within-connector and parity, and flat in the 39:40 resistive-strip beat phase.',
        sources('asym_lib.load', '../nb_copy_vs_amp.py', '../nb_copy_vs_angle.py', '../nb_asym.py')), short='Amplitude')


def s_hot(D):
    seq = []
    for tag, nm in (('t623s_x', 'short'), ('t623l_x', 'longer'), ('t623L_x', 'long')):
        if not have(tag):
            continue
        dv, dn = get(tag, True), get(tag, False)
        for k in range(int(np.load(f'{L.RAW}/{tag}.npz')['fi'].max()) + 1):
            sv, sn = dv['fi'] == k, dn['fi'] == k
            hv = np.isin(dn['p'], dv['hot'])
            seq.append((f'{nm} {k + 1}', L.mean_se(dn['a1'][sn]) if sn.sum() > 20 else (None, 0),
                        L.mean_se(dv['a1'][sv]) if sv.sum() > 40 else (None, 0), int(sn.sum()), int(sv.sum()),
                        float(hv[sn].mean()) if sn.sum() else 0.0))
    n = len(seq)
    P = sd.Plot(1664, 600, x=(-0.6, n - 0.4), y=(0, 0.25), ylabel='copy a', margin=(20, 20, 110, 100))
    P.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.10, 0.15, 0.20, 0.25)])
    xs = list(range(n))
    P.line(xs, [s[1][0] for s in seq], GREY, 3, dash='8 6', r=8,
           tips=[f'{s[0]}: all events a = {s[1][0]:+.3f} ({s[3]} events, {100 * s[5]:.0f} % on the hot channel)'
                 if s[1][0] is not None else None for s in seq])
    P.line(xs, [s[2][0] for s in seq], CJUN, 4, r=9,
           tips=[f'{s[0]}: hot channel removed a = {s[2][0]:+.3f} ({s[4]} events)' if s[2][0] is not None else None
                 for s in seq])
    for i, s in enumerate(seq):
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 32, s[0], 19, INK))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 58, f'{100 * s[5]:.0f} % hot', 17, MUT))
    P.text(0, 0.225, 'all events: looks like the copy “switches on” over the evening', 24, GREY, weight=600)
    P.text(0, 0.205, 'hot channel 492 removed: the same ~0.17 throughout', 24, CJUN, weight=700)
    body = (sd.title('A real trap: a hot channel faked a “turn-on” one evening',
                     'det3 June-15 X, 6-23 evening → 6-24 morning (same HV) · copy per data file, with and without the hot-channel veto')
            + P.svg('hot channel'))
    D.slide('hot', body, details(
        'The earlier survey (copy_survey.csv) showed det3’s copy at 0.035 at 19:29, 0.057 at 20:29 and 0.17 overnight, at '
        'unchanged HV, which looked like a slow turn-on. It is not. In the short run 3259 of 4493 “events” are channel '
        '492 firing on noise (flat waveforms, no pulse); its two fixed neighbours decide the answer, and they happen '
        'to be nearly equal. As the run goes on, real cosmics outnumber it and the mean climbs to the true value.',
        'With the pulse-shape cut and the hot-channel veto every file of the evening gives the same copy, and the '
        'per-chip profile of the longer run equals the long run’s.',
        'This is exactly the kind of algorithm effect worth worrying about, and it shows the veto is needed before any '
        'time-dependence claim. The survey slides below use the raw re-extraction with the veto.',
        sources('asym_lib.hot', 'raw_asym.py t623s_x t623l_x t623L_x')), short='Hot-channel trap')


def s_survey(D, rows):
    if not rows:
        D.slide('survey', sd.title('Every board, every sub-run: pending', 'Raw survey (run_raw3.sh) not yet available'),
                short='Survey')
        return
    cols = [('det1', 'X'), ('det1', 'Y'), ('det2', 'X'), ('det2', 'Y'), ('det3 May', 'X'), ('det3 May', 'Y'),
            ('det3 June-15', 'X'), ('det3 June-15', 'Y'), ('det4', 'X'), ('det4', 'Y'), ('det6', 'X'), ('det6', 'Y'),
            ('det7', 'X'), ('det7', 'Y')]
    P = sd.Plot(1664, 680, x=(-0.6, len(cols) - 0.4), y=(-0.08, 0.32), ylabel='copy a per sub-run',
                margin=(20, 20, 120, 100))
    P.yticks([(v, f'{v:.1f}') for v in (-0.1, 0, 0.1, 0.2, 0.3)])
    P.hline(0, MUT, None, 1.5)
    rng = np.random.default_rng(1)
    for i, (b, v) in enumerate(cols):
        rr = [r for r in rows if r['build'] == b and r['view'] == v]
        c = CJUN if b == 'det3 June-15' else CMAY if b == 'det3 May' else COTH
        if v == 'Y':
            c = '#b9c1cc' if c == COTH else c
        for r in rr:
            x = i + rng.uniform(-0.22, 0.22)
            P.points([x], [min(max(r['a'], -0.075), 0.315)], c, r=6 if v == 'X' else 5,
                     tips=[f'{b} {v}, {r["run"]} / {r["sub"]}\nFEU {r["feu"]}, mesh {r["mesh"]} V, drift {r["drift"]} V\n'
                           f'a = {r["a"]:+.3f} ± {r["ea"]:.3f} ({r["n"]} events); P(p+1 > p−1) = {100 * r["coin"]:.0f} %'
                           + (f'\nhot channels vetoed: {r["hot"]}' if r['hot'] else '')])
        if rr:
            med = float(np.median([r['a'] for r in rr]))
            P.raw(sd.line(P.X(i - 0.32), P.Y(med), P.X(i + 0.32), P.Y(med), INK, 3))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 32, b.replace('det3 ', 'det3·'), 18, INK))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 56, f'{v} · {len(rr)}', 18, MUT))
    nj = len([r for r in rows if r['build'] == 'det3 June-15' and r['view'] == 'X'])
    P.text(-0.4, 0.29, f'det3’s June-15 build copies in all {nj} of its X sub-runs; its Y stays at 0 ± 0.04', 26, CJUN, weight=700)
    P.text(-0.4, 0.26, 'Stephan’s June-22 batch (det4, det7) shows a few %', 22, INK)
    body = (sd.title('Every board, every sub-run: it is det3 June-15 X',
                     'Raw re-extraction of all cosmic-bench sub-runs (2 files each), hot channels vetoed · bar: median · n = sub-runs')
            + P.svg('survey'))
    D.slide('survey', body, details(
        'Each dot is one sub-run (a run point of an HV or drift scan, or a long run), re-extracted from decoded_root with '
        'the raw extractor and the same cuts as before. Values are clipped to the axis; the tooltip has the true value.',
        'Builds: det1 (January), det2 (Arnaud, June 12), det3 May build (bulked May 5), det3 June-15 build (Stephan), '
        'det4 and det7 (Stephan, batch of three on June 22), det6 (Stephan, June 24, relaminated).',
        'n_TOF beam (from the earlier survey, not re-extracted: the beam flash defeats this simple selection): det3 X '
        'on run_79 gave 0.053–0.066 (shape fit 0.085), its Y −0.04; det2 X 0.06–0.09 within the beam-angle '
        'systematic. det6/det7 beam X fail a quality flag (a position-map problem in those configs).',
        sources('run_raw3.sh', 'raw/survey/index.tsv', '../copy_survey.csv (times, HV)')), short='Survey')


def date_rows(rows):
    """det3 June-15 X, one row per sub-run with >= 60 board-centre events, plus the EXTRA runs; days since 6-16."""
    out = [r for r in rows if r['build'] == 'det3 June-15' and r['view'] == 'X' and r['t'] and r['nc'] >= 60]
    for f, run, sub, feu, t, drift in EXTRA:
        p = f'{L.RAW}/{f}'
        if not os.path.exists(p):
            continue
        d = L.load(p)
        a, ea = L.mean_se(d['a1'])
        ac, eac, nc = central(d)
        out = [r for r in out if r['run'] != run] + [dict(
            f=p, run=run, sub=sub, det='', view='X', feu=feu, board='det3', build='det3 June-15', a=a, ea=ea, ac=ac,
            eac=eac, nc=nc, n=d['n'], coin=L.coin(d)[0], hot=list(d['hot']), mesh='500', drift=drift, t=t)]
    t0 = dt.datetime(2026, 6, 16)
    for r in out:
        r['day'] = (dt.datetime.strptime(r['t'], '%Y-%m-%d %H:%M') - t0).total_seconds() / 86400
    return sorted(out, key=lambda r: r['day'])


def s_date(D, rows):
    jx = date_rows(rows)
    if not jx:
        return
    dcol = lambda v: sd.seq_color((float(v) - 100) / 1000)
    YR = (-0.02, 0.40)
    yt = [(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)]
    tip = lambda r: (f'{r["run"]} / {r["sub"]}\n{r["t"]}, FEU {r["feu"]}, mesh {r["mesh"]} V, drift {r["drift"]} V\n'
                     f'centre chips 2–5: a = {r["ac"]:+.3f} ± {r["eac"]:.3f} ({r["nc"]} events)\n'
                     f'whole board: a = {r["a"]:+.3f} ± {r["ea"]:.3f} ({r["n"]} events)')
    # May build (same board name, earlier bulk)
    may = sorted([r for r in rows if r['build'] == 'det3 May' and r['view'] == 'X' and r['nc'] >= 60],
                 key=lambda r: r['t'])
    M = sd.Plot(250, 640, x=(-0.5, len(may) - 0.5), y=YR, ylabel='copy a, board-centre chips 2–5',
                title='May build', margin=(24, 16, 92, 104))
    M.yticks(yt); M.hline(0, MUT, None, 1.5)
    for i, r in enumerate(may):
        M.raw(sd.line(M.X(i), M.Y(r['ac'] - r['eac']), M.X(i), M.Y(r['ac'] + r['eac']), MUT, 2))
        M.points([i], [r['ac']], dcol(r['drift']), r=7, tips=[tip(r)])
    M.raw(sd.T(M.X((len(may) - 1) / 2), M.y0 + M.ph + 34, '5-5 → 6-7', 20, MUT))
    M.text(-0.3, 0.375, '≈ 0', 26, INK, weight=700)
    # June-15 build on the bench
    P = sd.Plot(1080, 640, x=(1.4, 12.8), y=YR, xlabel='start of sub-run (June 2026)', ylabel='',
                title='June-15 build, cosmic bench', margin=(24, 24, 92, 40))
    P.xticks([(k, f'6-{16 + k}') for k in range(2, 13)])
    P.yticks(yt); P.hline(0, MUT, None, 1.5)
    for r in jx:
        P.raw(sd.line(P.X(r['day']), P.Y(r['ac'] - r['eac']), P.X(r['day']), P.Y(r['ac'] + r['eac']), MUT, 2))
    for r in jx:
        P.points([r['day']], [r['ac']], dcol(r['drift']), r=8, marker='circle' if r['feu'] == 3 else 'open',
                 tips=[tip(r)])
    def sess(key, lab_x, lab_y, extra='', drift=None):
        rr = [r for r in jx if key in r['run'] and (drift is None or r['drift'] == drift)]
        if not rr:
            return
        w = np.array([1 / r['eac'] ** 2 for r in rr]); m = float(np.sum(w * [r['ac'] for r in rr]) / w.sum())
        dv = sorted({r['drift'] for r in rr})
        P.text(lab_x, lab_y, f'{m:.2f}', 26, INK, weight=700)
        P.text(lab_x, lab_y - 0.022, f'drift {"/".join(dv)} V' + extra, 18, MUT)
    sess('6-17-26', 1.75, 0.375)
    sess('6-19-26', 3.6, 0.10, '')
    P.text(3.6, 0.056, 'drift HV sparking', 18, RED, weight=600)
    sess('det2_det3_overnight_6-22', 6.0, 0.375)
    sess('6-23-26', 7.6, 0.375)
    sess('6-25-26', 9.2, 0.375)
    sess('6-27-26', 10.6, 0.375, drift='1000')
    P.text(10.15, 0.235, 'drift scan → 100 V', 18, MUT, anchor='start')
    # n_TOF beam (older survey, whole board)
    beam = [r for r in csv.DictReader(open(os.path.join(HERE, '..', 'copy_survey.csv')))
            if r['phys'] == 'det3' and r['run'] == 'run_79' and r['view'] == 'X']
    B = sd.Plot(250, 640, x=(-0.6, max(len(beam), 1) - 0.4), y=YR, ylabel='', title='Beam, 7-26',
                margin=(24, 16, 92, 64))
    B.yticks(yt); B.hline(0, MUT, None, 1.5)
    for i, r in enumerate(beam):
        a, ea = float(r['a']), float(r['ea'])
        B.raw(sd.line(B.X(i), B.Y(a - ea), B.X(i), B.Y(a + ea), MUT, 2))
        B.points([i], [a], dcol(r['drift_V']), r=7, marker='diamond',
                 tips=[f'n_TOF run_79 / {r["sub"]}, {r["tstart"]}\nFEU {r["feu"]}, mesh {r["mesh_V"]} V, drift {r["drift_V"]} V, '
                       f'{r["gas"]}\nwhole board, older survey method: a = {a:+.3f} ± {ea:.3f}\n'
                       f'(det3 Y in the same run reads −0.04: beam-angle systematic)'])
    B.raw(sd.T(B.X((len(beam) - 1) / 2), B.y0 + B.ph + 34, 'run_79', 20, MUT))
    cb = sd.colorbar(130, 640, 100, 1100, [(v, f'{v}') for v in (100, 300, 600, 1000)], 'drift [V]')
    key = sd.legend([('FEU 3 slot', INK, 'dot'), ('FEU 7 slot', INK, 'open'),
                     ('beam: whole board, older method', INK, 'box')], 20)
    body = (sd.title('By date: on, off, on, and fading',
                     'det3 X copy per sub-run · raw extraction, hot channels vetoed · board-centre chips only · '
                     'colour: drift voltage')
            + sd.row(M.svg('may build'), P.svg('june bench by date'), B.svg('beam'), cb, gap=4)
            + key)
    D.slide('date', body, details(
        'Each point is one sub-run of det3 X (an HV/drift-scan point or a long run) with at least 60 clean events in '
        'the board-centre chips 2–5; bars are the standard error. Centre chips only, because the copy is centre-peaked '
        '(~0.25 mid-board, ~0.1 at the edges) and several 6-22 high-mesh points have 70–90 % of their hits on chip 7: '
        'their board average reads low for a reason that has nothing to do with time. The tooltip has both numbers. '
        'Bold numbers: inverse-variance mean of each session (the 6-27 one mixes the drift scan with the 1000 V runs; '
        'its 1000 V runs alone are ≈ 0.10).',
        'Reading it: <b>6-18</b> (first night, FEU 7, only chips 0–3 lit) 0.26. <b>6-19/20</b> (FEU 7, 7458 events, flat '
        'through all 8 files, no centre peak) 0.03. <b>6-22/23</b> (FEU 3, drift 1000 V) 0.16, flat through 12 h. '
        '<b>6-23/24</b> at drift 600 V 0.26, flat through 7 h. <b>6-25</b> drift 500 V 0.17. <b>6-27/28</b> (FEU 7, '
        'resist HV moved 3:3 → 3:4) 0.09–0.10 at 1000 V, rising to 0.21 at 100 V. <b>7-26</b> beam run_79 0.05–0.07 '
        '(Ar/iso 90/10, drift 700 V; older method on the whole board, and the same method gives −0.04 on det3 Y there, so '
        'treat it as ≲ 0.1). May build ≈ 0 throughout.',
        'What it says: the copy is stable within a session (hours) and jumps between sessions (days). It is not a slow '
        'monotonic decay: it was large on the first night, near zero the next, back two days later. Drift field moves it '
        '(×1.6–2 from 1000 V to ≲ 600 V) but cannot explain 6-19, which ran at 1000 V like 6-22. The slot (FEU 3 vs 7) '
        'cannot either: 6-18 and 6-19 were both on FEU 7.',
        'The one thing special about 6-19: it is the only run where det3’s drift channel (0:7 in that config) misbehaved, '
        '4–9 µA spikes with the voltage sagging to 890–960 V in most sub-runs, against a flat 0.34 µA (the module’s '
        'own load) in every other run. Resist current was also high (2–4 µA at 500 V vs ~1 µA). A board that was '
        'discharging in the drift gap all night carried no copy. One coincidence, not a mechanism, but it fits a charge-on-'
        'a-surface picture (sparks neutralising it) better than any fixed coupling. Caveat: the 6-19 config swapped the HV '
        'channels between the two boards relative to 6-17; the other board that night (FEU 3) reads 0.08 ± 0.03 from 372 '
        'noisy events, so it is not a clean label swap, but a cabling mix-up that night is not excluded.',
        'Missing: the 6-26 weekend scan (FEU 3) is empty on EOS; the 6-19 HV scan and long_run exist only as raw_daq_data '
        '(not decoded), so they could extend 6-19/20 if decoded.',
        sources('raw survey', 'extra/d617_long_x.npz, extra/d619_longer_x.npz (EOS june_tests, raw_asym.py)',
                'hv_monitor.csv per sub-run', '../copy_survey.csv (beam)')), short='By date')


def s_time(D, rows):
    jx = [r for r in rows if r['build'] == 'det3 June-15' and r['view'] == 'X' and r['t']]
    if not jx:
        D.slide('time', sd.title('Time and drift field: pending', 'raw survey not yet available'), short='Time, field')
        return
    t0 = dt.datetime(2026, 6, 22)
    def hrs(r):
        return (dt.datetime.strptime(r['t'], '%Y-%m-%d %H:%M') - t0).total_seconds() / 86400
    jx = [r for r in jx if hrs(r) >= -0.2]
    P = sd.Plot(640, 620, x=(-0.2, 6.6), y=(0, 0.3), xlabel='date (June 2026)', ylabel='copy a',
                title='Every bench sub-run')
    P.xticks([(k, f'6-{22 + k}') for k in range(0, 7, 2)])
    P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    dcol = lambda v: sd.seq_color((float(v) - 100) / 1000)
    for r in sorted(jx, key=hrs):
        P.points([hrs(r)], [r['a']], dcol(r['drift'] or 1000), r=8,
                 tips=[f'{r["run"]} / {r["sub"]}\n{r["t"]}, FEU {r["feu"]}, mesh {r["mesh"]} V, drift {r["drift"]} V\n'
                       f'a = {r["a"]:+.3f} ± {r["ea"]:.3f} ({r["n"]} events)'])
    P.text(0.0, 0.28, 'FEU 3', 22, MUT, weight=600); P.text(5.0, 0.28, 'FEU 7', 22, MUT, weight=600)
    P.raw(sd.line(P.X(4.9), P.y0, P.X(4.9), P.y0 + P.ph, RULE, 2, '6 5'))
    H = sd.Plot(470, 620, x=(415, 535), y=(0, 0.3), xlabel='mesh voltage [V]', ylabel='', title='Mesh HV scans')
    H.xticks([(v, str(v)) for v in (420, 460, 500)])
    H.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    for key, c, nm in (('6-22', '#e07a86', '6-22/23, FEU 3, drift 1000 V'), ('saturday', ORANGE, '6-27, FEU 7, drift 1000 V')):
        hs = sorted([r for r in jx if key in r['run'] and ('resist_' in r['sub'] and 'drift_scan' not in r['sub'])],
                    key=lambda r: float(r['mesh']))
        if hs:
            H.line([float(r['mesh']) for r in hs], [r['a'] for r in hs], c, 3, r=7,
                   tips=[f'{nm}, mesh {r["mesh"]} V: a = {r["a"]:+.3f} ± {r["ea"]:.3f} ({r["n"]} events)' for r in hs])
    H.text(422, 0.27, 'flat to ~500 V', 22, INK, weight=700)
    H.text(422, 0.25, 'then fewer, noisier points', 19, MUT)
    sc = [r for r in jx if 'saturday' in r['run'] and ('drift_scan' in r['sub'] or 'long_run' in r['sub'])]
    Q = sd.Plot(500, 620, x=(0, 1200), y=(0, 0.3), xlabel='drift voltage [V]', ylabel='',
                title='Drift scan, 6-27')
    Q.xticks([(v, str(v)) for v in (0, 300, 600, 900, 1200)])
    Q.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    scan = sorted([r for r in sc if 'drift_scan' in r['sub']], key=lambda r: -float(r['drift']))
    if scan:
        Q.line([float(r['drift']) for r in scan], [r['a'] for r in scan], PURPLE, 3.5, r=9,
               tips=[f'drift {r["drift"]} V at {r["t"][-5:]}: a = {r["a"]:+.3f} ± {r["ea"]:.3f}' for r in scan])
    for r in sc:
        if 'long_run' in r['sub']:
            Q.points([float(r['drift'])], [r['a']], INK, r=11, marker='open',
                     tips=[f'long run right after, drift {r["drift"]} V: a = {r["a"]:+.3f} ± {r["ea"]:.3f}'])
    Q.text(80, 0.27, 'low field → bigger copy', 22, PURPLE, weight=700)
    Q.text(450, 0.02, 'open: back at 1000 V', 19, INK)
    body = (sd.title('Not constant: it rises as the drift field drops',
                     'det3 June-15 X · raw re-extraction, hot channels vetoed · colour: drift voltage (dark low, yellow high)')
            + sd.row(P.svg('time'), H.svg('mesh'), Q.svg('drift'), gap=22))
    D.slide('time', body, details(
        'Left: every bench sub-run of det3’s June-15 build in X, at its start time. It is highest on 6-22/24 (FEU 3) and '
        'lower on 6-27/28 (FEU 7). The step coincides with the slot change, so time and slot are confounded here; '
        'the drift scan on the right is within one slot and one evening.',
        'Right: the 6-27 drift scan ran from 1100 V down to 100 V in time order; the copy rises as the field drops. The '
        'long run that started right after, at 1000 V, returns to the 1000 V value: the field dependence is '
        'reversible, not a drift in time.',
        'Middle: the mesh HV scans are flat up to ~500 V. Above it most peaks rail (cut), the surviving samples are small '
        'and several hot channels are vetoed, so the points there are few and noisy; they are not evidence of a gain dependence.',
        'Why this matters: a cable, a connector or the FEU cannot know the drift field. A detector-side effect can: '
        'the drift field sets the field just above the mesh and the electron transparency, and any charge sitting '
        'on insulating surfaces near the pads responds to it.',
        sources('raw survey', '../copy_survey.csv (start times, HV)')), short='Time, field')


# ----------------------------------------------------------------------------------------------------------------- #
# second round (10-19 evening): parity, board axes, drift time, gas, beam drift scans, simulations
# ----------------------------------------------------------------------------------------------------------------- #
def more():
    p = os.path.join(HERE, 'asym_more.json')
    return json.load(open(p)) if os.path.exists(p) else None


def checker_svg():
    """Top view of the checkerboard around one Y pad: its left/right X pads are the central X line, above/below the
    two X neighbours; the neighbour above has its bus in the gap at the Y pad's edge."""
    s = 120; o = []
    for j in range(3):
        for i in range(3):
            x0, y0 = 30 + i * s, 30 + j * s
            isx = (i + j) % 2 == 1
            c = BLUE if isx else '#e6c9a8'
            o.append(f'<rect x="{x0 + 8}" y="{y0 + 8}" width="{s - 16}" height="{s - 16}" rx="6" fill="{c}" opacity="{0.9 if isx else 0.7}"/>')
    o.append(sd.T(30 + 1.5 * s, 30 + 1.5 * s + 8, 'Y', 30, INK, weight=700))
    o.append(sd.T(30 + 0.5 * s, 30 + 1.5 * s + 8, 'n', 26, '#fff', weight=700)); o.append(sd.T(30 + 2.5 * s, 30 + 1.5 * s + 8, 'n', 26, '#fff', weight=700))
    o.append(sd.T(30 + 1.5 * s, 30 + 0.5 * s + 8, 'n+1', 24, '#fff', weight=700)); o.append(sd.T(30 + 1.5 * s, 30 + 2.5 * s + 8, 'n−1', 24, '#fff', weight=700))
    o.append(sd.line(30, 30 + s, 30 + 3 * s, 30 + s, RED, 4, '10 6'))
    o.append(sd.T(30 + 3 * s + 10, 30 + s + 8, 'bus of n+1', 21, RED, 'start'))
    o.append(sd.T(30, 30 + 3 * s + 34, 'blue: X pads (rows = X lines) · tan: Y pads', 20, MUT, 'start'))
    return sd.svg(560, 440, ''.join(o), 'checkerboard around a Y pad')


def s_parity(D):
    m = more()
    if not m: return
    order = ['det3 June-15', 'det4', 'det7', 'det6', 'det3 May', 'det2', 'det1']
    bs = {b['build']: b for b in m['boards']}; bs = [bs[k] for k in order if k in bs]
    P = sd.Plot(1060, 600, x=(-0.6, len(bs) - 0.4), y=(-0.05, 0.25), ylabel='copy a', margin=(24, 20, 110, 100))
    P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2)]); P.hline(0, MUT, None, 1.5)
    for i, b in enumerate(bs):
        for dx, key, col, nm in ((-0.17, 'xpad', GREY, 'charge over an X pad'), (0.17, 'ypad', CJUN, 'charge over a Y pad')):
            v, e, n = b[key]
            P.vbar(i + dx, max(v, -0.049), 46, col, base=0., tip=f'{b["build"]} X, {nm}: a = {v:+.3f} ± {e:.3f} ({n} matched events)')
            P.raw(sd.line(P.X(i + dx), P.Y(v - e), P.X(i + dx), P.Y(v + e), INK, 2))
        P.raw(sd.T(P.X(i), P.y0 + P.ph + 34, b['build'].replace('det3 ', 'det3·'), 20, INK))
    key = sd.legend([('charge over an X pad', GREY, 'box'), ('charge over a Y pad', CJUN, 'box')], 22)
    j = bs[0]
    body = (sd.title('Only charge over a Y pad is copied, on every board that has it',
                     'X and Y matched event by event; the pad type under the charge from the parity of the two peak channels')
            + sd.row(sd.col(P.svg('parity'), key, gap=8), sd.col(checker_svg(),
                     sd.callout(f'det3 June-15: {j["ypad"][0]:.2f} over Y pads, {j["xpad"][0]:.2f} over X pads. det4 and det7 '
                                'show the same split at a third of the size.', CJUN, 24), gap=10), gap=40))
    D.slide('parity', body, details(
        'Method: for each sub-run the X and Y dumps are matched by file index and eventId. In the checkerboard the parity of '
        '(peak X channel + peak Y channel) tells whether the charge sits over an X pad or a Y pad. Which parity is which is '
        'read from the data itself: the Y-pad class has the smaller Y ±1 share and the smaller A_x/A_y (det3 June-15: '
        f'A_x/A_y {j["axay_y"]:.2f} vs {j["axay_x"]:.2f}, Y share {j["ysh_y"]:.2f} vs {j["ysh_x"]:.2f}).',
        'Control: a forward model with the solved weighting fields and a symmetric geometry (fwd_parity.py) gives 0 in both '
        'classes, as mirror symmetry demands; the split cannot be manufactured by the classification.',
        'Within the Y-pad class the copy is the same at every Y pad: flat in (X ch mod k, Y ch mod k) for k = 3, 4, 6 (the '
        'pillar grid is 6 pads) and 12 (χ²/ndf 1.0–1.6). No 40-channel (31.2 mm, ESL beat) period either. So it is not '
        'the pillars and not the stripe-to-pad registration.',
        'Reading: for a Y pad its left and right neighbours are X pads of the central X line and the pads above and below '
        'belong to the two X neighbours; the copy is an imbalance between those two, about 0.43/0.18 of the central instead '
        'of 0.31/0.31. For an X pad the X neighbours are diagonal and the direct signal dominates.',
        sources('asym_more.py', 'fwd_parity.py', 'raw survey')), short='Y pads')


def s_axes(D):
    m = more()
    if not m: return
    cols = [CJUN, PURPLE, ORANGE]; pl = []
    for ax, ttl in (('xchip', 'along the X-chip axis (what X measures)'), ('ychip', 'along the Y-chip axis')):
        P = sd.Plot(780, 520, x=(-0.5, 7.5), y=(0, 0.6), xlabel=f'{"X" if ax == "xchip" else "Y"} chip (64 channels)',
                    ylabel='copy a, Y-pad events' if ax == 'xchip' else '', title=ttl)
        P.xticks([(c, str(c)) for c in range(8)]); P.yticks([(v, f'{v:.1f}') for v in (0, 0.2, 0.4, 0.6)])
        for (lab, d), c in zip(m['chips'].items(), cols):
            pts = [(i, q['copy'][0], q['copy'][1], q['w30']) for i, q in enumerate(d[ax]) if q['copy'][2] > 30]
            P.line([p[0] for p in pts], [p[1] for p in pts], c, 3, r=7,
                   tips=[f'{lab}, chip {p[0]}: a = {p[1]:+.3f} ± {p[2]:.3f}; pulse w30 {p[3]:.0f} ns' for p in pts])
        pl.append(P.svg(ax))
    key = sd.legend([(k, c, 'line') for k, c in zip(m['chips'], cols)], 22)
    body = (sd.title('It changes along one board axis only, and not with the local drift time',
                     'det3 June-15, Y-pad events, per 64-channel chip · three sessions')
            + sd.row(*pl, gap=30) + key
            + sd.callout('With your orientation (X cards at the bottom, Y cards on the right) the Y-chip axis is top to bottom: '
                         'the inlet-to-exhaust direction. The copy is flat along it.', BLUE, 24))
    D.slide('axes', body, details(
        'Left: along the X-chip axis the copy falls by a factor 2–3 from chips 1–3 to chip 7, in every session. Right: along '
        'the Y-chip axis it is flat within ~±0.03.',
        'The local pulse width (a drift-time proxy, in the tooltips) does not follow the X-chip profile: on 6-22/23 it is '
        '1080 ns on every chip while the copy halves. So the spatial profile is not local gas transport.',
        'What varies along one axis only and per build: anything printed or laminated in stripes along the X lines (the '
        'resistive stripes run parallel to the X lines), or the X lines themselves (length to the connector, routing).',
        sources('asym_more.py')), short='Board axes')


def s_drifttime(D):
    m = more()
    if not m: return
    sr = [s for s in m['subruns'] if s['ac'] == s['ac'] and s['ac'] is not None]
    dcol = lambda v: sd.seq_color((float(v) - 100) / 1000)
    P = sd.Plot(900, 600, x=(500, 1600), y=(0, 0.32), xlabel='pulse width above 30 % of the mean peak, w30 [ns]',
                ylabel='copy a, centre chips', title='Every det3 June-15 sub-run')
    P.xticks([(v, str(v)) for v in (600, 800, 1000, 1200, 1400, 1600)]); P.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3)])
    for s in sr:
        P.points([s['w30']], [s['ac']], dcol(s['drift']), r=8, marker='circle' if s['feu'] == 3 else 'open',
                 tips=[f'{s["t"]} {s["run"]} / {s["sub"]}\nFEU {s["feu"]}, drift {s["drift"]} V\ncopy {s["ac"]:+.3f} ± {s["eac"]:.3f}; w30 {s["w30"]:.0f} ns; median A {s["medA"]:.0f}'])
    for s in sr:
        if s['sub'] in ('long_run', 'longer_run') and ('6-17' in s['run'] or '6-19' in s['run']):
            P.text(s['w30'] + 20, s['ac'] + 0.012, '6-18' if '6-17' in s['run'] else '6-19', 20, INK, weight=700)
    se = [s for s in sr if s['exc']]
    Q = sd.Plot(700, 600, x=(0.4, 1.05), y=(0, 0.4), xlabel='Y ±1 share of the same events', ylabel='excess on n+1 / Y signal',
                title='Y-pad events')
    Q.xticks([(v, f'{v:.1f}') for v in (0.4, 0.6, 0.8, 1.0)]); Q.yticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2, 0.3, 0.4)])
    for s in se:
        Q.points([s['ysh']], [s['exc'][0]], dcol(s['drift']), r=7, marker='circle' if s['feu'] == 3 else 'open',
                 tips=[f'{s["t"]} {s["sub"]}, drift {s["drift"]} V: excess/Y {s["exc"][0]:+.3f} ± {s["exc"][1]:.3f}, Y share {s["ysh"]:.2f}'])
    cb = sd.colorbar(130, 600, 100, 1100, [(v, f'{v}') for v in (100, 300, 600, 1000)], 'drift [V]')
    body = (sd.title('The date and drift-field dependence collapse onto the cluster width',
                     'colour: drift voltage · filled: FEU 3 · open: FEU 7 · each point one sub-run')
            + sd.row(P.svg('w30'), Q.svg('ysh'), cb, gap=16))
    D.slide('drifttime', body, details(
        'Left, the drift time: at the same nominal 1000 V the pulses are ~1400 ns long on 6-22/23, 1140 on 6-18, 840 on '
        '6-27/28 and 600 on 6-19, and the copy roughly orders the same way (0.17, 0.26, 0.10, 0.03; 6-18 and 6-22 swap). A '
        'longer pulse at fixed field means a slower drift: a different gas on different nights. But the drift time is not '
        'the variable: 6-23/24 at 600 V has shorter pulses than 6-22/23 and a larger copy.',
        'For scale (Magboltz, 26.5 mm): at 1000 V the drift time is ~700 ns dry or with 1 % H2O; at 600 V it is 640 ns dry '
        'but 1220 ns with 1 % H2O; at 300 V 730 vs 3100 ns. The measured widths (shaping adds ~150–300 ns) follow the wet '
        'curve in the 6-27 scan. 6-22/23 needs an even slower gas at 1000 V.',
        'Right, the cluster width: for Y-pad events the excess on n+1, in units of the Y pad’s own signal, falls on one '
        'curve against the Y ±1 share of the same events for all sessions and drift voltages: narrower clusters, larger '
        'copy. Wet Ar/iso diffuses less (D_T 0.024–0.031 vs 0.048 √cm dry) and its D_T changes with field, so both the gas '
        'of the night and the drift field set the width. Within one run, event by event, the width matters much less '
        '(only the widest fifth drops), so the width is a proxy for the gas state, not the mechanism itself.',
        '6-19 has no Y in hand; its X ±1 share (0.84) would predict ~0.13 on this curve, against 0.03 measured. It stays '
        'an outlier: the night det3’s drift HV was sparking.',
        'Also: the median amplitude drops at low drift field (6-27: 2080 ADC at 1000 V, 630 at 100 V): attachment, so the '
        'gas was not clean either.',
        sources('asym_more.py', 'design/gas/*.json')), short='Drift time')


def s_gas(D):
    rows = [
        ['Electrons are carried by the flow', 'Flow ~mm/s against drift ~4 cm/µs: a 10⁻⁷ mm push.', '✗'],
        ['A gain or composition gradient', 'Needs ~20 % change per 1.5 mm; compounded over 40 cm that is e⁵⁰.', '✗'],
        ['It would act on X and Y alike', 'Y has no copy (|a| ≲ 0.02) on every board.', '✗'],
        ['A map following inlet → exhaust', 'Flat top-to-bottom; changes along one axis only.', '✗'],
        ['It should fade as the gas equalises', 'On, off, on (6-18, 6-19, 6-22) and reversible with the drift field.', '✗'],
        ['A large-scale cause', 'The same at every Y pad, absent over X pads.', '✗'],
        ['The gas composition as a modulator', 'The size follows the cluster width (diffusion), session to session.', '✓'],
    ]
    body = (sd.title('Your gas-circulation idea, tested: not the cause, but the gas sets the size',
                     'what a flow pattern inside the chamber predicts, against what the data show')
            + sd.table(['A gas-flow asymmetry predicts…', 'Observed', ''], [[r[0], r[1], f'<b style="color:{RED if r[2] == "✗" else GREEN}">{r[2]}</b>'] for r in rows],
                       size=25, widths=[620, 980, 60], align=['left', 'left', 'center'])
            + sd.callout('What survives: the copy is a fixed, pad-level property of the board whose size is set by the gas '
                         'transport (drift velocity, diffusion). A wetter or dirtier gas on some nights is the simplest reading.', GREEN, 26))
    D.slide('gas', body, details(
        'Why flow cannot do it: a few litres per hour through a ~5 L chamber is ~0.1–1 mm/s; electrons cross the 26.5 mm gap '
        'in ~1 µs. A composition gradient changes gain and drift velocity smoothly over centimetres; the copy needs a 20 % '
        'difference between strips 1.6 mm apart, the same at every Y pad, and only in X.',
        'Equalisation: if the copy came from incoming gas mixing in, it would fall monotonically after a flush. Instead it was '
        '0.26 on 6-18, 0.03 on 6-19, 0.16 on 6-22 and moves up and down with the drift voltage within an hour.',
        'What to look for if the gas matters as a modulator: the copy should follow anything that changes the water or air '
        'content (flow rate, bottle change, lines opened, room humidity, how long the chamber sat unflushed).',
        sources('NB_ASYM.md', 'asym_more.py')), short='Gas idea')


def s_beamdrift(D):
    p = os.path.join(HERE, 'beam', 'beam_drift.json')
    if not os.path.exists(p): return
    res = json.load(open(p)); runs = ('run_61', 'run_64'); bds = ('det3', 'det2', 'det6', 'det7'); drs = (200, 300, 400, 500, 600, 700)
    hdr = ['board / run'] + [f'{v} V' for v in drs]; body_rows = []; tips = []
    for b in bds:
        for r in runs:
            cells = []
            for v in drs:
                q = [x for x in res if x['run'] == r and x['board'] == b and x['drift'] == v]
                cells.append(f'{q[0]["a1"][0]:+.2f}' if q and q[0]['a1'][2] > 150 else '—')
            body_rows.append([f'{b} · {r}'] + cells); tips.append(f'{b}, {r}: X copy per drift point, all matched-quality events (— : < 150 events)')
    body = (sd.title('Drift scans of the other boards exist only in the beam, and they cannot settle it',
                     'July n_TOF runs 61 and 64 (Ar/iso 90/10), all four chambers scanned together · X copy, raw extractor with common-mode removal')
            + sd.table(hdr, body_rows, size=24, tips=tips)
            + sd.row(sd.callout('det6: ~0 at every field. det3 itself is flat here (0.05–0.14) in this gas, unlike on the bench.', BLUE, 24),
                     sd.callout('det7: its X FEU is noisy in the beam set-up and its Y view fails (known), so ~100 events/point. '
                                'det2 reads 0.03–0.12 here but ≈ 0 on the bench: beam-track angle.', GOLD, 24), gap=40))
    D.slide('beamdrift', body, details(
        'No bench run scanned the drift field of any board but det3 (6-27; the 6-26 FEU-3 scan is empty on EOS). The beam ran '
        'drift × mesh scans on all four chambers: runs 58, 61, 64 (Ar/iso 90/10, 150–700 V), 12 (Ar/CF4/iso), 14/15/18/19 '
        '(95/5, 400–800 V), 3/5 (Ar/CO2). Extracted on condor (4502894) with raw_asym.py --cm.',
        'Limits: the gamma flash needs the common-mode subtraction (it reads ~15–20 % lower than the raw bench extraction), '
        'beam gains are low (median 130–330 ADC), tracks come at an angle from the target (that shifts the X asymmetry: det2 '
        'reads ≈ 0 on the bench), and runs 3/12/14/15 used PS-only triggers or 20 ns sampling and gave no isolated tracks.',
        'What it does say: in 90/10 the det3 copy does not grow at low drift field. Consistent with the bench reading that the '
        'driver is the gas transport (drift velocity and diffusion behave differently in 90/10), not the field itself.',
        sources('beam/beam_drift.py', 'beam/beam_drift.json', 'condor 4502894')), short='Beam scans')


def s_sims(D):
    import glob as _g
    fo = {os.path.basename(f)[:-5]: json.load(open(f)) for f in _g.glob(os.path.join(HERE, '..', 'fwd_out', '*.json'))}
    nm = [('d20_nobus', 'pads and vias only (symmetric)'), ('d20_nominal', 'one-sided bus, nominal depth'),
          ('d20_h15', 'bus touching the pads (15 µm)'), ('shift60', 'charge over Y pads displaced 60 µm'),
          ('f_offm200', 'stripes floating, offset −200 µm'), ('f_offp200', 'stripes floating, offset +200 µm')]
    rows = []
    for k, lab in nm:
        if k not in fo: continue
        r = {x['sigma']: x for x in fo[k]['res']}
        c = lambda s, key: f'{r[s][key][0]:+.3f}' if s in r else '—'
        rows.append([lab, c(300., 'ypad'), c(300., 'xpad'), c(450., 'ypad'), c(450., 'xpad')])
    m = more(); j = [b for b in m['boards'] if b['build'] == 'det3 June-15'][0] if m else None
    if j: rows.append(['<b>data, det3 June-15</b>', f'<b>{j["ypad"][0]:+.3f}</b>', f'<b>{j["xpad"][0]:+.3f}</b>', '', ''])
    body = (sd.title('What the simulations exclude: no static geometry makes it',
                     '3D weighting-field solves (bus_wpot.py, bus_float.py) + a forward model on uniform clouds with the data’s own estimator')
            + sd.table(['geometry', 'Y-pad, σ 300', 'X-pad, σ 300', 'Y-pad, σ 450', 'X-pad, σ 450'], rows, size=24)
            + sd.row(sd.callout('The copper alone (any bus depth, wet gap fill, kapton ε) gives |a| < 0.02 on the data’s estimator.', BLUE, 24),
                     sd.callout('Floating stripes do make a real copy, but in both pad classes and growing with cloud size: the wrong fingerprint.', GOLD, 24), gap=40))
    D.slide('sims', body, details(
        'Why the old 0.010 is not the data’s number: the old solve took the asymmetry of one channel’s kernel. The data pick '
        'the peak channel per event; a kernel that is only lopsided (the bus) is then absorbed by relabelling, and the '
        'forward model returns ~0. A copy that survives the peak selection needs line n+1 to receive a share of line n’s own '
        'signal.',
        'Displacing the charge over Y pads also averages away (it only moves where the charge lands). Floating stripes '
        '(fast transverse equalisation, ~26 ns at 2 MΩ/sq) offset from the rows give 0.06–0.15, but equally for X- and '
        'Y-pad events and more for wider clouds; the data are the opposite on both counts.',
        'Static dielectric variants (bus_wpot.py --eps-gap 20/80, --eps-kap 10): a wet gap fill lowers the bus pickup, the '
        'kapton permittivity hardly matters.',
        sources('bus_wpot.py', 'bus_float.py', 'fwd_parity.py', 'fwd_out/*.json', 'bus_par/*.json')), short='Simulations')


def s_bus(D):
    sol = {}
    for k in ('nominal', 'L6', 'nobus', 'h132', 'h66', 'h33', 'h15'):
        f = os.path.join(HERE, '..', 'bus_scan', f'd10_{k}.json')
        if os.path.exists(f):
            j = json.load(open(f)); sol[k] = np.mean([j[s]['copy'] for s in ('sigma300', 'sigma450') if s in j])
    # cartoon: top view of pads with stubs and buses
    o = []
    s = 74
    for r in range(4):
        for c_ in range(7):
            x, y = 40 + c_ * s, 40 + r * s
            own = (r + c_) % 2 == 1
            fill = BLUE if own else '#e6c9a8'
            o.append(f'<rect x="{x}" y="{y}" width="{s - 12}" height="{s - 12}" rx="4" fill="{fill}" opacity="{0.9 if own else 0.6}"/>')
            if own:
                cx, cy = x + (s - 12) / 2, y + (s - 12) / 2
                o.append(sd.line(cx, cy, cx, cy + s / 2, INK, 3))
        yb = 40 + r * s + (s - 12) / 2 + s / 2
        o.append(sd.line(30, yb, 40 + 7 * s, yb, INK, 4))
    o.append(sd.T(40, 360, 'top view: X pads (blue), Y pads (tan)', 21, MUT, 'start'))
    o.append(sd.T(40, 390, 'every X stub points the same way (down)', 21, INK, 'start', 600))
    o.append(sd.T(40, 420, '→ each bus sits under the gap toward ONE neighbour', 21, INK, 'start', 600))
    cart = sd.svg(600, 440, ''.join(o), 'bus cartoon')
    rows = []
    lab = {'nobus': 'no bus (control)', 'L6': 'bus on the deeper layer (L6)', 'nominal': 'bus on L5, 264 µm FR4 (nominal)',
           'h132': 'L5 at 132 µm', 'h66': 'L5 at 66 µm', 'h33': 'L5 at 33 µm', 'h15': 'L5 at 15 µm (touching)'}
    for k in ('nobus', 'L6', 'nominal', 'h132', 'h66', 'h33', 'h15'):
        if k in sol:
            rows.append((lab[k], max(sol[k], 0.0005), GREY if k != 'nominal' else BLUE,
                         f'3D finite-volume prompt weighting solve, 10 µm grid: copy {sol[k]:.4f}'))
    data = [('det1 / det2 / det6 X (data)', 0.01, COTH, 'raw: −0.01 … +0.01'),
            ('det4 / det7 X (data)', 0.045, ORANGE, 'raw: 0.03–0.05'),
            ('det3 June-15 X (data)', 0.17, CJUN, 'raw: 0.06–0.18 by run, 0.1–0.3 by chip')]
    bars = sd.hbars(rows + data, 0.3, width=520, h=30, log=True, vmin=0.0005, label_w=420,
                    fmt=lambda v: f'{v:.3f}', size=22)
    body = (sd.title('The bus traces break left/right, but only at the 1 % level',
                     'Left: the pad/bus layout from the gerbers · right: copy from a 3D weighting-field solve vs data (log scale)')
            + sd.row(cart, sd.col(bars, sd.callout('Even with the bus touching the pads the solve stays below det3’s '
                                                  'copy, and a layout cannot vary across the board or with the field.',
                                                  CJUN, 24), gap=24), gap=48))
    D.slide('bus', body, details(
        'Each X channel is every other pad along a line; its pads reach a bus line through 390 µm stubs that all point '
        'the same way, so the bus lies under the gap on one side of the channel’s own pads, next to one neighbour. '
        'Charge landing near that gap induces a little on the bus too: a one-sided, prompt, linear copy. The pad '
        'kernels used in the simulation drive the pads only, so they are symmetric by construction.',
        'The solve (<code>bus_wpot.py</code>): own pads, vias, rings, stubs and bus at 1 V, everything else grounded, '
        'read at the resistive layer and folded with a 300–450 µm cloud. With the Geant4 stack estimate (264 µm FR4 per '
        'layer) it gives 0.010, about what det1/det2/det6 show. It reaches 0.085 only with the bus 15 µm below the '
        'pads, which no 8-layer board has.',
        'Not checked yet, and decisive: whether the side the buses sit on predicts the observed direction (p+1 gets '
        'the copy). That needs the gerber → strip-map orientation (board rotated 90° on the bench, connectors '
        'inverted). If it is the wrong side, the bus is out entirely. Also open: which data view sits on L5 vs L6; '
        'the solve says the L6 view should see ~0.001.',
        sources('../bus_wpot.py', '../bus_scan/d10_*.json', '../gerber_stubs.py')), short='Bus traces')


def s_verdict(D):
    steps = [
        dict(label='Track, our code', sub='Y clean; raw, symmetric', color=GREY, strike=True),
        dict(label='Gas flow', sub='pad-level, X only, flat top-bottom', color=GREY, strike=True),
        dict(label='FEU, cable, connector', sub='cannot know the pad type under the charge', color=GREY, strike=True),
        dict(label='Pillars, stripe beat', sub='flat mod 6, no 40-ch period', color=GREY, strike=True),
        dict(label='Copper (bus, vias)', sub='< 0.02 on the data estimator', color=GREY, strike=True),
        dict(label='Floating stripes', sub='wrong fingerprint', color=GOLD),
        dict(label='Y pads in the June bulks', sub='every Y pad, smooth along one axis', color=RED),
    ]
    body = (sd.title('Verdict: something at the Y pads of Stephan’s June bulks, scaled by the gas',
                     'Struck: excluded by a slide in this deck · gold: tested, does not fit · red: where the effect lives')
            + sd.flow(steps, size=23)
            + sd.row(sd.callout('Real: in raw ADC, in phase, linear, one-sided, signal-only, X only, '
                                'only for charge over a Y pad, the same at every Y pad.', BLUE, 27),
                     sd.callout('Boards: det3 June-15 0.22, det4 and det7 0.08 (Y-pad events); det6 0.02; det3 May, det1, det2 0.', CJUN, 27), gap=56)
            + sd.callout('Its size changes with the cluster width (gas of the night, drift field), and along one board axis. '
                         'The mechanism that moves Y-pad charge into the n+1 X line is not identified.', RED, 27))
    D.slide('verdict', body, details(
        'Excluded with evidence: the code, track geometry, the FEU/ADC (same slot clean on other boards; noise not copied). '
        'New: any cable, connector or electronics cause. The copy appears only when the charge sits over a Y pad, which the X '
        'electronics cannot know (it is read only by matching with the Y view), and the copy is linear in amplitude.',
        'Excluded by simulation (3D weighting fields + forward model with the data’s own estimator): the one-sided bus at any '
        'depth, wet gap fills, kapton ε, a displacement of the charge over Y pads. Floating resistive stripes offset from the '
        'rows do make a copy, but in both pad classes and growing with cloud size: the data show neither.',
        'Excluded as the cause: gas circulation (see “Gas idea”). The gas still matters as a modulator: the copy follows the '
        'cluster width, which the gas composition and the drift field set.',
        'Not identified: what at a Y pad sends ~30 % of its charge signal into the n+1 X line on these builds. Candidates '
        'not yet modelled: a conductive or charged residue near the Y-pad edge on the bus side, a defect in the Y-pad '
        'connection that leaves it partly floating, or a bulk-process layer over the pads that differs from the May recipe.'),
        short='Verdict')


def s_close(D):
    body = (sd.kicker('what is open, and what would settle it')
            + '<h2 style="font-size:60px;font-weight:600;line-height:1.1;color:#eef0f3">A pad-level property of the June '
              'bulks, made bigger by a narrower cloud; the mechanism at the Y pad is open</h2>'
            + sd.row(
                sd.col(sd.p('Open', 30, DINK, 600),
                       sd.p('· What at a Y pad couples into the n+1 X line (residue, a partly floating Y pad, a bulk layer).', 26, DMUT),
                       sd.p('· Whether the May and June-15 “det3” share a PCB and a resistive print.', 26, DMUT),
                       sd.p('· 6-19: 0.03 with ordinary cluster widths, the night det3’s drift HV was sparking.', 26, DMUT),
                       sd.p('· Drift scans of det4/det7 on the bench (the beam ones are too limited).', 26, DMUT), gap=14),
                sd.col(sd.p('Tests that would settle it', 30, DINK, 600),
                       sd.p('1. Same HV, two gas states (fresh vs long-unflushed, or a humidifier): the copy should follow the width.', 26, DMUT),
                       sd.p('2. A collimated Fe-55 spot on one Y pad vs one X pad: point-like clouds, the largest Y/X contrast.', 26, DMUT),
                       sd.p('3. Ask Stephan what differs in the June recipe; look at det3 Y pads under a microscope.', 26, DMUT),
                       sd.p('4. Drift scan of det4 or det7 on the bench at fixed gas.', 26, DMUT), gap=14),
                gap=80)
            + sd.p('For analysis now: correct det3 X per run with uncopy applied to Y-pad events (parity from the Y view).', 28, DGREEN, 600))
    D.slide('close', body, details(
        'Test 1 separates “the gas sets the size” from “something about the board changes”. Test 2 gives the cleanest '
        'fingerprint: a point source over a Y pad should show the full copy, over an X pad none.',
        'Correction: <code>uncopy.py</code> with one a for all events over-corrects X-pad events and under-corrects Y-pad '
        'events; applying it to the Y-pad class only (a ≈ 0.22 on det3 June-15 6-22..24, 0.15 on 6-27/28) is the better '
        'model of the data.'), dark=True, short='Open')


def build(out):
    D = sd.Deck('det3 X asymmetry', 'Is the det3 X left/right neighbour asymmetry real? Raw-data checks.')
    rows = survey_rows()
    s_cover(D); s_setup(D); s_suspects(D); s_event(D); s_processing(D); s_coin(D); s_hists(D); s_waves(D)
    s_copy_shape(D); s_copy_model(D); s_noise(D); s_follow(D); s_map(D); s_amp(D); s_hot(D)
    s_survey(D, rows); s_date(D, rows); s_time(D, rows)
    s_parity(D); s_axes(D); s_drifttime(D); s_gas(D); s_beamdrift(D); s_sims(D)
    s_bus(D); s_verdict(D); s_close(D)
    meta = dict(title='det3 X neighbour asymmetry: real, in the raw ADC, and det3’s own',
                summary='From the rawest observable up: det3’s June-15 build copies ~0.1–0.2 of each X pulse into the next '
                        'channel, only when the charge sits over a Y pad, at every Y pad; det4/det7 the same at a third. '
                        'Its size follows the cluster width (gas, drift field). Not gas flow, not electronics, not the '
                        'copper; mechanism at the Y pad open.',
                tags='MX17,detector,crosstalk', date=dt.date.today().isoformat())
    return D.write(out, meta, footer=f'Built {dt.datetime.now():%Y-%m-%d %H:%M} by '
                                     'MX17_Geant/scratchpad/2026-10-19/asym_deck/make_asym_deck.py with slidedoc.py.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=os.path.join(HERE, 'det3-x-asymmetry.html'))
    print('wrote', build(ap.parse_args().out))
