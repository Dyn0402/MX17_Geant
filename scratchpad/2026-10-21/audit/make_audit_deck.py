"""det3 X neighbour asymmetry: audit deck (2026-10-21). A clean remake of the det3-x-asymmetry note as an
up-to-date summary: observation -> fingerprint -> what is excluded (and how solid) -> leading candidates and their
problems -> way forward. Every number is read from the analysis outputs:
  this dir:   audit_data.json (audit_data.py + sessions.py), bands.json, share_vs_a1_nc8.json, late_vs_a1.json,
              toy_curves.json
  ../yview:   yview_class.json (main session: Y view by pad class, route 14 within board)
  ../../2026-10-20/r5: sim_delta.json, timing_width*.json, paired_odd*.json, probe_charge.json
  ../../2026-10-20:    edge_test.json, strip_orient2.json;  ../../2026-10-19/noise_corr_inject.log
Build:   ~/PycharmProjects/nTof_x17/.venv/bin/python make_audit_deck.py   -> det3-x-asymmetry.html (this dir)
Publish: python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py det3-x-asymmetry.html --slug det3-x-asymmetry
         --force --deploy"""
import os, sys, json, re, math
import numpy as np
sys.path.insert(0, os.path.expanduser(os.environ.get('SLIDEDOC_DIR', '~/PycharmProjects/dylan-cern-site/scripts')))
import slidedoc as sd
from slidedoc import (BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, INK, MUT, RULE, DINK, DMUT, DBLUE, DRED, DGREEN,
                      DGREY, term, Plot, legend, row, col, p, card, callout, table, T, line, svg, arrow, esc)

HERE = os.path.dirname(os.path.abspath(__file__))
R5 = os.path.join(HERE, '..', '..', '2026-10-20', 'r5')
S20 = os.path.join(HERE, '..', '..', '2026-10-20')
S19 = os.path.join(HERE, '..', '..', '2026-10-19')
def J(*a): return json.load(open(os.path.join(*a)))
AD = J(HERE, 'audit_data.json'); BANDS = J(HERE, 'bands.json'); SHARE = J(HERE, 'share_vs_a1_nc8.json')
LATE = J(HERE, 'late_vs_a1.json'); TOY = J(HERE, 'toy_curves.json')
YV = J(HERE, '..', 'yview', 'yview_class.json') if os.path.exists(os.path.join(HERE, '..', 'yview', 'yview_class.json')) else None
SIMD = J(R5, 'sim_delta.json')['col']; PROBE = J(R5, 'probe_charge.json'); EDGE = J(S20, 'edge_test.json')
TW = J(R5, 'timing_width.json'); TWA = J(R5, 'timing_width_glax10_glax10w690f1.json')

B = AD['build']
COL = {'det3 June-15': RED, 'det2': ORANGE, 'det4': PURPLE, 'det6': GOLD, 'det7': GREEN, 'det3 May': BLUE, 'det1': GREY}
JUNE = ['det3 June-15', 'det2', 'det4', 'det6', 'det7']; QUIET = ['det3 May', 'det1']
DESC = {'det3 June-15': 'det3, bulked by Stephan June 15', 'det2': 'det2, bulked by Arnaud June 12',
        'det4': 'det4, Stephan June-22 batch (board C)', 'det6': 'det6, Stephan June 24 (board D)',
        'det7': 'det7, Stephan June-22 batch (board B)', 'det3 May': 'det3, bulked May 5 (same FEU 3 as June-15)',
        'det1': 'det1, January build (strips run along x)'}
f3 = lambda v: f'{v:+.3f}'

def ebar(P, x, v, e, color, cap=8):
    X = P.X(x); a, b = P.Y(v - e), P.Y(v + e)
    P.raw(line(X, a, X, b, color, 2.5) + line(X - cap, a, X + cap, a, color, 2.5) + line(X - cap, b, X + cap, b, color, 2.5))

def field_rms(cells, key):
    """Typical local |a1|: rms of the cell means about zero, with the statistical noise subtracted."""
    v = np.array([c[key][0] for c in cells]); e = np.array([c[key][1] for c in cells])
    return float(np.sqrt(max((v ** 2).mean() - (e ** 2).mean(), 0))), len(v)

D = sd.Deck('det3 X asymmetry: audit',
            'An audit of the week-long det3 X neighbour-asymmetry study: what the data show, what is excluded, '
            'the leading candidates and their problems, and the way forward.')
sl = D.slide

# ------------------------------------------------------------------------------------------------ 1 cover
b3 = B['det3 June-15']
n_sub = len([s for s in AD['survey']])
cover = (sd.kicker('Audit · MX17 micromegas · June-2026 cosmic bench · raw DREAM ADC')
         + '<h1 style="font-size:70px;font-weight:600;line-height:1.08;letter-spacing:-2px;width:1660px">'
           'The June boards’ X view shifts charge from one neighbour to the other without changing how much they share. '
           'Resistive-strip registration, as simulated, can’t do that</h1>'
         + p('A fresh look at a week of work: the observation is solid and well mapped. Many causes are firmly '
             'excluded. The leading explanation explains the map but fails a test it was never scored on, so the '
             'way forward is to change what the models are scored on, and to look at the strips themselves.',
             30, DMUT, extra='width:1560px')
         + '<div style="flex:1"></div>'
         + row(sd.bignum('7', 'board builds surveyed', DRED, '5 June builds, 2 quiet controls (det3 May, det1)',
                         tip='det2, det3 June-15, det4, det6, det7 (June 2026 bulks); det3 May build and det1 (January)'),
               sd.bignum(f'{b3["n"]:,}', 'matched X/Y events on det3 June-15', DBLUE,
                         'single clusters, X and Y matched by event, edge channels cut',
                         tip=f'r5lib.ev: {n_sub} survey sub-run views in the raw-ADC survey'),
               sd.bignum('40', 'routes tried in a week', DGREEN, 'ledger in FINDINGS.md; re-graded here',
                         tip='scratchpad/2026-10-19/asym_deck/FINDINGS.md, routes 1–40; audit in scratchpad/2026-10-21/audit/AUDIT.md'),
               gap=56))
sl('cover', cover, 'This deck replaces the 42-slide history at the same address. The audit text, with every route re-graded and '
   'the new checks, is <code>scratchpad/2026-10-21/audit/AUDIT.md</code> (MX17_Geant). The previous deck is in git: '
   '<code>scratchpad/2026-10-19/asym_deck/det3-x-asymmetry.html</code>.', dark=True, short='Cover')

# ------------------------------------------------------------------------------------------------ 2 setup diagram
def setup_svg():
    W, H = 900, 640; o = []; P0 = 110; x0, y0 = 70, 40
    # strips (along y), drawn first, offset by delta
    dlt = 34; sw = 0.69 * P0
    for k in range(5):
        cx = x0 + (k + 0.5) * P0 + dlt
        o.append(f'<rect x="{cx - sw / 2:.1f}" y="{y0 - 18}" width="{sw:.1f}" height="{5 * P0 + 36}" fill="#8a8f99" fill-opacity="0.22" rx="6"/>')
    for i in range(5):
        for j in range(5):
            xp = (i + j) % 2 == 0
            c = BLUE if xp else ORANGE
            cx, cy = x0 + (i + 0.5) * P0, y0 + (j + 0.5) * P0; s = P0 * 0.36
            o.append(f'<polygon points="{cx},{cy - s} {cx + s},{cy} {cx},{cy + s} {cx - s},{cy}" fill="{c}" fill-opacity="0.85"/>')
    # cloud on the Y pad at column n (i=2), row j=1 (Y pad when (i+j) odd)
    cx, cy = x0 + 2.5 * P0, y0 + 1.5 * P0
    o.append(f'<circle cx="{cx}" cy="{cy}" r="58" fill="{RED}" fill-opacity="0.22" stroke="{RED}" stroke-width="3"/>')
    o.append(T(cx, cy + 8, 'e⁻', 26, RED, weight=600))
    for i, lab in ((1, 'n−1'), (2, 'n'), (3, 'n+1')):
        o.append(T(x0 + (i + 0.5) * P0, y0 + 5 * P0 + 52, lab, 26, INK, weight=600))
    o.append(T(x0 + 2.5 * P0, y0 + 5 * P0 + 90, 'X channels = pad columns (read the x position)', 22, MUT))
    o.append(arrow(x0 + 5 * P0 + 30, y0 + 4.5 * P0, x0 + 5 * P0 + 30, y0 + 0.5 * P0, MUT))
    o.append(T(x0 + 5 * P0 + 60, y0 + 2.5 * P0, 'strips run along y', 22, MUT, rot=90))
    o.append(T(x0 + 3.5 * P0 + dlt, y0 - 26, 'δ: strip centre − pad centre', 21, GREY))
    return svg(W, H, ''.join(o), 'checkerboard pads with resistive strips')
setup = sd.title('What is measured: on the X view, the two neighbours of the peak column are unequal',
                 'Checkerboard readout: X pads (blue) join into columns, Y pads (orange) into rows; resistive strips (grey) on top')
setup += row(setup_svg(), col(
    legend([('X pad', BLUE, 'box'), ('Y pad', ORANGE, 'box'), ('resistive strip', '#b9bcc3', 'box')]),
    callout(f'{term("a1", "a1 = α(n+1) − α(n−1), DAQ channel order. α = the neighbour’s in-phase amplitude divided by the central channel’s (projection of the neighbour waveform on the central pulse).")}'
            ' = α(n+1) − α(n−1): how much bigger the upper neighbour is', RED),
    callout(f'{term("s1", "s1 = α(n+1) + α(n−1): the ±1 sharing. s1X on the X view; s1Y (= the Y share) on the Y view measures the cloud width independently.")}'
            ' = α(n+1) + α(n−1): how much the two neighbours share in total', BLUE),
    callout(f'{term("Pad type", "Old label: parity of (X peak + Y peak). It uses the X peak, which the asymmetry moves, so it mislabels events (peak migration). X-independent label: the Y share s, small when the charge sits on one Y row (over a Y pad), large when it splits between two (over an X pad).")}'
            ' under the charge alternates along y; the Y share tells it without using X', ORANGE),
    p('Dotted words and the points on every plot carry tooltips: values, counts, sources.', 22, MUT), w=640, gap=26), gap=60)
sl('setup', setup, '<p>Observable: per event, the peak X channel n and its neighbours n±1. α is the in-phase content (each neighbour waveform '
   'projected on the central pulse), so it is a charge fraction, insensitive to small timing shifts. a1 > 0 means DAQ channel n+1 '
   'carries more than n−1. With uniform illumination and a symmetric response, a1 = 0 on average.</p>'
   '<p>Extractor: <code>raw_asym.py</code> (pedestal-subtracted DREAM ADC, DAQ channel numbers, no strip map, no common-mode subtraction, '
   'no clustering, no M3). Good-pulse cut 200 &lt; A_c &lt; 3500, central peak sample 4–22, quiet pre-samples. Matched X/Y events by '
   '(file, eventId) in <code>scratchpad/2026-10-20/r5/r5lib.py</code>; channels within 4 of a connector edge cut.</p>'
   '<p>The strips run along y on every May and June build (the slow, along-strip view is Y on all of them); det0/det1 run the other way.</p>',
   short='Setup')

# ------------------------------------------------------------------------------------------------ 3 board means vs field strength
order = JUNE + QUIET
P1 = Plot(810, 560, x=(-0.6, len(order) - 0.4), y=(-0.04, 0.16), ylabel='board mean a1', title='Board mean')
P1.yticks([(v, f'{v:+.2f}' if v else '0') for v in (-0.04, 0, 0.04, 0.08, 0.12, 0.16)])
P1.xticks([(i, b.replace('det3 June-15', 'det3 J15').replace('det3 May', 'det3 May')) for i, b in enumerate(order)])
P1.hline(0, MUT, None, 1.5)
P2 = Plot(810, 560, x=(-0.6, len(order) - 0.4), y=(0, 0.16), ylabel='typical local |a1| (8 × 8 cells)', title='Local size of the asymmetry')
P2.yticks([(v, f'{v:.2f}') for v in (0, 0.04, 0.08, 0.12, 0.16)])
P2.xticks([(i, b.replace('det3 June-15', 'det3 J15')) for i, b in enumerate(order)])
for i, b in enumerate(order):
    for dx, key, op in ((-0.17, 'a1x', 1.0), (0.17, 'a1y', 0.35)):
        v, e, n = B[b][key]
        tip = f'{DESC[b]}\n{key[-1].upper()} view: a1 = {v:+.3f} ± {e:.3f}, n = {n}'
        P1.vbar(i + dx, v, 30, COL[b], tip=tip, base=0, opacity=op); ebar(P1, i + dx, v, e, INK, 5)
    for dx, key, op in ((-0.17, 'ax', 1.0), (0.17, 'ay', 0.35)):
        r, nc = field_rms(AD['maps'][b], key)
        P2.vbar(i + dx, r, 30, COL[b], tip=f'{DESC[b]}\n{"X" if key == "ax" else "Y"} view: noise-subtracted rms about zero of {nc} cell means = {r:.3f}', base=0, opacity=op)
P1.text(1.2, 0.135, 'det2 ≈ 0 on average…', 24, INK, weight=600)
P2.text(1.4, 0.15, '…but as large as det3’s locally', 24, INK, weight=600)
s3 = sd.title('Only the June boards’ X view is lopsided, and det2’s board average hid it',
              'Left: board-mean a1; right: the typical |a1| in a 64 × 64 cell. Solid = X view, pale = Y view; matched events')
s3 += row(P1.svg('board means'), P2.svg('field strength'), gap=44)
s3 += legend([('X view', INK, 'box'), ('Y view (pale)', '#c9ccd2', 'box'), ('June builds: det3 J15 red, det2 orange, det4 purple, det6 gold, det7 green', MUT, 'dot'),
              ('quiet controls: det3 May blue, det1 grey', BLUE, 'dot')], 22)
sl('means', s3, '<p>Local size = sqrt(mean(cell mean²) − mean(cell error²)) over the 8 × 8 map (64 × 64 columns × rows per cell, ≥ 25 events), '
   'i.e. the typical |a1| in a cell beyond statistical noise, whatever its sign. det2’s field is as strong as det3’s but changes sign in bands, so its '
   'board mean is zero. That is why det2 served as a "copy-free control" until the maps were made (route 23), and why several early '
   'conclusions that used it are re-checked later in this deck.</p>'
   '<p>A small common a1Y ≈ −0.01 appears on quiet and June boards alike: an instrumental baseline.</p>'
   '<p>Source: <code>audit_data.py</code> (r5lib matched events, edge cut).</p>', foot='Raw-ADC survey, cosmic bench June–July 2026; ±1σ statistical errors.', short='Board means')

# ------------------------------------------------------------------------------------------------ 4 maps
def mapplot(b, key, ttl, xl=True):
    P = Plot(372, 340, x=(0, 512), y=(0, 512), xlabel='true column' if xl else '', ylabel='row', margin=(16, 12, 62 if xl else 30, 72), title=ttl)
    P.xticks([(0, '0'), (256, '256'), (512, '512')] if xl else []); P.yticks([(0, '0'), (256, '256'), (512, '512')])
    cells = []
    for c in AD['maps'][b]:
        v, e = c[key]; cells.append((c['i'] * 64, (c['i'] + 1) * 64, c['j'] * 64, (c['j'] + 1) * 64, v,
                                     f'{b}, {"X" if key == "ax" else "Y"} view, columns {c["i"] * 64}–{c["i"] * 64 + 63}, rows {c["j"] * 64}–{c["j"] * 64 + 63}\na1 = {v:+.3f} ± {e:.3f}, n = {c["n"]}'))
    P.cells(cells, -0.3, 0.3, sd.DIVERGE, gap=1.0)
    return P.svg(ttl)
mp = []
for key, lab in (('ax', 'X'), ('ay', 'Y')):
    mp.append(row(*[mapplot(b, key, f'{b.replace("June-15", "J15")} · {lab}', key == 'ay') for b in ('det3 June-15', 'det2', 'det4', 'det3 May')], gap=8))
cb = sd.colorbar(130, 620, -0.3, 0.3, [(-0.3, '−0.3'), (-0.15, '−0.15'), (0, '0'), (0.15, '+0.15'), (0.3, '+0.3')], 'a1', sd.DIVERGE)
s4 = sd.title('A smooth field over each June board, X view only', 'a1 per 64 × 64 cell (true columns × rows), all matched events; top X view, bottom Y view')
s4 += row(col(*mp, gap=4), cb, gap=24)
sl('maps', s4, '<p>Constant-model χ²/ndf of the X maps (8 × 8): det3 June-15 332/63, det2 736/61, det6 47/25, det7 21/8; controls det3 May 78/56, '
   'det1 34/27. Y maps: det3 June-15 74/63, det2 59/61. The main session’s 4 × 4 maps by pad class (Y-pad class): det3 X 438/15 vs Y 30/15; '
   'det2 X 971/15 vs Y 21/15; det2’s X-pad-class map is flat (15/15).</p>'
   '<p>det2’s two runs (6-17 FEU 3, 6-22 FEU 6) show the same bands; det3’s two FEU-3 sessions give the same map phase at different '
   'amplitude. det4 changes sign across x. No step at any of the 7 connector boundaries (det3: every step ≤ 0.02 ± 0.03).</p>'
   '<p>Caveat: these maps pool sessions. det3’s cloud width also varies across x (Y share 0.89 at chip 0 to 0.48 at chip 7), and the amplitude '
   'depends on the width, so the map shape is the field times a width pattern. A per-session, width-matched map is on the to-do list.</p>',
   short='Maps')

# ------------------------------------------------------------------------------------------------ 5 band profiles
P = Plot(1100, 640, x=(-0.3, 7.3), y=(-0.36, 0.32), xlabel='band (64 rows for det2, 64 columns = one DREAM chip for the others)', ylabel='a1, all matched events')
P.xticks([(i, str(i)) for i in range(8)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)])
P.hline(0, MUT, None, 1.5)
for key, b, lab in (('det2|y', 'det2', 'det2, by row band'), ('det3 June-15|x', 'det3 June-15', 'det3 June-15, by column'),
                    ('det4|x', 'det4', 'det4, by column'), ('det6|x', 'det6', 'det6, by column'), ('det3 May|x', 'det3 May', 'det3 May, by column')):
    rr = [r for r in BANDS[key] if r['cls'] == 'all']
    xs = [r['band'] for r in rr if r['a1'][0] is not None and r['a1'][0] == r['a1'][0]]
    ys = [r['a1'][0] for r in rr if r['a1'][0] is not None and r['a1'][0] == r['a1'][0]]
    tips = [f'{lab}, band {r["band"]}: a1 = {r["a1"][0]:+.3f} ± {r["a1"][1]:.3f}, n = {r["n"]}' for r in rr if r['a1'][0] == r['a1'][0]]
    P.line(xs, ys, COL[b], tips=tips, dash='8 6' if b == 'det3 May' else None)
P.text(0.2, 0.27, 'Two signs on one board, smooth from band to band', 25, INK, weight=600)
s5 = sd.title('The field changes sign smoothly: by rows on det2, across columns on det4 and det6',
              'a1 of all matched events per 64-row (det2) or 64-column (others) band')
s5 += row(P.svg('band profiles'), col(legend([('det2 (rows)', ORANGE), ('det3 June-15', RED), ('det4', PURPLE), ('det6', GOLD), ('det3 May (control)', BLUE, 'dash')], 24),
                                    p('det2’s positive band is a small det3: same coin flip, same timing. Its outer bands are the mirror image.', 24, MUT), w=520, gap=30), gap=44)
sl('bands', s5, '<p>det2 inside its positive band (rows 128–255): a1 +0.142, coin flip 0.71, narrow → wide X-pad class −0.18 → +0.12, '
   'excess timing +27 / −126 ns (det3: +20 / −145). Rows 0–127 and 384–511: a1 −0.058, coin 0.40, narrow X-pad class +0.31 (route 23).</p>'
   '<p>Inside that band, a1 is positive in every 8-column block across all 512 columns, with no 1/40 harmonic (p = 0.57): whatever causes it '
   'repeats at the 780 µm pad pitch, not at 800 µm (route 24).</p><p>Source: <code>bands.py</code> → <code>bands.json</code>.</p>', short='Bands')

# ------------------------------------------------------------------------------------------------ 6 same over X and Y pads
YS = J(HERE, '..', 'm3class', 'ysh_class.json')
def mid(r): return 0.5 * (r['lo'] + r['hi'])
def first(v): return v[0] if isinstance(v, list) else v
P = Plot(1080, 620, x=(0, 1.75), y=(-0.5, 0.35), xlabel='Y share s (measured on the Y view only)', ylabel='a1 on the X view', margin=(24, 120, 92, 104))
P.xticks([(v, f'{v:g}') for v in (0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75)])
P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)])
m_ = lambda v: -0.5 + v * 0.85
P._y2f = 1; P._y2t = [(m_(v), f'{v:.1f}') for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)]; P.ylabel2 = 'det3 May: fraction over Y pads'
P.hline(0, MUT, None, 1.5)
d3 = [r for r in YS['det3 June-15'] if r['a1']]; dm = [r for r in YS['det3 May'] if r.get('fY') is not None]
P.line([mid(r) for r in dm], [m_(r['fY']) for r in dm], GREY, dash='4 5', tips=[f'det3 May (quiet, labels can’t migrate), s {r["lo"]}–{r["hi"]}: fraction over Y pads {r["fY"]:.2f}, n = {r["n"]}' for r in dm])
a0 = max(d3[0]['a1'][0], d3[1]['a1'][0]); f0 = dm[0]['fY']
P.line([mid(r) for r in dm], [a0 * r['fY'] / f0 for r in dm], RED, dash='2 6', markers=False, w=3,
       tip='If the asymmetry lived only over the Y pads, it would fall with the Y-pad fraction (scaled to the narrowest bins)')
for key, c, lab in (('a1_Ylab', '#e8a0a8', 'old label “Y pad”'), ('a1_Xlab', '#7fa7d1', 'old label “X pad”')):
    q = [r for r in d3 if r.get(key) is not None]
    P.line([mid(r) for r in q], [first(r[key]) for r in q], c, dash='8 6', r=4,
           tips=[f'det3 J15, s {r["lo"]}–{r["hi"]}, {lab}: a1 = {first(r[key]):+.3f}' for r in q])
P.line([mid(r) for r in d3], [r['a1'][0] for r in d3], RED, tips=[f'det3 June-15, s {r["lo"]}–{r["hi"]}: a1 = {r["a1"][0]:+.3f} ± {r["a1"][1]:.3f}, n = {r["n"]}' for r in d3])
d2 = [r for r in YS['det2 band+'] if r['a1']]
P.line([mid(r) for r in d2], [r['a1'][0] for r in d2], ORANGE, tips=[f'det2 positive band (rows 128–255), s {r["lo"]}–{r["hi"]}: a1 = {r["a1"][0]:+.3f} ± {r["a1"][1]:.3f}, n = {r["n"]}' for r in d2])
P.text(0.05, 0.3, 'Flat: the same over Y pads and X pads', 25, INK, weight=600)
s6 = sd.title('The asymmetry is the same over X pads and Y pads; the old pad-class split was peak migration',
              'a1 vs the Y share s, which tells the pad type without using X (det3 May: 80 % over Y pads at small s, 20 % at large s)')
s6 += row(P.svg('a1 vs Y share'), col(legend([('det3 June-15, all events', RED), ('det2 positive band', ORANGE), ('expected if only over Y pads', RED, 'dot'),
                                             ('old X-peak labels “Y pad” (pink) / “X pad” (blue)', '#7fa7d1', 'dash'), ('det3 May Y-pad fraction (right axis)', GREY, 'dash')], 21),
                                     p('Within one bin the X-peak labels split the same events into +0.26 and −0.48: the asymmetry moves the X peak and with it the label.', 23, MUT),
                                     p('So class-split fingerprints measure size × migration, not where the effect sits.', 23, MUT), w=440, gap=24), gap=30)
yv = ''
if YV:
    m = YV['det3 June-15']
    yv = (f'<p>The Y view is flat in both label classes too (main session, <code>yview/yview_class.py</code>): det3 June-15 a1Y −0.002 / −0.003 ± 0.004; '
          f'Y-pad-label map X {m["Ypad"]["map_x"]["chi2"]:.0f}/{m["Ypad"]["map_x"]["ndf"]} vs Y 30/15.</p>')
sl('pads', s6, '<p>Main session, <code>scratchpad/2026-10-21/m3class/</code> (<code>ysh_class.py</code>, <code>ysh_class_sim.py</code>, <code>m3_align.py</code>). The old label '
   '(p_x + p_y) % 2 uses the X peak; when n+1 outgrows n the peak moves, the label flips, and the event lands in the other class with a large negative a1. '
   'M3 can’t give an X-independent column (pointing 0.9–2.0 mm at the detector vs 0.78 mm pads). The Y share s can: on det3 May the Y-pad fraction falls '
   '0.80 (s &lt; 0.3) → 0.35 (1.25–1.5) → 0.20 (&gt; 1.5).</p>'
   '<p>det3 June-15 over all events: +0.17 / +0.19 at s &lt; 0.3, +0.11 … +0.14 at 0.5–1.0, +0.15 / +0.18 at s &gt; 1.25. An effect only over the Y pads would fall 3–4×. '
   'det2’s positive band falls mildly at s &gt; 1.1 (+0.10 / +0.03 ± 0.04, few events).</p>'
   '<p>Consequences: anything tied to the Y pads (vias, Y-pad copper, the bus beside the Y pads, Y→X transfer) is excluded; the "lives at the Y pads" reading (route 11) '
   'is overturned; every model score built on A_Y / A_X by X-peak class (routes 12, 13, 15–17, 26 and the floating-strip rejection, route 9) tested size convolved with '
   'migration, not location. The sim arms are also flat in s (±0.03 per bin) and reproduce the label artefact, so this does not discriminate between them.</p>' + yv,
   short='Both pads')

# ------------------------------------------------------------------------------------------------ 7 charge moved
cats = [('det3 June-15', 'Ypad', 'det3 J15, Y-pad'), ('det3 May', 'Ypad', 'det3 May, Y-pad'), ('det3 June-15', 'Xpad', 'det3 J15, X-pad'), ('det3 May', 'Xpad', 'det3 May, X-pad')]
P = Plot(980, 580, x=(-0.6, 3.6), y=(0, 0.26), ylabel='neighbour charge / ΣY')
P.xticks([(i, c[2]) for i, c in enumerate(cats)]); P.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25)])
for i, (b, c, lab) in enumerate(cats):
    q = PROBE[b][c]
    for dx, k, col_, nm in ((-0.16, 'Xm1', BLUE, 'n−1'), (0.16, 'Xp1', RED, 'n+1')):
        P.vbar(i + dx, q[k], 60, col_, tip=f'{lab}: X {nm} integral / ΣY = {q[k]:.3f}  (n = {q["n"]})', base=0, label=f'{q[k]:.2f}')
P.text(-0.4, 0.245, 'n+1 gains what n−1 loses', 25, INK, weight=600)
P2 = Plot(640, 580, x=(-0.6, 2.6), y=(0.8, 1.05), ylabel='ΣX / ΣY', title='Total X charge vs Y')
bl = ['det3 June-15', 'det3 May', 'det6']
P2.xticks([(i, b.replace('det3 June-15', 'det3 J15')) for i, b in enumerate(bl)]); P2.yticks([(v, f'{v:.2f}') for v in (0.8, 0.85, 0.9, 0.95, 1.0, 1.05)])
for i, b in enumerate(bl):
    for dx, c, op in ((-0.16, 'Ypad', 1.0), (0.16, 'Xpad', 0.4)):
        v, e = PROBE[b][c]['SXSY']
        P2.vbar(i + dx, v, 56, COL.get(b, GOLD), tip=f'{b}, {c}: ΣX/ΣY = {v:.3f} ± {e:.3f}', base=0.8, opacity=op); ebar(P2, i + dx, v, e, INK, 5)
yp, xp = PROBE['det3 June-15']['Ypad']['SXSY'][0], PROBE['det3 June-15']['Xpad']['SXSY'][0]
P2.raw(arrow(P2.X(-0.16), P2.Y(yp), P2.X(-0.16), P2.Y(yp + 0.10) - 4, GREY, 2.5))
P2.text(-0.5, yp + 0.115, 'Y→X transfer would add +0.10', 21, GREY)
s7 = sd.title('Charge is moved within X: not added, not taken from Y',
              'Left: X n−1 and n+1 integrals relative to ΣY; right: ΣX/ΣY by class (solid Y-pad, pale X-pad)')
s7 += row(P.svg('bookkeeping'), P2.svg('sum x over sum y'), gap=44)
s7 += legend([('n−1', BLUE, 'box'), ('n+1', RED, 'box')], 22)
sl('moved', s7, '<p>Integrals of positive samples per channel (<code>r5/probe_charge.py</code>). On det3 June-15 Y-pad events n+1 carries 0.216 of ΣY and n−1 only 0.116; '
   'det3 May has 0.160 / 0.163. A plain copy x(n+1) += a·x(n) would raise n+1 without lowering n−1 (and raise the central), so this is not a copy (route 10).</p>'
   '<p>Route 14 (Y→X transfer: a floating or cracked-via Y pad, Y-line crosstalk, bus pickup) was first checked against det2, which carries the field and has a '
   'different gain ratio (ΣX/ΣY spans 0.86–1.01 across boards). Redone within one board: a transfer of ≈ the n+1 − n−1 excess (0.10 of ΣY) would raise the '
   'Y-pad class above the X-pad class. Instead it is lower (0.907 vs 0.925). Excluded without any reference board.</p>', short='Moved')

# ------------------------------------------------------------------------------------------------ 8 on the detector
def parse_noise():
    out = {}
    for ln in open(os.path.join(S19, 'noise_corr_inject.log')):
        m = re.match(r'(\S+) X injected A = ([\d.]+):.*d1 excess ([+-][\d.]+)', ln)
        if m: out.setdefault(m.group(1), []).append((float(m.group(2)), float(m.group(3))))
    return out
NZ = parse_noise()
P = Plot(810, 560, x=(-0.02, 0.22), y=(-0.06, 0.26), xlabel='copy injected into the real noise, A', ylabel='noise correlation excess at distance 1', title='Noise is not copied')
P.xticks([(v, f'{v:.1f}') for v in (0, 0.1, 0.2)]); P.yticks([(v, f'{v:+.2f}' if v else '0') for v in (-0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25)])
for nm, c in (('g_det3', RED), ('sat_det3', '#e07a86'), ('g_det3_wknd', '#8f2433'), ('g_det2', ORANGE), ('g_det4', PURPLE)):
    if nm in NZ:
        xs, ys = zip(*NZ[nm]); P.line(xs, ys, c, tips=[f'{nm}: injected A = {a:.1f} → excess {v:+.3f}' for a, v in NZ[nm]])
P.band([-0.02, 0.22], [-0.04, -0.04], [0.04, 0.04], GREY, 0.12)
P.text(0.0, 0.225, 'det3 measures ≈ 0 (A = 0)', 24, INK, weight=600)
e3 = EDGE['e_d3jun23_x']; ein = e3['col%64==1..3']['phys']; eedge = e3['col%64==0 (col-1 on OTHER connector)']
P2 = Plot(810, 560, x=(-0.6, 3.6), y=(-0.05, 0.3), ylabel='a1 (physical neighbours)', title='The asymmetry crosses connector edges')
lab2 = [('inside a connector', ein, RED), ('across the edge', eedge['phys'], RED), ('det3 May at the edge', EDGE['e_d3may_x']['col%64==0 (col-1 on OTHER connector)']['phys'], BLUE),
        ('det2 at the edge', EDGE['e_d2_x']['col%64==0 (col-1 on OTHER connector)']['phys'], ORANGE)]
P2.xticks([(i, t) for i, t in enumerate(['inside', 'across edge', 'det3 May', 'det2'])]); P2.yticks([(v, f'{v:.1f}' if v else '0') for v in (0, 0.1, 0.2, 0.3)])
P2.hline(0, MUT, None, 1.5)
for i, (l, (v, e), c) in enumerate(lab2):
    P2.vbar(i, v, 80, c, tip=f'{l}: a1 = {v:+.3f} ± {e:.3f}', base=0); ebar(P2, i, v, e, INK)
s8 = sd.title('It is induced on the detector: the noise isn’t copied, and it follows the pads across connectors',
              'Left: channel-to-channel noise correlation vs an injected copy; right: physical-neighbour a1 at a connector’s first column, det3 June-15 6-23')
s8 += row(P.svg('noise'), P2.svg('connector edge'), gap=44)
s8 += legend([('det3 J15 runs (3 sessions)', RED), ('det2', ORANGE), ('det4', PURPLE), ('±0.04 noise floor', GREY, 'box')], 22)
sl('detector', s8, '<p>Noise test (<code>noise_corr_inject.py</code>): quiet events, common mode removed; the correlation between channels at distance 1 minus distance 2. '
   'An electronic copy (FEU, ADC mux, cable) copies noise too: injecting x(n) += A·x(n−1) into the real noise raises the excess by +0.06…0.12 at A = 0.1. '
   'det3 X measures +0.03 / −0.03 / −0.02, the same as boards without the effect.</p>'
   '<p>Connector edge (<code>raw_asym.py --edge</code>, <code>edge_test.py</code>): at a connector’s lowest column, the physical neighbour col−1 is on another connector, '
   'cable and DREAM chip. a1 there (+0.21 ± 0.03) equals the inside value; the DAQ neighbour on the same chip (127 columns away) gets nothing. On 6-28 (FEU 7) '
   'it is +0.03 ± 0.03 vs +0.09 inside, a 2σ difference.</p><p>With the two-signed field, these exclude the electronics, connectors, M1 card and any fixed '
   'copper layout (bus traces, stubs, vias: one sign everywhere).</p>', short='On the detector')

# ------------------------------------------------------------------------------------------------ 9 sessions + width
SS = AD['sessions']
P = Plot(810, 580, x=(-0.6, len(SS) - 0.4), y=(0, 0.3), ylabel='a1, centre chips 2–5', title='det3 June-15 by session')
P.xticks([(i, (s['t'] or '')[5:10].replace('-', '/')) for i, s in enumerate(SS)]); P.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)])
for i, s in enumerate(SS):
    c = GREY if 'other board' in s['run'] else (DRED if '6-19/20 longer' in s['run'] else RED)
    P.vbar(i, s['a'], 46, c, tip=f'{s["run"]}\na1 = {s["a"]:+.3f} ± {s["e"]:.3f}, n = {s["n"]}\nFEU {s["feu"] or "?"}, drift {", ".join(s["drift"]) or "?"} V', base=0)
    ebar(P, i, s['a'], s['e'], INK, 5)
i19 = [i for i, s in enumerate(SS) if '6-19/20 longer' in s['run']][0]
P.text(i19 - 0.3, 0.12, '6-19: drift-HV trouble', 21, INK)
P2 = Plot(810, 580, x=(0.2, 1.45), y=(0, 0.4), xlabel='Y share (cloud width, measured on the other view)', ylabel='a1, all matched events', title='Amplitude vs cloud width')
P2.xticks([(v, f'{v:.1f}') for v in (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)]); P2.yticks([(v, f'{v:.1f}' if v else '0') for v in (0, 0.1, 0.2, 0.3, 0.4)])
WLEG = []
WC = [RED, '#e07a86', '#8f2433', '#d9822b', '#b04a5a', '#f0a0a8']
SW = J(HERE, 'sess_width.json')
for k, (sess, pts) in enumerate(sorted(AD['width']['det3 June-15'].items())):
    lab = re.sub(r'mx17_|_overnight|_6-\d\d-26', '', sess).replace('_', ' ')
    WLEG.append((lab, WC[k % len(WC)]))
    P2.line([q['ysh'] for q in pts], [q['a1'][0] for q in pts], WC[k % len(WC)], w=2, r=3, dash='3 5',
            tips=[f'det3 J15 {sess}, within the night: Y share {q["lo"]}–{q["hi"]}, a1 = {q["a1"][0]:+.3f} ± {q["a1"][1]:.3f}, n = {q["n"]}' for q in pts])
    if sess in SW:
        w = SW[sess]
        P2.points([w['ysh_med']], [w['a1']], WC[k % len(WC)], r=13,
                  tips=[f'{sess}: whole night a1 = {w["a1"]:+.3f} ± {w["a1e"]:.3f}, median Y share {w["ysh_med"]:.2f}, n = {w["n"]}'])
ss = sorted((SW[k]['ysh_med'], SW[k]['a1']) for k in SW)
P2.line([q[0] for q in ss], [q[1] for q in ss], INK, w=2.5, markers=False, tip='night-to-night trend (whole-night a1 vs median Y share)')
P2.text(0.62, 0.375, 'Between nights: narrower, bigger', 24, INK, weight=600)
s9 = sd.title('The size changes between nights and tracks each night’s cloud width, with one night off the curve',
              'Left: det3 June-15 per session (board centre); right: whole-night a1 vs median Y share (big dots), and within each night (dotted)')
s9 += row(P.svg('sessions'), P2.svg('width'), gap=44)
s9 += legend([('det3 J15 session', RED, 'box'), ('6-19/20 (FEU 7)', DRED, 'box'), ('other board on FEU 3, 6-19', GREY, 'box')] + WLEG + [('night-to-night trend', INK)], 20)
sl('sessions', s9, '<p>Sessions: 6-18 long run and 6-19/20 longer run are EOS extractions (FEU 7, X only); 6-22/23 and 6-23/24 FEU 3 (drift 1000 / 600 V); 6-25 drift 500 V; '
   '6-27 drift scan (1100 → 100 V: a1 rose 0.05 → 0.16, back to the 1000 V value right after, i.e. reversible); 6-27/28 FEU 7. The n_TOF beam run_79 (July 26, '
   'Ar/iso 90/10) gives 0.05–0.07 with the older whole-board method.</p>'
   '<p>A fixed map × one scale per session fits det3 (0.96 / 1.57 / 0.44; χ² 47/60, 126/61, 47/59), and det2’s two runs differ by ×2.5 with the same bands '
   '(median Y share 0.35 vs 1.07).</p>'
   '<p><b>Within a night the trend is not there</b> (dotted lines, all events): the per-event Y share mostly says which pad type the charge sits on, '
   'and a1 is the same over both. The older "one curve" (a1 of Y-pad-labelled events vs Y share, all sessions on one curve) used the X-peak labels, whose migration '
   'grows as the Y share falls, so it overstated the per-event width dependence. What survives is night-to-night: median Y share 0.52 / 0.77 / 0.93 with a1 '
   '+0.18 / +0.14 / +0.10 (<code>sess_width.py</code>), plus the drift scan and det2’s two runs.</p>'
   '<p><b>6-19/20 is the exception</b>: 0.03 where its short pulses (narrow clouds) would predict ~0.13. It is also the only det3 night with drift-current '
   'spikes (4–9 µA) and a sagging drift HV (890–960 V), and the HV channels were swapped between the two boards relative to 6-17. Under "fixed map × width" '
   'it is a falsifier unless it is a different board or a different board state. Its raw_daq_data (HV scan + long run) can still be decoded.</p>',
   short='Sessions')

# ------------------------------------------------------------------------------------------------ 10 timing by width
P = Plot(1120, 620, x=(250, 800), y=(-130, 100), xlabel='central pulse width (samples above half max × 60 ns)', ylabel='excess delay vs central (ns)')
P.xticks([(v, str(v)) for v in (300, 400, 500, 600, 700, 800)]); P.yticks([(v, f'{v:+d}' if v else '0') for v in (-120, -80, -40, 0, 40, 80)])
P.hline(0, MUT, None, 1.5)
for key, src, c, dash, lab in (('data det3 June-15 X', TW, RED, None, 'data, det3 June-15'), ('sim xo0 |delta|>=150', TW, BLUE, '8 6', 'sim, isotropic strips'),
                               ('sim glax10 |delta|>=150', TWA, PURPLE, '8 6', 'sim, anisotropic σxx = 10 σyy'),
                               ('sim glax10w690f1 |delta|>=150', TWA, GREEN, '8 6', 'sim, anisotropic + 700 µm + funnel')):
    pts = src.get(key)
    if not pts: continue
    P.line([q['wid'] * 60 for q in pts], [q['delay'] for q in pts], c, dash=dash,
           tips=[f'{lab}, width quartile {q["q"]}: {q["wid"]:.1f} samples, a = {q["a"]:+.3f}, delay {q["delay"]:+.0f} ns, excess peak {q["epk"]:+.0f} ns, n = {q["n"]}' for q in pts])
P.text(270, 85, 'Short pulses: excess late. Long: early', 25, INK, weight=600)
s10 = sd.title('The excess is late on short pulses and early on long ones', 'Delay of the n+1 − n−1 excess relative to the central, by central-pulse-width quartile (fit D = a·W₀ + b·W₀′, delay = −b/a)')
s10 += row(P.svg('timing'), col(legend([('data, det3 June-15', RED), ('sim, isotropic strips', BLUE, 'dash'), ('sim, anisotropic', PURPLE, 'dash'),
                                        ('sim, anisotropic + 700 µm + funnel', GREEN, 'dash')], 22),
                                 p('The all-angle average looks prompt (−10 ns); the two ends cancel.', 24, MUT),
                                 p('The sim’s pulses are longer: its shortest quartile is the data’s q2, so there is no like-for-like comparison for vertical tracks.', 22, MUT), w=500, gap=26), gap=44)
sl('timing', s10, '<p><code>r5/timing_width.py</code>: the cloud_basics §42 estimator (nTof_x17) on the same events. Pulse width mixes track angle (inclined → longer), '
   'drift speed and session, so the quartiles are not pure angle bins. Route 21’s "prompt" and §42’s "delayed 114 ns" are the same data under different angle '
   'mixes (route 39).</p><p>Only the anisotropic-strip sim (σxx = 10 σyy, one tuned parameter) reproduces the sign swing; its Y view was not computed (the hybrids keep the '
   'production isotropic Y kernels).</p>', short='Timing')

# ------------------------------------------------------------------------------------------------ 11 what's excluded
def chip(t, c):
    return f'<span style="display:inline-block;padding:2px 12px;border-radius:8px;background:{c};color:white;font-weight:600">{t}</span>'
SOL, MOD, WK = chip('solid', GREEN), chip('moderate', GOLD), chip('weak', GREY)
ex = [('Code / clustering / maps', 'independent raw-ADC extractor; 0.50 on every other board and view', SOL),
      ('Electronics after the input (FEU, mux, cable)', 'noise not copied (det3’s own channels); crosses connector edges; two signs', SOL),
      ('Connector, M1 card, any linear channel mixing', 'smooth through all 7 edges; mixing fits χ² 410–750/28 (det3 May, det1 hosts)', SOL),
      ('Copper: bus traces, stubs, vias', 'any fixed layout has one sign; the field has two', SOL),
      ('A plain copy n → n+1', 'n−1 depleted, not raised; ΣX not raised', SOL),
      ('Y → X transfer (floating / cracked Y pad, crosstalk)', 'ΣX/ΣY lower in the Y-pad class (within board, no reference)', SOL),
      ('Anything tied to the Y pads (vias, Y-pad copper, bus beside them)', 'a1 is the same over X and Y pads (X-independent label)', SOL),
      ('Track angle', 'maps hold for vertical tracks; an x tilt can’t make y bands', SOL),
      ('An 800 µm strip beat (31.2 mm)', 'per-column periodogram, permutation null, every board', SOL),
      ('Uniformly low strip resistivity (2, 0.5 MΩ/sq)', 'sim keeps 0.12 / 0.001 of the Y events vs 0.9 in data', MOD),
      ('Gap conduction as the local cause', 'uniform over each board while the field changes sign (next slide)', MOD),
      ('Gas flow pattern', 'flat inlet → exhaust, X only', MOD),
      ('Extra sensitivity regions, skip transfer, static relocation', 'fast static model, which misses the full sim’s own baseline', WK)]
s11 = sd.title('Excluded, and what each exclusion rests on', 'Re-graded in the audit; "weak" means the test can’t carry the verdict')
s11 += table(['hypothesis', 'what kills it', 'how solid'], [[a, b, c] for a, b, c in ex], 24, widths=[600, 860, 160], align=['left', 'left', 'center'])
sl('excluded', s11, '<p>Route numbers (FINDINGS.md): code 1–2, electronics 4, connectors/mixing 5, copper 8 + 22, copy 10 + 20, Y→X 14, angle 6, beat 3 + 24, ρ_s 36, '
   'gap conduction 29 + 34 (+ this audit), gas 7, fast-model routes 15, 16, 31.</p><p>"Weak": routes 15–17 and 31 were scored with <code>fwd_region.py</code> / '
   '<code>esl_fast.py</code>, a static model that gives 0 where the full sim gives ±0.05. It is not a calibrator, so its exclusions are not evidence either way. In '
   'particular the skip-like transfer n−1 → n+1 (route 16) is close to what the data’s x-conditional shows (n+1 takes a share of what n−1 should get) and deserves '
   'a proper test.</p>', short='Excluded')

# ------------------------------------------------------------------------------------------------ 12 controls
ok, part, bad = chip('holds', GREEN), chip('weakened', GOLD), chip('not evidence', GREY)
ct = [('Y→X transfer: ΣX/ΣY vs det2 (route 14)', 'det2 carries the field; gains differ across boards', 'redone within det3: still excluded', ok),
      ('Moved, not added: integrals vs det2 (route 20)', 'det2 board average', 'det3 May and within-board: holds', ok),
      ('Same FEU clean on det6 (route 4)', 'det6 carries a field', 'det3 May (same FEU 3) and the noise test suffice', ok),
      ('Parity injection on det2 (route 11)', 'det2 host', 'det1, det3 May hosts give the same migration', ok),
      ('Pad-class fingerprints (A_Y, A_X by X-peak parity) used to score models', 'labels use the X peak, which the effect moves', 'X-independent label: same over both pads; void as location tests', bad),
      ('Mixed responses scored as det3 − det2 (route 13)', 'det2’s average carries class structure', 'marginal-a1 argument survives; class scores don’t', part),
      ('"Prompt" excess (route 21)', 'all-angle average', 'late / early by pulse width (route 39)', part),
      ('v5a / v4c tuned on the det3 "+ side" as copy-free', 'that side is depleted, not clean', 'retarget on det3 May / det1', part),
      ('Bus-side regions, skip transfer (routes 15–17)', 'unvalidated fast model', 'untested in effect', bad)]
s12 = sd.title('Controls that were not clean, and what survives once they are redone', 'det2 and det6 carry the field; board averages and fast models hid it')
s12 += table(['conclusion', 'problem with the control', 'redone / status', ''], [list(r) for r in ct], 23, widths=[560, 470, 500, 150], align=['left', 'left', 'left', 'center'])
sl('controls', s12, '<p>The pattern: every conclusion that matters for the exclusions survives a within-board or det3-May redo. What doesn’t survive are fine-grained '
   'model scores that subtracted det2 as a baseline and the fast-model exclusions.</p><p>One leftover from the copy model: NB_ASYM §2 used the det3 X side with share '
   '0.30 as "copy-free" and concluded the sims overshare (v4c 0.39). Under "charge moved" that side is depleted; the symmetric value is ~0.36, which v5a (0.35) '
   'already matches. Sim targets for the X ±1 share and width should come from det3 May / det1. <code>uncopy.py</code> assumes the copy model and should be retired '
   'in favour of <code>wft.xasym</code>.</p>', short='Controls')

# ------------------------------------------------------------------------------------------------ 13 gap conduction
def latepanel(key, lab, xl):
    d = LATE[key]; rr = d['rows']
    P = Plot(810, 560, x=(-0.3, 7.3), y=(-0.3, 0.3), xlabel=xl, ylabel='a1 (all events)', title=lab)
    P.xticks([(i, str(i)) for i in range(8)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)])
    P.y2(0.3 / 40, [(0, '0'), (10, '10'), (20, '20'), (30, '30'), (40, '40')], 'X ±2 late content (×10⁻³)')
    P.hline(0, MUT, None, 1.5)
    P.line([r['band'] for r in rr], [r['a1'] for r in rr], RED, tips=[f'{lab}, band {r["band"]}: a1 = {r["a1"]:+.3f} ± {r["a1e"]:.3f}, n = {r["n"]}' for r in rr])
    P.line([r['band'] for r in rr], [P.Y2(r['late']) - 0.3 + 0 for r in rr], GREEN,
           tips=[f'{lab}, band {r["band"]}: late ±2 = {r["late"]:.1f} ± {r["latee"]:.1f} ×10⁻³' for r in rr])
    return P
# the right axis maps 0..40 onto -0.3..+0.3: value v -> -0.3 + v * 0.6/40
def latepanel2(key, lab, xl):
    d = LATE[key]; rr = d['rows']
    P = Plot(810, 560, x=(-0.3, 7.3), y=(-0.3, 0.3), xlabel=xl, ylabel='a1 (all events)', title=lab, margin=(24, 110, 92, 104))
    P.xticks([(i, str(i)) for i in range(8)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)])
    m = lambda v: -0.3 + v * 0.6 / 40
    P._y2f = 1; P._y2t = [(m(v), str(v)) for v in (0, 10, 20, 30, 40)]; P.ylabel2 = 'X ±2 late content (×10⁻³)'
    P.hline(0, MUT, None, 1.5)
    P.line([r['band'] for r in rr], [r['a1'] for r in rr], RED, tips=[f'{lab}, band {r["band"]}: a1 = {r["a1"]:+.3f} ± {r["a1e"]:.3f}, n = {r["n"]}' for r in rr])
    P.line([r['band'] for r in rr], [m(r['late']) for r in rr], GREEN, tips=[f'{lab}, band {r["band"]}: late ±2 = {r["late"]:.1f} ± {r["latee"]:.1f} ×10⁻³' for r in rr])
    return P
PA = latepanel2('det2', 'det2, by 64-row band', 'band'); PB = latepanel2('det3 June-15', 'det3 June-15, by chip', 'band')
PB.text(0, 0.26, 'Leak flat, field not', 25, INK, weight=600)
so = J(S20, 'strip_orient2.json')
s13 = sd.title('The June boards conduct across their strip gaps, but uniformly: the field doesn’t follow',
               'a1 (red, left axis) and the across-strip late charge (green, right axis) per band')
s13 += row(PA.svg('det2 leak'), PB.svg('det3 leak'), gap=40)
s13 += legend([('a1', RED), ('X ±2 late content, 600–960 ns after the central', GREEN)], 22)
sl('leak', s13, f'<p>The late ±2 content in the across-strip view (X) is the only physical June-vs-May difference found in the data: June boards 12–31 ×10⁻³ '
   f'(det2 {1e3 * so["det2"]["X"]["n2_late"]:.0f}, det3 June-15 {1e3 * so["det3 June-15"]["X"]["n2_late"]:.0f}, det6 {1e3 * so["det6"]["X"]["n2_late"]:.0f}), '
   f'May and earlier 1–9 (det3 May {1e3 * so["det3 May"]["X"]["n2_late"]:.1f}, det4 May {1e3 * so["det4 May (5-6 bulk)"]["X"]["n2_late"]:.1f}). The sim reproduces it with a gap '
   'conductance of 5–10 % of the strips’ (route 34).</p><p>New in the audit (<code>late_vs_a1.py</code>): inside det2 it stays 8–13 in every band while a1 runs −0.21 → '
   '+0.16 → −0.11; on det3 27–32 while a1 runs 0.09 → 0.22. Across boards det6 has the most leak and a weak field, det2 the least leak and a strong field. So conduction '
   'is a June-build marker, not what sets the local size or sign. This now rests on data, not only on the sim.</p>', short='Gap leak')

# ------------------------------------------------------------------------------------------------ 14 registration: what it predicts
def reg_svg():
    W, H = 760, 560; o = []; P0 = 200; x0, y0 = 60, 330
    for k in range(-1, 2):
        cx = x0 + (k + 1.5) * P0
        o.append(f'<rect x="{cx - 0.45 * P0:.0f}" y="{y0}" width="{0.9 * P0:.0f}" height="46" fill="{BLUE}" fill-opacity="0.8" rx="4"/>')
        o.append(T(cx, y0 + 84, ['n−1', 'n', 'n+1'][k + 1], 26, INK, weight=600))
    d = 60
    for k in range(-1, 2):
        cx = x0 + (k + 1.5) * P0 + d
        o.append(f'<rect x="{cx - 0.35 * P0:.0f}" y="{y0 - 52}" width="{0.7 * P0:.0f}" height="26" fill="#8a8f99" fill-opacity="0.55" rx="4"/>')
    cx = x0 + 1.5 * P0
    xs = np.linspace(cx - 2 * P0, cx + 2 * P0, 120); ys = y0 - 70 - 170 * np.exp(-0.5 * ((xs - cx + 20) / 95) ** 2)
    o.append(sd.poly(list(xs), list(ys), RED, 3, dash='8 6'))
    ys2 = y0 - 70 - 200 * np.exp(-0.5 * ((xs - cx - d) / 70) ** 2)
    o.append(sd.poly(list(xs), list(ys2), RED, 4))
    o.append(arrow(cx - 10, y0 - 250, cx + d - 4, y0 - 250, INK, 2.5))
    o.append(T(cx + 120, y0 - 262, 'pulled to the strip centre', 22, INK, 'start'))
    o.append(line(cx, y0 - 20, cx, y0 + 4, INK, 2)); o.append(line(cx + d, y0 - 20, cx + d, y0 + 4, INK, 2))
    o.append(T(cx + d / 2, y0 + 120, 'δ', 28, INK, weight=600))
    o.append(T(W / 2, y0 + 170, 'odd effect: n+1 gains (a1 ∝ δ)', 23, RED))
    o.append(T(W / 2, y0 + 204, 'even effect: the cloud is concentrated → s1 changes', 23, BLUE))
    return svg(W, H, ''.join(o), 'registration cartoon')
sd_rows = [r for r in SIMD if r]
s0 = [r for r in sd_rows if abs(r['d']) < 30]; s0v = np.mean([r['s1'] for r in s0])
P = Plot(860, 600, x=(-360, 360), y=(-0.2, 0.55), xlabel='true strip offset δ in the sim (µm)', ylabel='change vs δ = 0')
P.xticks([(v, f'{v:+d}' if v else '0') for v in (-300, -150, 0, 150, 300)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)])
P.hline(0, MUT, None, 1.5)
P.line([r['d'] for r in sd_rows], [r['a1'] for r in sd_rows], RED, tips=[f'sim δ = {r["d"]:+.0f} µm: a1 = {r["a1"]:+.3f} ± {r["e"]:.3f}, n = {r["n"]}' for r in sd_rows])
P.line([r['d'] for r in sd_rows], [r['s1'] - s0v for r in sd_rows], BLUE, tips=[f'sim δ = {r["d"]:+.0f} µm: s1 = {r["s1"]:.3f} (s1 − s1(0) = {r["s1"] - s0v:+.3f})' for r in sd_rows])
P.text(-250, 0.48, 'Sharing moves as much as a1', 25, INK, weight=600)
s14 = sd.title('The leading candidate, strip registration, predicts an even partner to a1',
               'Left: strips offset by δ pull the cloud; right: the full sim (strips at 800 µm over 780 µm pads) per δ bin')
s14 += row(reg_svg(), col(P.svg('sim delta'), legend([('a1 (odd in δ)', RED), ('s1 − s1(δ=0) (even in δ)', BLUE)], 22), gap=8), gap=40)
sl('registration', s14, '<p>Registration: if the strip centres sit δ off the pad columns and charge evens out across a strip quickly, narrow clouds are pulled toward the strip '
   'centre. a1 ∝ sin(2πδ/P) and the field geometry follows from the print’s scale and rotation errors (det3 ~0.07 % scale + ~0.25 mrad; det2 ~2 mrad '
   'rotation → bands ~P/θ ≈ 390 mm apart). It explains X only, charge moved, pad-pitch locking and the width dependence.</p>'
   '<p>Requires a strip pitch commensurate with 780 µm; the quoted pitch is 800 µm ("borrowed from another detector"). With 800 µm the sim predicts a 31.2 mm beat '
   'of amplitude 0.12 on every board, which no board shows (≤ 0.03).</p>'
   '<p>Sim table: <code>r5/sim_delta.py</code> (36k X-only events, xo0 th00 + thy*, seeds 5–16). s1 is 0.63 at δ ≈ 0 and 0.77–1.14 near the gaps; '
   'a fit over |δ| ≤ 275 µm gives s1 = 0.64 + 2.2 |a1|.</p>', short='Registration')

# ------------------------------------------------------------------------------------------------ 15 sim levers
LEV = [('isotropic, 540 µm strips (grid control)', 1.0, 1.0, '✓', 0.760), ('420 µm strips', 0.769, 0.687, '✓', 0.738), ('leaky gaps f = 0.1', 0.902, 1.041, '✓', 0.828),
       ('anisotropic σxx = 10 σyy', 1.063, 1.007, '? (Y not computed)', 0.761), ('gap funnel', 1.397, 1.535, '✓', 0.765), ('700 µm strips', 1.448, 1.535, '✓', 0.808),
       ('700 µm + funnel', 1.583, 1.776, '✓', 0.812), ('anisotropic + 700 µm + funnel', 1.792, 2.217, '? (Y not computed)', 0.824),
       ('2 MΩ/sq + 700 µm + funnel', 2.263, 2.198, '✗ (Y/X 0.12)', 0.869), ('0.5 MΩ/sq + 700 µm + funnel', 2.382, 2.115, '✗ (Y/X 0.001)', 0.893)]
PO = {}
for f in ('paired_odd_aniso.json', 'paired_odd_funnel.json', 'paired_odd_gapleak.json', 'paired_odd_width.json', 'paired_odd_w690f1.json',
          'paired_odd_r2_vs_g0.json', 'paired_odd_r05f1_vs_g0.json'):
    PO.update(J(R5, f))
KEYS = [None, 'glw410', 'glg0.1', 'glax10', 'fu1', 'glw690', 'glw690f1', 'glax10w690f1', 'glr2w690f1', 'glr05w690f1']
P = Plot(1180, 640, x=(0, 5.2), y=(-0.5, len(LEV) - 0.5), xlabel='odd amplitude relative to the grid control (paired arms, same angles and seeds)', margin=(20, 40, 92, 440))
P.xticks([(v, f'×{v:g}') for v in (0, 1, 2, 3, 4, 5)])
P.raw(f'<rect x="{P.X(2.5):.1f}" y="{P.y0}" width="{P.X(5.0) - P.X(2.5):.1f}" height="{P.ph}" fill="{RED}" fill-opacity="0.10"/>', back=True)
P.text(2.6, len(LEV) - 0.9, 'data need ×2.5–5', 23, RED, weight=600)
for i, ((lab, r, rn, yok, s1), k) in enumerate(zip(LEV, KEYS)):
    yy = len(LEV) - 1 - i
    if k and k in PO: r, rn, s1 = PO[k]['all']['ratio'], PO[k]['nar']['ratio'], PO[k]['s1']['var']
    c = RED if yok.startswith('✗') else (GOLD if yok.startswith('?') else BLUE)
    P.raw(f'<rect x="{P.X(0):.1f}" y="{P.Y(yy) - 16:.1f}" width="{P.X(r) - P.X(0):.1f}" height="32" fill="{c}" rx="4"{sd.tipattr(f"{lab}: ×{r:.2f} all, ×{rn:.2f} narrow clouds; board s1 {s1:.3f}; Y view {yok}")}/>')
    P.raw(T(P.x0 - 14, P.Y(yy) + 8, lab, 22, INK, 'end'))
    P.raw(T(P.X(r) + 10, P.Y(yy) + 8, f'×{r:.2f}', 21, INK, 'start'))
s15 = sd.title('Every size lever stays inside the same family, and none reaches the data with the Y view intact',
               'Full sim (v5a recipe), paired arms vs the grid control; colour = Y view kept (blue), not computed (gold), destroyed (red)')
s15 += row(P.svg('levers'), col(legend([('Y view kept', BLUE, 'box'), ('Y view not computed', GOLD, 'box'), ('Y view destroyed', RED, 'box')], 22),
                                 p('Each lever also raises the board s1 (0.76 → 0.82–0.89): it concentrates or spreads charge, it doesn’t only shift it.', 24, MUT), w=440, gap=24), gap=40)
sl('levers', s15, '<p>Odd amplitude O = ⟨sgn δ · a1⟩ over |δ| ≥ 50 µm (<code>r5/paired_odd.py</code>), each arm paired event-for-event with the same (angle, seed) control. '
   'Grid control glg0 = ×0.95 of the production arm. Leaky gaps: a kernel-level ×1.9 became ×0.9 in the full sim (the gain arrives late). Low ρ_s: ×2.1–2.4 and the data’s '
   'timing trend, but the charge runs along the strips over many rows and the Y view loses its compact clusters.</p>'
   '<p>Anisotropic strips (A = 10) give the data’s timing swing; the Y-box solve (64 families) that would check the Y view has not been run.</p>', short='Levers')

# ------------------------------------------------------------------------------------------------ 16 the even-term test (det2 bands)
rr = [r for r in BANDS['det2|y'] if r['cls'] == 'all']
P = Plot(1000, 620, x=(-0.3, 7.3), y=(-0.4, 0.6), xlabel='det2 row band (64 rows each)', ylabel='per band, all matched events')
P.xticks([(i, str(i)) for i in range(8)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.4, -0.2, 0, 0.2, 0.4, 0.6)])
P.hline(0, MUT, None, 1.5)
ks = SHARE['sim']['k_abs']; base = float(np.mean([r['s1x'][0] - r['s1y'][0] for r in rr if abs(r['a1'][0]) < 0.05]))
P.line([r['band'] for r in rr], [r['a1'][0] for r in rr], RED, tips=[f'det2 band {r["band"]}: a1 = {r["a1"][0]:+.3f} ± {r["a1"][1]:.3f}, n = {r["n"]}' for r in rr])
P.line([r['band'] for r in rr], [r['s1x'][0] - r['s1y'][0] for r in rr], BLUE,
       tips=[f'det2 band {r["band"]}: s1X = {r["s1x"][0]:.3f} ± {r["s1x"][1]:.3f}, s1Y = {r["s1y"][0]:.3f}  →  s1X − s1Y = {r["s1x"][0] - r["s1y"][0]:+.3f}' for r in rr])
P.line([r['band'] for r in rr], [base + ks * abs(r['a1'][0]) for r in rr], BLUE, dash='8 6', markers=False,
       tip=f'Registration as simulated: s1X − s1Y = baseline {base:+.3f} + {ks:.1f} |a1| (slope from r5/sim_delta.json)')
P.text(0.0, 0.54, 'a1 swings, sharing doesn’t', 25, INK, weight=600)
sh = SHARE['det2']['ypad']; sa = SHARE['det2']['all']
cats = [('sim registration', ks, 0.15, BLUE), ('det2, all cells', sa['k_abs'], sa['k_abs_e'], ORANGE), ('det2, Y-pad-label cells', sh['k_abs'], sh['k_abs_e'], ORANGE),
        ('det3 May, all cells', SHARE['det3 May']['all']['k_abs'], SHARE['det3 May']['all']['k_abs_e'], GREY)]
P2 = Plot(620, 620, x=(-0.6, len(cats) - 0.4), y=(-1.4, 2.6), ylabel='Δs1X per unit |a1| (cell fit)', title='Slope over 8 × 8 cells')
P2.xticks([(i, c[0].replace(', ', ', ')) for i, c in enumerate(cats)]); P2.yticks([(v, f'{v:+.0f}' if v else '0') for v in (-1, 0, 1, 2)])
P2.xticks([(i, ['sim', 'det2', 'det2 Y-lab', 'det3 May'][i]) for i in range(len(cats))])
P2.y0 += 0
P2.hline(0, MUT, None, 1.5)
for i, (lab, v, e, c) in enumerate(cats):
    P2.vbar(i, v, 70, c, tip=f'{lab}: slope {v:+.2f} ± {e:.2f}' + (' (from the δ table, no cell fit)' if i == 0 else ''), base=0); ebar(P2, i, v, e, INK)
s16 = sd.title('The data fail that prediction: on det2 the sharing ignores a field of ±0.2',
               'Left: det2 per row band: a1 (red), s1X − s1Y (blue; Y share removes the cloud width), registration’s prediction (dashed). Right: fitted slopes')
s16 += row(P.svg('even term'), P2.svg('slopes'), gap=44)
s16 += legend([('a1', RED), ('s1X − s1Y, data', BLUE), ('s1X − s1Y predicted by the sim’s registration', BLUE, 'dash')], 22)
sl('even', s16, f'<p>New in the audit (<code>bands.py</code>, <code>share_vs_a1.py</code>). det2 is the test board: its field takes both signs, so an even term (∝ |a1|) is '
   f'separable from an odd one (∝ a1). Per cell (8 × 8, ≥ 30 events) the fit s1X = c0 + c1·s1Y + k|a1| + k′a1 gives, for det2 (all events), '
   f'k = {sa["k_abs"]:+.2f} ± {sa["k_abs_e"]:.2f} (Y-pad-labelled events {sh["k_abs"]:+.2f} ± {sh["k_abs_e"]:.2f}), against {ks:.1f} for the sim’s registration. det3 June-15’s field is '
   f'mostly one-signed there, so k and k′ are degenerate (their sum ≈ −0.3); the quiet boards give k ≈ 0 at |a1| rms 0.03.</p>'
   '<p>The earlier quadrature test (<code>r5/quad_test.py</code>, 2026-10-20) looked for the same thing as a cos φ term on the fitted map phase and found −0.007 ± 0.009 on det3. It was '
   'listed "against" and then dropped from the summary. Together with route 20 ("n+1 gains what n−1 loses, the ±1 sum unchanged") it says the data move charge sideways '
   'without concentrating it.</p><p>Caveat: s1Y is a proxy for the cloud width; det2’s bands also differ in width (s1Y 0.33–0.97), so a residual width term could hide a '
   'small even term. It cannot hide one of the predicted size (+0.3 to +0.7 at |a1| = 0.15–0.33).</p>', short='Even term')

# ------------------------------------------------------------------------------------------------ 17 any attractor
P = Plot(1060, 640, x=(0, 0.25), y=(-0.1, 0.4), xlabel='|a1|', ylabel='Δs1 that comes with it')
P.xticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.1, 0, 0.1, 0.2, 0.3, 0.4)])
P.hline(0, MUT, None, 1.5)
TC = {(250., 50.): '#2c5f9e', (250., 150.): BLUE, (250., 300.): '#7fa7d1', (450., 50.): '#5b3a8c', (450., 150.): PURPLE, (450., 300.): '#b79be0'}
for t in TOY:
    k = int(np.argmax([q[1] for q in t['pts']])); pts = [q for q in t['pts'][:k + 1] if q[1] >= 0]
    P.line([q[1] for q in pts], [q[2] for q in pts], TC[(t['sigma'], t['sig_ind'])], w=3, r=4,
           tips=[f'toy: cloud σ {t["sigma"]:.0f} µm, induction width {t["sig_ind"]:.0f} µm, 700 µm strips, δ = {q[0]:.0f} µm → a1 {q[1]:+.3f}, Δs1 {q[2]:+.3f}' for q in pts])
dx = [(abs(r['a1'][0]), (r['s1x'][0] - r['s1y'][0]) - base, r) for r in rr]
P.points([d[0] for d in dx], [d[1] for d in dx], ORANGE, r=9,
         tips=[f'det2 band {d[2]["band"]}: |a1| = {d[0]:.3f}, s1X − s1Y relative to the quiet bands = {d[1]:+.3f}' for d in dx])
P.text(0.005, 0.37, 'Any strip pull that makes a1 also makes Δs1', 25, INK, weight=600)
s17 = sd.title('It isn’t a detail of the sim: any pull toward the strips concentrates the charge',
               '1-D toy: cloud landing on 700 µm strips at the pad pitch, spread across each strip; curves for cloud σ and induction width. Orange: det2 bands')
s17 += row(P.svg('toy'), col(legend([('cloud σ 250 µm (blues)', BLUE), ('cloud σ 450 µm (purples)', PURPLE), ('det2 bands, data', ORANGE, 'dot')], 22),
                              p('Darker = narrower induction (50, 150, 300 µm).', 22, MUT),
                              p('With σ = 450 µm clouds the toy can’t exceed |a1| ≈ 0.06 at all, even with full spreading across 700 µm strips.', 24, MUT), w=500, gap=24), gap=40)
sl('toy', s17, '<p><code>toy_strips.py</code> / <code>toy_curves.py</code>: 1-D, pads of pitch 780 µm read through box ⊗ Gaussian (induction width); a Gaussian cloud at uniform '
   'positions; the part landing on a strip is spread uniformly across that strip’s width (full equalisation), gap charge stays; peak-pad selection. Δs1/a1 is 0.2–4 '
   'depending on σ, induction width and δ; at the δ needed for a1 ≥ 0.1, Δs1 ≈ +0.05 … +0.25.</p>'
   '<p><code>toy_attractor.py</code>: the most forgiving case, a fraction f = 0.2 of the cloud moved fully to a line at the attractor, gives a1 = 0.2 with Δs1 ≈ +0.05. '
   'Even that is above the det2 bands (≈ 0 ± 0.05), and it is not strip physics.</p>'
   '<p>What would escape: a response that is skewed rather than attracting. If X pad n’s weighting field leans toward n+1 (for example steered by a partly '
   'conducting, offset strip), the centroid shifts with no first-order change of the sum, and the pad type under the charge doesn’t matter. Only a static version '
   'of this was ever modelled (route 9), and it was rejected on a pad-class criterion that is now void.</p>', short='Any attractor')

# ------------------------------------------------------------------------------------------------ 18 fingerprint matrix
Y_, N_, Q_ = chip('✓', GREEN), chip('✗', RED), chip('?', GREY)
FP = [('X view only; Y flat in both pad classes', Y_, Y_, Q_),
      ('Odd n−1 → n+1 transfer, ±1 sum unchanged', N_, '✓ (first order)', Q_),
      ('Smooth two-signed field, fixed per board', Y_, Q_, Q_),
      ('Pad-pitch locked, no 31.2 mm beat', '✓ only with a 780 µm print', Y_, Q_),
      ('Same over X pads and Y pads', Y_, Y_, Q_),
      ('Size: board mean 0.14, cells ±0.2', '✗ (≤ ×1.8 of 1/3)', Q_, Q_),
      ('Bigger on narrow-cloud nights and at low drift field', Y_, Q_, 'plausible (field)'),
      ('Late on short pulses, early on long', '✓ anisotropic only', Q_, Q_),
      ('June builds only', 'needs a print difference', Q_, 'June material'),
      ('6-19/20 off the width curve', N_, Q_, 'plausible (HV trouble)')]
s18 = sd.title('The fingerprint any explanation must reproduce', 'Rows: data facts from the previous slides. Columns: the candidates still standing')
s18 += table(['fact', 'strip registration (as simulated)', 'skewed X response (untested)', 'surface charging (untested)'], [list(r) for r in FP], 23,
             widths=[640, 380, 330, 320], align=['left', 'center', 'center', 'center'])
sl('fingerprint', s18, '<p><b>Skewed X response</b>: with uniform illumination, a1 ≠ 0 needs a response that is not mirror-symmetric about the pad. Registration gets there by '
   'attracting the charge (and pays with the even term). The other way is a weighting field that leans toward one neighbour, without moving the charge: the centroid '
   'shifts, the ±1 sum is unchanged at first order, and the pad type doesn’t matter. Offset strips acting as partly conducting electrodes in the induction are one '
   'way to get it, and a smooth two-signed lean would follow from the print registration. Only its static version was modelled (route 9, rejected on a class '
   'criterion now void, and because it grew with cloud size).</p>'
   '<p><b>Surface charging</b>: insulating surfaces (coverlay, gap kapton) on fresh June bulks charging up and adding a lateral field. It would naturally give time and '
   'field dependence, and an HV incident could reset it (6-19). It has no worked prediction for the map geometry or pad-pitch locking yet.</p>'
   '<p>These are candidates to test, not results.</p>', short='Fingerprint')

# ------------------------------------------------------------------------------------------------ 19 way forward
WF = [('1', 'Score the even term for every model arm (Δs1 vs a1 per δ; s1X − s1Y per cell in data); drop X-peak-class scores', 'archive + existing sim', 'retires or rescues the lever family'),
      ('2', 'a1 vs time inside long runs and around every HV change (6-27/28 resist 3:3 → 3:4, drift-scan steps, 6-19 spikes); beam vs bench at matched width', 'archive', 'charging vs pure width'),
      ('3', 'Decode 6-19/20 raw_daq_data and map it', 'archive', 'fixed map at low scale, or a different board/state'),
      ('4', 'Per-session, width-matched maps for det3', 'archive', 'is the geometry really fixed'),
      ('5', 'Split the existing sim arms’ a1 into attraction and induction skew (even vs odd per δ); redo route 9 with resistive offset strips in the weighting field', 'sim, small', 'whether any strip model passes the even term'),
      ('6', 'Photograph the strips (or get the mask) on det2 rows 0–63 vs 128–255, and on det3 May', 'needs the board', 'registration: opposite δ there, 780 µm pitch'),
      ('7', 'Collimated source across one pad; HV ramps and dwell; dry-gas flush', 'needs the detector', 'attractor vs displacement; charging; humidity')]
s19 = sd.title('Way forward, ranked: cheapest decisive tests first', 'Steps 1–5 need nothing new; 6 needs a photo; 7 needs the detector')
s19 += table(['#', 'test', 'needs', 'decides'], [list(r) for r in WF], 23, widths=[50, 980, 250, 380], align=['center', 'left', 'left', 'left'])
sl('forward', s19, '<p>Step 1 is the one that changes the direction of the work: until now every lever was scored on the odd amplitude and the timing only. A sim arm passes only '
   'if its even term per unit |a1| is ≲ 0.1–0.2 (det2’s Y-pad slope +0.11 ± 0.05).</p>'
   '<p>Step 6 is the most decisive physical measurement: registration predicts strip-centre offsets of opposite sign (±100–300 µm) in det2’s positive and negative bands, '
   'and a 780 µm pitch. A phone photo through a microscope eyepiece at a few points would do. Read strip width and gap too.</p>'
   '<p>Analysis does not wait for the mechanism: <code>wft.xasym</code> (nTof_x17, commit 30198e4) corrects det2 and det3 June-15 positions with a measured map × '
   'per-run scale (K = −0.77 ± 0.13 vs M3). It removes a positional systematic of up to ±70 µm and leaves the resolution unchanged.</p>', short='Way forward')

# ------------------------------------------------------------------------------------------------ 20 close
cl = (sd.kicker('Where it stands')
      + '<h2 style="font-size:58px;font-weight:600;line-height:1.1;width:1660px">The effect is well measured and the analysis can be corrected. Its cause is open: '
        'registration explains the map but not the sharing</h2>'
      + row(card(p('What this does not rule out', 30, DINK, 600)
                 + p('• Strips that skew the induced signal without concentrating the charge (registration through induction).', 26, DMUT)
                 + p('• A difference in the June print (pitch, width, fill) that the mask or a photo would show.', 26, DMUT)
                 + p('• A surface or charging effect specific to fresh June bulks.', 26, DMUT)
                 + p('• Something in the June bulk process nobody has named yet (Stephan: "nothing different").', 26, DMUT), bg='#1f2533'),
            card(p('Open decisions', 30, DINK, 600)
                 + p('• Stop tuning attractor-type levers until the even-term score is in place?', 26, DMUT)
                 + p('• Can someone photograph a June board’s strips (det2 ideally) and det3 May?', 26, DMUT)
                 + p('• Adopt wft.xasym in production (needs a re-alignment and per-run scales)?', 26, DMUT)
                 + p('• Was det3’s board or cabling different on 6-19/20?', 26, DMUT), bg='#1f2533'), gap=40))
sl('close', cl, '<p>Audit text with every route re-graded (excluded / weakened / supported / untested, flagged for contaminated controls, board averages, selection effects '
   'and unvalidated models): <code>scratchpad/2026-10-21/audit/AUDIT.md</code> in MX17_Geant. Builder: <code>scratchpad/2026-10-21/audit/make_audit_deck.py</code>.</p>',
   dark=True, short='Open')

OUTP = os.path.join(HERE, 'det3-x-asymmetry.html')
D.write(OUTP, note_meta=dict(title='det3 X neighbour asymmetry: an audit',
                             summary='Fresh-eyes summary of the week-long study. Every June board carries a smooth, two-signed, X-only field of unequal ±1 neighbours that moves charge '
                                     'sideways without changing the ±1 sharing. Electronics, connectors, copper, copies, Y→X transfer and gap conduction are excluded on '
                                     'within-board tests. Strip registration explains the map, but it predicts a sharing change the data don’t show (det2: ×20 too large). '
                                     'The deck lays out the fingerprint and a ranked plan.',
                             tags='MX17,detector,crosstalk', date='2026-10-21'))
print('wrote', OUTP)
