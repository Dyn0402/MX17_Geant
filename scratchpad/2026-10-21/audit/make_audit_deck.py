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
           'The June boards’ X view moves charge straight from one neighbour to the other. Its cause is open: '
           'a pull onto narrow resistive strips did it in a 1-D toy, but not in the full sim</h1>'
         + p('A fresh look at a week of work, updated the same evening: the observation is solid and well mapped, and many causes are firmly '
             'excluded. The size is set by the cloud width alone (not gain, not time). The sharing test first thought to exclude registration '
             'does not decide it. Narrow strips with a charged-gap funnel, which reach the data’s size in a 1-D toy, give a smaller asymmetry than '
             'the production strips in the full sim.',
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

# ------------------------------------------------------------------------------------------------ 7b shape of the transfer
DP = J(HERE, '..', 'even', 'det2_profile.json'); AS = J(HERE, '..', 'even', 'antisym.json')
ks_ = [-3, -2, -1, 1, 2, 3]
fd = {r['k']: r for r in DP['det2']['fit']}
mb = [r for r in AS if r['lo'] >= 0.6 and r['hi'] <= 1.05]
def avg_norm(getter):
    out = []
    for i, k in enumerate((-3, -2, -1, 0, 1, 2, 3)):
        if k == 0: continue
        out.append(float(np.mean([getter(r)[i] / (getter(r)[4] - getter(r)[2]) for r in mb])))
    return out
meas = avg_norm(lambda r: r['d']); meas = [0.5 * (meas[i] - meas[5 - i]) for i in range(6)]; anti = avg_norm(lambda r: r['anti']['pred']); copy = avg_norm(lambda r: r['copy']['pred'])
P = Plot(1640, 570, x=(-3.6, 3.6), y=(-0.65, 0.95), xlabel='neighbour k (DAQ order; n+1 is the excess side)', ylabel='change of α(k) per unit change of a1')
P.xticks([(k, f'{k:+d}') for k in ks_]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)]); P.hline(0, MUT, None, 1.5)
SER = [('det2 (cell slope)', [fd[k]['b_a1'] for k in ks_], [fd[k]['b_a1_e'] for k in ks_], ORANGE),
       ('det3 June-15 − May (odd part)', meas, None, RED),
       ('linear antisymmetric coupling', anti, None, GREY), ('one-sided copy', copy, None, BLUE)]
TR = J(HERE, '..', 'loop', 'toy_shape_rep.json')
SER.insert(2, (f'cloud-level pull toy ({TR["f"]:.0%}, δ {TR["delta"]:.0f} µm)', [0.0 if k in (-3, 3) else TR['per_a1'][str(k)] for k in ks_], None, PURPLE))
for j, (lab, v, e, c) in enumerate(SER):
    for i, k in enumerate(ks_):
        x = k + (j - 2) * 0.15
        P.vbar(x, v[i], 26, c, tip=f'{lab}: k = {k:+d}: {v[i]:+.3f}' + (f' ± {e[i]:.3f}' if e else ''), base=0)
        if e: ebar(P, x, v[i], e[i], INK, 4)
P.text(-3.4, 0.85, 'Data: n−1 → n+1 directly; ±2 carries 1/8 (det2) – 1/4 (det3)', 24, INK, weight=600)
s7b = sd.title('The transfer goes straight from n−1 to n+1: not a copy or a fixed coupling between channels',
               'Change of each neighbour’s in-phase content per unit a1: data (det2 cells, det3 June-15 vs May) against a strip-pull toy and two channel-level operators')
s7b += P.svg('transfer shape') + legend([(s[0], s[3], 'box') for s in SER], 20)
sl('transfer', s7b, '<p>det2 (<code>even/det2_profile.py</code>, 32 cells, a1 −0.30 … +0.26, Y share and Y ±2 as width controls): per unit a1, −0.54 / +0.47 at ∓1, −0.08 / +0.04 at ∓2, '
   '~0 at ∓3. det3 June-15 − May (<code>even/antisym.py</code>, Y-share bins 0.6–1.05, odd part): the same shape, with ±2 at ~1/4 of ±1. Operators applied to May’s mean profile and scaled to the measured Δa1: '
   'a nearest-neighbour antisymmetric coupling dq_n = (ε/2)(q_{n−1} − q_{n+1}) puts ~40 % of the ±1 change at ±2; a copy dq_n = a·q_{n−1} raises n+1 and gives n−1 nothing.</p>'
   '<p>So it isn’t a fixed linear channel operator (consistent with the failed mixing fits, route 5). It depends on where the charge sits, i.e. a position-level mechanism; a pull '
   'toward the strips is one. The May → June-15 ±2/±3 even drop (−0.05 / −0.02) is mostly a bulk difference (e.g. the photon-feedback halo): det2’s within-board even slope '
   'predicts only −0.01 at a1 = 0.14.</p>'
   f'<p>Strip-pull toy (<code>loop/toy_shape.py</code>, the attractor toy of slide "Even term"): at the data’s operating point (a1 0.10–0.20 with |Δs1| ≤ 0.03 vs no pull) the ±2 '
   f'odd part is 0.01–0.14 of the ±1 one (σ 450 µm: 0.10–0.14) with a small negative ±2 even part, like det2. Shown: {TR["f"]:.0%} of the charge pulled onto lines {TR["delta"]:.0f} µm '
   f'off the pad centres, σ {TR["sigma"]:.0f} µm clouds, 200 µm induction: a1 {TR["a1"]:+.3f}, Δs1 {TR["ds1"]:+.3f}. ±3 is not computed in the toy (drawn at 0). This toy pulls part of the whole cloud to the line nearest its centre, which a linear detector can’t do; it is a shape check, and slide “Physical pull” has the per-electron version.</p>', short='Transfer shape')


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
    c = GREY if '6-19/20 longer' in s['run'] else (DRED if 'other board' in s['run'] else RED)
    P.vbar(i, s['a'], 46, c, tip=f'{s["run"]}\na1 = {s["a"]:+.3f} ± {s["e"]:.3f}, n = {s["n"]}\nFEU {s["feu"] or "?"}, drift {", ".join(s["drift"]) or "?"} V', base=0)
    ebar(P, i, s['a'], s['e'], INK, 5)
i19 = [i for i, s in enumerate(SS) if '6-19/20 longer' in s['run']][0]
P.text(i19 - 0.45, 0.042, '= det2', 20, INK)
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
s9 = sd.title('The size changes between nights and tracks each night’s cloud width; the one outlier was det2',
              'Left: det3 June-15 per session (board centre); right: whole-night a1 vs median Y share (big dots), and within each night (dotted)')
s9 += row(P.svg('sessions'), P2.svg('width'), gap=44)
s9 += legend([('det3 J15 session', RED, 'box'), ('6-19/20 FEU 7: det2 (mislabelled)', GREY, 'box'), ('6-19/20 FEU 3: the real det3 (372 events)', DRED, 'box')] + WLEG + [('night-to-night trend', INK)], 20)
sl('sessions', s9, '<p>Sessions: 6-18 long run and 6-19/20 longer run are EOS extractions (FEU 7, X only); 6-22/23 and 6-23/24 FEU 3 (drift 1000 / 600 V); 6-25 drift 500 V; '
   '6-27 drift scan (1100 → 100 V: a1 rose 0.05 → 0.16, back to the 1000 V value right after, i.e. reversible); 6-27/28 FEU 7. The n_TOF beam run_79 (July 26, '
   'Ar/iso 90/10) gives 0.05–0.07 with the older whole-board method.</p>'
   '<p>A fixed map × one scale per session fits det3 (0.96 / 1.57 / 0.44; χ² 47/60, 126/61, 47/59), and det2’s two runs differ by ×2.5 with the same bands '
   '(median Y share 0.35 vs 1.07).</p>'
   '<p><b>Within a night the trend is not there</b> (dotted lines, all events): the per-event Y share mostly says which pad type the charge sits on, '
   'and a1 is the same over both. The older "one curve" (a1 of Y-pad-labelled events vs Y share, all sessions on one curve) used the X-peak labels, whose migration '
   'grows as the Y share falls, so it overstated the per-event width dependence. What survives is night-to-night: median Y share 0.52 / 0.77 / 0.93 with a1 '
   '+0.18 / +0.14 / +0.10 (<code>sess_width.py</code>), plus the drift scan and det2’s two runs.</p>'
   '<p><b>6-19/20 was det2</b> (<code>s619/</code>): the "det3 FEU 7" data (a1 0.03, all 8 chips lit) has det2’s readout-chain noise pattern (r = 0.957 with det2 X on 6-22) '
   'and not det3’s on the same FEU 7 slot on 6-17; the pattern follows chain + board, not the FEU slot. det3’s X readout lit only chips 0–3 on 6-17, and on 6-19 that '
   '4-chip readout is the FEU 3 file: 372 events, a1 +0.077 ± ~0.03 (≈ 2σ low, few events). The raw_daq_data HV scan + long run could add more.</p>',
   short='Sessions')

# ------------------------------------------------------------------------------------------------ 9b gain and time
TA = J(HERE, '..', 'time', 'time_a1_det3.json')
P = Plot(810, 560, x=(440, 535), y=(0, 0.3), xlabel='resist (mesh) voltage, V', ylabel='a1 (all events)', title='Resist-voltage scans')
P.xticks([(v, str(v)) for v in (450, 470, 490, 510, 530)]); P.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)])
for run, c, lab in (('mx17_det2_det3_overnight_6-22-26', RED, '6-22/23, drift 1000 V'), ('mx17_det3_det4_overnight_6-23-26', '#8f2433', '6-23/24, drift 600 V'),
                    ('mx17_det3_saturday_scan_6-27-26', '#d9822b', '6-27, drift 1000 V')):
    pts = [r for r in TA if r['run'] == run and r['sub'].startswith(('resist_', 'hv_scan_')) and r['chunks'][0]['amp'] > 400
           and not ('6-22' in run and int(r['mesh']) > 500) and sum(c_['n'] for c_ in r['chunks']) >= 300]
    pts.sort(key=lambda r: int(r['mesh']))
    vals = [(int(r['mesh']), float(np.mean([c_['a1'] for c_ in r['chunks']])), float(np.sqrt(np.mean([c_['a1_e'] ** 2 for c_ in r['chunks']]) / len(r['chunks']))),
             float(np.median([c_['amp'] for c_ in r['chunks']]))) for r in pts]
    P.line([v[0] for v in vals], [v[1] for v in vals], c, tips=[f'{lab}, {v[0]} V: a1 = {v[1]:+.3f} ± {v[2]:.3f}, median central amplitude {v[3]:.0f} ADC' for v in vals])
    for v in vals: ebar(P, v[0], v[1], v[2], c, 4)
P.text(445, 0.28, 'Gain ×6 across a scan: no gain trend', 23, INK, weight=600)
P2 = Plot(810, 560, x=(0.4, 1.15), y=(0, 0.3), xlabel='median Y share of the sub-run (cloud width)', ylabel='a1 (all events)', title='Every det3 June-15 sub-run')
P2.xticks([(v, f'{v:.1f}') for v in (0.4, 0.6, 0.8, 1.0)]); P2.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)])
SC = {'mx17_det2_det3_overnight_6-22-26': RED, 'mx17_det3_det4_overnight_6-23-26': '#8f2433', 'mx17_det3_day_6-25-26': '#e07a86',
      'mx17_det3_saturday_scan_6-27-26': '#d9822b', 'mx17_det3_p2_det1_overnight_6-27-26': '#b04a5a'}
for run, c in SC.items():
    pts = []
    for r in TA:
        if r['run'] != run: continue
        ys = [c_['ysh'] for c_ in r['chunks'] if c_['ysh']]
        if not ys or sum(c_['n'] for c_ in r['chunks']) < 400: continue
        a = float(np.mean([c_['a1'] for c_ in r['chunks']])); e = float(np.sqrt(np.mean([c_['a1_e'] ** 2 for c_ in r['chunks']]) / len(r['chunks'])))
        pts.append((float(np.mean(ys)), a, e, r['sub'], r['t']))
    P2.points([q[0] for q in pts], [q[1] for q in pts], c, r=8, tips=[f'{run} {q[3]} ({q[4]}): a1 = {q[1]:+.3f} ± {q[2]:.3f}, median Y share {q[0]:.2f}' for q in pts])
P2.text(0.62, 0.28, 'One curve in width, all nights', 23, INK, weight=600)
s9b = sd.title('At fixed cloud width the size ignores the gain and hours of running: not charging-up',
               'Direct a1 per sub-run from the raw survey (≤ 2 files each); left: resist scans, right: all sub-runs vs their median Y share')
s9b += row(P.svg('gain'), P2.svg('width all'), gap=44)
s9b += legend([('6-22/23', RED), ('6-23/24', '#8f2433'), ('6-25', '#e07a86'), ('6-27 scans', '#d9822b'), ('6-27/28 overnight', '#b04a5a')], 21)
sl('gain', s9b, '<p><code>time/time_a1.py</code>. Charging-up of insulating surfaces scales with the deposited charge (gain × rate) and builds over minutes to hours. At fixed width the '
   'resist scans change the central amplitude ~×6 (6-23: 525 → 465 V) with a1 flat at 0.15–0.21; 6-22’s 8-hour scan and 6-27 afternoon → overnight show no drift. '
   'Runs where most pulses rail (6-22 at 505–525 V) or sit below ~400 ADC are dropped from the left panel. The 6-23/24 step at 490 V coincides with a width change '
   '(median Y share 0.49 → 0.57), as the right panel’s single curve expects.</p>'
   '<p>Retraction: a "5× jump at fixed HV on 6-23" seen briefly in the survey table came from its copy-model estimator, which misreads low-amplitude runs; on the direct a1 '
   'the three 495 V / 600 V runs are +0.06 ± 0.05, +0.16 ± 0.02, +0.19 ± 0.01.</p>', short='Gain and time')


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

# ------------------------------------------------------------------------------------------------ 10b early sharing by pulse length (2026-10-22)
C22 = os.path.join(HERE, '..', '..', '2026-10-22', 'cloud')
FB = J(C22, 'fwhm_bins_dtn80.json'); SP = J(C22, 'short_pulses.json')
P = Plot(1080, 600, x=(250, 750), y=(0, 2.1), xlabel='central pulse FWHM (ns), |Δt(n+1 − n−1)| < 80 ns', ylabel='early / peak sharing  r(5 % rise) / r(peak)')
P.xticks([(v, str(v)) for v in (300, 400, 500, 600, 700)]); P.yticks([(v, f'{v:.1f}') for v in (0, .5, 1, 1.5, 2)])
P.hline(1, MUT, None, 1.5)
for key, c, dash, lab in (('det3 May X', BLUE, None, 'data det3 May X'), ('det3 May Y', GREEN, None, 'data det3 May Y'),
                          ('sim glg0 X', GREY, '8 6', 'sim baseline X'), ('sim glg0 Y', PURPLE, '8 6', 'sim baseline Y'),
                          ('sim cdry X', ORANGE, '3 5', 'sim, dry-gas D_T, X')):
    pts = [q for q in FB[key] if q['s1'] == [0.45, 0.7]]
    xs = [min(.5 * (q['fw'][0] + q['fw'][1]), 700) if q['fw'][0] else 300 for q in pts]
    ys = [q['r05'][0] / q['peak'][0] for q in pts]; es = [ys[i] * math.hypot(q['r05'][1] / q['r05'][0], q['peak'][1] / q['peak'][0]) for i, q in enumerate(pts)]
    P.line(xs, ys, c, dash=dash, tips=[f'{lab}, FWHM {q["fw"][0]}–{q["fw"][1]} ns, s1 0.45–0.7: r05 {q["r05"][0]:.3f}, peak {q["peak"][0]:.3f}, '
                                        f'f50 {q["f50"][0]:.3f}, n = {q["n"]}' for q in pts])
    if dash is None:
        for x, y, e in zip(xs, ys, es):                       # clipped to the axis (det3 May X at 300 ns: 1.98 ± 0.7)
            lo, hi = max(y - e, 0), min(y + e, 2.1); ebar(P, x, .5 * (lo + hi), .5 * (hi - lo), c)
P.text(430, 1.95, '← short pulses: data only, prompt-wide', 24, INK, weight=600)
sp = lambda k, i: SP[k][i]
rows = [[k, f'{100 * sp(k, 0)["frac"]:.0f} %', f'{sp(k, 0)["onset"][1]:.1f} / {sp(k, -1)["onset"][1]:.1f}', f'{sp(k, 0)["f50"]:.1f} / {sp(k, -1)["f50"]:.1f}']
        for k in ('det3 May X', 'det3 May Y', 'det3 June-15 X', 'sim glg0 X')]
s10b = sd.title('The data’s “flat early sharing” is a population of short pulses the sim does not make',
                'Early/peak sharing vs central pulse length (s1 0.45–0.7). At the sim’s own pulse lengths det3 May X agrees with the sim; the cut of 10-21 was dominated by short pulses')
s10b += row(P.svg('earlyfw'), col(legend([('data det3 May X', BLUE), ('data det3 May Y', GREEN), ('sim baseline X', GREY, 'dash'), ('sim baseline Y', PURPLE, 'dash'),
                                          ('sim dry-gas D_T X', ORANGE, 'dash')], 22),
                                  table(['set', 'FWHM<350', 'onset short/long', 'fall short/long'], rows, 20, widths=[190, 110, 150, 150], align=['left', 'center', 'center', 'center']),
                                  p('Samples of 60 ns from the trigger. det3 May’s short pulses start ~4 samples late (no charge from the top ~9 mm); det3 June-15 (600 V drift) ones end early.', 21, MUT),
                                  w=620, gap=22), gap=36)
sl('earlyfw', s10b, '<p><code>2026-10-22/cloud/</code> (MX17_Geant). The 10-21 near-vertical cut keeps central FWHM ≥ the set’s own 25th percentile. The data’s pulses are short '
   '(det3 May X |Δt| < 30 ns: median 423 ns; sim 590 ns), so that cut compared different populations. In absolute FWHM bins det3 May X matches the baseline sim '
   '(shape χ² 23.9 / 18 on the cells every arm fills). Dry-gas D_T is excluded by the time profile (86 / 18); depth-independent clouds (no diffusion, σ0 250 / 400 µm) '
   'fail (565 / 57); σ0 400 µm (17.5) and feedback spread 400 µm (19.0) are no better than the baseline.</p>'
   '<p>The short pulses are not small x-tilts: sheared th00 tracks at 1–3° and Stage A thx05/thx10 almost never pass |Δt| < 80 ns, and those that do have <i>lower</i> '
   'early sharing. They are partial drift columns (late onset or early end), 29–43 % of det3 May’s sample, 62–84 % of det3 June-15’s, 1–4 % of the sim’s. '
   'What makes them: next slide (§12).</p>', short='Early sharing')

# ------------------------------------------------------------------------------------------------ 10c partial columns: what makes them (2026-10-22, §12)
P22 = os.path.join(HERE, '..', '..', '2026-10-22', 'partial')
EE = J(P22, 'event_ends.json'); RS = J(P22, 'rescore_chi2.json')
P = Plot(1060, 600, x=(6, 31), y=(0, 0.38), xlabel='end of the column (last 30 % crossing of the 9-pad sum, 60 ns samples)', ylabel='fraction of events')
P.xticks([(v, str(v)) for v in (6, 10, 14, 18, 22, 26, 30)]); P.yticks([(v, f'{v:.1f}') for v in (0, .1, .2, .3)])
for key, c, dash, lab in (('d3may_x', BLUE, None, 'det3 May, 900 V'), ('s051_x', GREEN, None, 'det3 June, 6-27, 1000 V'),
                          ('s024_x', GOLD, None, 'det3 June, 6-25, 500 V'), ('s052_x', RED, None, 'det3 June, 6-23, 600 V (the §8–10 reference)')):
    h = np.array(EE[key]['hist_end9'], float); h = h / h.sum(); xs = list(range(6, 32))
    P.line(xs, [h[i] for i in xs], c, dash=dash, tips=[f'{lab}: end at sample {i}: {h[i]:.3f} (n = {EE[key]["n"]} M3-vertical tracks, |θ| < 2°)' for i in xs])
P.text(7, 0.36, '6-23: no column end, a ~0.5 µs burst at a random depth', 23, INK, weight=600)
rows = [[lab, f'{RS[k]["chi2"]:.0f} / {RS[k]["ndof"]}', f'{RS[k]["early_pk"]:.2f}', f'{RS[k]["late_pk"]:.2f}']
        for lab, k in (('June 6-27, 1000 V', 'det3 June 6-27 1000 V (s051+s025) [dtn]'), ('June 6-23, 600 V', 'det3 June 6-23 600 V (d3jun) [dtn]'),
                       ('May, 900 V', 'det3 May 900 V [dtn]'), ('June 6-27, M3 vert.', 'det3 June 6-27 1000 V (s051+s025) [m3v]'), ('May, M3 vert.', 'det3 May 900 V [m3v]'))]
rows.append(['sim glg0', '—', f'{RS["det3 June 6-27 1000 V (s051+s025) [dtn]"]["sim_early_pk"]:.2f}', f'{RS["det3 June 6-27 1000 V (s051+s025) [dtn]"]["sim_late_pk"]:.2f}'])
s10c = sd.title('The partial columns are inclined tracks (May) and a degraded session (June 6-23); a clean det3 June X matches the sim',
                'Left: where each M3-vertical track’s drift column ends. Right: X timing at matched pulse length vs the sim (χ² over r(5 %, 25 % rise, peak, 50 % fall))')
s10c += row(P.svg('colend'), col(legend([('det3 May, 900 V', BLUE), ('det3 June 6-27, 1000 V', GREEN), ('det3 June 6-25, 500 V', GOLD), ('det3 June 6-23, 600 V', RED)], 22),
                                  table(['set', 'χ² vs sim', 'early/pk', 'late/pk'], rows, 20, widths=[210, 120, 110, 110], align=['left', 'center', 'center', 'center']),
                                  p('det3 May: the short pulses have M3 |θ| 6° across the pads (long ones 1.3°); the |Δt| < 80 ns cut lets them in.', 21, MUT),
                                  w=620, gap=22), gap=36)
sl('colend', s10c, '<p><code>2026-10-22/partial/</code> (MX17_Geant), M3 tracks joined by eventId (Chi2X, Chi2Y < 5). The view X ↔ M3 y mapping holds in May too (chip medians). '
   '<b>det3 May</b>: the |Δt(n+1 − n−1)| < 80 ns sample keeps inclined tracks; FWHM < 350 ns pulses have |θ_meas| 6.3° (long 1.3°). With M3 |θ| < 2° only 11 % are short.</p>'
   '<p><b>det3 June 6-23</b> (the d3jun reference, every 600 V sub-run): per vertical muon the column arrives as a ~0.5 µs burst at a random depth (onset p10–p90 4–15 '
   'samples, end flat 9–31); X and Y (FEU 3 and 4) agree event by event (onset r = 0.83); det4 in the same events is uncorrelated (r = 0.04), so not the trigger; 95 % '
   'have the M3 track on the hit chip, so not other particles; uniform over angle, board position, M3 multiplicity and time in the run; no dropped-sample floors. '
   '9-pad amplitude 4–5× below the other sessions and nearly flat over the 465–525 V resist scan (×1.5 vs ×8 on 6-27); α+1 0.08–0.2 vs 0.3–0.47. HV readback normal '
   '(600 V, 495 V; resist current 0.62 vs ~1.1 µA). Not the drift field: 6-25 at 500 V and the 6-27 scan at 500–1100 V have sharp column ends. Mechanism open '
   '(a resistive layer not following its supply, i.e. charging, fits the flat gain; electron loss in the gas does not explain random depths).</p>'
   '<p><b>Step 2</b>: clean det3 June X (6-27 long_run + overnight, 1000 V = sim) at matched absolute FWHM: χ² 34 / 16 vs the sim (det3 May 121 / 36, 6-23 639 / 24); '
   'late growth not suppressed (2.02 vs sim 1.87). Early sharing sits above the sim on both boards alike (May 0.83, June 0.78, sim 0.58). So §8–10’s "det3 June '
   'prompt-wide, flat, suppressed late growth" fingerprint is the 6-23 session, not the June bulk. The asymmetry itself is in the clean sessions (survey a ≈ 0.07–0.11).</p>',
   short='Partial columns')

# ------------------------------------------------------------------------------------------------ 10d degraded nights + det2 re-score (2026-10-22, §13)
STB = J(P22, 'slot_table.json'); D2 = J(P22, 'det2_rescore.json'); BB = J(P22, 'board_vs_board.json'); SC17 = J(P22, 'scan617.json')
PHYS = {('mx17_det1_det2_overnight_6-17-26', 'mx17_1'): 'det2', ('mx17_det2_det3_overnight_6-22-26', 'mx17_2'): 'det2', ('mx17_det3_long_run_5-6-26', 'mx17_1'): 'det3 May'}
trows = []
for r in sorted(STB, key=lambda r: (r['run'][-7:], r['z'])):
    if r['n_vert'] < 50: continue
    bad = r['early_end'] > 0.1
    nm = PHYS.get((r['run'], r['det']), r['det'].replace('mx17_', 'det'))
    cell = (lambda t: f'<b style="color:{RED}">{t}</b>') if bad else (lambda t: t)
    trows.append([cell(re.search(r'(\d+-\d+)-26', r['run']).group(1)), cell(r['sub'].split('_resist')[0].replace('_p2_det1_sanity_check', '')[:10]), cell(nm), cell('bottom' if r['z'] < 500 else 'top'), cell(str(r['drift'])),
                  cell(f"{r['A9']:.0f}"), cell(f"{100 * r['early_end']:.0f} %"), cell(f"{r['s1']:.2f}")])
d2 = lambda sel: D2[f'det2 6-22 (s004+s006) [{sel}]']
s10d = sd.title('Two degraded nights, both in the bottom slot; det2’s clean timing is sim-like too',
                'Left: every detector-session with M3 tracks (M3-vertical tracks; red = columns cut short). Right: det2 6-22 (top slot, 1000 V) vs the sim and vs det3 June')
s10d += row(table(['night', 'run', 'board', 'slot', 'drift V', 'A9 (ADC)', 'cut short', 's1'], trows, 19, widths=[80, 140, 100, 90, 90, 110, 100, 70],
                  align=['left', 'left', 'left', 'left', 'center', 'center', 'center', 'center']),
            col(table(['det2 6-22', 'χ² vs sim', 'r25/pk', 'late/pk'],
                      [[f'|Δt| < 80', f"{d2('dtn')['chi2']:.0f} / {d2('dtn')['ndof']}", f"{d2('dtn')['ratios'][1]:.2f} ({d2('dtn')['ratios'][4]:.2f})", f"{d2('dtn')['ratios'][2]:.2f} ({d2('dtn')['ratios'][5]:.2f})"],
                       [f'M3 |θ| < 2°', f"{d2('m3v')['chi2']:.0f} / {d2('m3v')['ndof']}", f"{d2('m3v')['ratios'][1]:.2f} ({d2('m3v')['ratios'][4]:.2f})", f"{d2('m3v')['ratios'][2]:.2f} ({d2('m3v')['ratios'][5]:.2f})"]],
                      20, widths=[150, 120, 140, 140], align=['left', 'center', 'center', 'center']),
                p('(sim in brackets). det2 − det3 June: r25/pk ' + f"{BB['dtn']['r25/pk diff'][0]:+.2f} ± {BB['dtn']['r25/pk diff'][1]:.2f}, late/pk {BB['m3v']['f50/pk diff'][0]:+.2f} ± {BB['m3v']['f50/pk diff'][1]:.2f} (M3): "
                  'a ≤ 2σ lean to “narrow-prompt, late-wide”, not the clean contrast of §8.', 21, MUT),
                p('Degraded = 9-pad amplitude ~4× low, columns cut at random depths, gain hardly following the resist HV (6-23: ×1.5, 6-17: ×2.6 over 60 V; healthy ≈ ×6).', 21, MUT),
                w=620, gap=22), gap=40)
sl('degraded', s10d, '<p><code>2026-10-22/partial/slot_table.py</code>, <code>scan617.py</code>, <code>det2_rescore.py</code>, <code>board_vs_board.py</code>. Board names: physical boards as in '
   'copy_survey (the 6-17 run_config’s descriptions are swapped). Both degraded nights put the board in the bottom slot (z 232 mm) on drift 0:7 and resist 3:3, but the same '
   'slot and channels were healthy on 6-22 (det3, slower columns), 6-25 (det3) and 6-26 (det6), so it is a per-night condition, not the slot. HV readback normal where logged. '
   'A resistive layer not following its supply (bad contact, floating) fits the flat gain; nothing here tests it. The bench logbook for 6-17 and 6-23 is the next check.</p>'
   '<p>det2’s §8 timing used 6-22 (clean). Re-scored with the same cells as det3: near the sim, prompt sharing not narrower than the sim’s, late growth slightly above it. '
   'The §8 statement “det2 narrow-prompt / late-wide, det3 June the opposite” is reduced to a ≤ 2σ difference between two sim-like boards: timing no longer constrains '
   'the strip routes either way.</p>', short='Degraded nights')

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
      ('Charging-up of insulating surfaces', 'a1 flat with gain (×6) and over hours at fixed width', MOD),
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
    o.append(T(W / 2, y0 + 204, 'even effect: the cloud is concentrated → s1 changes with δ', 23, BLUE))
    return svg(W, H, ''.join(o), 'registration cartoon')
sd_rows = [r for r in SIMD if r]
s0 = [r for r in sd_rows if abs(r['d']) < 30]; s0v = np.mean([r['s1'] for r in s0])
P = Plot(860, 600, x=(-360, 360), y=(-0.2, 0.55), xlabel='true strip offset δ in the sim (µm)', ylabel='change vs δ = 0')
P.xticks([(v, f'{v:+d}' if v else '0') for v in (-300, -150, 0, 150, 300)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)])
P.hline(0, MUT, None, 1.5)
P.line([r['d'] for r in sd_rows], [r['a1'] for r in sd_rows], RED, tips=[f'sim δ = {r["d"]:+.0f} µm: a1 = {r["a1"]:+.3f} ± {r["e"]:.3f}, n = {r["n"]}' for r in sd_rows])
P.line([r['d'] for r in sd_rows], [r['s1'] - s0v for r in sd_rows], BLUE, tips=[f'sim δ = {r["d"]:+.0f} µm: s1 = {r["s1"]:.3f} (s1 − s1(0) = {r["s1"] - s0v:+.3f})' for r in sd_rows])
P.text(-250, 0.48, 'Sharing moves with δ (vs δ = 0)', 25, INK, weight=600)
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

# ------------------------------------------------------------------------------------------------ 16 the even term doesn't decide
LOOP = J(HERE, '..', 'loop', 'toy_loop.json'); SP = J(HERE, '..', 'even', 'same_pcb.json')
P = Plot(860, 600, x=(-0.05, 0.32), y=(-0.55, 0.25), xlabel='a1', ylabel='Δs1 vs an unpulled board')
P.xticks([(v, f'{v:.2f}') for v in (0, 0.1, 0.2, 0.3)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2)])
P.hline(0, MUT, None, 1.5); P.vline(0, MUT, None, 1.5)
LC = {(450., 0.2): PURPLE, (300., 0.2): BLUE, (600., 0.2): GREEN}
for t in LOOP:
    key = (t['sigma'], t['f'])
    if t['sig_ind'] != 200. or t['sig_line'] != 200. or key not in LC: continue
    half = [(d, a, s) for d, a, s in zip(t['delta'], t['a1'], t['ds1']) if d <= 390 and a <= 0.31]
    P.line([h[1] for h in half], [h[2] for h in half], LC[key], w=3, r=4,
           tips=[f'toy: cloud σ {t["sigma"]:.0f} µm, pulled fraction f {t["f"]}, δ = {d:.0f} µm → a1 {a:+.3f}, Δs1 {s:+.3f} (vs f = 0)' for d, a, s in half])
dm, dj = SP['builds']['det3 May'], SP['builds']['det3 June-15']
diffs = [(j[0] - m[0], np.hypot(j[1], m[1])) for m, j in zip(dm, dj) if m and j][1:-1]
w = np.array([1 / e ** 2 for _, e in diffs]); dmean = float(sum(d * wi for (d, _), wi in zip(diffs, w)) / w.sum()); derr = float(1 / np.sqrt(w.sum()))
a1J = B['det3 June-15']['a1x'][0] if 'a1x' in B['det3 June-15'] else 0.137
P.points([a1J], [dmean], GOLD, r=14, tips=[f'det3 June-15 vs det3 May (same PCB) at matched Y share: Δs1X = {dmean:+.3f} ± {derr:.3f} (7 inner bins), a1 = {a1J:+.3f}'])
P.text(0.005, 0.12, 'δ → P/2 (strip over the boundary)', 21, MUT); P.text(0.005, -0.45, 'δ = 0 (strip over the pad centre)', 21, MUT)
P.text(0.005, 0.2, 'A pull can leave the sharing unchanged', 23, INK, weight=600)
P2 = Plot(760, 600, x=(0.0, 1.45), y=(0.3, 1.0), xlabel='Y share (cloud width, other view)', ylabel='X sharing s1X', title='Same PCB, two bulks')
P2.xticks([(v, f'{v:.1f}') for v in (0, 0.3, 0.6, 0.9, 1.2)]); P2.yticks([(v, f'{v:.1f}') for v in (0.3, 0.5, 0.7, 0.9)])
E_ = SP['edges']
for nm, c in (('det3 May', BLUE), ('det3 June-15', RED)):
    rows_ = SP['builds'][nm]; xs = [0.5 * (E_[k] + E_[k + 1]) for k in range(len(rows_)) if rows_[k]]; ys = [r[0] for r in rows_ if r]
    P2.line(xs, ys, c, tips=[f'{nm}, Y share {E_[k]:.2f}–{E_[k+1]:.2f}: s1X {r[0]:.3f} ± {r[1]:.3f}, a1 {r[2]:+.3f}, n = {r[3]}' for k, r in enumerate(rows_) if r])
P2.text(0.08, 0.95, 'a1 0 → 0.14, sharing unchanged', 23, INK, weight=600)
s16 = sd.title('Against an unpulled board the strip pull’s sharing change can be zero: the even term doesn’t decide',
               'Left: 1-D attractor toy over a half period of δ, Δs1 relative to no pull; gold = det3 June-15 vs May. Right: same-PCB X sharing at matched width')
s16 += row(col(P.svg('loop'), legend([('cloud σ 300, f 0.2', BLUE), ('σ 450, f 0.2', PURPLE), ('σ 600, f 0.2', GREEN),
                                     ('det3 June-15 − May', GOLD, 'dot')], 20), gap=8),
           col(P2.svg('same pcb'), legend([('det3 May (quiet)', BLUE), ('det3 June-15', RED)], 21), gap=8), gap=40)
sl('even', s16, f'<p><b>Correction (2026-10-21 evening).</b> The first version of this deck said the sharing (the even term) excludes registration "×20". That compared the sharing '
   'at large |a1| with the sharing at δ = 0 inside the same attractor (the sim’s s1 = 0.64 + 2.2 |a1| fit, and the det2 cell regression). An attractor narrows the cluster most '
   'at δ = 0 and leaves it about unchanged near δ = P/2, while a1 is zero at both. Measured against an unpulled board, a pull at offsets between the a1 maximum and P/2 gives '
   'a1 0.1–0.3 with Δs1 ≈ 0.</p>'
   f'<p>Same PCB (<code>even/same_pcb.py</code>): det3 May and June-15 X sharing at matched Y share differ by {dmean:+.3f} ± {derr:.3f} over the inner bins while a1 goes 0 → 0.14: '
   'consistent with no pull, or with a pull sampled around the a1 maximum. Toy: <code>loop/toy_loop.py</code> (pad response box ⊗ 200 µm Gaussian, pulled charge smeared 200 µm, '
   'peak-pad selection, uniform positions). The 25 sim arms’ E/O = 1.2–2.5 (<code>even/sim_even.py</code>) is the same within-attractor contrast, not an exclusion.</p>', short='Even term')

# ------------------------------------------------------------------------------------------------ 17 the quadrature loop in data
CR = J(HERE, '..', 'loop', 'det2_crossings_np.json')
def looppanel(key, ttl, xl):
    rs = CR[key]; P = Plot(810, 560, x=(-0.4, len(rs) - 0.6), y=(-0.25, 0.25), xlabel=xl, ylabel='per band, all events', title=ttl)
    P.xticks([(i, str(r['band'])) for i, r in enumerate(rs)]); P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.2, -0.1, 0, 0.1, 0.2)]); P.hline(0, MUT, None, 1.5)
    P.line(list(range(len(rs))), [r['a1'] for r in rs], RED, tips=[f'{key} band {r["band"]}: a1 = {r["a1"]:+.3f} ± {r["a1_e"]:.3f}, n = {r["n"]}' for r in rs])
    P.line(list(range(len(rs))), [r['r'] for r in rs], BLUE, tips=[f'{key} band {r["band"]}: s1X − F(s1Y) = {r["r"]:+.3f} ± {r["r_e"]:.3f}, mean Y share {r["s1Y"]:.2f}' for r in rs])
    for i, r in enumerate(rs): ebar(P, i, r['r'], r['r_e'], BLUE, 5)
    return P
PA = looppanel('det3 June-15 by x', 'det3 June-15 (asymmetric)', 'x band (64 true columns)')
PA.text(-0.2, 0.22, 'Sharing crosses zero at the a1 peak', 22, INK, weight=600)
PB = looppanel('det3 May by x', 'det3 May, same PCB (quiet)', 'x band (64 true columns)')
PB.text(-0.2, 0.22, 'Same-size sharing structure, no a1', 22, INK, weight=600)
s17 = sd.title('The June boards trace the loop a sweeping offset predicts, but the quiet board’s sharing varies as much',
               'Per 64-column band: a1 (red) and the X sharing residual s1X − F(s1Y) (blue; F = the board’s own binned X-vs-Y relation)')
s17 += row(PA.svg('loop J15'), PB.svg('loop May'), gap=44)
s17 += legend([('a1', RED), ('sharing residual (width removed)', BLUE)], 22)
d2 = CR['det2 by y']
sl('loop', s17, '<p>An attractor whose offset δ sweeps across a board (print scale or rotation) gives a1 ∝ sin φ and a sharing change ∝ −cos φ: the sharing residual should '
   'cross zero where a1 peaks and be extreme where a1 crosses zero. det3 June-15 does exactly that across x (+0.08 → 0.01 at the a1 maximum → −0.11). det2 by y '
   f'(16 bands): +{max(r["r"] for r in d2[:5]):.2f} at its rising a1 zero crossing, about 0 near the extrema, −0.04 … −0.08 at the falling crossing. det4 by x: +0.08 / +0.05 / −0.19 '
   'with a1 0.03 / 0.11 / 0.03.</p>'
   '<p>But det3 May, the same PCB with no asymmetry, has residuals of −0.15 and +0.13 in adjacent bands, so bulk-local sharing structure (gap, gain) is as large as the loop. '
   'June-15 doesn’t repeat May’s pattern, so it isn’t the PCB or the electronics. The loop is suggestive, not evidence (<code>loop/det2_crossings_np.py</code>).</p>', short='Loop')

# ------------------------------------------------------------------------------------------------ 17b what a physical pull needs
FR = J(HERE, '..', 'loop', 'toy_funnel_regime.json'); TL = J(HERE, '..', 'loop', 'toy_local.json')
P = Plot(810, 600, x=(60, 560), y=(-0.02, 0.3), xlabel='strip width (µm) at a 780 µm pitch', ylabel='largest a1 with unchanged sharing (|Δs1| ≤ 0.03)')
P.xticks([(v, str(v)) for v in (100, 200, 300, 400, 500)]); P.yticks([(v, f'{v:.2f}') for v in (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)])
P.raw(f'<rect x="{P.x0}" y="{P.Y(0.20):.1f}" width="{P.pw}" height="{P.Y(0.10) - P.Y(0.20):.1f}" fill="{RED}" fill-opacity="0.10"/>', back=True)
P.text(360, 0.185, 'det3 June-15 range of a1', 21, RED)
FC = {100.: '#5b3a8c', 150.: PURPLE, 200.: '#b79be0'}
for sg in (100., 150., 200.):
    rs = sorted([r for r in FR if r['sig_ind'] == 350. and r['sigma'] == sg], key=lambda r: r['w'])
    xs = [r['w'] for r in rs] + [540]; ys = [r['best'] for r in rs] + [0.0]
    P.line(xs, ys, FC[sg], tips=[f'saturated-gap funnel, cloud σ {sg:.0f} µm, induction 350 µm (no-pull s1 {rs[0]["s1_nopull"]:.2f}), {w:.0f} µm strips: a1 {v:+.3f}' for w, v in zip(xs, ys)][:-1] + ['540 µm strips: funnel alone gives a1 ≈ 0 (toy_funnel.py)'])
loc = [r for r in TL if r['kind'] == 'local' and r['sigma'] == 450. and np.isfinite(r['best'])]
P.hline(max(r['best'] for r in loc), GREY, '6 5', 2, tip='toy_local.py: per-electron pulls of 20–100 % onto 50/200 µm lines, σ 450 µm Gaussian cloud, 100/200 µm induction')
P.text(150, 0.28, 'Only narrow strips + a narrow cloud reach it', 23, INK, weight=600)
P2 = Plot(810, 600, x=(50, 1150), y=(0, 1.1), xlabel='drift voltage (V, 6-27 drift scan)', title='Low drift field: narrower cloud, bigger a1')
P2.xticks([(v, str(v)) for v in (100, 300, 500, 700, 900, 1100)]); P2.yticks([(v, f'{v:.1f}') for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)])
ds = sorted([r for r in TA if r['sub'].startswith('drift_scan')], key=lambda r: int(r['drift']))
dv = [int(r['drift']) for r in ds]; da = [float(np.mean([c_['a1'] for c_ in r['chunks']])) for r in ds]; dy = [r['chunks'][0]['ysh'] for r in ds]
P2.line(dv, dy, BLUE, tips=[f'{v} V: median Y share {y:.2f}' for v, y in zip(dv, dy)])
P2.line(dv, [5 * a for a in da], RED, tips=[f'{v} V: a1 = {a:+.3f} (drawn ×5)' for v, a in zip(dv, da)])
s17b = sd.title('In a 1-D toy a physical pull reaches the data’s size onto narrow strips; the full sim does not',
                'Left: 1-D toy, gaps charged to saturation funnel every gap electron onto the nearer strip edge; right: the drift scan')
s17b += row(col(P.svg('funnel width'), legend([('cloud σ 100 µm', '#5b3a8c'), ('σ 150 µm', PURPLE), ('σ 200 µm', '#b79be0'), ('per-electron pull on a σ 450 µm cloud', GREY, 'dash')], 20), gap=8),
            col(P2.svg('drift scan'), legend([('median Y share (cloud width)', BLUE), ('a1 × 5', RED)], 20), gap=8), gap=40)
sl('funnel', s17b, '<p>A detector responds electron by electron, so a real pull moves each electron by where it lands (the "cloud-level" toy of the transfer slide does not; '
   'it is a shape check only). Per-electron pulls on σ 300–450 µm Gaussian clouds stay at a1 ≤ 0.04–0.07 with unchanged sharing (<code>loop/toy_local.py</code>).</p>'
   '<p>Charged insulating gaps saturate at σ = ε₀E₀ ≈ 3.5 nC/cm² (39 kV/cm): then no field line ends on a gap and every gap electron lands on the nearer strip edge. '
   'With the quoted 540 µm strips that funnel gives a1 ≈ 0 (each gap feeds both edges). The sharing fixes only cloud width ⊗ induction: with a narrow electron '
   'cloud (σ 100–200 µm) and wide induction (350 µm; no-pull s1 0.59–0.67, the data’s range) the funnel onto narrow strips reaches a1 0.1–0.24 '
   '(<code>loop/toy_funnel_regime.py</code>).</p>'
   '<p>A narrow cloud is plausible: a vertical track’s electrons are a depth mixture (30 mm gap), whose pitch-scale content is ~10× that of a Gaussian of the same '
   'rms, and attachment in the ~1 % H₂O gas at low drift field removes the deep electrons; the drift scan shows the Y share narrowing (0.98 → 0.56) as a1 rises.</p>'
   '<p>Prediction tested: narrow strips (≲ 300 µm) with wide gaps on the June prints. <b>The full sim says no.</b> Paired arms vs the 550 µm control '
   '(36 angle/seed pairs each, 190 / 290 µm X kernels, matched funnel; <code>paired_odd_narrow.json</code>): odd amplitude ×0.21 (190 µm), ×0.59 (190 µm + funnel), '
   '×0.52 (290 µm), ×0.92 (290 µm + funnel), against the ×2.5–5 needed; the funnel arms also narrow the sharing (Δs1 −0.05 … −0.10). In the full sim the asymmetry '
   'grows with strip width, opposite to the toy: the sim’s real depth-mixture clouds and resistive induction are not the toy’s narrow-cloud regime.</p>', short='Physical pull')


# ------------------------------------------------------------------------------------------------ 18 fingerprint matrix
Y_, N_, Q_ = chip('✓', GREEN), chip('✗', RED), chip('?', GREY)
FP = [('X view only; Y flat in both pad classes', Y_, Y_, Q_),
      ('Odd n−1 → n+1 transfer, ±1 sum unchanged', '✓ (allowed vs no pull)', '✓ (first order)', Q_),
      ('Direct n−1 → n+1: ±2 carries ~1/8 of ±1', 'plausible (position-level)', Q_, Q_),
      ('Smooth two-signed field, fixed per board', Y_, Q_, Q_),
      ('Pad-pitch locked, no 31.2 mm beat', '✓ only with a 780 µm print', Y_, Q_),
      ('Same over X pads and Y pads', Y_, Y_, Q_),
      ('Size: board mean 0.14, cells ±0.2', '✗ narrow strips: toy yes, full sim ×0.2–0.9', Q_, Q_),
      ('Bigger on narrow-cloud nights and at low drift field', Y_, Q_, 'plausible (field)'),
      ('Independent of gain (×6) and hours at fixed width', Y_, Y_, N_),
      ('Late on short pulses, early on long', '✓ anisotropic only', Q_, Q_),
      ('June builds only', 'needs a print difference', Q_, 'June material')]
s18 = sd.title('The fingerprint any explanation must reproduce', 'Rows: data facts from the previous slides. Columns: the candidates still standing')
s18 += table(['fact', 'strip registration (as simulated)', 'skewed X response (in the sim; ion paths, fast strips, gradients: ✗)', 'surface charging (untested)'], [list(r) for r in FP], 23,
             widths=[640, 380, 330, 320], align=['left', 'center', 'center', 'center'])
sl('fingerprint', s18, '<p><b>Skewed X response</b>: with uniform illumination, a1 ≠ 0 needs a response that is not mirror-symmetric about the pad. Registration gets there by '
   'attracting the charge (and pays with the even term). The other way is a weighting field that leans toward one neighbour, without moving the charge: the centroid '
   'shifts, the ±1 sum is unchanged at first order, and the pad type doesn’t matter. Offset strips acting as partly conducting electrodes in the induction are one '
   'way to get it, and a smooth two-signed lean would follow from the print registration. Only its static version was modelled (route 9, rejected on a class '
   'criterion now void, and because it grew with cloud size).</p>'
   '<p><b>Surface charging</b>: insulating surfaces (coverlay, gap kapton) on fresh June bulks charging up and adding a lateral field. It would naturally give time and '
   'field dependence, and an HV incident could reset it (6-19). It has no worked prediction for the map geometry or pad-pitch locking yet.</p>'
   '<p>These are candidates to test, not results.</p>'
   '<p><b>Update (FOLLOWUP §7):</b> the skewed response is not untested after all. Every full-sim arm runs on W2 kernels whose X box is the 31.2 mm '
   'superperiod, so offset resistive strips at every δ are in them, and a 2-D floating-strip toy at the sim’s sharing gives the sim’s size (O ≈ 0.04). '
   'The one piece the sim lacked, avalanche ions returning along the funnel’s field lines, reduces the funnel’s asymmetry by 5–10 % (<code>skew/ion_path_toy.py</code>).</p>', short='Fingerprint')

# ------------------------------------------------------------------------------------------------ 19 way forward
WF = [('done', 'Even term vs a no-pull baseline; transfer shape; a1 vs gain and time; 6-19/20 identified as det2', 'archive + sim', 'even term doesn’t decide; not charging; no outlier'),
      ('1', 'Score models on the transfer shape (±1 vs ±2) and the size, not on pad classes or the within-attractor even contrast', 'existing sim', 'which position-level mechanisms survive'),
      ('done', 'Narrow strips (190 / 290 µm) + saturated-gap funnel in the full sim: odd amplitude ×0.2–0.9 of the 550 µm control, not ×2.5–5', 'sim', 'the toy’s size does not survive the real kernels'),
      ('done', 'Skewed induction: the W2 kernels already carry offset resistive strips at every δ (31.2 mm box) and sit at the conducting-strip limit (2-D toy: O ≈ 0.04 at the sim’s sharing); funnel ion return paths subtract 5–10 %', 'sim + toy', 'no strip-induction route reaches ×2.5–5 at the data’s sharing'),
      ('done', 'Prompt sharing, data vs sim at matched s1 (near-vertical tracks): det3 May = sim; det2 narrow-prompt / late-wide (the strip regime); det3 June-15 the opposite, wide before the peak and flat in time, which no sim arm makes', 'archive + sim', 'no strip route covers both asymmetric boards (det3 June side withdrawn by §12: the 6-23 session)'),
      ('done', 'det3 June’s flat X sharing is response, not coupling: flat from 5 % of the rise, s1 floor allows ≤ 0.06 per neighbour, same in every session; it looks like strips that equalise across their width before the first sample', 'archive', 'fast-equalising strips proposed; refuted by the full sim (next row)'),
      ('done', 'Fast-equalising strips (σ_xx = 100 σ_yy) in the full sim: ×1.04 (×1.41 with funnel), r(f) still rising; a resistivity gradient across the strip (19:1) changes a1 by ≤ 0.012', 'sim + toy', 'neither; det3 June X instead lacks LATE X spreading'),
      ('done', 'Sim early-cloud width: at matched pulse length det3 May X = sim (dry D_T and depth-independent clouds excluded); the “flat early sharing” is a data-only population of short, partial-column pulses (29–84 % of the cut), not tilt', 'sim arms + archive', 'compare timing only at matched FWHM; what are the partial columns?'),
      ('done', 'Partial columns (§12, M3 tracks): det3 May’s are tracks inclined across the pads (6° vs 1.3°) that pass the |Δt| cut; det3 June’s come from the 6-23 session, where det3 collected each muon column as a ~0.5 µs burst at a random depth (not trigger, tracks, DAQ, position or field). Clean det3 June X (6-27, 1000 V) vs sim at matched FWHM: χ² 34 / 16, late growth normal', 'archive + M3', 'det3 June X timing has no June-specific feature; 6-23 cause open'),
      ('done', 'det2’s regime re-scored with M3 angles on its clean 6-22 sessions: near the sim (χ² 56 / 20), a ≤ 2σ lean to narrow-prompt / late-wide vs det3 June. The degraded state recurs once more: det2 on 6-17, also bottom slot, gain hardly following the resist HV', 'archive + M3', 'timing no longer constrains the strip routes'),
      ('3', 'What degraded the bottom-slot board on 6-17 and 6-23 (bench logbook: resist HV contact / cable at the bottom slot)', 'logbook', 'whether a floating resistive layer explains it'),
      ('4', 'Per-session, width-matched maps for det3; decode 6-19/20 FEU 3 for more events', 'archive', 'is the geometry really fixed'),
      ('5', 'Photograph the strips (or get the mask) on det2 rows 0–63 vs 128–255, and on det3 May', 'needs the board', 'strip width and gap; opposite δ there, 780 µm pitch'),
      ('6', 'Collimated source across one pad; HV ramps and dwell; dry-gas flush', 'needs the detector', 'attractor vs displacement; charging; humidity')]
s19 = sd.title('Way forward, ranked: cheapest decisive tests first', 'Steps 1–4 need nothing new; 5 needs a photo; 6 needs the detector')
s19 += table(['#', 'test', 'needs', 'decides'], [list(r) for r in WF], 23, widths=[50, 980, 250, 380], align=['center', 'left', 'left', 'left'])
sl('forward', s19, '<p>The even-term score proposed earlier in this audit turned out not to discriminate (slide "Even term"): the sharing change of a pull runs both ways relative '
   'to an unpulled board, and bulk-local sharing structure is ±0.1. The shape of the transfer (±1 vs ±2) and the size are the scores that remain.</p>'
   '<p>Step 6 is the most decisive physical measurement: registration predicts strip-centre offsets of opposite sign (±100–300 µm) in det2’s positive and negative bands, '
   'and a 780 µm pitch. A phone photo through a microscope eyepiece at a few points would do. Read strip width and gap too.</p>'
   '<p>Analysis does not wait for the mechanism: <code>wft.xasym</code> (nTof_x17, commit 30198e4) corrects det2 and det3 June-15 positions with a measured map × '
   'per-run scale (K = −0.77 ± 0.13 vs M3). It removes a positional systematic of up to ±70 µm and leaves the resolution unchanged.</p>', short='Way forward')

# ------------------------------------------------------------------------------------------------ 20 close
cl = (sd.kicker('Where it stands')
      + '<h2 style="font-size:58px;font-weight:600;line-height:1.1;width:1660px">The effect is well measured and the analysis can be corrected. Its cause is open: '
        'a charged-gap funnel onto narrow strips fits in a toy; the strips’ width decides</h2>'
      + row(card(p('What this does not rule out', 30, DINK, 600)
                 + p('• Strips that skew the induced signal without concentrating the charge (registration through induction).', 26, DMUT)
                 + p('• A difference in the June print (pitch, width, fill) that the mask or a photo would show.', 26, DMUT)
                 + p('• A static (saturated) surface charge on fresh June bulks: gain-independent, unlike charging-up.', 26, DMUT)
                 + p('• Something in the June bulk process nobody has named yet (Stephan: "nothing different").', 26, DMUT), bg='#1f2533'),
            card(p('Open decisions', 30, DINK, 600)
                 + p('• Next model: an electrostatic pull from charged gaps, or skewed induction?', 26, DMUT)
                 + p('• Can someone photograph a June board’s strips (det2 ideally) and det3 May?', 26, DMUT)
                 + p('• Adopt wft.xasym in production (needs a re-alignment and per-run scales)?', 26, DMUT)
                 + p('• What differed in the June print (strip width, pitch, gap)?', 26, DMUT), bg='#1f2533'), gap=40))
sl('close', cl, '<p>Audit text with every route re-graded (excluded / weakened / supported / untested, flagged for contaminated controls, board averages, selection effects '
   'and unvalidated models): <code>scratchpad/2026-10-21/audit/AUDIT.md</code> in MX17_Geant. Builder: <code>scratchpad/2026-10-21/audit/make_audit_deck.py</code>.</p>',
   dark=True, short='Open')

OUTP = os.path.join(HERE, 'det3-x-asymmetry.html')
D.write(OUTP, note_meta=dict(title='det3 X neighbour asymmetry: an audit',
                             summary='Fresh-eyes summary of the week-long study, updated the same evening. Every June board carries a smooth, two-signed, X-only field that moves '
                                     'charge straight from neighbour n−1 to n+1, with a size set by the cloud width alone (not gain, not time). Electronics, connectors, copper, '
                                     'copies, linear coupling, Y→X transfer, charging and gap conduction are excluded; the 6-19/20 outlier was det2. A pull toward the strips '
                                     'is still possible (the sharing test does not decide it) but no simulated version is big enough. Fingerprint and ranked plan.',
                             tags='MX17,detector,crosstalk', date='2026-10-21'))
print('wrote', OUTP)
