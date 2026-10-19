"""Third-round slides for the det3 X-asymmetry deck (2026-10-20), imported by ../2026-10-19/asym_deck/make_asym_deck.py.
Inputs (all in this directory): esl_bins.json (esl_bins.py), inject_parity.json (inject_parity.py),
uncopy_parity.json (uncopy_parity.py), fwd_variants.json (fwd_variants.py), ypad_probe.json (ypad_probe.py),
migrate.py output quoted in the text."""
import os, json, numpy as np
import slidedoc as sd
from slidedoc import INK, MUT, BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN
HERE = os.path.dirname(os.path.abspath(__file__))
CJUN = RED


def _j(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def _det(*paras):
    return ''.join(f'<p>{p}</p>' for p in paras)


def _src(*items):
    return '<p>Sources: ' + ', '.join(f'<code>{i}</code>' for i in items) + '</p>'


def _err(P, x, v, e, c=INK):
    P.raw(sd.line(P.X(x), P.Y(v - e), P.X(x), P.Y(v + e), c, 2))


def s_stripes(D):
    e = _j('esl_bins.json')
    if not e: return
    j = e['det3 June-15']; offs = np.arange(10) * 80 + 40
    P = sd.Plot(800, 540, x=(0, 800), y=(0, 0.4), xlabel='stripe offset under the pad column (µm, 4 columns per bin)',
                ylabel='copy, Y-pad events', title='The copy against the stripe position')
    P.xticks([(v, str(v)) for v in (0, 200, 400, 600, 800)]); P.yticks([(v, f'{v:.1f}') for v in (0, .1, .2, .3, .4)])
    xs = np.linspace(0, 800, 81); m = j['ypad']['mean']
    P.band(list(xs), [m - 0.1] * len(xs), [m + 0.1] * len(xs), GOLD, 0.12,
           tip='Swing a floating-stripe coupling would give over one beat (±0.1, fwd_out/f_off±200)')
    P.line(list(xs), list(m + 0.1 * np.sin(2 * np.pi * xs / 800)), GOLD, 2.5, dash='8 6', markers=False,
           tip='Floating-stripe model, ±200 µm offset gives ±0.10 (phase arbitrary): what a stripe-mediated copy would look like')
    v = [c[0] for c in j['ypad']['copy']]; s = [c[1] for c in j['ypad']['copy']]
    for x_, a, b in zip(offs, v, s): _err(P, x_, a, b, CJUN)
    P.line(list(offs), v, CJUN, 3, r=7, tips=[f'columns {4*k}–{4*k+3} mod 40: copy {a:+.3f} ± {b:.3f}' for k, (a, b) in enumerate(zip(v, s))])
    P2 = sd.Plot(640, 540, x=(0, 800), y=(0.4, 1.1), xlabel='stripe offset (µm)', ylabel='Y ±1 share', title='Charge sharing too')
    P2.xticks([(v_, str(v_)) for v_ in (0, 400, 800)]); P2.yticks([(v_, f'{v_:.1f}') for v_ in (.4, .6, .8, 1.0)])
    for cls, c, nm in (('ypad', CJUN, 'Y-pad events'), ('xpad', GREY, 'X-pad events')):
        v = [q[0] for q in j[cls]['ysh']]; s = [q[1] for q in j[cls]['ysh']]
        for x_, a, b in zip(offs, v, s): _err(P2, x_, a, b, c)
        P2.line(list(offs), v, c, 3, r=6, tips=[f'{nm}, bin {k}: Y share {a:.3f} ± {b:.3f}' for k, (a, b) in enumerate(zip(v, s))])
    amp, ea = j['ypad']['sine_copy']
    body = (sd.title('Strip or kapton gap? The resistive stripes leave no mark',
                     'det3 June-15 · true pad columns (each 64-ch connector reversed) · 39 stripes on 40 columns: the stripe '
                     'slides 20 µm per column')
            + sd.row(P.svg('stripe copy'), P2.svg('stripe share'), gap=30)
            + sd.row(sd.callout(f'Over the full beat the copy is flat: sine amplitude {amp:.3f} ± {ea:.3f} on {m:.2f}. A stripe-mediated '
                                'coupling (gold) would swing by ±0.1.', CJUN, 24),
                     sd.callout('Which one a given event hit, stripe or gap, can’t be told: the cloud (σ ≈ 450 µm) covers more than '
                                'one 800 µm period, and the peak pad fixes x only to ±390 µm.', BLUE, 24), gap=36))
    D.slide('stripes', body, _det(
        'Geometry: ESL stripes 550 µm wide with 250 µm kapton gaps, 800 µm pitch, running along y; pads 680 µm on 780 µm. '
        'Pad column k sits at stripe phase (k·780 + x0) mod 800, so 40 columns sweep the stripe once under every pad. The '
        'column phase is exact per pad, so this tests every mechanism tied to how the stripes sit over the pads (a stripe '
        'over the Y pad vs over the gap next to it, a stripe over the bus gap, a floating stripe).',
        f'The earlier “no 40-channel period” check binned DAQ channel numbers, which are reversed inside each 64-channel '
        f'connector and so mix the phases between connectors. Redone with true columns: Y-pad class sine amplitude '
        f'{amp:.3f} ± {ea:.3f}, X-pad class {j["xpad"]["sine_copy"][0]:.3f} ± {j["xpad"]["sine_copy"][1]:.3f}; Y share '
        f'{j["ypad"]["sine_ysh"][0]:.3f} ± {j["ypad"]["sine_ysh"][1]:.3f}. det3 May and det2 are flat as well. Not even '
        'the sharing shows the stripes.',
        'Why the per-event question has no handle: the charge cloud of a cosmic track has σ ≈ 450 µm (the forward model '
        'matches the data’s Y share 0.65 at σ 450), and inside the peak pad the landing point is uniform over ±390 µm, '
        'about one stripe period. The 800 µm harmonic is attenuated by exp(−2π²σ²/P²) ≈ 0.002 and again by the pad '
        'window. Every event covers stripes and gaps in about the 550:250 ratio. A per-event test needs a point source '
        'and short drift (collimated Fe-55, or X-rays converting near the mesh).',
        _src('../../2026-10-20/esl_cache.py', '../../2026-10-20/esl_bins.py', '../../2026-10-20/esl_share.py',
             '../../2026-10-20/esl_phase.py', 'fwd_out/f_off*.json')), short='Stripes')


def s_inject(D):
    inj = _j('inject_parity.json')
    if not inj: return
    g = {(o['build'], o['a']): o for o in inj}
    tb = [('det2', 'det2'), ('det3 May', 'det3 May'), ('det1', 'det1')]
    P = sd.Plot(1100, 560, x=(-0.6, 3 * 5 - 0.4 + 1.6), y=(-0.05, 0.26), ylabel='copy in the class', margin=(24, 20, 110, 100))
    P.yticks([(v, f'{v:.1f}') for v in (0, .1, .2)]); P.hline(0, MUT, None, 1.5)
    xi = 0
    for b, lab in tb:
        for a in (0, .1, .2, .3):
            o = g.get((b, a))
            if not o: continue
            for dx, key, c in ((-0.18, 'ypad', CJUN), (0.18, 'xpad', GREY)):
                v, e, n = o[key]
                P.vbar(xi + dx, max(v, -0.049), 18, c, base=0., tip=f'{b} + injected copy {a}: {key[0].upper()}-pad class {v:+.3f} ± {e:.3f}')
            P.raw(sd.T(P.X(xi), P.y0 + P.ph + 30, f'+{a:.1f}', 18, MUT))
            xi += 1
        P.raw(sd.T(P.X(xi - 2.5), P.y0 + P.ph + 62, lab, 21, INK, weight=600))
        xi += 0.8
    o = g.get(('det3 June-15', 0))
    if o:
        for dx, key, c in ((-0.18, 'ypad', CJUN), (0.18, 'xpad', GREY)):
            v, e, n = o[key]; P.vbar(xi + dx, v, 18, c, base=0., tip=f'det3 June-15 data: {key[0].upper()}-pad class {v:+.3f} ± {e:.3f}')
        P.raw(sd.T(P.X(xi), P.y0 + P.ph + 62, 'det3 June-15', 21, CJUN, weight=600))
    key = sd.legend([('Y-pad class', CJUN, 'box'), ('X-pad class', GREY, 'box')], 22)
    up = _j('uncopy_parity.json'); clo = ''
    if up:
        c = [o_ for o_ in up if o_['build'] == 'det2' and o_.get('gen') == 0.3 and abs(o_['a'] - 0.3) < 1e-9]
        if c: clo = f' Closure: inject 0.3 into det2, undo 0.3 → {c[0]["ypad"][0]:+.3f} / {c[0]["xpad"][0]:+.3f}.'
    body = (sd.title('Caution: a plain channel copy also makes a parity split',
                     'Copy-free boards with x(n+1) += a·x(n) injected into the raw windows, then the identical peak + parity analysis')
            + P.svg('inject') + key
            + sd.row(sd.callout('Inject 0.3 into det2: 0.20 for Y-pad events, 0.05 for X-pad ones. Most of the 0.22 / 0.02 contrast '
                                'is what peak selection does to any copy.', GOLD, 24),
                     sd.callout('So the split alone does not show that the effect knows the pad type, or rule out electronics. '
                                'The noise test still does.' + clo, BLUE, 24), gap=36))
    D.slide('inject', body, _det(
        'Why: when the copy lifts channel n+1 above n, the peak moves to n+1, the parity flips, and the event lands in the '
        'other class with a large negative asymmetry (its old centre is now at n−1). Y-pad events have an X central that '
        'is barely larger than its neighbours (the X pads all sit 780 µm from the charge), so they migrate most. The X-pad '
        'class collects the migrants and its mean is pulled down; the Y-pad class keeps the events where the copy stayed '
        'below the central.',
        'The earlier control (forward model, symmetric geometry → 0 in both classes) only showed that the classification '
        'creates no asymmetry from nothing. It did not test what it does to a copy that is there. This injection does.',
        'det4 and det7 (0.08 / 0.01 in the two classes) are consistent with a plain copy of a ≈ 0.15–0.20: undoing 0.20 '
        'leaves −0.01 / −0.06 and −0.02 / −0.04 (±0.02). For them the “Y-pad only” reading is not established.',
        _src('../../2026-10-20/inject_parity.py', '../../2026-10-20/uncopy_parity.py')), short='Caution')


def s_flip(D):
    up, inj, fw = _j('uncopy_parity.json'), _j('inject_parity.json'), _j('fwd_variants.json')
    if not (up and inj): return
    d3 = [o for o in up if o['build'] == 'det3 June-15' and o.get('gen') is None and o['a'] == 0][0]
    d3u = [o for o in up if o['build'] == 'det3 June-15' and o.get('gen') is None and abs(o['a'] - 0.35) < 1e-9][0]
    i2 = [o for o in inj if o['build'] == 'det2' and o['a'] == 0.3][0]
    cen = [0.3, 0.475, 0.625, 0.775, 0.925, 1.15, 1.45]
    P = sd.Plot(1080, 600, x=(0.2, 1.55), y=(-0.5, 0.35), xlabel='Y ±1 share (copy-free view): narrow cluster → wide',
                ylabel='X asymmetry  α(n+1) − α(n−1)', margin=(24, 30, 92, 110))
    P.xticks([(v, f'{v:.1f}') for v in (0.4, 0.6, 0.8, 1.0, 1.2, 1.4)]); P.yticks([(v, f'{v:+.1f}') for v in (-.4, -.2, 0, .2)])
    P.hline(0, MUT, None, 1.5)
    def bins_of(o, key, form):
        out = []
        for c_, b in zip(cen, o[key]):
            if form == 'up': out.append((c_, b[3], b[4]) if b[3] is not None else (c_, None, None))
            else: out.append((c_, (b[3] - b[4]) if b[3] is not None else None, None))
        return out
    series = [
        (bins_of(d3, 'ypad_bins', 'up'), CJUN, None, 'det3 June-15, Y-pad class'),
        (bins_of(d3, 'xpad_bins', 'up'), CJUN, '10 7', 'det3 June-15, X-pad class'),
        (bins_of(i2, 'ypad_bins', 'inj'), GREY, None, 'det2 + plain copy 0.3, Y-pad'),
        (bins_of(i2, 'xpad_bins', 'inj'), GREY, '10 7', 'det2 + plain copy 0.3, X-pad'),
    ]
    if fw:
        cb = [o for o in fw if o['kind'] == 'combo' and o['v'] == [-200, 0.1]]
        if cb:
            series.append(([(c_, v, None) for c_, v in zip(cen, cb[0]['xpad'])], GOLD, '3 6', 'forward model, X offset 200 µm + copy 0.1, X-pad'))
    for pts, c, dash, lab in series:
        xs = [p[0] for p in pts if p[1] is not None and p[1] == p[1]]; ys = [p[1] for p in pts if p[1] is not None and p[1] == p[1]]
        P.line(xs, ys, c, 3.5 if c == CJUN else 2.5, dash=dash, r=6 if c == CJUN else 4,
               tips=[f'{lab}: Ysh ≈ {x_:.2f}, {y_:+.2f}' for x_, y_ in zip(xs, ys)])
        for p in pts:
            if p[1] is not None and p[2] is not None and p[1] == p[1]: _err(P, p[0], p[1], p[2], c)
    key = sd.legend([(s_[3], s_[1], 'line') for s_ in series], 20)
    nb = d3u['xpad_bins'][0][3]
    body = (sd.title('det3 is not a plain copy: its X-pad events flip sign with cluster width',
                     'det3 June-15 by Y share, against a copy-free board with a plain copy injected')
            + sd.row(sd.col(P.svg('flip'), key, gap=6),
                     sd.col(sd.callout('Narrow clusters filed as X-pad: −0.40. A plain copy gives ≈ 0 there at any a, and undoing one '
                                       f'(a = 0.35) leaves {nb:+.2f}.', CJUN, 23),
                            sd.callout('These are Y-pad charges whose X peak sits one column toward n+1: 13 % of narrow clusters on '
                                       'det3, 4–5 % on det2/det1 even with a 0.3 copy.', BLUE, 23),
                            sd.callout('Pattern fits better as the X view’s sensitivity shifted ~200 µm toward n+1 relative to the Y '
                                       'view (gold), plus a small copy.', GOLD, 23), gap=14, w=640), gap=30))
    D.slide('flip', body, _det(
        'Inverse test (<code>uncopy_parity.py</code>): undo a single copy a on the raw windows, re-find the peak, redo the '
        'split. On injected controls this closes exactly (both classes and every bin → 0 at the injected a). On det3 '
        'June-15, a = 0.33 flattens the Y-pad class in every bin, but the narrow X-pad bins stay at −0.43, −0.32, −0.19 '
        'for every a from 0 to 0.35, and the wide X-pad bins go negative. No single copy describes both.',
        'Two-copy fit (<code>fit_two_copy.py</code>, separate a for Y-pad and X-pad charge, templates det2/det3 May/det1 '
        'reweighted to det3’s Y-share spectrum): best a_Y 0.30–0.35, a_X 0.05–0.15, but χ² 230–380 for 21 bins. The '
        'narrow X-pad bins stay at −0.07 … +0.05 against −0.40, whatever a_X is.',
        'Migration (<code>migrate.py</code>): among clusters with Y share < 0.4, the X peak sits on the X-pad parity in '
        '13.4 % of det3 June-15 events, 77 % of them with the larger neighbour at n−1 (the true Y-pad column). det2 and det1 '
        'give 4.0–4.5 % at a = 0–0.3 and need a ≈ 0.8 to reach 13 %. Caveat: det3 May (same PCB design) has a 20 % '
        'baseline here, so this number alone is not decisive. The −0.40 bin is.',
        'Forward model (<code>fwd_variants.py</code>, solved weighting fields): a plain copy gives X-pad bins ≈ 0 → +0.2; '
        'an X-sensitivity offset of 200 µm toward n+1 gives −0.39 → +0.21, as in the data; offset + copy 0.1 also lifts '
        'the Y-pad class (+0.25 → −0.07, data +0.26 → +0.09). The model’s pad-type purity vs Y share is lower than the data’s '
        '(0.53–0.65 vs 0.87 in the narrowest bin), so only the shape is compared, not a χ².',
        _src('../../2026-10-20/uncopy_parity.py', '../../2026-10-20/inject_parity.py', '../../2026-10-20/fit_two_copy.py',
             '../../2026-10-20/migrate.py', '../../2026-10-20/fwd_variants.py')), short='Not a copy')


def s_inside(D):
    yp = _j('ypad_probe.json')
    if not yp: return
    j = yp['det3 June-15']
    P = sd.Plot(640, 520, x=(-1.3, 1.3), y=(0, 0.35), xlabel='y position inside the Y pad (η_y, from the Y view)',
                ylabel='copy, Y-pad events', title='Where along y')
    P.xticks([(v, f'{v:+.0f}' if v else '0') for v in (-1, 0, 1)]); P.yticks([(v, f'{v:.1f}') for v in (0, .1, .2, .3)])
    pts = [b for b in j['vs_eta_y'] if b[1] is not None]
    for b in pts: _err(P, b[0], b[1], b[2], CJUN)
    P.line([b[0] for b in pts], [b[1] for b in pts], CJUN, 3, r=7, tips=[f'η_y ≈ {b[0]:+.2f}: {b[1]:+.3f} ± {b[2]:.3f} ({b[3]} ev)' for b in pts])
    sh = j['shape']; t = np.arange(len(sh['xc'])) * 60
    P2 = sd.Plot(760, 520, x=(0, 1500), y=(-0.05, 1.05), xlabel='time (ns)', ylabel='mean pulse / its own peak', title='When in time')
    P2.xticks([(v, str(v)) for v in (0, 500, 1000, 1500)]); P2.yticks([(v, f'{v:.1f}') for v in (0, .5, 1)])
    for k, c, lab in (('yc', BLUE, 'Y central (the charge’s own pad)'), ('xc', INK, 'X central'), ('x1', GREY, 'X ±1, own signal'),
                      ('ex', CJUN, 'excess on n+1')):
        w = np.array(sh[k]); w = w / w.max()
        P2.line(list(t), list(w), c, 3.5 if k == 'ex' else 2.5, markers=False, tip=f'{lab}: peak at {np.argmax(w) * 60} ns')
    key = sd.legend([('Y central', BLUE, 'line'), ('X central', INK, 'line'), ('X ±1 own', GREY, 'line'), ('excess on n+1', CJUN, 'line')], 20)
    f = j['fit_Xc+Yc']
    body = (sd.title('Inside the Y pad: flat along y, prompt in time',
                     'det3 June-15, Y-pad events · y position from the copy-free Y view · waveforms normalised to their own peak')
            + sd.row(P.svg('eta'), sd.col(P2.svg('shape'), key, gap=6), gap=30)
            + sd.row(sd.callout('Same copy at both y edges and the centre: not a spot (via, stub, defect) at one side of the Y pad.', CJUN, 24),
                     sd.callout(f'The excess peaks with the Y central, ~180 ns before the neighbours’ own signal: it rides on the '
                                f'directly collected charge. Fit {f["coef"][0]:.2f}·X_c + {f["coef"][1]:.2f}·Y_c leaves {100*f["rms"]:.0f} %.', BLUE, 24), gap=36))
    D.slide('inside', body, _det(
        'Left: η_y = α_y(+1) − α_y(−1) places the charge above or below the Y pad’s centre (Y has no copy, so this is '
        'clean). The copy is 0.21–0.25 in every bin. Copy vs A_x/A_y falls from 0.27 to 0.17, a term that does not scale '
        'with the X central.',
        f'Right: peak sample of the mean waveforms: Y central {sh["peak"]["yc"]*60:.0f} ns, excess {sh["peak"]["ex"]*60:.0f} ns, '
        f'X central {sh["peak"]["xc"]*60:.0f} ns, X ±1 own {sh["peak"]["x1"]*60:.0f} ns. One-template fits leave X_c 13 %, '
        'Y_c 31 %; the derivative of Y_c adds nothing. So the extra charge on n+1 is induced while the avalanche charge '
        'arrives, not later by spreading.',
        _src('../../2026-10-20/ypad_probe.py', 'ypad_probe.json')), short='Inside')


def s_readout(D):
    pr, mx = _j('m1_profile.json'), _j('mix_reach2.json')
    if not (pr and mx): return
    j = pr['det3 June-15']
    P = sd.Plot(1000, 520, x=(0, 512), y=(0, 0.4), xlabel='pad column (true order) · dashed: 64-channel connector edges',
                ylabel='copy, Y-pad events', title='Across the board: smooth through every connector edge')
    P.xticks([(v, str(v)) for v in range(0, 513, 128)]); P.yticks([(v, f'{v:.1f}') for v in (0, .1, .2, .3, .4)])
    for k in range(1, 8): P.vline(64 * k, MUT, '6 6', 1.5)
    pts = [p for p in j['profile_ypad'] if p[1] is not None]
    for p in pts: _err(P, p[0], p[1], p[2], CJUN)
    P.line([p[0] for p in pts], [p[1] for p in pts], CJUN, 3, r=6, tips=[f'columns {p[0]-8}–{p[0]+7}: {p[1]:+.3f} ± {p[2]:.3f} ({p[3]} ev)' for p in pts])
    jr = np.sqrt(np.mean([s[1] ** 2 for s in j['jumps_at_connector']])); ir = np.sqrt(np.mean([s[1] ** 2 for s in j['steps_inside']]))
    yc = [v[1] for v in j['vs_ychip_ypad'] if v[1] is not None]
    best = {k: v for k, v in mx.items() if k != 'data'}
    rows = [[k, f'{v["chi2_zero"]:.0f}', ', '.join(f'{m:+.1f}' for m in v['best']), f'{v["chi2"]:.0f} / {v["ndf"]}',
             f'{v["model"][2][0]:+.2f}  {v["model"][7][0]:+.2f}'] for k, v in best.items()]
    rows.append(['<b>det3 June-15 data</b>', '', '', '', f'<b>{mx["data"][2][0]:+.2f}  {mx["data"][7][0]:+.2f}</b>'])
    body = (sd.title('Not the one-sided readout (M1 cards, connectors, cables)',
                     'Anything after the pads is a fixed channel-to-channel mixing: it cannot know where the hit is, and it '
                     'steps where the cards do')
            + sd.row(P.svg('readout profile'),
                     sd.col(sd.table(['template', 'χ² no mix', 'best m(−1,+1,−2,+2)', 'χ² best', 'narrow X-pad'], rows, size=19),
                            sd.callout('No linear mixing of neighbours, of any sign, makes the narrow X-pad bins negative while '
                                       'keeping the shares: −0.10 … +0.23 against −0.40 / −0.29.', CJUN, 22), gap=14, w=760), gap=30)
            + sd.row(sd.callout(f'Steps at connector edges rms {jr:.3f}, steps between bins inside a connector {ir:.3f}: no card '
                                'boundaries in the profile.', BLUE, 23),
                     sd.callout(f'Along the X line (distance from the readout edge, by Y chip): {min(yc):.2f}–{max(yc):.2f}, '
                                'flat. Coupling along the line or in the fan-out would grow with it.', GOLD, 23), gap=36))
    D.slide('readout', body, _det(
        'Every X line is read out at the same board edge, so the readout is not handed in x: line n and its two neighbours are '
        'read the same way. A passive coupling between neighbouring traces (fan-out, M1 card, cable) is reciprocal: it '
        'raises both neighbours and creates no left/right asymmetry, unless the routing order differs from the pad order. Any '
        'such effect, and any electronics effect, is a fixed linear map from channels to channels.',
        'Test (<code>mix_fit.py</code>, <code>mix_reach2.py</code>): inject x′(n) = x(n) + Σ m_k x(n+k), k = ±1, ±2, '
        'negative values allowed (a derivative-like mixing that the noise test is blind to at first order), into copy-free '
        'boards, run det3’s full analysis, and score class asymmetries and ±1 sums in seven Y-share bins (the class fraction is '
        'left out: it differs between boards with no mixing). Grid m₋₁ 0–1, m₊₁ −0.5–0.5, m±₂ −0.2–0.2. The best mixing still '
        'leaves χ² 410–750 / 28, and the narrow X-pad bins never go below −0.19. The mixings that get there (m₋₁ ≈ 1) wreck '
        'the shares.',
        'The reason is structural: det3 behaves as if the X view sensed the charge displaced relative to the Y pads. That '
        'is a change in the X–Y registration, and nothing acting on the X channels alone can make it.',
        'The DAQ-neighbour pairs across a connector edge are tested directly on the next slide (re-extraction with chip-edge '
        'peaks kept). det4/det7 are consistent with a plain copy, so for them the readout question rests on that slide, the '
        'smooth profile and the noise test.',
        _src('../../2026-10-20/m1_profile.py', '../../2026-10-20/mix_fit.py', '../../2026-10-20/mix_reach2.py')), short='Readout?')


def s_sessions(D):
    f = _j('fwd_session_fit.json')
    if not f: return
    rows = []
    for lab, o in f.items():
        rows.append([lab, f'{o["n"]}', f'{o["median_ysh"]:.2f}', f'{o["offset"]:+.0f} µm', f'{o["copy"]:.2f}',
                     f'{o["chi2"]:.0f}', f'{o["pure_copy"]["chi2"]:.0f} (a {o["pure_copy"]["a"]:.2f})'])
    pts = [(o['median_ysh'], -o['offset'], lab) for lab, o in f.items() if not lab.startswith('all')]
    P = sd.Plot(640, 480, x=(0.4, 1.1), y=(0, 300), xlabel='session median Y share (wider →)', ylabel='fitted X offset (µm)',
                title='Offset follows the cluster width')
    P.xticks([(v, f'{v:.1f}') for v in (.5, .7, .9, 1.1)]); P.yticks([(v, str(v)) for v in (0, 100, 200, 300)])
    P.points([p[0] for p in pts], [p[1] for p in pts], CJUN, r=9, tips=[f'{p[2]}: {p[1]:.0f} µm at Ysh {p[0]:.2f}' for p in pts])
    body = (sd.title('Per session: an offset plus a small copy fits every night far better than a copy',
                     'Forward model (solved weighting fields) fitted to the 12 class asymmetries in Y-share sextiles')
            + sd.row(sd.table(['session', 'events', 'median Ysh', 'X offset', 'copy', 'χ² /12', 'χ² pure copy'], rows, size=21),
                     P.svg('sessions'), gap=30)
            + sd.callout('The copy stays at 0.10–0.15; the offset carries the night-to-night and drift-field change, and it '
                         'grows as the clusters narrow. Whether one fixed geometry does that through the width alone needs '
                         'a model whose Y-share scale matches the data (its medians are 0.37–0.49).', BLUE, 24))
    D.slide('sessions', body, _det(
        'Bins are Y-share quantiles in data and model separately, because the model’s Y-share scale and pad-type purity '
        'differ from the data’s (the static weighting fields have no resistive spreading and no shaping). Grid: offset '
        '0 to −300 µm, copy 0–0.3, five cloud-width ranges; a 0.02 floor is added to each bin’s error.',
        'Same night, 6-27: 1000 V drift gives −100 µm, 100–500 V gives −250 µm, reversibly. With quantile bins the absolute '
        'width is not used, so a fixed sensitivity extension whose effect grows for narrower clouds would show up here as an '
        'offset that grows with narrowness. That is what is seen; it is not yet evidence of a field-dependent geometry.',
        _src('../../2026-10-20/fwd_session_fit.py', 'fwd_session_fit.json')), short='Sessions')


def s_edge(D):
    e = _j('edge_test.json')
    if not e: return
    runs = [('e_d3jun23_x', 'det3 June-15, 6-23 (FEU 3)', CJUN), ('e_d3wk_x', 'det3 June-15, 6-28 (FEU 7)', ORANGE),
            ('e_d2_x', 'det2', GREY), ('e_d3may_x', 'det3 May', GREEN), ('e_d4_x', 'det4', PURPLE)]
    cats = [('col%64==1..3', 'inside, low end'), ('col%64==0 (col-1 on OTHER connector)', 'lowest column:\ncol−1 on the other connector'),
            ('col%64==60..62', 'inside, high end')]
    P = sd.Plot(1100, 540, x=(-0.6, 2.6), y=(-0.1, 0.3), ylabel='α(col−1) − α(col+1), physical neighbours', margin=(24, 20, 120, 110))
    P.yticks([(v, f'{v:+.1f}') for v in (-.1, 0, .1, .2, .3)]); P.hline(0, MUT, None, 1.5)
    for i, (k, lab) in enumerate(cats):
        for j, (tag, nm, c) in enumerate(runs):
            r = e.get(tag, {}).get(k)
            if not r: continue
            x = i + (j - 2) * 0.16; v, s = r['phys']
            P.vbar(x, v, 22, c, base=0., tip=f'{nm}, {lab.replace(chr(10), " ")}: {v:+.3f} ± {s:.3f} ({r["n"]} ev)')
            _err(P, x, v, s, INK)
        for li, line in enumerate(lab.split('\n')):
            P.raw(sd.T(P.X(i), P.y0 + P.ph + 34 + 26 * li, line, 21, INK, weight=600 if i == 1 else 400))
    key = sd.legend([(nm, c, 'box') for _, nm, c in runs], 20)
    j = e['e_d3jun23_x']['col%64==0 (col-1 on OTHER connector)']
    body = (sd.title('The copy crosses connector boundaries with the pads, not with the electronics',
                     'Peaks at the edge of a 64-channel connector · re-extracted in true column order (raw_asym.py --edge)')
            + P.svg('edge') + key
            + sd.row(sd.callout(f'det3 6-23: when col−1 is on the other connector, cable and DREAM chip, it still gets the excess: '
                                f'{j["phys"][0]:+.2f} ± {j["phys"][1]:.2f}, as inside a connector.', CJUN, 23),
                     sd.callout(f'Its DAQ n+1 there (same FEU, next chip, 127 columns away) gets α = {j["a_daqp1"]:.2f}. '
                                'A copy made after the pads would land on DAQ n+1.', BLUE, 23), gap=36))
    D.slide('edge', body, _det(
        'Wiring: inside each 64-channel connector the DAQ order runs opposite to the columns, so DAQ n+1 is the next-lower '
        'column, except at the connector’s lowest column, where DAQ n+1 is the far (highest) column of the next connector '
        'and the physical neighbour col−1 is on the previous connector. The two orders separate there.',
        'Extraction: peak search and the isolation cut (all channels over 5σ within ±3 columns) in column order, peaks with '
        'col mod 64 in 0..3 or 60..63; both the column-order and the DAQ-order 9-channel windows are stored. The usual good-pulse '
        'cut and hot-channel veto.',
        'The 6-28 weekend run (FEU 7, copy ~0.09 inside) reads +0.03 ± 0.03 at the lowest column, 2σ below its inside value: '
        'not decisive alone. At col mod 64 = 63 (col+1 on the other connector) every board, controls included, reads '
        '+0.04 … +0.14: the cross-connector neighbour comes out low on all of them, a chip-to-chip or edge-pad effect common '
        'to the readout, which is why the clean comparison is the col mod 64 = 0 case.',
        _src('asym_deck/raw_asym.py --edge', '../../2026-10-20/edge_test.py', 'edge_test.json')), short='Connector edge')


BOARDS = [  # (build, run_config description, Y-pad copy (asym_more), stripe orientation from data)
    ('det3 June-15', '“Bulked by Stephan June 15” (beam: mx17_A)', 'Y slow (+420 ns at ±2)'),
    ('det4', '“Stephan batch of 3, June 22. Was board C. A few bubbles … just no caps”', 'Y slow'),
    ('det7', '“Stephan batch of 3, June 22. Was board B. One or two bubbles”', 'Y slow (+240 ns)'),
    ('det6', '“Stephan June 24 (?). Was board D. Relaminated … column of waves”', 'Y slow (+240 ns)'),
    ('det2', '“Bulked by Arnaud June 12. Giant pillars on parts” (beam: mx17_B)', 'Y slow (+240 ns)'),
    ('det3 May', '“Bulked 5-5-26. Strip resist vertical. ESL no silver paste”', 'Y slower than X'),
    ('det1', 'Jan build (no description); April re-bulk “Strip resist horizontal … strips broken”', 'X slower than Y (reversed)'),
]


def s_boards(D, more_json):
    m = more_json; fp = _j('pcb_fingerprint.json') or {}
    cp = {b['build']: b['ypad'] for b in m['boards']} if m else {}
    def dead(bu):
        ks = [k for k in fp if k.startswith(bu + '|')]
        bo = {}
        for k in ks:
            _, v, feu = k.split('|'); bo.setdefault(v, []).append(set(fp[k]))
        out = []
        for v, ss in sorted(bo.items()):
            com = set.intersection(*ss) if len(ss) > 1 else ss[0]
            if com: out.append(f'{v} {",".join(map(str, sorted(com)))}' + ('' if len(ss) > 1 else ' (1 slot)'))
        return '; '.join(out) or '—'
    rows = [[bu, desc, f'{cp[bu][0]:+.2f}' if bu in cp else '—', ori, dead(bu)] for bu, desc, ori in BOARDS]
    body = (sd.title('The boards: what the run configs say, and what the copy does',
                     'Descriptions verbatim from run_config.json · Y-pad-class copy · stripe direction read from the slow sharing · '
                     'dead lines (true order) on the board’s own FEU slot(s)')
            + sd.table(['board', 'run_config.json', 'Y-pad copy', 'stripes (from data)', 'dead lines'], rows, size=19)
            + sd.row(sd.callout('Stephan: nothing different about the June bulks. The copy is per board: 0.22 on det3 June-15, '
                                '0.07–0.08 on the June-22 batch, ~0 on det2, det6, det3 May and det1.', CJUN, 23),
                     sd.callout('No shared dead line between det3 May and det3 June-15 (its Y row 23 is dead on both of its FEU '
                                'slots, fine in May): no sign they are the same PCB.', BLUE, 23), gap=36))
    D.slide('boards', body, _det(
        'Bulk dates for the July/August beam: mx17_A = det3, “Bulked by Stephan June 15”; mx17_B = det2, “Bulked by Arnaud '
        'June 12. Giant pillars on parts of the detector.” (all 162 beam run configs). mx17_C = det6, mx17_D = det7. The May det3 '
        '(bench 5-5 to 5-6, and “Det3 from May beam test” in the 6-6 to 6-17 runs) was “Bulked 5-5-26”.',
        'Stripe direction (<code>stripe_orient.py</code>): charge spreads along the stripes, so the view measuring along them '
        'gets the slow, late sharing; at ±2, the late content is larger in Y than in X on every board but det1, with Y peaking '
        '240–420 ns late on the June boards. det3 June-15 looks like the others: orientation is not what singles it out.',
        'Dead lines (<code>pcb_fingerprint.py</code>): a line that is the peak < 15 % as often as its neighbours. X 112 is dead '
        'on every board read through FEU 3 (an electronics channel), so only lines that follow a board across slots, or a '
        'board’s only slot, are listed. A PCB-level break would survive a re-bulk; bulk-level ones would not.',
        _src('run_config.json (local + EOS)', '../../2026-10-20/stripe_orient.py', '../../2026-10-20/pcb_fingerprint.py')), short='Boards')
