"""Rounds 4-8 slides for the det3 X-asymmetry deck (2026-10-21), imported by ../2026-10-19/asym_deck/make_asym_deck.py.
Inputs (read, never typed in): ../2026-10-20/r5/{corr_maps, sim_delta, kernel_delta_prod, kernel_delta_gap1560,
kernel_delta_width, paired_odd_*, sim_late2, timing_width}.json and ./{reco_resid_map, reco_correct_test}.json.
A slide whose input is missing is skipped."""
import os, json, numpy as np
import slidedoc as sd
from slidedoc import INK, MUT, BLUE, ORANGE, RED, GOLD, PURPLE, GREY, GREEN, DINK, DMUT, DGREEN, DRED, DBLUE
HERE = os.path.dirname(os.path.abspath(__file__))
R5 = os.path.join(HERE, '..', '2026-10-20', 'r5')
CJUN = RED


def _j(name, d=R5):
    p = os.path.join(d, name)
    return json.load(open(p)) if os.path.exists(p) else None


def _det(*paras):
    return ''.join(f'<p>{p}</p>' for p in paras)


def _src(*items):
    return '<p>Sources: ' + ', '.join(f'<code>{i}</code>' for i in items) + '</p>'


def _err(P, x, v, e, c=INK, w=2):
    P.raw(sd.line(P.X(x), P.Y(v - e), P.X(x), P.Y(v + e), c, w))


# --------------------------------------------------------------------------- maps
def s_fields(D):
    cm = _j('corr_maps.json')
    if not cm: return
    VM = 160.
    panels = []
    for b, lab in (('det3 June-15', 'det3 June-15 build'), ('det2', 'det2 (June build)'), ('det3 May', 'det3 May build (control)')):
        o = cm['builds'][b]; G = o['dx_um']; N = o['n']; E = o['err_um']
        P = sd.Plot(520, 560, x=(0, 8), y=(0, 8), xlabel='x (64-strip blocks, 49.9 mm)', ylabel='y (blocks)' if b == 'det3 June-15' else '',
                    title=lab, margin=(70, 20, 120, 100) if b == 'det3 June-15' else (70, 20, 120, 40))
        P.xticks([(v, str(v)) for v in (0, 2, 4, 6, 8)]); P.yticks([(v, str(v)) for v in (0, 2, 4, 6, 8)] if b == 'det3 June-15' else [])
        cells = []
        for i in range(8):
            for j in range(8):
                v = G[i][j]
                tp = (f'{lab}, cell x {i} y {j}\nX-centroid bias {v:+.0f} ± {E[i][j]:.0f} µm, n {N[i][j]}' if v is not None
                      else f'{lab}, cell x {i} y {j}: n {N[i][j]} (< 25, not shown)')
                cells.append((i, i + 1, j, j + 1, v, tp))
        P.cells(cells, -VM, VM, sd.DIVERGE, missing='#e9e9ec', gap=1.5)
        P.raw(sd.T(P.x0 + P.pw / 2, P.y0 + P.ph + 104, f'mean {o["mean_um"]:+.0f} µm · cell rms {o["rms_um"]:.0f} µm', 22, INK, 'middle', 600))
        panels.append(P.svg(lab))
    cb = sd.colorbar(150, 560, -VM, VM, [(-150, '−150'), (-75, '−75'), (0, '0'), (75, '+75'), (150, '+150')],
                     label='bias (µm)', cmap=sd.DIVERGE)
    d3 = cm['builds']['det3 June-15']; d2 = cm['builds']['det2']; ss = cm.get('det3_sessions', {})
    stxt = ', '.join(f'{k} ×{v["scale"]:.2f} (χ² {v["chi2"]:.0f}/{v["ncells"] - 1})' for k, v in ss.items())
    body = (sd.title(f'Every June board carries a smooth X field; det2’s bands cancel, det3’s does not ({d3["mean_um"]:+.0f} µm)',
                     'Where the 3-strip X centroid is pulled by the n±1 asymmetry, per 8×8 cell of the board, matched X/Y cosmics')
            + sd.row(*panels, cb, gap=28))
    D.slide('fields', body, _det(
        'dx = −(α₊₁ − α₋₁)/(1 + α₊₁ + α₋₁) × 780 µm along the true column axis (= +x in the wft frame: pos = 0.78 mm × '
        'true column on both views), α the in-phase content of the DAQ ±1 channels. Minus the symmetric reference '
        f'(det3 May + det1, {cm["ref_um"]:+.1f} µm). Cells with fewer than 25 events are grey.',
        'Fixed geometry, variable size: one det3 map times one scale per session fits every session: ' + stxt + '. '
        'The size follows the cloud width (gas, drift field).',
        f'det2: mean {d2["mean_um"]:+.1f} µm but ±70 µm locally, in bands along y; that is why det2 looked like a clean control. '
        'det4/det6/det7 show fields too but have too few bench events for a cell map (13, 26 and 9 cells).',
        _src('r5/corr_maps.py', 'r5/corr_maps.json', 'r5/map_fit.py (the registration-field fits)')),
        foot='Raw survey, all bench sessions per board; colour clipped at ±160 µm', short='Fields')


# --------------------------------------------------------------------------- sim vs stripe offset
def s_simdelta(D):
    sd_ = _j('sim_delta.json')
    if not sd_: return
    rows = [r for r in sd_['col'] if r]
    xs = [r['d'] for r in rows]
    P = sd.Plot(1080, 640, x=(-400, 400), y=(-0.3, 0.3), xlabel='stripe offset δ under the peak column (µm)',
                ylabel='a1 = α(n+1) − α(n−1), sim', title='v5a full sim, X-only events, D = 0 arms')
    P.xticks([(v, f'{v:+d}' if v else '0') for v in (-400, -200, 0, 200, 400)])
    P.yticks([(v, f'{v:+.1f}' if v else '0') for v in (-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)])
    P.band([-400, 400], [0.14, 0.14], [0.29, 0.29], CJUN, 0.10, tip='Data: det3 June-15 cell / session a1, +0.14 … +0.29')
    P.band([-400, 400], [-0.29, -0.29], [-0.14, -0.14], CJUN, 0.10, tip='Data: det2 negative bands, down to −0.29')
    P.hline(0, MUT, None, 1.5)
    tips = [f'δ {r["d"]:+.0f} µm: a1 {r["a1"]:+.3f} ± {r["e"]:.3f}, n {r["n"]}\nnarrow half {r["narrow"]:+.3f}, wide half {r["wide"]:+.3f}' for r in rows]
    for r in rows: _err(P, r['d'], r['a1'], r['e'], BLUE)
    P.line(xs, [r['a1'] for r in rows], BLUE, 3, tips=tips)
    P.line(xs, [r['narrow'] for r in rows], PURPLE, 2.5, dash='8 6', markers=False, tip='Narrow-cloud half (X sharing below the median)')
    P.text(-380, 0.24, 'data, det3 June-15', 22, CJUN, weight=600)
    P.text(-380, -0.27, 'data, det2 negative bands', 22, CJUN, weight=600)
    P.text(40, -0.17, 'odd in δ, both signs: the shape is right,', 24, INK, weight=600)
    P.text(40, -0.21, 'the size is ~1/3 of the data', 24, INK, weight=600)
    body = (sd.title('Stripe registration makes the asymmetry in the sim, with both signs, at ~1/3 of the size',
                     'The sim’s 800 µm stripes slide 20 µm per 780 µm column, so each column has its own offset δ')
            + sd.row(P.svg('sim a1 vs delta'),
                     sd.col(sd.legend([('sim, all X-only events', BLUE), ('sim, narrow-cloud half', PURPLE, 'dash'),
                                       ('data range', CJUN, 'box')], 24),
                            sd.callout('A stripe sitting off the pad column’s centre registers the landed charge toward one '
                                       'neighbour. The data then need a print whose offset drifts across the board.', BLUE, 25),
                            sd.callout('That needs a 780 µm-commensurate stripe pitch, against the 800 µm quote: '
                                       'the mask will tell.', GOLD, 25), w=520), gap=48))
    D.slide('simdelta', body, _det(
        'δ = signed distance from the peak column’s centre to the nearest stripe centre in the sim geometry (stripes at k·800 µm, '
        'columns at PAD_ORIGIN + c·780 µm). X-only events of the xo0 arms (th00 + thy*, seeds 5–16).',
        'At δ = +150…+350 the sim also reproduces det3’s pad-class pattern (narrow-cloud X-pad events negative, wide positive); '
        'δ = −350…−150 gives the mirror, like det2’s negative bands. A rotation alone cannot make this: it only shortens the moiré.',
        _src('r5/sim_delta.py', 'r5/sim_delta.json', 'FINDINGS.md sixth round')), short='Sim vs δ')


# --------------------------------------------------------------------------- kernel time build-up
def s_kernel(D):
    kp = _j('kernel_delta_prod.json'); kg = _j('kernel_delta_gap1560.json'); kw = _j('kernel_delta_width.json')
    if not (kp and kg and kw): return
    pr = kp['gx_slice_rho5M.npz']
    P = sd.Plot(800, 600, x=(0, 1000), y=(0, 0.06), xlabel='time after the charge lands (ns)', ylabel='odd asymmetry |⟨sgn δ · A⟩|',
                title='Production kernel (5 MΩ/sq)')
    P.xticks([(v, str(v)) for v in (0, 250, 500, 750, 1000)]); P.yticks([(v, f'{v:.2f}') for v in (0, .02, .04, .06)])
    for sg, c in (('250', BLUE), ('450', GREY)):
        t = pr[sg]['t_ns']; v = [abs(a) for a in pr[sg]['odd_amp']]
        P.line(t, v, c, 3.5, tips=[f'σ {sg} µm, t {a:.0f} ns: {b:.4f}' for a, b in zip(t, v)])
    P.text(300, 0.040, 'zero at t = 0: built by sheet transport', 23, INK, weight=600)
    # variants at 300 ns and 1 us, sigma 250
    def odd(dct, key, tt):
        r = dct[key]['250']; i = int(np.argmin(np.abs(np.array(r['t_ns']) - tt))); return abs(r['odd_amp'][i])
    V = [('f = 0', kg, 'gx_w2_rho5M_gap0_nx1560.npz', GREY), ('gap leak f = 0.03', kg, 'gx_w2_rho5M_gap0.03_nx1560.npz', GOLD),
         ('gap leak f = 0.1', kg, 'gx_w2_rho5M_gap0.1_nx1560.npz', ORANGE), ('gap leak f = 0.3', kg, 'gx_w2_rho5M_gap0.3_nx1560.npz', RED),
         ('stripes 420 µm', kw, 'gx_w2_rho5M_gap0_w410_nx1560.npz', PURPLE), ('stripes 700 µm', kw, 'gx_w2_rho5M_gap0_w690_nx1560.npz', GREEN)]
    Q = sd.Plot(800, 600, x=(0, len(V)), y=(0, 0.06), xlabel='', ylabel='odd asymmetry, σ 250 µm', title='Kernel variants (nx = 1560)')
    Q.yticks([(v, f'{v:.2f}') for v in (0, .02, .04, .06)])
    for i, (lab, dct, key, c) in enumerate(V):
        a3, a1 = odd(dct, key, 300), odd(dct, key, 1000)
        Q.vbar(i + 0.32, a3, 44, c, tip=f'{lab}: {a3:.4f} at 300 ns')
        Q.vbar(i + 0.68, a1, 44, c, opacity=0.45, tip=f'{lab}: {a1:.4f} at 1 µs')
        Q.text(i + 0.5, -0.004, lab.replace('gap leak ', 'leak ').replace('stripes ', ''), 18, INK, 'middle')
    body = (sd.title('In the kernel the asymmetry is sheet transport: a leak doubles it late, wider stripes early',
                     'Odd part of A = (G[c+1] − G[c−1])/G[c] over the 40 column phases, gaussian clouds, peak-column selection')
            + sd.row(P.svg('kernel time'), sd.col(Q.svg('kernel variants'),
                     sd.legend([('σ 250 µm', BLUE), ('σ 450 µm', GREY), ('bar: at 300 ns (solid) / 1 µs (pale)', INK, 'box')], 22)), gap=44))
    D.slide('kernel', body, _det(
        'Production slice: G_X of greens_comb_w2_rho5M_dk50um_g19um.npz at 9 times (sliced on lxplus). Variants: X-box W2 solves '
        'with conducting gaps (σ_gap = f·σ_strip, <code>--gap-frac</code>) or another strip width (<code>--esl-width-um</code>; '
        '410/690 µm are 420/700 µm effective on the 20 µm grid), insulating otherwise. The nx = 1560 f = 0 control reproduces '
        'the production slice (0.024 at 300 ns).',
        'Wide clouds (σ 450 µm) stay below 0.004 in every variant: the kernel effect is a narrow-cloud one, as in the data. '
        'Lower ρ_s only makes it arrive sooner (saturated by 92 ns at 0.5 MΩ/sq), never larger. Too much leak (f = 0.3) '
        'makes the sheet uniform and kills it.',
        _src('r5/kernel_delta.py', 'response/solver/w2_production.py (family --nx --gap-frac --esl-width-um, combine-x)',
             '2026-10-21/run_gx_gap.sh', 'condor 4520685, 4520736')), foot='X kernels only; the full-sim test of each lever is the next slide', short='Kernel')


# --------------------------------------------------------------------------- full-sim levers
def s_levers(D):
    R = {}
    for f in ('paired_odd_gapleak.json', 'paired_odd_width.json', 'paired_odd_w690f1.json'):
        j = _j(f)
        if j: R.update(j)
    fu = _j('paired_odd_funnel.json') or {}
    R.update({k: v for k, v in fu.items()})
    for f in ('paired_odd_r2_vs_g0.json', 'paired_odd_r05f1_vs_g0.json'):
        j = _j(f)
        if j: R.update({k + '@g0': v for k, v in j.items()})
    if not R: return
    rows_def = [('fu0.5', 'funnel 50 % of gap electrons', GOLD), ('fu1', 'funnel all gap electrons', ORANGE),
                ('glg0.03', 'gap leak f = 0.03', GREY), ('glg0.1', 'gap leak f = 0.1', GREY),
                ('glw410', 'stripes 420 µm', PURPLE), ('glw690', 'stripes 700 µm', GREEN),
                ('glg0.1f1', 'leak 0.1 + funnel all', ORANGE), ('glw690f1', 'stripes 700 µm + funnel all', RED),
                ('glr2w690f1@g0', '… at 2 MΩ/sq (Y fails)', '#c9a5ab'), ('glr05w690f1@g0', '… at 0.5 MΩ/sq (Y fails)', '#c9a5ab')]
    P = sd.Plot(1180, 680, x=(0, 3.2), y=(0, len(rows_def)), xlabel='odd asymmetry, variant / its paired baseline (narrow-cloud half)',
                ylabel='', title='', margin=(30, 30, 90, 340))
    P.xticks([(v, f'×{v:g}') for v in (0, 0.5, 1, 1.5, 2, 2.5, 3)])
    P.band([2.5, 3.2], [0, 0], [len(rows_def), len(rows_def)], CJUN, 0.10, tip='What the data need: ×2.5–5 over the sim baseline')
    P.vline(1, MUT, '6 6')
    for i, (k, lab, c) in enumerate(rows_def):
        y = len(rows_def) - i - 0.5
        P.raw(sd.T(P.x0 - 16, P.Y(y) + 8, lab, 23, INK, 'end'))
        if k not in R or 'nar' not in R[k]:
            P.raw(sd.T(P.X(0.05), P.Y(y) + 8, 'pending', 22, MUT, 'start')); continue
        o, a = R[k]['nar'], R[k]['all']
        rel = o['err'] / max(o['base'], 1e-9)
        tp = (f'{lab}: narrow x{o["ratio"]:.2f} (diff {o["diff"]:+.4f} ± {o["err"]:.4f}), all x{a["ratio"]:.2f}; '
              f'{o["pairs"]} paired (angle, seed) sets')
        P.raw(f'<rect x="{P.X(0):.1f}" y="{P.Y(y) - 15:.1f}" width="{P.X(o["ratio"]) - P.X(0):.1f}" height="30" fill="{c}" rx="4"'
              + sd.tipattr(tp) + '/>')
        P.raw(sd.line(P.X(o['ratio'] - rel), P.Y(y), P.X(o['ratio'] + rel), P.Y(y), INK, 2.5))
        P.raw(sd.T(P.X(o['ratio']) + 14, P.Y(y) + 8, f'×{o["ratio"]:.2f}', 22, INK, 'start', 600))
    P.text(2.56, len(rows_def) - 0.55, 'data need', 22, CJUN, weight=600); P.text(2.56, len(rows_def) - 0.95, '×2.5–5', 22, CJUN, weight=600)
    body = (sd.title('In the full sim, wider stripes and funnelling raise it up to ×1.8; a leaky gap adds nothing',
                     'Paired arms: same angles, seeds, noise and Y as the baseline; odd amplitude ⟨sgn δ · a1⟩, |δ| ≥ 50 µm')
            + sd.row(P.svg('levers'),
                     sd.col(sd.callout('The kernel’s late gain from a leak never reaches the peak-sample a1.', GREY, 25),
                            sd.callout('A faster sheet (2 or 0.5 MΩ/sq) reaches ×2.1–2.4, but its charge runs along the stripes '
                                       'over many Y rows: Y/X kept 0.12 / 0.001 against 0.9 in the data.', '#c9a5ab', 25),
                            sd.callout('Width and funnelling act at the peak, but they don’t multiply: narrow gaps leave '
                                       'little to funnel.', ORANGE, 25), w=420), gap=40))
    D.slide('levers', body, _det(
        'Baselines: funnel arms against the production-kernel D = 0 arms (condor 4507242); hybrid-kernel arms (production Y + nx = '
        '1560 X kernels, <code>2026-10-21/build_hybrid.py</code>) against the f = 0 hybrid control glg0 (×0.95 of the production '
        'arms: the grid). Error bars from the spread of the per-(angle, seed) paired differences.',
        'Funnel: <code>--gap-funnel F</code> moves a fraction F of the avalanches that land in an ESL gap onto the nearer strip edge '
        '(20 µm inside) before the kernel lookup, with its own rng stream so the pairing is exact.',
        'All arms are at 5 MΩ/sq except the two pale rows, which are against the 5 MΩ/sq control: 2 MΩ/sq (condor 4520945) '
        'and 0.5 MΩ/sq (4520942). Both fail the Y view: for vertical tracks the fraction of Y events with a compact '
        'cluster (≤ 12 channels over 5σ, within ±3) relative to X is 0.12 and 0.001, against 0.92 for the 5 MΩ/sq sim '
        'and ≈ 0.9 in the data (det3 June-15, det2, det3 May). At 0.5 MΩ/sq the timing also follows the data’s trend '
        '(+28 → −32 ns over the width quartiles).',
        _src('r5/paired_odd.py', 'r5/paired_odd_*.json', 'condor 4520688, 4520810, 4520870, 4520871')), short='Levers')


# --------------------------------------------------------------------------- late cross-stripe content
def s_leak(D):
    L = _j('sim_late2.json'); so = _j('stripe_orient2.json', os.path.join(HERE, '..', '2026-10-20'))
    if not (L and so): return
    jun = [v['X']['n2_late'] * 1e3 for b, v in so.items() if b in ('det2', 'det3 June-15', 'det4', 'det6', 'det7')]
    may = [v['X']['n2_late'] * 1e3 for b, v in so.items() if b in ('det3 May', 'det4 May (5-6 bulk)', 'Apr-24 bulk "horizontal"')]
    arms = [('xo0', 'sim, production', GREY), ('glg0', 'sim, f = 0 control', GREY), ('glg0.03', 'sim, gap leak 0.03', GOLD),
            ('glg0.1', 'sim, gap leak 0.1', ORANGE), ('glw410', 'sim, stripes 420', PURPLE), ('glw690', 'sim, stripes 700', GREEN)]
    arms = [a for a in arms if a[0] in L]
    P = sd.Plot(1100, 620, x=(0, len(arms)), y=(-5, 40), xlabel='', ylabel='late ±2 content across the stripes (×10⁻³)',
                title='X view, mean ±2 waveform 600–960 ns after the central peak')
    P.yticks([(v, str(v)) for v in (0, 10, 20, 30, 40)])
    P.band([0, len(arms)], [min(jun)] * 2, [max(jun)] * 2, CJUN, 0.12, tip='June builds: ' + ', '.join(f'{v:.1f}' for v in jun))
    P.band([0, len(arms)], [min(may)] * 2, [max(may)] * 2, BLUE, 0.12, tip='May and earlier builds: ' + ', '.join(f'{v:.1f}' for v in may))
    P.hline(0, MUT, None, 1.5)
    for i, (k, lab, c) in enumerate(arms):
        v = L[k]['late2'] * 1e3; e = L[k]['err'] * 1e3
        P.vbar(i + 0.5, v, 70, c, base=0, tip=f'{lab}: {v:+.1f} ± {e:.1f} ×10⁻³, {L[k]["n"]} events')
        P.text(i + 0.5, -8.5, lab.replace('sim, ', ''), 19, INK, 'middle')
    P.text(0.1, max(jun) + 1.5, 'June builds', 22, CJUN, weight=600)
    P.text(0.1, max(may) + 1.5, 'May / earlier', 22, BLUE, weight=600)
    body = (sd.title('The June gaps conduct: 5–10 % of the stripe conductance gives their late cross-stripe charge',
                     'Same metric on data and sim (stripe_orient2.py); the sim had none until the gaps were made to conduct')
            + sd.row(P.svg('late cross-stripe'),
                     sd.col(sd.legend([('June builds (data)', CJUN, 'box'), ('May and earlier (data)', BLUE, 'box')], 23),
                            sd.callout('This closes “why the June boards leak across the stripes”. It does not make the '
                                       'asymmetry: the same arms leave a1 at ×0.9–1.0.', GOLD, 25), w=480), gap=44))
    D.slide('leak', body, _det(
        'Data per build from <code>2026-10-20/stripe_orient2.json</code> (X view = across the stripes on every board with '
        'stripes along y). Sim arms scored with <code>r5/sim_late2.py</code> on the same estimator.',
        'Stripe width does not move this number (420 / 700 µm: −2.8 / +0.2), so width and gap conduction are separate, '
        'measurable properties of a board.',
        _src('r5/sim_late2.py', 'r5/sim_late2.json', 'condor 4520810, 4520870')), short='Gap leak')


# --------------------------------------------------------------------------- timing by width
def s_timing(D):
    tw = _j('timing_width.json')
    if not tw: return
    d = tw['data det3 June-15 X']; s = tw['sim xo0 |delta|>=150']
    P = sd.Plot(1080, 620, x=(4, 13), y=(-140, 80), xlabel='central pulse width (samples ≥ half max, 60 ns each)',
                ylabel='excess vs central: first-order delay (ns)', title='Width quartiles; short pulses ≈ near-vertical tracks')
    P.xticks([(v, str(v)) for v in (4, 6, 8, 10, 12)]); P.yticks([(v, f'{v:+d}' if v else '0') for v in (-120, -80, -40, 0, 40, 80)])
    P.hline(0, MUT, None, 1.5)
    P.line([r['wid'] for r in d], [r['delay'] for r in d], CJUN, 3.5,
           tips=[f'data q{r["q"]}: width {r["wid"]:.1f}, delay {r["delay"]:+.0f} ns, a {r["a"]:+.3f}, peak {r["epk"]:+.0f} ns, n {r["n"]}' for r in d])
    P.line([r['wid'] for r in s], [r['delay'] for r in s], BLUE, 3.5,
           tips=[f'sim q{r["q"]}: width {r["wid"]:.1f}, delay {r["delay"]:+.0f} ns, a {r["a"]:+.3f}, peak {r["epk"]:+.0f} ns, n {r["n"]}' for r in s])
    P.text(4.3, 62, 'late for short pulses (the beam study’s +114 ns)', 22, CJUN, weight=600)
    P.text(8.2, -125, 'the sim is early everywhere', 23, BLUE, weight=600)
    body = (sd.title('The excess is late for near-vertical tracks and early for inclined ones; the sim is always early',
                     'Fit D = a·W0 + b·W0′ to the mean (n+1 − n−1) waveform, delay = −b/a (the cloud_basics §42 estimator)')
            + sd.row(P.svg('timing by width'),
                     sd.col(sd.legend([('det3 June-15 X (data)', CJUN), ('sim, |δ| ≥ 150 µm', BLUE)], 24),
                            sd.callout('“Prompt” (this deck) and “delayed 114 ns” (beam study) are one result under '
                                       'different angle mixes.', GREY, 25),
                            sd.callout('The sim’s pulses are longer too: its shortest quartile is the data’s third. '
                                       'At matched width it leads by ~50 ns more.', BLUE, 25), w=520), gap=48))
    D.slide('timing', body, _det(
        'Over all angles the same estimator gives −14 ns on det3. The nTof_x17 study '
        '(<code>sps_beam_test_26/analysis/cloud_basics/FINDINGS.md</code> §42) selects tan ≈ −0.005 and reads +114 ns.',
        'A late excess for vertical tracks is what sheet transport predicts (charge spreading across its stripe onto the '
        'neighbouring pad). The sim’s registration excess never comes out late, at any width or in any variant '
        '(funnel −110…−130 ns, wide stripes −101…−129 ns).',
        _src('r5/timing_width.py', 'r5/timing_width.json', 'r5/delta_timing.py')), short='Timing')


# --------------------------------------------------------------------------- the reconstruction
def s_reco(D):
    rr = _j('reco_resid_map.json', HERE); ct = _j('reco_correct_test.json', HERE)
    if not (rr and ct): return
    runs = [k for k in rr if 'vs_raw_map' in rr[k]]
    ks = [(k, rr[k]['vs_raw_map']['fx0fy0']['k'], rr[k]['vs_raw_map']['fx0fy0']['k_err']) for k in runs]
    w = np.array([1 / e ** 2 for _, _, e in ks]); K = float(np.sum(w * np.array([k for _, k, _ in ks])) / w.sum()); EK = float(1 / np.sqrt(w.sum()))
    P = sd.Plot(800, 620, x=(-3, 3), y=(0, len(ks)), xlabel='K: M3 residual per unit raw centroid bias', ylabel='',
                title='Transfer into the production x_p0', margin=(60, 20, 90, 200))
    P.xticks([(v, f'{v:+d}' if v else '0') for v in (-3, -2, -1, 0, 1, 2, 3)])
    P.band([K - EK, K + EK], [0, 0], [len(ks), len(ks)], CJUN, 0.15, tip=f'weighted mean K = {K:+.2f} ± {EK:.2f}')
    P.vline(0, MUT, '6 6')
    for i, (k, v, e) in enumerate(ks):
        y = len(ks) - i - 0.5; c = CJUN if k.startswith('det3') else (GREEN if k.startswith('det2') else GREY)
        P.raw(sd.line(P.X(max(v - e, -3)), P.Y(y), P.X(min(v + e, 3)), P.Y(y), c, 3))
        P.points([v], [y], c, 8, tips=[f'{k}: K {v:+.2f} ± {e:.2f}, {rr[k]["vs_raw_map"]["fx0fy0"]["n"]} cells, {rr[k]["n"]} events'])
        P.raw(sd.T(P.x0 - 14, P.Y(y) + 8, k, 21, INK, 'end'))
    P.text(K + 0.1, len(ks) - 0.25, f'K = {K:+.2f} ± {EK:.2f}', 23, CJUN, weight=600)
    Q = sd.Plot(800, 620, x=(0, len(ct)), y=(0, 4.0), xlabel='', ylabel='M3 residual-map χ² / cell', title='Leave-one-out correction')
    Q.yticks([(v, f'{v:.0f}') for v in (0, 1, 2, 3, 4)]); Q.hline(1, MUT, '6 6', 1.5)
    for i, (k, o) in enumerate(ct.items()):
        c = CJUN if k.startswith('det3') else (GREEN if k.startswith('det2') else GREY)
        b0, b1 = o['chi2_before'] / o['ncell'], o['chi2_after'] / o['ncell']
        Q.vbar(i + 0.32, min(b0, 4.0), 34, c, opacity=0.4, tip=f'{k}: before {o["chi2_before"]:.0f}/{o["ncell"]}, rms {o["rms_before"]:.0f} µm')
        Q.vbar(i + 0.68, min(b1, 4.0), 34, c, tip=f'{k}: after {o["chi2_after"]:.0f}/{o["ncell"]}, rms {o["rms_after"]:.0f} µm (K {o["K_loo"]:+.2f}, scale {o["scale"]:.2f})')
        Q.text(i + 0.5, -0.16, k.split()[0], 18, INK, 'middle'); Q.text(i + 0.5, -0.32, ' '.join(k.split()[1:]), 16, MUT, 'middle')
    body = (sd.title(f'The production X position carries the field (K = {K:+.2f} ± {EK:.2f}); correcting it helps det2 and det3',
                     'Bench runs with wft bench caches: x_p0 − M3 per 49.9 mm cell vs the raw-ADC bias map, alignment terms removed')
            + sd.row(P.svg('transfer'), sd.col(Q.svg('correction'),
                     sd.legend([('det3', CJUN, 'box'), ('det2', GREEN, 'box'), ('det4/6/7', GREY, 'box'), ('pale: before, solid: after', INK, 'box'), ('χ²/cell = 1', MUT, 'dash')], 21)), gap=44))
    D.slide('reco', body, _det(
        'Frame: pos_mm = 0.78 × true column on both views of every board (read from the bench-cache windows), so the raw '
        'maps’ cells are the reco frame’s 49.92 mm cells with no flips. K is negative: the forward-model fit is pulled the '
        'opposite way to the 3-strip centroid, at ~80 % of its size.',
        'Correction x − K·dx·scale with K from the other runs (leave-one-out); offset and rotation removed, which is what the '
        'alignment absorbs (it has no x-scale term). det2 136 → 108 and 34 → 23; det3 111 → 96, 99 → 94, 97 → 92; det4 '
        'worse (162 → 174) and det6/7 unchanged: their maps are too sparse, so they are excluded.',
        'Shipped opt-in in nTof_x17: <code>wft/xasym.py</code> (<code>correct_x(x_p0, y_p0, build, scale)</code>, a new column, '
        'no product changes), maps in <code>wft/xasym_maps.json</code>, tests <code>wft/tests/test_xasym.py</code>. Using it '
        'needs a re-alignment on corrected positions and the run’s scale.',
        _src('2026-10-21/reco_resid_map.py', 'reco_transfer.py', 'reco_correct_test.py', 'chpos_lookup.json')), short='Reco')


# --------------------------------------------------------------------------- verdict and close
def s_verdict2(D):
    steps = [
        dict(label='Code, track, electronics', sub='raw ADC; noise not copied; follows the board', color=GREY, strike=True),
        dict(label='Connector, copper bus', sub='smooth across connectors; ≤ 0.01', color=GREY, strike=True),
        dict(label='Copy, offset, mixtures', sub='no fit to the class pattern', color=GREY, strike=True),
        dict(label='Gap leak', sub='explains late cross-stripe charge, not a1', color=GOLD),
        dict(label='Stripe registration', sub='right shape, both signs; ~1/3 size', color=RED),
        dict(label='+ wide stripes, funnelling', sub='×1.5–1.8 in the full sim', color=RED),
    ]
    body = (sd.title('Where it stands: stripe registration on the June print, still ~1.5–3× short in size',
                     'Struck: excluded · gold: real but not the asymmetry · red: what makes it, so far')
            + sd.flow(steps, size=23)
            + sd.row(sd.callout('What it is: a fixed, smooth, two-signed X-only field on every June board, scaled by the cloud '
                                'width. It is charge moved n−1 → n+1, locked to the pad pitch.', BLUE, 26),
                     sd.callout('What is missing: the last factor in size, and the timing: late for vertical tracks in the data, '
                                'always early in the sim.', CJUN, 26), gap=56)
            + sd.callout('For the analysis now: the production X position carries the field (K ≈ −0.8 of the centroid bias); '
                         '<code>wft.xasym</code> corrects det2 and det3 June-15 (opt-in, new column, re-align after).', GREEN, 26))
    D.slide('verdict2', body, _det(
        'Ledger of all 37 routes with verdicts: <code>FINDINGS.md</code> top. Rounds 4–8 here: X-offset sims, stripe direction '
        'per board, bookkeeping (moved, not added), the board fields, sim a1(δ), the kernel’s time build-up, the gap-leak, '
        'width and funnel arms, the timing split by width, and the transfer into the reconstruction.'), short='Verdict')


def s_close2(D):
    body = (sd.kicker('what is open, and what would settle it')
            + '<h2 style="font-size:60px;font-weight:600;line-height:1.1;color:#eef0f3">The June print sets it; the mask '
              'decides the pitch, the width and the offsets</h2>'
            + sd.row(
                sd.col(sd.p('Open', 30, DINK, 600),
                       sd.p('· Size: the best sim is ×1.6–1.8; the data need ×2.5–5.', 26, DMUT),
                       sd.p('· Timing: late for vertical tracks in data, early in every sim variant.', 26, DMUT),
                       sd.p('· The pitch: registration needs a 780 µm-commensurate print, the quote says 800.', 26, DMUT),
                       sd.p('· Why only the June builds: their gaps conduct (shown), and their print?', 26, DMUT), gap=14),
                sd.col(sd.p('Next', 30, DINK, 600),
                       sd.p('1. The mask: stripe pitch, width and alignment, June vs May. Predict δ(x, y) per board.', 26, DMUT),
                       sd.p('2. A faster sheet gives the size and the timing, but fails Y (2 and 0.5 MΩ/sq). Fast across a stripe, slow along it: what print does that?', 26, DMUT),
                       sd.p('3. For the analysis: wft.xasym on det2 and det3 June-15, with a re-alignment.', 26, DMUT),
                       sd.p('4. The det3 residual maps keep an x-scale term the alignment doesn’t fit.', 26, DMUT), gap=14),
                gap=80))
    D.slide('close2', body, _det(
        'No detector access since 10-13: everything here is archive data and simulation. The mask is the one external input '
        'that could confirm or kill the registration picture.'), dark=True, short='Open')
