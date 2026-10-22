"""Consequences of the photographed strip film being the true pattern (2026-10-22, FOLLOWUP §15).
Film (§14): pitch 760.7 µm (±2.3 spread, ±1 % systematic), fill 0.694 (dark 528 µm). Pads 780 µm.
 A. Geometry: registration offset delta(c) of each pad column to the nearest strip centre, its advance per column, the beat
    period (and its range for ±1 %), strips across the 399.36 mm active width, the exact-box pitch 780 x 39 / 40 = 760.5 µm.
 B. Sim prediction on that geometry: the full sim's a1(delta) and s1(delta) (r5/sim_delta.json, 800 µm kernels; the curves
    depend on delta / pitch, rescaled to 760.5 µm) laid onto the 512 columns; (i) the first-harmonic amplitude the per-column
    beat search sees, (ii) the a1 and s1 means of 64-column cells (the 8 x 8 maps), (iii) the board-wide mean a1.
 C. Data at the 761 µm beat: per-column first-harmonic amplitude of s1 and d1 (= a1) over f in the ±1 % band, with a
    permutation null (r5lib survey loads, as r5/beat_col.py / beat_sig.py), per build and view.
Usage: consequences.py   (-> consequences.json)"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R5 = os.path.join(HERE, '..', '..', '2026-10-20', 'r5'); sys.path.insert(0, R5)
PS, PP, NPAD = 760.7, 780.0, 512; SYS = 0.01; OUT = {}

# ---------------------------------------------------------------- A. geometry
def beat(ps): return abs(1 / ps - 1 / PP) ** -1                     # µm
lo, hi = PS * (1 - SYS), PS * (1 + SYS)
A = dict(pitch=PS, advance_um_per_col=PP - PS, beat_mm=beat(PS) / 1e3, beat_cols=beat(PS) / PP,
         beat_mm_range=[beat(hi) / 1e3, beat(lo) / 1e3], f_per_col=PP / beat(PS), f_range=[PP / beat(hi), PP / beat(lo)],
         strips_active=NPAD * PP / PS, exact_box_pitch=PP * 39 / 40, n_beats_on_board=NPAD * PP / beat(PS),
         design800=dict(advance=-20.0, beat_mm=31.2, f=1 / 40, strips_active=NPAD * PP / 800))
OUT['geometry'] = A
print(f"A. pitch {PS} µm vs pads {PP}: offset advances {PP - PS:+.1f} µm per column (800 µm design: -20); beat {A['beat_mm']:.1f} mm = "
      f"{A['beat_cols']:.1f} columns (±1 % pitch: {A['beat_mm_range'][0]:.1f}-{A['beat_mm_range'][1]:.1f} mm), f = {A['f_per_col']:.4f}/column "
      f"({A['f_range'][0]:.4f}-{A['f_range'][1]:.4f}); {A['n_beats_on_board']:.1f} beats across the board; {A['strips_active']:.0f} strips over the active width "
      f"(800 µm: {A['design800']['strips_active']:.0f}); exact-box pitch 780 x 39/40 = {A['exact_box_pitch']:.1f} µm")

# ---------------------------------------------------------------- B. sim prediction on the 761 geometry
SD = json.load(open(os.path.join(R5, 'sim_delta.json')))['col']; SD = [r for r in SD if r]
dd = np.array([r['d'] for r in SD]) / 800.0; a1d = np.array([r['a1'] for r in SD]); s1d = np.array([r['s1'] for r in SD]); nd = np.array([r['n'] for r in SD])
fa = lambda u: np.interp(((u + 0.5) % 1) - 0.5, dd, a1d, period=1.0); fs = lambda u: np.interp(((u + 0.5) % 1) - 0.5, dd, s1d, period=1.0)
print(f"\nB. sim a1(delta) (r5/sim_delta.json, xo0 arms): range {a1d.min():+.3f}..{a1d.max():+.3f}, event-weighted mean {np.average(a1d, weights=nd):+.4f}, "
      f"uniform-delta mean {a1d.mean():+.4f}; s1 {s1d.min():.2f}..{s1d.max():.2f}")
col = np.arange(NPAD); res = {}
for nm, ps in (('761', PS), ('800 (sim design)', 800.0), ('780 (pad-locked)', 780.0)):
    rows = []
    for x0 in np.linspace(0, 1, 40, endpoint=False):
        u = (col * PP / ps + x0)                                         # strip-phase of each column centre (cycles)
        a, s = fa(u), fs(u); f = PP / beat(ps) if ps != PP else 0.0
        h = (lambda v: 2 * abs(np.mean((v - v.mean()) * np.exp(2j * np.pi * f * col)))) if f else (lambda v: 0.0)
        cells_a = a.reshape(8, 64).mean(1); cells_s = s.reshape(8, 64).mean(1)
        rows.append((h(a), h(s), np.abs(cells_a).max(), cells_a.std(), np.ptp(cells_s), a.mean()))
    r = np.array(rows); res[nm] = dict(harm_a1=float(r[:, 0].mean()), harm_s1=float(r[:, 1].mean()), cell_max_abs_a1=float(r[:, 2].max()),
                                       cell_rms_a1=float(r[:, 3].mean()), cell_ptp_s1=float(r[:, 4].mean()), board_mean_a1=[float(r[:, 5].min()), float(r[:, 5].max())])
    q = res[nm]; print(f"   {nm:17s}: beat-harmonic amplitude a1 {q['harm_a1']:.3f}, s1 {q['harm_s1']:.3f};  64-column cells: max |mean a1| {q['cell_max_abs_a1']:.3f}, "
                      f"rms {q['cell_rms_a1']:.3f}, s1 spread {q['cell_ptp_s1']:.3f};  board mean a1 {q['board_mean_a1'][0]:+.3f}..{q['board_mean_a1'][1]:+.3f} (over offsets x0)")
OUT['sim_on_geometry'] = res

# ---------------------------------------------------------------- C. data at the 761 beat
import r5lib as R
F = np.arange(A['f_range'][0] - 0.002, A['f_range'][1] + 0.002, 0.00025); rng = np.random.default_rng(11); D = {}
def detrend(c, o, sig=30.):
    n = np.bincount(c, minlength=512).astype(float); s = np.bincount(c, o, minlength=512); k = np.exp(-0.5 * (np.arange(-120, 121) / sig) ** 2)
    return o - (np.convolve(s, k, 'same') / np.maximum(np.convolve(n, k, 'same'), 1e-9))[c]
def colamp(cm, f):                                                       # per-column means, first harmonic over columns
    ok = np.isfinite(cm); c = np.flatnonzero(ok); v = cm[ok] - cm[ok].mean()
    return np.array([2 * abs(np.mean(v * np.exp(2j * np.pi * ff * c))) for ff in np.atleast_1d(f)])
print(f"\nC. data, per-column first harmonic over f = {F[0]:.4f}..{F[-1]:.4f} (761 µm ±1 % + margin), detrended (sigma 30 columns), "
      f"permutation null (300) of the column means")
for b in ('det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1'):
    for view in ('X', 'Y'):
        p, s1, d1 = [], [], []
        for r in R.rows():
            if r['build'] != b or r['view'] != view: continue
            d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
            p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
        if not p: continue
        p, s1, d1 = map(np.concatenate, (p, s1, d1)); c = 64 * (p // 64) + 63 - p % 64
        m = (np.abs(s1) < 3) & (np.abs(d1) < 3); c, s1, d1 = c[m], s1[m], d1[m]; line = f'   {b:13s} {view} n {len(c):6d}'
        for nm, v in (('s1', detrend(c, s1)), ('d1', detrend(c, d1))):
            n = np.bincount(c, minlength=512); cm = np.where(n >= 20, np.bincount(c, v, minlength=512) / np.maximum(n, 1), np.nan)
            As = colamp(cm, F); i = int(As.argmax()); okc = np.isfinite(cm)
            null = np.array([colamp(np.where(okc, np.concatenate([rng.permutation(cm[okc]), np.full((~okc).sum(), np.nan)])[np.argsort(np.argsort(~okc, kind='stable'), kind='stable')], np.nan), F).max() for _ in range(300)])
            D[f'{b}|{view}|{nm}'] = dict(ncol=int(okc.sum()), Amax=float(As[i]), fbest=float(F[i]), p=float((null >= As[i]).mean()), null95=float(np.percentile(null, 95)),
                                        A761=float(colamp(cm, A['f_per_col'])[0]))
            q = D[f'{b}|{view}|{nm}']; line += f" | {nm}: max {q['Amax']:.4f} at f {q['fbest']:.4f} (p {q['p']:.2f}, null95 {q['null95']:.4f}), at 761 µm {q['A761']:.4f}"
        print(line, flush=True)
OUT['data_beat'] = D
json.dump(OUT, open(os.path.join(HERE, 'consequences.json'), 'w'), indent=1)
