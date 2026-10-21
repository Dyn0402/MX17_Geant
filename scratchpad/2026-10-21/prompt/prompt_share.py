"""Prompt vs in-phase X sharing, data vs sim (2026-10-21 night; FOLLOWUP 'What this leaves').
Every strip route needs a cloud the strip pins while it is narrow, with the sharing widened afterwards. If so, the data's
neighbour/central ratio should start lower on the central's rise than the sim's at the same in-phase sharing s1, and
catch up later. Per event (asym_lib.load: good + hot veto; columns 4..59 of each connector as paired_odd.py):
s1 = alpha_+1 + alpha_-1 (in-phase); r(f) = (W_+1 + W_-1) / W_c at the central's rising-edge crossings of f = 25/50/75 %
of its peak (linear interpolation between samples, so the different data/sim pulse lengths compare at the same phase of
the pulse), at the peak, and on the tail at 50 % falling. Ratio of sums over events in s1 bins (noise-safe).
Also the rise time 25 -> 75 % (ns, 60 ns samples) per bin.
SEL=vert: |t(n+1) - t(n-1)| < DTN_MAX ns (40) and central FWHM >= its FW_PCT percentile (50) per set.
SIMS=tag,tag adds sim arms (any ~/x17/nb_asym_raw/sim*/ dir); ONLY=substr,... keeps matching sets.
Usage: [SEL=all|vert DTN_MAX=40 FW_PCT=50 SIMS=... ONLY=...] prompt_share.py [out.json]   (data: ~/x17/nb_asym_raw/{d3jun,d2,d3may}_x.npz; sim: simgapleak/glg0_*_x.npz)"""
import os, sys, glob, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L
RAW = os.path.expanduser('~/x17/nb_asym_raw'); HW = 4; DT = 60.
SETS = {'det3 June-15': [f'{RAW}/d3jun_x.npz'], 'det2': [f'{RAW}/d2_x.npz'], 'det3 May': [f'{RAW}/d3may_x.npz'],
        'sim glg0': sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_x.npz')),
        'sim glw190f1': sorted(glob.glob(f'{RAW}/simnarrow/glw190f1_th*_x.npz'))}
for tag in filter(None, os.environ.get('SIMS', '').split(',')):     # extra sim arms by tag, from any sim dir
    SETS[f'sim {tag}'] = sorted(glob.glob(f'{RAW}/sim*/{tag}_th*_x.npz'))
if os.environ.get('ONLY'): SETS = {k: v for k, v in SETS.items() if any(o in k for o in os.environ['ONLY'].split(','))}
BINS = np.array([.3, .45, .6, .75, .9, 1.1]); FR = (.25, .5, .75)
RNG = np.random.default_rng(3); SEL = os.environ.get('SEL', 'all'); DTN_MAX = float(os.environ.get('DTN_MAX', 40.)); FW_PCT = float(os.environ.get('FW_PCT', 50.))

def cross(c, pk, f, rising=True):
    """Fractional sample index where c crosses f * c[pk] before (rising) or after the peak; nan if not found."""
    n, ns = c.shape; thr = f * c[np.arange(n), pk]; out = np.full(n, np.nan)
    for i in range(n):
        p = pk[i]
        if rising:
            for s in range(p, 0, -1):
                if c[i, s - 1] < thr[i] <= c[i, s]:
                    out[i] = s - 1 + (thr[i] - c[i, s - 1]) / (c[i, s] - c[i, s - 1]); break
        else:
            for s in range(p, ns - 1):
                if c[i, s] >= thr[i] > c[i, s + 1]:
                    out[i] = s + (c[i, s] - thr[i]) / (c[i, s] - c[i, s + 1]); break
    return out

def at(W, tf):
    """W (n, ns) linearly interpolated at fractional sample tf (n,)."""
    ok = np.isfinite(tf); t = np.where(ok, tf, 0); i0 = np.clip(np.floor(t).astype(int), 0, W.shape[1] - 2); u = t - i0
    v = (1 - u) * W[np.arange(len(t)), i0] + u * W[np.arange(len(t)), i0 + 1]
    return np.where(ok, v, np.nan)

def one(files):
    W, al1, alm = [], [], []
    for f in files:
        d = L.load(f); p = d['p'] % 512; m = (p % 64 >= 4) & (p % 64 <= 59)
        W.append(d['W'][m]); al1.append(d['al'][1][m]); alm.append(d['al'][-1][m])
    W = np.concatenate(W).astype(float); s1 = np.concatenate(al1) + np.concatenate(alm)
    c = W[:, HW]; nb = W[:, HW - 1] + W[:, HW + 1]; pk = c.argmax(1)
    pts = {f'r{int(f * 100)}': cross(c, pk, f) for f in FR}; pts['peak'] = pk.astype(float); pts['fall50'] = cross(c, pk, .5, False)
    num = {k: at(nb, t) for k, t in pts.items()}; den = {k: at(c, t) for k, t in pts.items()}
    rise = (pts['r75'] - pts['r25']) * DT
    # x-inclination proxies: neighbour centroid-time difference (vertical: n+1 and n-1 see the same drift column) and the
    # central FWHM (vertical: the whole column sits over the central pad, longest pulse)
    ts = np.arange(W.shape[2]); pos = lambda w: np.clip(w, 0, None)
    tc = lambda w: (pos(w) * ts).sum(1) / np.maximum(pos(w).sum(1), 1e-9)
    dtn = (tc(W[:, HW + 1]) - tc(W[:, HW - 1])) * DT; fwhm = (pts['fall50'] - pts['r50']) * DT
    out = dict(n=int(len(s1)), sel=SEL, bins=[], dtn_q=[float(v) for v in np.nanpercentile(np.abs(dtn), [25, 50, 75])],
               fwhm_q=[float(v) for v in np.nanpercentile(fwhm, [25, 50, 75])])
    vert = (np.abs(dtn) < DTN_MAX) & (fwhm >= np.nanpercentile(fwhm, FW_PCT)) if SEL == 'vert' else np.ones(len(s1), bool)
    out['n_sel'] = int(vert.sum())
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        m = vert & (s1 >= lo) & (s1 < hi); row = dict(s1=[float(lo), float(hi)], n=int(m.sum()), s1_mean=float(s1[m].mean()) if m.any() else None)
        for k in pts:
            g = m & np.isfinite(num[k]) & np.isfinite(den[k])
            row[k] = float(num[k][g].sum() / den[k][g].sum()) if g.sum() > 20 else None
            if g.sum() > 20:                                   # bootstrap error of the ratio of sums
                ii = np.flatnonzero(g); bs = RNG.integers(0, len(ii), (200, len(ii)))
                row[k + '_err'] = float(np.std(num[k][ii][bs].sum(1) / den[k][ii][bs].sum(1)))
        g = m & np.isfinite(rise); row['rise_ns'] = float(np.median(rise[g])) if g.sum() > 20 else None
        out['bins'].append(row)
    return out

if __name__ == '__main__':
    R = {}
    for name, files in SETS.items():
        if not files: continue
        R[name] = r = one(files); print(f'\n{name}  ({r["n"]} events, {r["n_sel"]} selected [{SEL}]; |dt_nb| q25/50/75 '
              f'{np.round(r["dtn_q"])} ns, FWHM q25/50/75 {np.round(r["fwhm_q"])} ns)')
        print('  s1 bin      n     <s1>   r@25%         r@50%  r@75%  peak  fall50        rise25-75(ns)')
        for b in r['bins']:
            f = lambda v: f'{v:6.3f}' if v is not None else '   -  '
            e = lambda k: f'±{b[k + "_err"]:.3f}' if b.get(k + '_err') is not None else '      '
            print(f'  {b["s1"][0]:.2f}-{b["s1"][1]:.2f} {b["n"]:6d}  {f(b["s1_mean"])} {f(b["r25"])}{e("r25")} {f(b["r50"])} {f(b["r75"])} '
                  f'{f(b["peak"])} {f(b["fall50"])}{e("fall50")}  {b["rise_ns"] if b["rise_ns"] is None else round(b["rise_ns"])}')
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'prompt_share.json')
    json.dump(R, open(out, 'w'), indent=1)
