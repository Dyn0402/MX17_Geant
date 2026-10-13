#!/usr/bin/env python3
"""Peak-strip pulse shape per chamber x gas x track inclination (FRESH_LOOK §13).
Same metric as t14_compare / scratchpad 2026-10-12 wf_matched.table: 10-90 % rise of the
peak strip by linear interpolation on 60 ns samples, unsaturated (max < 3500), here also
peak >= 300 ADC (well above noise). Inclination = M3 slope ALONG the view's coordinate
(x-view strips measure ref_y, so tan_y; y-view: tan_x). DREAM settings identical in all
runs (180 ns, 60 ns; decode_dream_cfg --sweep, 2026-10-13).
Per dataset, all clean voltages pooled (rise is gain-independent to first order; checked
by splitting low/high voltage halves)."""
import glob, numpy as np, pandas as pd
C = '/home/dylan/.claude/jobs/913ff27a/tmp/c/'
DT = 60.0
def cross(w, frac, side):
    ipk = int(np.argmax(w)); lvl = frac * w[ipk]
    if side < 0:
        b = np.flatnonzero(w[:ipk + 1] < lvl)
        if not len(b): return np.nan
        i = b[-1]; j = i + 1
    else:
        b = np.flatnonzero(w[ipk:] < lvl)
        if not len(b): return np.nan
        j = ipk + b[0]; i = j - 1
    return i + (lvl - w[i]) / (w[j] - w[i]) if w[j] != w[i] else float(i)
def metrics(W):
    r = []
    for w in W:
        w = np.nan_to_num(w.astype(float))
        l10, l50, l90, r50 = cross(w, .1, -1), cross(w, .5, -1), cross(w, .9, -1), cross(w, .5, 1)
        r.append(((l90 - l10) * DT, (r50 - l50) * DT, (np.argmax(w) - l10) * DT))
    return np.array(r)
if __name__ == '__main__':
    SETS = [('det0 Ar/iso', 'det0__*Ar_Iso*', (50, 150, -200, 150), (440, 480)),
            ('det0 He/C2H6', 'det0__*He_HV_Scan_4-2*', (50, 150, -200, 150), (495, 515)),
            ('det0 He/C2H6 (4-1)', 'det0__*He_HV_Scan_4-1*', (50, 150, -200, 150), (480, 495)),
            ('det3o Ar/iso', 'det3old__*HV_Scan_5-5*', (-200, 100, -200, 50), (460, 490)),
            ('det3o Ar/CF4', 'det3old__*ArCF4*', (-200, 100, -200, 50), (590, 630))]
    BINS = [(0.0, 0.05, 'vert'), (0.15, 0.22, '~10deg'), (0.33, 0.40, '~20deg')]
    rows = []
    for name, pat, (x0, x1, y0, y1), (vlo, vhi) in SETS:
        fs = sorted(glob.glob(C + pat + '.parquet'))
        if not fs: continue
        for f in fs:
            d = pd.read_parquet(f); W = np.load(f.replace('.parquet', '_wf.npy'))[:, 2, :].astype(np.float32)
            d['i'] = np.arange(len(d)); d = d[(d.volt >= vlo) & (d.volt <= vhi) & ~d.subrun.str.startswith('final')]
            d = d[(d.ref_x > x0) & (d.ref_x < x1) & (d.ref_y > y0) & (d.ref_y < y1) & (d.peak_amp >= 300) & (d.peak_amp < 3500)]
            if not len(d): continue
            d['tan'] = np.where(d.view == 'x', d.tan_y.abs(), d.tan_x.abs())
            for lo, hi, lab in BINS:
                g = d[(d.tan >= lo) & (d.tan < hi)]
                if not len(g): continue
                m = metrics(W[g.i.values])
                for (view, volt), k in g.groupby(['view', 'volt']).indices.items():
                    for row in m[k]:
                        rows.append(dict(set=name, view=view, volt=volt, bin=lab, rise=row[0], fwhm=row[1], delay=row[2]))
    R = pd.DataFrame(rows); R.to_parquet(C + 'rise_all.parquet')
    pd.set_option('display.width', 250)
    t = R.groupby(['set', 'view', 'bin']).agg(n=('rise', 'size'), rise=('rise', 'median'), fwhm=('fwhm', 'median'), delay=('delay', 'median'))
    print(t.round(0).to_string())
    # gain dependence check: lower vs upper half of the voltage range
    R['half'] = R.groupby('set').volt.transform(lambda v: np.where(v <= v.median(), 'lowV', 'highV'))
    print(R.groupby(['set', 'bin', 'half']).rise.median().unstack().round(0).to_string())
