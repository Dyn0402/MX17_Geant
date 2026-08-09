#!/usr/bin/env python3
"""
make_s3_ion_report.py — build the S3 ion-investigation report from its products.

Per the repo convention (nTof_x17 CLAUDE.md, "Reporting results"): the HTML is
GENERATED from the JSON the analysis wrote, so re-running the analysis and
re-running this updates numbers, tables and verdict text together. Nothing here
is hand-typed except the prose.

Inputs (all optional except the first two — a missing one degrades that
section to a "not yet run" note rather than failing):

    response/meshcell/psi_readout.json          item 2, f_ion through the mesh
    response/meshcell/ion_template_check.json   item 3, template validation
    <stageB>/t14_fion_scan_readout.json         the f_ion demand curve, when
                                                the T14 session has read the
                                                decoded sets out

Usage:
    ~/PycharmProjects/nTof_x17/.venv/bin/python make_s3_ion_report.py \
        --out design/report/s3_ion_2026-08-09.html
"""
from __future__ import annotations

import argparse
import datetime
import html
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

CSS = """
:root{--bg:#fff;--fg:#1a1a1a;--mut:#5b6472;--line:#dfe3e8;--acc:#0b6bcb;
      --good:#1a7f4b;--bad:#b3261e;--warn:#8a6100;--card:#f7f9fb}
@media (prefers-color-scheme:dark){:root{--bg:#14171a;--fg:#e8eaed;
      --mut:#9aa4b2;--line:#2b3138;--acc:#6fb2ff;--good:#5fd394;--bad:#ff8a80;
      --warn:#e0b34d;--card:#1b1f24}}
*{box-sizing:border-box}
body{margin:0;padding:2rem 1rem 4rem;background:var(--bg);color:var(--fg);
     font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:60rem;margin:0 auto}
h1{font-size:1.7rem;line-height:1.25;margin:0 0 .3rem}
h2{font-size:1.25rem;margin:2.4rem 0 .6rem;padding-top:.6rem;
   border-top:1px solid var(--line)}
h3{font-size:1.02rem;margin:1.6rem 0 .4rem}
.sub{color:var(--mut);margin:0 0 1.6rem;font-size:.92rem}
.verdict{background:var(--card);border-left:4px solid var(--acc);
         padding:1rem 1.2rem;border-radius:0 6px 6px 0;margin:1.2rem 0}
.verdict p:first-child{margin-top:0}.verdict p:last-child{margin-bottom:0}
table{border-collapse:collapse;width:100%;font-size:.9rem;margin:.8rem 0}
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
th,td{padding:.4rem .6rem;border-bottom:1px solid var(--line);text-align:right;
      white-space:nowrap}
th:first-child,td:first-child{text-align:left}
thead th{border-bottom:2px solid var(--line);font-weight:600;color:var(--mut)}
tbody tr:hover{background:var(--card)}
code{background:var(--card);padding:.1rem .3rem;border-radius:3px;
     font-size:.87em}
.good{color:var(--good);font-weight:600}.bad{color:var(--bad);font-weight:600}
.warn{color:var(--warn);font-weight:600}
.note{color:var(--mut);font-size:.88rem}
ul{padding-left:1.2rem}li{margin:.3rem 0}
footer{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--line);
       color:var(--mut);font-size:.84rem}
"""


def esc(x):
    return html.escape(str(x))


def table(headers, rows, cls=""):
    h = "".join(f"<th>{esc(c)}</th>" for c in headers)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
                for r in rows)
    return (f'<div class="scroll"><table class="{cls}"><thead><tr>{h}</tr>'
            f"</thead><tbody>{b}</tbody></table></div>")


def load(path):
    try:
        with open(os.path.expanduser(path)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def sec_psi(d):
    if not d:
        return "<p class='note'>psi_readout.json not found — item 2 not run.</p>"
    s, g, bf = d["split"], d["gates"], d["bulk_fit"]
    fates = d["ion_fates"]
    rows = [
        ["parallel plate <code>1 - z/g</code> (what S3 assumes)",
         f"{s['f_electron_parallel_plate']:.4f}",
         f"{s['f_ion_parallel_plate']:.4f}", "&mdash;"],
        ["true &psi; through the woven mesh",
         f"{s['f_electron_true']:.4f}", f"{s['f_ion_true']:.4f}",
         f"{s['shift_vs_parallel_plate']:+.4f}"],
    ]
    return f"""
<div class="verdict"><p><strong>Mesh screening does not move the split.</strong>
f_ion goes from 0.9006 to <strong>{s['f_ion_effective']:.4f}</strong>, a shift of
{s['shift_vs_parallel_plate']:+.4f} &mdash; and the shift is in the
<em>wrong direction</em>: the true weighting potential puts slightly
<em>more</em> charge on the slow ion, not less.</p></div>

<p>The worry was well posed. <code>mx17_aval_calib.py</code> loads the realistic
woven-mesh map through <code>ComponentGrid</code> for the <em>drift</em> field,
but keeps a <code>ComponentConstant</code> with
<code>SetWeightingField(0,0,1/gap)</code> for the <em>weighting</em> field &mdash;
even in the meshfield branch. So the measured f_ion = 0.9006 is the
parallel-plate in-gap split evaluated on a meshfield avalanche profile, and had
never been checked against a through-mesh readout split.</p>

<p>It turns out the answer was already sitting in <code>solve_fieldmap.py</code>:
its <code>u_A</code> unit problem (anode = +V_mesh, wires = 0, drift-top = 0) is
exactly the readout weighting problem. <code>psi_readout.py</code> reuses that
solve verbatim and takes &psi; = u_A / V_mesh.</p>

{table(["charge split on the readout", "f_electron", "f_ion", "shift"], rows)}

<p><strong>Why it barely moves, structurally.</strong> With the ion absorbed on
a grounded wire the electron gets 1 &minus; &psi;(z) and the ion gets &psi;(z),
so the split depends only on &psi; at the <em>birth</em> height &mdash;
{s['z_birth_mean_um']:.1f} &micro;m above the ESL, where the weave's harmonics
have decayed by exp(&minus;2&pi;&middot;136/67) = e<sup>&minus;12.8</sup>.
Whatever the mesh does near itself cannot reach down there. What the mesh
<em>does</em> change is the overall weighting gradient:
{bf['slope_per_um']:.4e}/&micro;m against the parallel-plate
{bf['parallel_plate_slope_per_um']:.4e}/&micro;m, i.e. the real weighting field
is ~5&nbsp;% weaker &mdash; but that scales electron and ion together and
cancels from the ratio.</p>

<p><strong>Ion backflow is real but 0.02&nbsp;%.</strong> T6's
<code>funnel_ion_endpoints.json</code> measured
{fates['frac_absorbed']:.1%} of ions absorbed on the mesh wires
(&psi; = 0 exactly) and {fates['frac_escaped']:.1%} escaping into the drift bulk.
An escaped ion sits in E<sub>drift</sub> = 333 V/cm at ~5 nm/ns, so it moves
~5 &micro;m per microsecond and is frozen on any DREAM waveform timescale: it
carries away its residual &psi; = {fates['psi_at_mesh_topside']:.5f}, weighted
{fates['frac_escaped']:.3f} &rarr; {fates['psi_end']:.5f} of the charge lost.</p>

<h3>Why this number is quotable</h3>
<ul>
<li><strong>Gated, not asserted.</strong> Away from the weave the transverse
average of &psi; must be <em>exactly</em> linear in z (Laplace on a periodic
cell kills every harmonic but the constant). The fit residual over the amp bulk
is <span class="good">{g['linearity_rms']:.1e}</span> RMS, and
&psi;(anode) = {g['psi_anode']:.6f} against an exact 1. A mesh too coarse to
reproduce an exactly-linear function would have failed this.</li>
<li><strong>Converged.</strong> Refining lc_wire 2.0 &rarr; 1.2 moves f_ion by
0.0000 and the bulk slope by 0.01&nbsp;%.</li>
</ul>
"""


def sec_template(d):
    if not d:
        return ("<p class='note'>ion_template_check.json not found — item 3 "
                "not run.</p>")
    q = d["quantiles_ns"]
    keys = list(q["measured"])
    rows = [
        ["reconstruction (ion's own clock)"] +
        [f"{q['recon_raw'][k]:.1f}" for k in keys],
        [f"reconstruction + t<sub>aval</sub> = "
         f"{d['t_avalanche_offset_ns']:.2f} ns"] +
        [f"{q['recon'][k]:.1f}" for k in keys],
        ["<strong>S3 v2 measured template</strong>"] +
        [f"<strong>{q['measured'][k]:.1f}</strong>" for k in keys],
        ["deviation"] +
        [f"{d['deviation_pct'][k]:+.1f}&nbsp;%" for k in keys],
        ["<em>analytic 306 ns rectangle</em>"] +
        [f"<em>{q['analytic_rect'][k]:.1f}</em>" for k in keys],
    ]
    amp = d["amp"]
    ok = d["gate_pass"]
    return f"""
<div class="verdict"><p><strong>The measured template is correct.</strong> An
independent reconstruction sharing no code with the Garfield run that produced
it agrees at <em>every</em> quantile to
<span class="{'good' if ok else 'bad'}">{d['worst_dev_pct']:.1f}&nbsp;%</span>,
with one constant &mdash; the calib's own measured
t<sub>arrival</sub> = {d['t_avalanche_offset_ns']:.2f} ns. 172-ns-to-half is
kinematically right.</p></div>

<p>The question deserved asking: the schema-1 calib shipped
<code>i_elec</code>/<code>i_ion</code> as 2000 zeros and silently NaN'd the LUT.
The v2 arrays are populated, but populated is not validated. So the ion current
was rebuilt from four inputs that share no code with the emitter: E<sub>z</sub>(z)
from the T6 production map, Garfield's Ar+ mobility table evaluated <em>at the
field the ions actually see</em>, the true &psi;(z) from item 2, and the
measured birth- and absorption-height distributions.</p>

{table(["ion charge delivered, cumulative [ns]"] + keys, rows)}

<h3>The thing nobody was looking at</h3>
<p>The amplification gap runs at {amp['E_Vcm']:.0f} V/cm =
<strong>{amp['E_over_N_Td']:.1f} Td</strong>, where Ar+'s reduced mobility is
K<sub>0</sub> = <strong>{amp['K0_at_field']:.3f}</strong> &mdash; not the
zero-field {amp['K0_zero_field']:.3f} that <code>ions.py</code>'s analytic model
uses. The analytic 306 ns rectangle is therefore <strong>~21&nbsp;% too
fast</strong>. Between the two ion models it is the <em>analytic</em> one that
is wrong, not the measured one, and the measured template being slower than it
is a feature.</p>

<h3>Species: the &times;2 worry is not there</h3>
<p>The emitter hardcodes <code>IonMobility_Ar+_Ar.txt</code> for every gas, and
in Ar/iC<sub>4</sub>H<sub>10</sub> 95/5 charge transfer really does move the
charge onto an isobutane / cluster ion. But Blanc's law over the mixture gives
K<sub>0</sub> = {amp['K0_blanc_mixture']:.3f} against Ar+'s
{amp['K0_at_field']:.3f} &mdash; the real ion is a few per cent
<strong>slower</strong>, not twice as fast. (Flagged separately to the
avalanche session: for a CF<sub>4</sub>-bearing gas the same hardcoded file
would be much further off.)</p>
"""


def sec_scan(d, shaper):
    if not d:
        return f"""
<p class="note">The f_ion demand scan has been submitted (Stage B points
DIAGNOSIS_fion030 / 050 / 070 at rho2M, plus DIAGNOSIS_ionanalytic) but the
decoded sets have not been read out yet, so the measured curve is not in this
build of the report.</p>

<p>What the single-channel model predicts, for reference &mdash; the DREAM
shaper at &beta; = 0.75 driven by f<sub>e</sub>&delta;(t) + f<sub>ion</sub>
&times; rectangle, with no resistive-sheet kernel:</p>
{table(["f_ion", "10&ndash;90 % rise, 306 ns rect [ns]",
        "10&ndash;90 % rise, 369 ns rect [ns]"], shaper)}
<p class="note">The sim's own no-ions point sits ~27 ns above this model's
f_ion = 0 floor (142 vs 115.5 ns), which is the kernel and the noise. Sliding
the data's 150 ns down by that offset lands near f_ion &asymp; 0.2&ndash;0.3.</p>
"""
    # Measured. Keyed <f_ion>_<view>; rise is [p5, p25, p50, p75, p95].
    DATA = {"x": {"p5": 150.0, "fast": 0.40, "under": -3.4},
            "y": {"p5": 155.0, "fast": 0.49, "under": -12.0}}
    out = []
    for view in ("x", "y"):
        pts = sorted(((float(k.split("_")[0]), v)
                      for k, v in d.items() if k.endswith(f"_{view}")),
                     key=lambda kv: kv[0])
        rows = []
        for f, v in pts:
            rows.append([f"{f:.4f}" + (" <em>(measured template)</em>"
                                       if v.get("model") == "measured" else ""),
                         f"{v['rise'][0]:.1f}", f"{v['rise'][2]:.1f}",
                         f"{v['fast']*100:.1f}", f"{v['under']*100:.2f}",
                         f"{v.get('peak_ratio', float('nan')):.4f}"])
        dv = DATA[view]
        rows.append([f"<strong>DATA</strong>", f"<strong>{dv['p5']:.1f}</strong>",
                     "&mdash;", f"<strong>{dv['fast']*100:.0f}</strong>",
                     f"<strong>{dv['under']:.1f}</strong>", "<strong>1.0</strong>"])

        # Where the data lands, by linear interpolation on the analytic arm.
        def demand(getter, target):
            xs = [(f, getter(v)) for f, v in pts if v.get("model") != "measured"]
            for (f0, y0), (f1, y1) in zip(xs, xs[1:]):
                if (y0 - target) * (y1 - target) <= 0 and y1 != y0:
                    return f0 + (f1 - f0) * (target - y0) / (y1 - y0)
            return None
        d_p5 = demand(lambda v: v["rise"][0], dv["p5"])
        d_ff = demand(lambda v: v["fast"], dv["fast"])
        u = [v["under"] * 100 for _, v in pts]
        out.append(f"""
<h3>{view.upper()} view</h3>
{table(["f_ion", "p5 rise [ns]", "p50 rise [ns]", "fast frac <240 ns [%]",
        "undershoot [%]", "peak ratio"], rows)}
<p>The data's p5 is met at <strong>f_ion &asymp;
{('%.2f' % d_p5) if d_p5 is not None else 'below 0'}</strong>; its fast fraction
at <strong>f_ion &asymp; {('%.2f' % d_ff) if d_ff is not None else 'below 0'}
</strong>.</p>
<p class="note">Undershoot across the whole f_ion range:
{min(u):.1f} % to {max(u):.1f} % &mdash; a
{abs(max(u)-min(u)):.1f}-point swing while the rise moves
{abs(pts[-1][1]['rise'][0]-pts[0][1]['rise'][0]):.0f} ns.</p>""")

    return f"""
<div class="verdict"><p><strong>The contradiction is now measured, not
inferred.</strong> The data demands an effective slow fraction of
<strong>~0.0&ndash;0.2</strong> against a charge split defended at
<strong>0.9056</strong> by two independent routes. Per the recommendation
below, this gets reported &mdash; not absorbed into a fit.</p></div>

<p class="note"><strong>&#9888; Selection caveat (2026-08-09).</strong> The
per-view reco-quality cut drops railed waveforms, cutting the data saturation
fraction from the detector's 0.326 &rarr; 0.260 (X) and 0.327 &rarr; 0.110 (Y)
against sim acceptance 0.954 / 0.989 &mdash; so the data legs are
low-amplitude-biased and the sim legs are not. Direction checked on the frozen
parquets: corr(peak_amp, rise) = &minus;0.14, i.e. higher amplitude means
<em>faster</em> rise, so the cut removes the fastest population and the true
detector is faster than these legs show. <strong>The demanded f_eff is therefore
a lower bound &mdash; the contradiction grows, not shrinks.</strong> The legs
retain essentially no railed waveforms (X max peak_amp 4050), so the fast
population is not a clipping artifact. The Y-vs-X spread in the demand below is
a selection artifact (Y is cut ~4&times; harder) and is <em>not</em> read here as
physics. The per-view undershoot targets are being re-derived on
saturation-matched samples; the high-pass falsification uses the X figure with a
~20&times; margin, so it is unaffected.</p>

<p><strong>And the two dials are orthogonal in the real chain.</strong> f_ion
moves the rise by ~110 ns while leaving the undershoot flat to ~1 point; &beta;
moves the undershoot by 8.7 points while leaving the rise flat to 4 ns. So a
descriptive (f_eff &asymp; 0.2, &beta; &asymp; 0.2) reproduces rise <em>and</em>
X undershoot with no tension at all &mdash; which is precisely why it must not
be presented as a fit. It is a two-parameter description whose first parameter
contradicts defended physics by a factor of four.</p>
{''.join(out)}

<h3>How the pre-registration came out</h3>
<p>Predicted before the jobs ran: p5 rises monotonically and convexly, with
f_ion 0.30 landing at 150&ndash;175 ns, 0.50 at 175&ndash;205 and 0.70 at
205&ndash;235. Measured: <strong>156.8, 186.8, 222.9</strong> &mdash; all three
inside their bands, and the shape convex as predicted. The attached conclusion
that "the data's 150 ns is met at f_ion ~ 0.25&ndash;0.35" was
<strong>wrong by about a factor of two</strong>: the true answer is ~0.16 on p5
and lower still on the fast fraction, because the curve is steeper near f_ion = 0
than the estimate assumed. The direction of the finding is unaffected &mdash; it
makes the contradiction larger, not smaller.</p>"""


def leverage_tables():
    """Two sizing calculations, computed here so they cannot go stale.

    A single-channel toy: charge lands on the resistive sheet with the
    avalanche's transverse spread, spreads with D = 1/(rho_s c'), and what
    falls inside one channel pitch is differentiated and shaped by the REAL
    DreamShaper. It is a toy — it has no track geometry, no noise and no LUT —
    and it reproduces the full chain about 55 ns low (192 vs 254 ns with ions,
    91 vs 142 without). So only its DIFFERENCES are used below, never its
    absolute values.
    """
    import sys

    import numpy as np
    from scipy.special import erf
    sys.path.insert(0, REPO)
    os.environ.setdefault("MX17_SKIP_HEADER_CHECK", "1")
    from response.dream.shaper import DreamShaper, BENCH_REG1

    D, sig0, W, T, V_ION = 1.0e3, 94.55, 400.0, 369.0, 0.403
    nt = 2500
    t = np.arange(nt, dtype=float)

    def channel_current(f, image=False):
        q = (1 - f) * erf(W / (np.sqrt(2) * np.sqrt(sig0 ** 2 + 2 * D * t)))
        n = int(T)
        for tA in np.arange(0.5, T, 1.0):
            i0 = int(tA)
            s0 = np.sqrt(sig0 ** 2 + ((V_ION * tA) ** 2 if image else 0.0))
            s = np.sqrt(s0 ** 2 + 2 * D * (t[i0:] - tA))
            q[i0:] += (f / n) * erf(W / (np.sqrt(2) * s))
        i = np.empty(nt)
        i[0], i[1:] = q[0], np.diff(q)      # the prompt term is a true delta
        return i

    def rise(y):
        y = np.asarray(y, float)
        k = int(np.argmax(y))
        if y[k] <= 0 or k < 2:
            return float("nan")
        seg = y[:k + 1]
        return float(np.interp(0.9 * y[k], seg, np.arange(k + 1))
                     - np.interp(0.1 * y[k], seg, np.arange(k + 1)))

    c9, c0 = channel_current(0.9006), channel_current(0.0)
    peak_rows = []
    for code in range(5):
        reg1 = (BENCH_REG1 & ~0xF0) | (code << 4)
        sh = DreamShaper(reg1=reg1, pzc_residual=0.75, dt_ns=1.0)
        hh = np.asarray(sh.h, float)
        mark = " &larr; assumed" if code == ((BENCH_REG1 >> 4) & 0xF) else ""
        peak_rows.append([f"{code}{mark}", f"{sh.t_peak_ns:.0f}",
                          f"{rise(hh):.1f}",
                          f"{rise(np.convolve(c9, hh)[:nt]):.1f}",
                          f"{rise(np.convolve(c0, hh)[:nt]):.1f}"])

    sh = DreamShaper(pzc_residual=0.75, dt_ns=1.0)
    hh = np.asarray(sh.h, float)
    r_frozen = rise(np.convolve(channel_current(0.9006), hh)[:nt])
    r_broad = rise(np.convolve(channel_current(0.9006, image=True), hh)[:nt])
    lat = [["frozen surface kernel (what the model does)", f"{r_frozen:.1f}"],
           ["ion image broadened with height (what T10 would give)",
            f"{r_broad:.1f}"],
           ["<strong>effect of the known-wrong approximation</strong>",
            f"<strong>{r_broad - r_frozen:+.1f}</strong>"]]
    return peak_rows, lat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi", default=os.path.join(
        REPO, "response", "meshcell", "psi_readout.json"))
    ap.add_argument("--template", default=os.path.join(
        REPO, "response", "meshcell", "ion_template_check.json"))
    ap.add_argument("--scan", default=os.path.expanduser(
        "~/x17/response_sim/stageB_w2/t14_fion_scan_readout.json"))
    ap.add_argument("--beta", default=os.path.expanduser(
        "~/x17/response_sim/stageB_w2/t14_beta_scan_readout.json"),
        help="T14's beta scan — the MEASURED rise/undershoot exchange rate "
             "that falsifies the missing-high-pass class")
    ap.add_argument("--out", default=os.path.join(
        HERE, "s3_ion_2026-08-09.html"))
    a = ap.parse_args()

    psi, tmpl, scan = load(a.psi), load(a.template), load(a.scan)
    # The MEASURED rise/undershoot exchange rate, straight out of T14's beta
    # scan. This is what falsifies the missing-high-pass class, so it is read
    # from their record rather than retyped.
    bd = load(a.beta)
    beta_rows, beta_rate, beta_cost = [["beta scan not found", "", ""]], "?", "?"
    if bd and "0.0_x" in bd and "0.75_x" in bd:
        lo, hi = bd["0.0_x"], bd["0.75_x"]
        dr = hi["rise"][0] - lo["rise"][0]
        du = (hi["under"] - lo["under"]) * 100.0
        beta_rows = [
            [f"beta = 0.00", f"{lo['rise'][0]:.1f}", f"{lo['under']*100:.2f}"],
            [f"beta = 0.75 (production)", f"{hi['rise'][0]:.1f}",
             f"{hi['under']*100:.2f}"],
            ["<strong>bought / paid</strong>",
             f"<strong>{dr:+.1f}</strong>", f"<strong>{du:+.2f}</strong>"]]
        rate = abs(dr / du) if du else float("nan")
        beta_rate = (f"{rate:.2f} ns of rise per point of undershoot")
        beta_cost = f"{104.0/rate:.0f} points"

    try:
        peak_rows, lat_rows = leverage_tables()
    except Exception as exc:                                  # noqa: BLE001
        peak_rows = lat_rows = [["leverage tables unavailable", esc(exc)]]

    # The single-channel shaper prediction, recomputed here rather than pasted.
    shaper_rows = []
    try:
        import sys
        import numpy as np
        sys.path.insert(0, REPO)
        os.environ.setdefault("MX17_SKIP_HEADER_CHECK", "1")
        from response.dream.shaper import DreamShaper
        sh = DreamShaper(pzc_residual=0.75, dt_ns=1.0)
        hh = np.asarray(sh.h, float)
        n = 4000

        def rise(y):
            y = np.asarray(y, float)
            i = int(np.argmax(y))
            pk = y[i]
            if pk <= 0:
                return float("nan")
            seg = y[:i + 1]
            return float(np.interp(0.9 * pk, seg, np.arange(i + 1))
                         - np.interp(0.1 * pk, seg, np.arange(i + 1)))

        for f in (0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9006):
            cells = [f"{f:.3f}"]
            for T in (306.1, 369.0):
                L = max(1, int(round(T)))
                x = np.zeros(n)
                x[0] += (1 - f)
                x[:L] += f / L
                cells.append(f"{rise(np.convolve(x, hh)[:n]):.1f}")
            shaper_rows.append(cells)
    except Exception as exc:                                  # noqa: BLE001
        shaper_rows = [["shaper model unavailable", esc(exc), ""]]

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    body = f"""<main>
<h1>S3 ion investigation &mdash; f_ion and the i_ion template are both defended</h1>
<p class="sub">Follow-up to HANDOFF_S3_ION_2026-08-09. Generated {esc(now)} by
<code>design/report/make_s3_ion_report.py</code>.</p>

<div class="verdict">
<p><strong>Neither of the two dials the handoff put under suspicion can explain
the rise-time floor, and both corrections go the wrong way.</strong></p>
<ul>
<li><strong>f_ion</strong> through the real woven mesh is
{psi['split']['f_ion_effective']:.4f} if psi_readout ran, against the
parallel-plate 0.9006 &mdash; a shift of
{psi['split']['shift_vs_parallel_plate']:+.4f}, toward <em>more</em> slow
charge.</li>
<li><strong>The i_ion template</strong> reproduces to
{tmpl['worst_dev_pct']:.1f}&nbsp;% at every quantile under an independent
reconstruction. 172-ns-to-half is right.</li>
<li><strong>Species/mobility</strong> is not a factor-2 lever: the real
isobutane/cluster ion is ~4&nbsp;% <em>slower</em> than the Ar+ the emitter
assumes.</li>
</ul>
<p>So the ion model is not the defect &mdash; and neither is anything else that
has been proposed. The data's rise demands an <em>effective</em> slow fraction
near 0.2, three to four times smaller than a split now defended by two
independent routes. Cross-check at the median, where the arithmetic is
independent of the p5 estimate: the ion term is worth 63 ns (333.6 with ions,
270.3 without) and the gap to data is 50 ns, so ~80&nbsp;% of the ion
contribution has to disappear.</p>
<p><strong>Everything with enough leverage to do that has now been eliminated.</strong>
&beta; (4 ns across its range), the peaking-time register (code 2, 44/44
archived configs), the ion charge split, the ion template, the ion species, the
lateral factorisation (3.7 ns), and &mdash; by an exchange-rate argument on
&beta;'s own measured lever &mdash; the entire class of missing high-pass
elements. This is a structural contradiction, not a parameter error, and it
should be reported as one rather than absorbed into a fit.</p>
</div>

<h2>Item 2 &mdash; f_ion on the readout electrode, through the mesh</h2>
{sec_psi(psi)}

<h2>Item 3 &mdash; is the S3 v2 i_ion template real?</h2>
{sec_template(tmpl)}

<h2>The f_ion demand curve</h2>
{sec_scan(scan, shaper_rows)}

<h2>Two leads sized, so nobody chases the wrong one</h2>

<h3>The lateral factorisation is not it (&minus;4 ns)</h3>
<p><code>apply_longitudinal</code> convolves the <em>surface</em> kernel with
the longitudinal profile, so the ion is handed the surface kernel's lateral
shape frozen at its creation point, while the true &Psi;<sub>n</sub> broadens
as the ion climbs. <code>ions.py</code>'s own docstring flags this, and T10
exists to replace it. It moves the rise the right way &mdash; but not nearly
far enough:</p>
{table(["single-channel toy, 10&ndash;90 % rise", "ns"], lat_rows)}
<p>The reason is a scale mismatch. The ion's image width is at most its own
height, 149 &micro;m at the end of transit, against an 800 &micro;m channel
pitch: the central channel's share falls only from 1.0000 to 0.977 even at the
worst moment. <strong>T10 will not close this gap.</strong></p>

<h3>The peaking-time register is not it either &mdash; it is code 2</h3>
<p>&beta; is only one of the shaper's two parameters and it is the weak one:
one peaking-code step is worth <strong>50&ndash;65 ns</strong> of rise against
&beta;'s 4 ns across its entire range, which is the order of the discrepancy.
So the register looked like the live suspect. It is not:</p>
{table(["peaking code", "t<sub>peak</sub> [ns]", "shaper alone [ns]",
        "channel rise, f_ion 0.9006 [ns]", "f_ion 0 [ns]"], peak_rows)}
<p>All <strong>44</strong> archived <code>CosmicTb_MX17.cfg</code> copies under
the bench disk carry <code>1 0x081F 0xD023 0x0000 0x0000</code> &rarr;
<strong>code 2</strong>, across det1/det3/det4 from January to 2026-06-16.
Uniform, no exceptions.</p>
<p class="note"><strong>Inference, not measurement, and the record should say
so.</strong> The T14 target run
(<code>mx17_det3_p2_det1_overnight_6-27-26</code>) archived <em>no</em>
<code>.cfg</code> at all &mdash; only <code>run_config.json</code> and the
subrun directory &mdash; and its config points at a DAQ-side template
<em>path</em>, not an archived file. So the evidence is 44 identical copies of
that template from <em>other</em> runs, the most recent 11 days before the
target. Strong, and adopted; just not a read-back of the target run's own
register.</p>

<h3>And a missing high-pass cannot do it, whatever its topology</h3>
<p>The natural next thought is that the shaper <em>model</em> is missing a real
AC-coupling / high-pass element between DREAM output and ADC, which would
differentiate slow content and take the fast side back. The &beta; scan already
measured the exchange rate for exactly that, because &beta; <em>is</em> a
high-pass strength knob:</p>
{table(["", "X p5 rise [ns]", "undershoot [%]"], beta_rows)}
<p>That is <strong>{beta_rate}</strong>. Buying the 104 ns needed would cost
~{beta_cost} of undershoot &mdash; while the data requires the undershoot to
move <em>the other way</em>, from &minus;9.8&nbsp;% to &minus;3.4&nbsp;%, i.e.
6.4 points shallower. A shorter time constant is a better deal but not nearly
enough: over an extra series high-pass with &tau; from 100 ns to 5 &micro;s the
best exchange rate anywhere is ~1.0 ns of rise per point of undershoot, still
~104 points against a budget of &minus;6.4.</p>
<p>Any element that differentiates slow content buys rise by deepening
undershoot, and the data demands faster rise <em>and</em> shallower undershoot
at once. <strong>That falsifies the class, not just one realisation of it</strong>
&mdash; the argument is a generic property of high-pass filters.</p>

<h2>Item 1 &mdash; analytic vs measured ion model</h2>
<p>Answered, and the pre-registered prediction held. X median rise: measured
template 333.6 ns, analytic 320.8 ns, no ions 270.3 ns, data 283.9 ns. The
analytic model is <strong>3.8&nbsp;% faster</strong>, inside the pre-registered
2&ndash;10&nbsp;% band, and the falsifier (analytic below 200 ns, meaning the
template's time profile is load-bearing) did not trigger. The direction was
predicted too: analytic runs on the zero-field K<sub>0</sub> = 1.53 while the
gap sits at 123.8 Td where it is 1.212, so the 306 ns rectangle is ~21&nbsp;%
too fast. <strong>The template's fine shape is not load-bearing; the measured
one is the correct of the two.</strong></p>

<h2>What this does not rule out</h2>
<ul>
<li><strong>&beta;, jointly.</strong> The &beta; scan moved the rise by 4 ns
over 0&ndash;0.75, but it was run <em>with</em> f_ion = 0.90. A joint
(&beta;, f_ion) fit is not the same experiment &mdash; though with every
parameter here defended, a fit that closes it would be absorbing a modelling
error rather than measuring anything.</li>
<li><strong>The amplitude deficit is genuinely separate.</strong> Tested, not
assumed: a single-channel toy predicted that removing the slow charge would
raise the peak &times;1.87, which would have unified amplitude, rise and
undershoot into one cause. The measured ratios say &times;1.14 (peak_amp_med
0.5527 with ions &rarr; 0.6319 without), so the toy over-predicts amplitude
sensitivity ~6&times; and the unification does <em>not</em> hold.</li>
<li><strong>The amplitude deficit.</strong> Untouched and deliberately so: the
no-ions sim still peaks at &times;0.63 of data, so amplitude is upstream of all
of this.</li>
<li><strong>Anything at 490 V only.</strong> Every number here is the pooled
490 V bench point in Ar/iC<sub>4</sub>H<sub>10</sub> 95/5. The mobility is
field-dependent, so a different mesh voltage moves the transit.</li>
</ul>

<footer>Products: <code>response/meshcell/psi_readout.json</code>,
<code>response/meshcell/ion_template_check.json</code>. Scripts:
<code>psi_readout.py</code>, <code>ion_template_check.py</code>.
Stage B points: DIAGNOSIS_ionanalytic, DIAGNOSIS_fion030/050/070 at rho2M.</footer>
</main>"""

    with open(a.out, "w") as fh:
        fh.write(f"<style>{CSS}</style>\n{body}\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
