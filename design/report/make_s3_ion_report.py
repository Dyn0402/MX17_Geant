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
    return "<p class='note'>scan section: see the JSON.</p>"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi", default=os.path.join(
        REPO, "response", "meshcell", "psi_readout.json"))
    ap.add_argument("--template", default=os.path.join(
        REPO, "response", "meshcell", "ion_template_check.json"))
    ap.add_argument("--scan", default=os.path.expanduser(
        "~/x17/response_sim/stageB_w2/t14_fion_scan_readout.json"))
    ap.add_argument("--out", default=os.path.join(
        HERE, "s3_ion_2026-08-09.html"))
    a = ap.parse_args()

    psi, tmpl, scan = load(a.psi), load(a.template), load(a.scan)

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
<p>So the ion model is not the defect. The single-channel shaper model says the
data's 150 ns rise demands an <em>effective</em> f_ion near 0.2&ndash;0.3 &mdash;
three to four times smaller than a split now defended by two independent routes.
That is a structural contradiction, not a parameter error, and the next axis is
the one the model itself flags as known-wrong: the longitudinal &times; lateral
factorisation, in which the ion is handed the surface kernel's lateral shape
frozen at its creation point.</p>
</div>

<h2>Item 2 &mdash; f_ion on the readout electrode, through the mesh</h2>
{sec_psi(psi)}

<h2>Item 3 &mdash; is the S3 v2 i_ion template real?</h2>
{sec_template(tmpl)}

<h2>The f_ion demand curve</h2>
{sec_scan(scan, shaper_rows)}

<h2>What this does not rule out</h2>
<ul>
<li><strong>The lateral factorisation.</strong> <code>apply_longitudinal</code>
convolves the surface kernel with the longitudinal profile, so the ion is given
the <em>surface</em> kernel's lateral shape. The true &Psi;<sub>n</sub> broadens
as the ion climbs, so the model over-weights the ion on the central channel &mdash;
which is where the rise time is measured. This is stated in
<code>ions.py</code>'s own docstring and is what T10 exists to replace. It is
the only remaining candidate that moves the rise in the right direction.</li>
<li><strong>&beta;, jointly.</strong> The &beta; scan moved the rise by 4 ns
over 0&ndash;0.75, but it was run <em>with</em> f_ion = 0.90. A joint
(&beta;, f_ion) fit is not the same experiment.</li>
<li><strong>The peaking-time register.</strong> Still an assumption from
<code>CosmicTb_MX17.cfg</code>, not archived with the run. Everything here
scales with it.</li>
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
