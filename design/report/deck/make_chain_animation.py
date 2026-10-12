#!/usr/bin/env python3
"""Animated walk-through of the MX17 response simulation chain, for the site.

One muon crosses det3; the page follows it through every stage of the chain
(Geant4 → drift → mesh → avalanche → ions → induction → DREAM → FEU → wft),
says which code does each part and how far it is trusted, and ends by
overlaying the simulated peak-strip pulse on the real average det3 pulse at the
same angle. Toggling the ion signal shows the T14 rise-time problem directly.

The physics in the page is the independent toy (scratchpad/2026-10-12/toy_pulse.py):
its 180 ns DREAM shaper is computed here and embedded; the electron cloud,
drift and arrival times are generated in the browser with the toy's numbers.
The data and chain average pulses come from the angled waveform dump, as in
make_sim_status_deck.py.

    python3 make_chain_animation.py --wf-dir <dir with waveforms_angled_pm3.npz> --out mx17-response-chain.html
    python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py mx17-response-chain.html \
        --slug mx17-response-chain --force --deploy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from make_sim_status_deck import GRID, DT, avg_shape, TOY

EPS_MESH = 0.955       # response/meshcell/transparency_curve.csv, 333 V/cm, production_accepted
GAIN_SIM = 24094       # FRESH_LOOK §1


def payload(wf_dir):
    za = np.load(Path(wf_dir) / 'waveforms_angled_pm3.npz')
    h = TOY.shaper('manual_1', 180.0)[:1600]
    shapes = {}
    for tag, th in (('th00', 0), ('thx10', 10), ('thx20', 20)):
        d, nd = avg_shape(za[f'{tag}_data_wf'][:, 3])
        s, ns = avg_shape(za[f'{tag}_sim_wf'][:, 3])
        shapes[th] = dict(data=[round(x, 4) for x in d], chain=[round(x, 4) for x in s], nd=nd, ns=ns)
    return dict(h=[round(float(x), 5) for x in h[::2]], h_dt=2.0,
                grid_ns=[float(g * DT) for g in GRID], shapes=shapes,
                v=TOY.BASE['v'], dT=0.04766, dL=0.02425, eps=EPS_MESH, gain=GAIN_SIM,
                f_ion=TOY.BASE['f_ion'], t_ion=TOY.BASE['t_ion'])


PAGE = r'''<!--note
title: MX17 response simulation, animated
summary: One muon through the whole simulation chain (Geant4 → Magboltz drift → Garfield avalanche and ions → resistive induction → DREAM → FEU → wft), which code does each part, and how its pulse compares with the real det3 average. Toggle the ion signal to see the T14 rise-time problem.
tags: MX17,simulation,animation
date: 2026-10-12
-->
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>MX17 simulation chain</title>
<meta name="description" content="Animated walk-through of the MX17 Micromegas response simulation: what each stage does, which code does it, and where the simulation disagrees with det3 data.">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {
  --bg:#0f1219; --panel:#151a24; --panel2:#1b2130; --rule:#283044; --ink:#e7e9f0; --muted:#8b93a8; --soft:#c4c9d6;
  --mu:#9ad17a; --ele:#52b6ff; --ion:#ff9b54; --aval:#f2c94c; --strip:#45c2ae; --data:#f1f2f6; --chain:#ff5d7a;
  --ok:#5cc28f; --warn:#e7b43c; --bad:#ff5d7a;
  --f-body:"IBM Plex Sans","Segoe UI",system-ui,sans-serif; --f-mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
  color-scheme: dark;
}
html,body{background:var(--bg);color:var(--ink)}
body{margin:0;font:400 15px/1.5 var(--f-body)}
.wrap{max-width:1240px;margin:0 auto;padding:22px 16px 40px;display:grid;gap:14px}
.crumb{font:500 12px var(--f-mono)}.crumb a{color:var(--muted);text-decoration:none}.crumb a:hover{color:var(--ink)}
.eyebrow{font:500 11.5px var(--f-mono);letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
h1{font:700 clamp(24px,3.6vw,36px)/1.1 var(--f-body);margin:4px 0 6px;letter-spacing:-.01em;text-wrap:balance}
.lede{margin:0;max-width:80ch;color:var(--soft)}
.bar{display:flex;flex-wrap:wrap;gap:10px 18px;align-items:center;padding:10px 12px;background:var(--panel);border:1px solid var(--rule);border-radius:8px;position:sticky;top:0;z-index:5}
.grp{display:flex;gap:6px;align-items:center;flex-wrap:wrap}
.grp>span{font:500 11px var(--f-mono);color:var(--muted);text-transform:uppercase;letter-spacing:.06em;margin-right:2px}
button{font:600 13px var(--f-body);color:var(--ink);background:var(--panel2);border:1px solid var(--rule);border-radius:6px;padding:6px 11px;cursor:pointer}
button:hover{border-color:#3d4864}
button[aria-pressed="true"]{background:#1f2c3f;border-color:var(--ele);color:#e6f3ff}
button:focus-visible{outline:2px solid var(--ele);outline-offset:2px}
#play{min-width:82px}
input[type=range]{accent-color:var(--ele);width:110px}
.cols{display:grid;grid-template-columns:minmax(0,7fr) minmax(0,4fr);gap:14px}
@media (max-width:900px){.cols{grid-template-columns:minmax(0,1fr)}}
.panel{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:12px;display:grid;gap:8px;min-width:0}
.ph{display:flex;justify-content:space-between;align-items:baseline;gap:10px;flex-wrap:wrap}
.ph h2{font:700 15px var(--f-body);margin:0}
.ph .note{font:400 11.5px var(--f-mono);color:var(--muted)}
canvas{display:block;width:100%;border-radius:5px;background:#0b0e14}
#cDet{aspect-ratio:1000/480}
#cWf{aspect-ratio:1000/330}
.stages{display:grid;gap:6px;margin:0;padding:0;list-style:none}
.stg{display:grid;grid-template-columns:22px 1fr auto;gap:8px;align-items:center;padding:7px 9px;border-radius:6px;background:var(--panel2);border:1px solid transparent;cursor:pointer;text-align:left;font:inherit;color:inherit;width:100%}
.stg .n{font:500 11px var(--f-mono);color:var(--muted)}
.stg b{font-weight:600;font-size:14px}
.stg small{display:block;color:var(--muted);font:400 11.5px var(--f-mono)}
.stg .dot{width:10px;height:10px;border-radius:50%}
.stg.on{border-color:var(--ele);background:#18263a}
.card{background:var(--panel2);border-radius:8px;padding:12px 14px;display:grid;gap:6px;border-left:4px solid var(--muted)}
.card h3{margin:0;font:700 17px var(--f-body)}
.card .code{font:400 12px var(--f-mono);color:var(--soft);word-break:break-word}
.card p{margin:0;color:var(--soft);font-size:14px}
.card .st{font:600 12.5px var(--f-mono)}
.legend{display:flex;flex-wrap:wrap;gap:4px 14px;font:400 11.5px var(--f-mono);color:var(--muted)}
.legend i{display:inline-block;width:12px;height:4px;border-radius:2px;margin-right:5px;vertical-align:3px}
.tiles{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px}
@media (max-width:560px){.tiles{grid-template-columns:repeat(2,minmax(0,1fr))}}
.tile{background:var(--panel2);border-radius:6px;padding:7px 10px;display:grid;gap:1px}
.tile .k{font:500 10.5px var(--f-mono);text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
.tile b{font:500 19px var(--f-mono);font-variant-numeric:tabular-nums}
.verdict{font-size:14px;color:var(--soft);min-height:2.6em}
.verdict strong{color:var(--ink)}
.foot{font:400 12px var(--f-mono);color:var(--muted)}
.foot a{color:var(--ele)}
</style>
</head>
<body>
<div class="wrap">
  <div class="crumb"><a href="index.html">← notes</a> · <a href="mx17-response-sim-status.html">state-of-play slides</a></div>
  <header>
    <div class="eyebrow">MX17 response simulation · det3 cosmic bench · Ar/iC₄H₁₀ 95/5 · 490 V mesh, 1000 V drift</div>
    <h1>One muon through the whole simulation chain</h1>
    <p class="lede">A cosmic muon crosses the 30 mm drift gap. Each stage of the chain is a separate code, and the list on the right
      says which one and how far it is trusted. At the end, the simulated peak-strip pulse is laid over the <em>real</em> average det3 pulse
      at the same angle. Switch the ion signal off to see what the data look like.</p>
  </header>

  <div class="bar" role="toolbar" aria-label="Animation controls">
    <div class="grp"><button id="play" aria-pressed="false">▶ Play</button><button id="again">New muon</button></div>
    <div class="grp"><span>Angle</span><button data-ang="0">0°</button><button data-ang="10">10°</button><button data-ang="20" aria-pressed="true">20°</button></div>
    <div class="grp"><span>Ion signal</span><button id="ionOn" aria-pressed="true">Ar⁺ ions (as simulated)</button><button id="ionOff" aria-pressed="false">none</button></div>
    <div class="grp"><span>Speed</span><input id="speed" type="range" min="40" max="500" value="160" aria-label="ns per second"><span id="spd" class="note" style="font:500 11px var(--f-mono);color:var(--muted)">160 ns/s</span></div>
  </div>

  <div class="cols">
    <div class="panel">
      <div class="ph"><h2>Detector cross-section</h2><span class="note" id="clock">t = 0 ns</span></div>
      <canvas id="cDet" width="1000" height="480" aria-label="Animated detector cross-section"></canvas>
      <div class="legend"><span><i style="background:var(--mu)"></i>muon</span><span><i style="background:var(--ele)"></i>drift electrons</span>
        <span><i style="background:var(--aval)"></i>avalanche</span><span><i style="background:var(--ion)"></i>ions</span><span><i style="background:var(--strip)"></i>strip signal</span></div>
    </div>
    <div class="panel">
      <div class="ph"><h2>The chain, stage by stage</h2><span class="note">click a stage to read it</span></div>
      <ol class="stages" id="stages"></ol>
    </div>
  </div>

  <div class="cols">
    <div class="panel">
      <div class="ph"><h2>Peak-strip signal</h2><span class="note">current on the strip → DREAM shaper → FEU samples every 60 ns</span></div>
      <canvas id="cWf" width="1000" height="330" aria-label="Induced current and shaped waveform"></canvas>
      <div class="legend"><span><i style="background:#6b7590"></i>induced current</span><span><i style="background:var(--strip)"></i>after DREAM shaper (this muon)</span>
        <span><i style="background:var(--data)"></i>real det3 average</span><span><i style="background:var(--chain)"></i>full chain average</span></div>
      <div class="verdict" id="verdict"></div>
    </div>
    <div class="panel">
      <div class="card" id="card"></div>
      <div class="tiles">
        <div class="tile"><span class="k">clusters</span><b id="tCl">–</b></div>
        <div class="tile"><span class="k">electrons</span><b id="tEl">–</b></div>
        <div class="tile"><span class="k">arrived</span><b id="tAr">–</b></div>
        <div class="tile"><span class="k">peak strip</span><b id="tPk">–</b></div>
      </div>
    </div>
  </div>
  <p class="foot">Physics in the page: the independent toy model (<code>MX17_Geant/scratchpad/2026-10-12/toy_pulse.py</code>): Magboltz v = 39.1 µm/ns and diffusion,
    91 e⁻/cm in 3.8-e⁻ clusters, mesh transparency 0.955, parallel-plate ion current (f = 0.91, T = 340 ns, Ar⁺), DREAM PZC + Sallen-Key shaper at 180 ns.
    The resistive spreading is drawn but not computed. Real and chain averages: peak strip, unsaturated, no quality cut, aligned at the leading 50 %.
    Built by <code>MX17_Geant/design/report/deck/make_chain_animation.py</code>. Status and numbers: <a href="mx17-response-sim-status.html">slides</a>.</p>
</div>
<script>
const P = __PAYLOAD__;
const STAGES = [
  {id:'g4', name:'Geant4', sub:'primary ionisation', st:'ok', code:'MX17_Geant mm_sim (Geant4, PAI) → ClusterTree',
   text:'The muon is tracked through the gas and the energy it loses is split into primary ionisation clusters (on average 3.8 electrons each).',
   check:'Checked: 91 e⁻/cm, the measured value.'},
  {id:'drift', name:'Drift', sub:'Magboltz gas table', st:'warn', code:'response/digitizer/gas.py · design/gas/ (Magboltz, dry 95/5)',
   text:'Each electron drifts 0–30 mm down to the mesh at 39 µm/ns (333 V/cm), spreading sideways and along the drift by diffusion.',
   check:'Nominal: the data-fitted drift speed is 36.6 µm/ns (6.9 % lower). Humidity was never measured on any run.'},
  {id:'mesh', name:'Mesh', sub:'transparency', st:'ok', code:'response/meshcell/ (Garfield++ field map of the woven mesh)',
   text:'At the mesh the field jumps ~100×. Most electrons are funnelled through the holes; a few end on the wires.',
   check:'Checked: transparency 0.955 at 333 V/cm, now counted once (double count fixed 2026-08-11).'},
  {id:'aval', name:'Avalanche', sub:'Garfield++ S3', st:'bad', code:'response/avalanche/mx17_aval_calib.py (Garfield++ microscopic, condor)',
   text:'In the 150 µm amplification gap each electron multiplies into an avalanche of ~24 000 electrons in about a nanosecond.',
   check:'Mismatch: data need ≈ 1.6× more charge (X), and the gain grows 1.4–1.6× faster with voltage than Garfield predicts. Penning, gap, water and field shape are excluded.'},
  {id:'ion', name:'Ions', sub:'drift back to the mesh', st:'bad', code:'Garfield++ AvalancheMC, Ar⁺-in-Ar mobility · response/digitizer/ions.py',
   text:'The avalanche ions crawl 150 µm back to the mesh in ~340 ns. In a parallel-plate gap they induce ~91 % of the signal, slowly.',
   check:'Mismatch: at 10° and 20° the real pulses are as narrow as if the ion signal were almost absent (f ≈ 0.2). No mechanism found.'},
  {id:'ind', name:'Induction', sub:'resistive strips → readout', st:'warn', code:'response/solver/ (S1/W2 semi-spectral weighting potentials) · digitizer/kernel_lut.py',
   text:'The moving charges induce current through the resistive ESL strips and the kapton onto the X and Y readout strips; the resistive layer spreads it to neighbours over hundreds of ns.',
   check:'Partly off: Y/X ratio 18 % low, ±1-neighbour delay about twice the data. Spreading is drawn here but not computed.'},
  {id:'dream', name:'DREAM', sub:'CSA + 180 ns shaper', st:'warn', code:'response/dream/shaper.py (PZC + Sallen-Key, code 2)',
   text:'The front-end chip integrates the current and shapes it with a 180 ns peaking time.',
   check:'Setting confirmed from the run\'s own config (2026-10-12). The actual analog response has never been measured: one FEU pulser shot would.'},
  {id:'feu', name:'FEU', sub:'60 ns samples, ADC, noise', st:'ok', code:'response/dream/daq.py, noise.py → sim_decoded_*.root',
   text:'The FEU samples each strip every 60 ns (32 samples), adds the measured pedestal noise and digitises.',
   check:'Checked: ADC scale end to end within 1.2 % of the datasheet chain; noise round-trip 0.2 %.'},
  {id:'wft', name:'wft', sub:'same reco as data', st:'ok', code:'nTof_x17 wft (FeuReader, t13_reco) · mx17_sim_wft/t14_compare.py',
   text:'The simulated file is read and reconstructed by exactly the code used on data, and compared waveform by waveform.',
   check:'Both legs identical. Compare the green pulse with the white data average below.'},
];
const STC = {ok:'var(--ok)', warn:'var(--warn)', bad:'var(--bad)'};
const STL = {ok:'checked against something outside the code', warn:'nominal / not measured', bad:'disagrees with data'};

// ------------------------------------------------------------------ state
const cv = document.getElementById('cDet'), cx = cv.getContext('2d');
const wv = document.getElementById('cWf'), wx = wv.getContext('2d');
let ang = 20, ions = true, playing = false, speed = 160, t = -250, pinned = null, ev = null, last = null;
const T_END = 1900, PITCH = 0.78, GAP = 30, XW = 13;   // lateral half-window, mm
let seed = 7;
function rnd(){ seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; }
function gauss(){ let u = 0; while (!u) u = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rnd()); }
function poisson(l){ let L = Math.exp(-l), k = 0, p = 1; do { k++; p *= rnd(); } while (p > L); return k - 1; }
function geom(m){ return 1 + Math.floor(Math.log(1 - rnd()) / Math.log(1 - 1 / m)); }

// single-electron strip response: (1-f) h + f/T ∫_0^T h(t-u) du, on a 2 ns grid
function kernel(f, T){
  const h = P.h, n = h.length, dt = P.h_dt, out = new Float64Array(n);
  const c = new Float64Array(n + 1); for (let i = 0; i < n; i++) c[i + 1] = c[i] + h[i];
  const m = Math.max(1, Math.round(T / dt));
  for (let i = 0; i < n; i++){ const lo = Math.max(0, i - m + 1); out[i] = (1 - f) * h[i] + (f > 0 ? f * (c[i + 1] - c[lo]) / m : 0); }
  return out;
}
let KER = kernel(P.f_ion, P.t_ion);
function erf(x){ const s = Math.sign(x); x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x);
  return s * (1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x)); }

function newEvent(){
  const tan = Math.tan(ang * Math.PI / 180), v = P.v; // mm/ns
  const x0 = (rnd() - 0.5) * PITCH - 0.5 * GAP * tan;
  const ncl = poisson(91.2 / 3.83 * GAP / 10), els = [], cls = [];
  for (let i = 0; i < ncl; i++){
    const z = rnd() * GAP, n = geom(3.83); cls.push({z, x: x0 + z * tan, n});
    for (let k = 0; k < n; k++){
      const s = Math.sqrt(z / 10) * 10;
      const x = x0 + z * tan + gauss() * P.dT * s;
      const ta = z / v + gauss() * P.dL * s / v;
      const lost = rnd() > P.eps;
      els.push({z, x0: x0 + z * tan, x, ta, lost, g: -Math.log(1 - rnd()) * 0.8 + 0.2, jx: (rnd() - 0.5) * 0.6});
    }
  }
  const nS = Math.ceil(XW / PITCH) + 1, strips = [];
  for (let s = -nS; s <= nS; s++) strips.push(s);
  const N = 1100, dt = P.h_dt, wf = strips.map(() => new Float64Array(N)), cur = strips.map(() => new Float64Array(N));
  const sp = 0.222;
  for (const e of els){
    if (e.lost) continue;
    const i0 = Math.round(e.ta / dt);
    strips.forEach((s, j) => {
      const w = 0.5 * (erf(((s + 0.5) * PITCH - e.x) / (sp * Math.SQRT2)) - erf(((s - 0.5) * PITCH - e.x) / (sp * Math.SQRT2)));
      if (w < 1e-4) return;
      const a = w * e.g;
      for (let k = 0; k < KER.length && i0 + k < N; k++) wf[j][i0 + k] += a * KER[k];
      // current: prompt electron spike + ion rectangle
      if (i0 >= 0 && i0 < N) cur[j][i0] += a * (ions ? 1 - P.f_ion : 1) * 6;
      if (ions) { const m = Math.round(P.t_ion / dt); for (let k = 0; k < m && i0 + k < N; k++) cur[j][i0 + k] += a * P.f_ion / m * 6 * 2.4; }
    });
  }
  let pk = 0, best = -1;
  wf.forEach((w, j) => { const m = Math.max(...w); if (m > best){ best = m; pk = j; } });
  const W = wf[pk], C = cur[pk], mx = best;
  // leading-50 % crossing of the shaped peak strip (2 ns grid)
  let ip = W.indexOf(mx), l50 = 0; for (let i = ip; i > 0; i--) if (W[i] < 0.5 * mx){ l50 = (i + (0.5 * mx - W[i]) / (W[i + 1] - W[i])) * dt; break; }
  const phase = rnd() * 60;
  const tot = strips.map((s, j) => Math.max(...wf[j]) / mx);
  els.sort((a, b) => a.ta - b.ta);
  const firstA = els.find(e => !e.lost)?.ta ?? 0;
  ev = {els, cls, strips, pk, W, C, mx, l50, phase, tot, x0, tan, firstA, cmax: Math.max(...C)};
  document.getElementById('tCl').textContent = ncl;
  document.getElementById('tEl').textContent = els.length;
  document.getElementById('tPk').textContent = (strips[pk] >= 0 ? '+' : '') + strips[pk];
  t = -250; render();
}

// ------------------------------------------------------------------ stage logic
function autoStage(){
  if (!ev) return 0;
  if (t < 0) return 0;
  const a = ev.firstA;
  if (t < a + 20) return 1;
  if (t < a + 90) return 2;
  if (t < a + 180) return 3;
  if (t < a + 320) return 4;
  if (t < a + 520) return 5;
  if (t < a + 760) return 6;
  if (t < T_END - 200) return 7;
  return 8;
}
const olStages = document.getElementById('stages');
STAGES.forEach((s, i) => {
  const li = document.createElement('li');
  li.innerHTML = `<button class="stg" data-i="${i}"><span class="n">${i + 1}</span><span><b>${s.name}</b><small>${s.sub}</small></span><span class="dot" style="background:${STC[s.st]}" title="${STL[s.st]}"></span></button>`;
  olStages.appendChild(li);
});
olStages.addEventListener('click', e => { const b = e.target.closest('.stg'); if (!b) return; const i = +b.dataset.i; pinned = pinned === i ? null : i; render(); });
let shownStage = -1;
function showStage(i){
  document.querySelectorAll('.stg').forEach((b, j) => b.classList.toggle('on', j === i));
  if (i === shownStage) return; shownStage = i;
  const s = STAGES[i], card = document.getElementById('card');
  card.style.borderLeftColor = STC[s.st];
  card.innerHTML = `<div class="eyebrow">Stage ${i + 1} of 9${pinned !== null ? ' · pinned (click again to follow the muon)' : ''}</div><h3>${s.name}</h3>
    <div class="code">${s.code}</div><p>${s.text}</p><p class="st" style="color:${STC[s.st]}">${STL[s.st]}</p><p>${s.check}</p>`;
}

// ------------------------------------------------------------------ drawing
const DX0 = 30, DX1 = 640, DY0 = 40, DY1 = 400, AMP = 34, RES = 8;
const xs = x => DX0 + (x + XW) / (2 * XW) * (DX1 - DX0);
const zs = z => DY1 - z / GAP * (DY1 - DY0);
function glowDot(ctx, x, y, r, col, a){ ctx.globalAlpha = a; ctx.fillStyle = col; ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.fill(); ctx.globalAlpha = 1; }
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

function drawDet(stage){
  const c = cx; c.clearRect(0, 0, 1000, 480);
  const col = {mu: css('--mu'), ele: css('--ele'), ion: css('--ion'), av: css('--aval'), st: css('--strip'), mut: '#5c6479', rule: '#283044'};
  // regions
  c.fillStyle = '#0f1520'; c.fillRect(DX0, DY0, DX1 - DX0, DY1 - DY0);
  c.fillStyle = '#16140f'; c.fillRect(DX0, DY1, DX1 - DX0, AMP);
  c.strokeStyle = '#7d8796'; c.lineWidth = 2; c.beginPath(); c.moveTo(DX0, DY0); c.lineTo(DX1, DY0); c.stroke();
  c.setLineDash([3, 3]); c.beginPath(); c.moveTo(DX0, DY1); c.lineTo(DX1, DY1); c.stroke(); c.setLineDash([]);
  c.font = '12px IBM Plex Mono, monospace'; c.fillStyle = col.mut;
  c.fillText('drift cathode, −1000 V', DX0 + 6, DY0 - 8);
  c.fillText('30 mm drift, 333 V/cm', DX1 - 168, DY0 + 18);
  c.fillText('mesh', DX1 + 6, DY1 + 4);
  c.fillText('150 µm gap (drawn ×19)', DX1 + 6, DY1 + 24);
  // resistive strips (800 µm pitch) + readout strips (780 µm)
  const yR = DY1 + AMP, yK = yR + RES;
  for (let x = -XW; x < XW; x += 0.8){ c.fillStyle = '#3a3346'; c.fillRect(xs(x), yR, (DX1 - DX0) / (2 * XW) * 0.55, 4); }
  c.fillStyle = '#20232c'; c.fillRect(DX0, yK, DX1 - DX0, 10);
  c.fillStyle = col.mut; c.fillText('resistive ESL', DX1 + 6, yR + 6); c.fillText('kapton', DX1 + 6, yK + 10);
  const yS = yK + 12;
  if (ev){
    ev.strips.forEach((s, j) => {
      const xl = xs((s - 0.5) * PITCH) + 1, w = (DX1 - DX0) / (2 * XW) * 0.68;
      if (xl < DX0 || xl + w > DX1) return;
      const lit = stage >= 5 ? Math.min(1, ev.tot[j] * revealFrac()) : 0;
      c.fillStyle = '#26303d'; c.fillRect(xl, yS, w, 16);
      if (lit > 0.01){ c.globalAlpha = Math.min(1, 0.15 + lit); c.fillStyle = col.st; c.fillRect(xl, yS, w, 16 * Math.max(0.15, lit)); c.globalAlpha = 1; }
      if (j === ev.pk){ c.strokeStyle = '#fff'; c.lineWidth = 1.5; c.strokeRect(xl - 1, yS - 1, w + 2, 18); }
    });
    c.fillStyle = col.mut; c.fillText('readout strips, 780 µm', DX1 + 6, yS + 12);
    if (stage >= 5){ c.fillStyle = css('--soft'); c.fillText('peak strip ▲', xs((ev.strips[ev.pk]) * PITCH) - 30, yS + 34); }
    // muon track
    const k = Math.min(1, (t + 250) / 160);
    const zb = -(AMP + 40) / ((DY1 - DY0) / GAP), xa = xs(ev.x0 + (GAP + 2) * ev.tan), ya = zs(GAP + 2), xb = xs(ev.x0 + zb * ev.tan), yb = zs(zb);
    c.strokeStyle = col.mu; c.lineWidth = 2; c.globalAlpha = t > 300 ? 0.35 : 0.9;
    c.beginPath(); c.moveTo(xa, ya); c.lineTo(xa + (xb - xa) * k, ya + (yb - ya) * k); c.stroke(); c.globalAlpha = 1;
    // clusters (pop as the muon passes)
    for (const cl of ev.cls){
      const passed = (GAP - cl.z) / GAP < k * 0.9;
      if (passed && t < 60) glowDot(c, xs(cl.x), zs(cl.z), 2 + Math.min(cl.n, 6), col.mu, 0.35);
    }
    // electrons
    let arrived = 0;
    for (const e of ev.els){
      if (t < 0){ if ((GAP - e.z) / GAP < k * 0.9) glowDot(c, xs(e.x0) + e.jx * 4, zs(e.z), 1.8, col.ele, 0.9); continue; }
      if (t < e.ta){
        const f = Math.max(0, t / e.ta), z = e.z * (1 - f), x = e.x0 + (e.x - e.x0) * Math.sqrt(f);
        glowDot(c, xs(x), zs(z), 1.8, col.ele, 0.95);
      } else if (e.lost){
        if (t - e.ta < 250) glowDot(c, xs(e.x), DY1, 2.5, '#888', 1 - (t - e.ta) / 250);
      } else {
        arrived++;
        const dt = t - e.ta;
        if (dt < 40){ // avalanche flash in the gap
          const r = 3 + dt / 4; c.fillStyle = col.av; c.globalAlpha = 0.8 * (1 - dt / 40);
          c.beginPath(); c.moveTo(xs(e.x), DY1 + 2); c.lineTo(xs(e.x) - r, DY1 + AMP); c.lineTo(xs(e.x) + r, DY1 + AMP); c.fill(); c.globalAlpha = 1;
        }
        if (ions && dt < P.t_ion) glowDot(c, xs(e.x), DY1 + AMP - (AMP - 3) * dt / P.t_ion, 1.6, col.ion, 0.85);
      }
    }
    document.getElementById('tAr').textContent = t < 0 ? '0' : arrived;
  }
  // zoom inset: one avalanche in the 150 µm gap
  drawInset(c, stage, col);
}

function drawInset(c, stage, col){
  const X0 = 700, Y0 = 60, Wd = 270, Hd = 300;
  c.fillStyle = '#121722'; c.fillRect(X0, Y0, Wd, Hd); c.strokeStyle = '#283044'; c.strokeRect(X0, Y0, Wd, Hd);
  c.fillStyle = css('--muted'); c.font = '12px IBM Plex Mono, monospace';
  c.fillText('zoom: one hole of the mesh', X0 + 8, Y0 - 10);
  const yM = Y0 + 60, yA = Y0 + Hd - 50;
  // mesh wires
  for (let i = 0; i < 6; i++){ c.fillStyle = '#9aa3b5'; c.beginPath(); c.arc(X0 + 30 + i * 46, yM, 6, 0, 7); c.fill(); }
  c.fillStyle = '#3a3346'; c.fillRect(X0 + 10, yA, Wd - 20, 6);
  c.fillStyle = css('--muted');
  c.fillText('mesh (wire Ø 19 µm)', X0 + 8, yM - 18); c.fillText('resistive anode', X0 + 8, yA + 22);
  c.fillText('150 µm, ~31 kV/cm', X0 + 8, (yM + yA) / 2 + 4);
  if (!ev) return;
  // the representative avalanche follows the first arrival in each 120 ns window
  const a = ev.firstA; if (t < a) { // electron coming in
    const f = Math.max(0, Math.min(1, (t - a + 120) / 120)); glowDot(c, X0 + 145, Y0 + 10 + f * (yM - Y0 - 10), 3, col.ele, f > 0 ? 1 : 0); return; }
  const dt = (t - a) % 520;
  const cxA = X0 + 145;
  // electrons: rapid multiplication down to the anode (few ns, shown slowed)
  const fe = Math.min(1, dt / 60);
  for (let i = 0; i < 70; i++){
    const u = (i * 0.618) % 1, yy = yM + 4 + (yA - yM - 4) * Math.min(1, fe * (0.5 + u));
    const spread = 4 + 46 * Math.pow((yy - yM) / (yA - yM), 1.3) * (u - 0.5) * 2;
    if (fe < 1) glowDot(c, cxA + spread, yy, 1.5, col.av, 0.9);
  }
  if (fe >= 1 && dt < 90) { c.fillStyle = col.av; c.globalAlpha = 0.5 * (1 - (dt - 60) / 30); c.fillRect(cxA - 50, yA - 6, 100, 6); c.globalAlpha = 1; }
  // ions: slow drift back up to the mesh
  if (ions){
    const fi = Math.max(0, Math.min(1, (dt - 20) / P.t_ion));
    if (fi < 1) for (let i = 0; i < 60; i++){
      const u = (i * 0.618) % 1, y0 = yA - 2 - 30 * u;
      const yy = y0 - (y0 - yM - 8) * fi, sx = (u - 0.5) * 90 * (1 - fi * 0.6);
      glowDot(c, cxA + sx, yy, 1.8, col.ion, 0.9);
    }
    c.fillStyle = col.ion; c.fillText(fi < 1 ? `ions → mesh: ${Math.round(fi * P.t_ion)} / ${P.t_ion} ns` : 'ions collected', X0 + 8, Y0 + Hd - 8);
  } else { c.fillStyle = col.ion; c.fillText('ion signal switched off', X0 + 8, Y0 + Hd - 8); }
  c.fillStyle = col.av; c.fillText(`gain ≈ ${P.gain.toLocaleString('en')} (sim)`, X0 + 8, Y0 + 20);
}

function revealFrac(){ if (!ev) return 0; const i = Math.min(ev.W.length - 1, Math.max(0, Math.floor(t / P.h_dt))); return Math.max(0, ev.W[i]) / ev.mx; }

function drawWf(stage){
  const c = wx; c.clearRect(0, 0, 1000, 330);
  const L = 60, R = 980, Tp = 20, B = 290, t0 = -100, t1 = T_END;
  const X = s => L + (s - t0) / (t1 - t0) * (R - L), Y = v => B - v * (B - Tp - 30);
  c.strokeStyle = '#283044'; c.lineWidth = 1;
  for (let v = 0; v <= 1; v += 0.25){ c.beginPath(); c.moveTo(L, Y(v)); c.lineTo(R, Y(v)); c.stroke(); }
  c.fillStyle = '#8b93a8'; c.font = '12px IBM Plex Mono, monospace';
  for (let s = 0; s <= t1; s += 300){ c.fillText(s + ' ns', X(s) - 14, B + 20); }
  c.fillText('1', L - 16, Y(1) + 4); c.fillText('0', L - 16, Y(0) + 4);
  if (!ev) return;
  const tt = Math.max(-100, Math.min(t, T_END)), dt = P.h_dt, n = Math.floor(tt / dt);
  // current
  c.strokeStyle = '#6b7590'; c.lineWidth = 1.2; c.beginPath();
  for (let i = 0; i <= n && i < ev.C.length; i++){ const v = ev.C[i] / ev.cmax * 0.55; i ? c.lineTo(X(i * dt), Y(v)) : c.moveTo(X(i * dt), Y(v)); }
  c.stroke();
  // shaped
  if (stage >= 6 || t > ev.firstA + 520){
    c.strokeStyle = css('--strip'); c.lineWidth = 2.2; c.beginPath();
    for (let i = 0; i <= n && i < ev.W.length; i++){ const v = ev.W[i] / ev.mx; i ? c.lineTo(X(i * dt), Y(v)) : c.moveTo(X(i * dt), Y(v)); }
    c.stroke();
  }
  // FEU samples
  if (stage >= 7){
    for (let s = ev.phase - 120; s <= tt; s += 60){
      if (s < 0) continue; const i = Math.round(s / dt); if (i >= ev.W.length) break;
      c.fillStyle = '#fff'; c.beginPath(); c.arc(X(s), Y(ev.W[i] / ev.mx), 3.2, 0, 7); c.fill();
    }
  }
  // comparison overlay
  if (t >= T_END - 200){
    const sh = P.shapes[ang], g = P.grid_ns;
    const line = (arr, color, w, dash) => { c.strokeStyle = color; c.lineWidth = w; c.setLineDash(dash || []); c.beginPath();
      arr.forEach((v, i) => { const x = X(g[i] + ev.l50); i ? c.lineTo(x, Y(v)) : c.moveTo(x, Y(v)); }); c.stroke(); c.setLineDash([]); };
    line(sh.chain, css('--chain'), 1.8, [6, 4]);
    line(sh.data, css('--data'), 2.6);
  }
  c.fillStyle = '#8b93a8'; c.fillText('normalised to peak', L + 6, Tp + 4);
}

function verdict(){
  const v = document.getElementById('verdict');
  if (!ev || t < T_END - 200){ v.innerHTML = 'The comparison with the real det3 average appears when the pulse is complete.'; return; }
  const sh = P.shapes[ang];
  if (ang === 0) v.innerHTML = `<strong>0°:</strong> a vertical track gives a ~770 ns flat-topped input, so the rise mostly measures where the peak lands. This angle can't separate the ion question; use 10° or 20°. (data n = ${sh.nd})`;
  else if (ions) v.innerHTML = `<strong>${ang}°, ions as simulated:</strong> the pulse peaks later and decays slower than the real one (white). This is the unresolved T14 problem: the 180 ns shaper is confirmed, so no plausible ion makes this fit. (data n = ${sh.nd})`;
  else v.innerHTML = `<strong>${ang}°, ion signal off:</strong> now the pulse matches the real data closely. The data behave as if the ions induce almost nothing, but nobody knows why. (data n = ${sh.nd})`;
}

function render(){
  const stage = pinned !== null ? pinned : autoStage();
  showStage(stage);
  const auto = autoStage();
  drawDet(auto); drawWf(auto); verdict();
  document.getElementById('clock').textContent = `t = ${Math.round(Math.max(t, 0))} ns${t < 0 ? ' (muon crossing)' : ''}`;
}

function frame(ts){
  if (!playing){ last = null; return; }
  if (last !== null) t += (ts - last) / 1000 * (t < 0 ? 160 : speed);
  last = ts;
  if (t > T_END){ t = T_END; setPlay(false); }
  render();
  requestAnimationFrame(frame);
}
function setPlay(p){ playing = p; const b = document.getElementById('play'); b.textContent = p ? '❚❚ Pause' : '▶ Play'; b.setAttribute('aria-pressed', p); if (p) requestAnimationFrame(frame); }
document.getElementById('play').onclick = () => { if (t >= T_END) t = -250; setPlay(!playing); };
document.getElementById('again').onclick = () => { newEvent(); setPlay(true); };
document.querySelectorAll('[data-ang]').forEach(b => b.onclick = () => {
  ang = +b.dataset.ang; document.querySelectorAll('[data-ang]').forEach(x => x.setAttribute('aria-pressed', x === b)); newEvent(); setPlay(true); });
function setIons(on){ ions = on; KER = kernel(on ? P.f_ion : 0, P.t_ion);
  document.getElementById('ionOn').setAttribute('aria-pressed', on); document.getElementById('ionOff').setAttribute('aria-pressed', !on);
  const keep = t; const s = seed; seed = evSeed; newEvent(); seed = s; t = keep; render(); }
let evSeed = seed;
const _ne = newEvent; newEvent = function(){ evSeed = seed; _ne(); };
document.getElementById('ionOn').onclick = () => setIons(true);
document.getElementById('ionOff').onclick = () => setIons(false);
const sp = document.getElementById('speed'); sp.oninput = () => { speed = +sp.value; document.getElementById('spd').textContent = speed + ' ns/s'; };
// deep link: #ang=10&ions=0&t=1700 opens paused at that moment
const H = new URLSearchParams(location.hash.slice(1));
if (H.has('ang')){ ang = +H.get('ang'); document.querySelectorAll('[data-ang]').forEach(x => x.setAttribute('aria-pressed', +x.dataset.ang === ang)); }
if (H.get('ions') === '0'){ ions = false; KER = kernel(0, P.t_ion); document.getElementById('ionOn').setAttribute('aria-pressed', false); document.getElementById('ionOff').setAttribute('aria-pressed', true); }
newEvent();
if (H.has('t')){ t = +H.get('t'); render(); }
else if (matchMedia('(prefers-reduced-motion: reduce)').matches){ t = T_END; render(); } else setPlay(true);
</script>
</body>
</html>
'''


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--wf-dir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(PAGE.replace('__PAYLOAD__', json.dumps(payload(a.wf_dir), separators=(',', ':'))))
    print('wrote', out)
