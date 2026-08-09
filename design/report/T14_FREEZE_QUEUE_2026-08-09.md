# T14 freeze queue — two inputs to fix BEFORE the blind comparison

**For Dylan. Two decisions, same kind, same discipline.** Both are simulation
*inputs* that the chain does not determine, and for both there is a value that
would obviously make T14 agree better. That is exactly why they have to be
frozen on independent evidence and written down **before** anyone looks at the
comparison — otherwise T14's answer has been laundered into its own input.

Neither is blocking T13, which runs on the data's own calibration bundle
whatever you choose. Both change what T14 means.

---

## 1. ρ_s — the ESL sheet resistivity

**Recommendation: freeze at ρ_s = 2 MΩ/sq.** rho1M sits below the allowed band.

ρ_s has never been measured on our modules; the plan scans {0.5, 1, 2, 5}
MΩ/sq for exactly that reason. The one independent constraint is T2b's charge
spread, which bounds the *product* ρ_s·c′. Re-quoting it for the production
insulator (50 µm kapton + 18.76 µm glue in series, `d_eff` = 70.5 µm) scales it
by c′(75)/c′(70.5) = 0.947:

| stack | band |
|---|---|
| bare 75 µm kapton (as the plan quotes T2b) | 1.5 – 2.7 MΩ/sq |
| bare 50 µm (plan's own re-quote) | 1.1 – 1.9 |
| **production, d_eff = 70.5 µm** | **1.42 – 2.56** |

The mapping is validated against the plan's own numbers: applying it to bare
50 µm reproduces 1.04–1.88 against the 1.1–1.9 the plan states.

⚠️ **The "ρ_s = 1 MΩ/sq nominal" line elsewhere in the plan is a pre-glue
artifact** (2026-08-07, bare d_k = 75 µm) and is not a T14 decision.

**Stakes:** ρ_s alone moved the T10 slow-path residual by 0.65 pp — comparable
to the entire W1→W2 boundary-model change. Stage B/C has been produced at
**both** rho2M and rho1M, so nothing is forced by what happens to exist.

## 2. Gas — dry vs wet Ar/iso 95/5 (plan P1 sub-item)

**No recommendation; this one is genuinely yours.** Both options are defensible
and the distinction between them is the whole point.

Stage B's Magboltz table gives **v_drift = 39.14 µm/ns** at 333 V/cm. The det3
bench bundle for the T14 target run measures **36.60 µm/ns** — the simulation
drifts **6.9 % fast**. The plan already flags this ("+~1 % H2O on the bench —
matters, measured v_drift 36.6 µm/ns is far below dry Magboltz"; water fraction
listed as an open P1 sub-item). This note only makes it a number.

Consequence if left as is: reconstructed depth is scaled ~6.5 % against the
sim's own truth, because reco converts time→depth with the data's v_drift.

- **(a) Keep the dry table**, and quote 6.9 % as an input systematic on every
  depth-dependent comparison. This is what T13 does today.
- **(b) Use a wet (~1 % H₂O) 95/5 table for Stage B.** This is **legitimate
  physics input if and only if** the water fraction comes from the independent
  June bench record — the measured Ar/iso + ~1 % H₂O, det3 drying 3 % → 1 %
  over a week, and the waveform-first result that v(E) matches wet-Magboltz
  95/5. It is **laundering** if the water fraction is chosen because it closes
  the 6.9 %.

The difference between (b)-as-physics and (b)-as-tuning is not in the table
used; it is in what the number was derived from and when it was written down.

---

## What has already been decided, for contrast

These were frozen before looking and are recorded with their evidence, which is
the pattern the two above should follow:

- **Avalanche calib** — pooled 490 V meshfield point, chosen on the grounds
  that all 56 slices used the same field map regardless of their voltage label
  (not on downstream agreement). Reconfirmed 2026-08-08.
- **Drift gap** — 30 mm, decided 2026-08-08, closed.
- **T14 target** — det3 `long_run_resist_490V_drift_1000V`, per P1.

## What is NOT on this list

The W2 kernel itself. It is not a free input: V6 showed W1 grounds copper that
does not exist, and W2 fixes an error of known sign and size (+25.66 % prompt
capture, pre-registered and reproduced). Choosing W1 over W2 would be choosing
a known-wrong boundary condition, not exercising a scan axis.
