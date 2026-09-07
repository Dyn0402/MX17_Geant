# Contaminant diagnosis grid — READ THIS BEFORE QUOTING ANY NUMBER HERE

**Everything under `results_diagnosis_grid/` and any file with `diagnosis_grid`
in its name is DIAGNOSIS-GRID / unconstrained-contaminant-search, NOT a gas
assay.** No humidity was ever measured on the det3 bench. Every water figure —
including the 1 % point, which is the existing best fit — is a Magboltz fit to
a slow measured drift velocity, not a direct measurement. This grid exists so
that *if* the frozen T14 default (dry Ar/iC4H10 95/5, `aval_calib_meshfield_
pooled.json`) disagrees with det3 data, there is already more than one
candidate composition to interpolate from — it does not pick or bless a
winner, and none of these points should be quoted as "the" gas.

## What's in the grid

| gas | H2O | N2 | Penning | voltages |
|---|---|---|---|---|
| Ar/iC4H10 94.5/5/0.5 | 0.5 % | — | manual rP=0.40 | 490 V |
| Ar/iC4H10 94/5/1 | 1.0 % | — | manual rP=0.40 | 480, 490, 500 V |
| Ar/iC4H10 93.5/5/1.5 | 1.5 % | — | manual rP=0.40 | 490 V |
| Ar/iC4H10 93/5/1/1 | 1.0 % | 1.0 % | manual rP=0.40 | 490 V |

1 % H2O is the leading candidate (the existing det3 v(E) fit), hence the
3-point mini-scan there; the others are single points at the bench voltage.
The +N2 point is the June best fit from the water_grid family
(`nTof_x17/garfield_sim/results/water_grid.json`, mixture `Ar_iso5_H2O1_N2_1`).

## Penning: why manual rP=0.40 and not auto

Garfield has no Ar/iC4H10/H2O(/N2) parameterisation, so `--penning auto` would
silently run at rP=0 while the dry Ar/iC4H10 95/5 reference these get compared
against runs at 0.40 — biasing every comparison. Both H2O (IP 12.62 eV) and N2
(IP 15.58 eV) sit above both Ar metastables (11.55/11.72 eV), so neither opens
a *new* Penning channel on energetic grounds; at these concentrations they
mainly steal metastables into non-ionising channels. rP=0.40 (the Ar/iC4H10
value) is therefore the upper bracket, not the central value — carry 0.30-0.40
as the systematic until there is enough gain data across this grid to say
otherwise. Full reasoning: `nTof_x17/garfield_sim/mm_config.py`, the
`Ar_iC4H10_H2O_94_5_1` entry and the three added alongside it.

## Ion mobility

`mx17_aval_calib.py` hardcodes Ar+/Ar ion mobility for the ion tail shape.
Fine here — every mixture in this grid is Ar-dominant (93-94.5 % Ar), so Ar+
stays the majority drifting ion; this is the same approximation already used
for the dry 95/5 and 90/10 campaigns, not a new one introduced for this grid.
Recorded in each point's `config.ion_mobility_file`, not just here.

## Field maps

Same gas-agnostic ladder as the dry HV scan (`s2/meshfield_ladder/` on EOS,
300-700V/10V) — the FEM solve is pure electrostatics, so no new maps were
needed for this grid.

## Machine-readable provenance

Every point's JSON carries `provenance.campaign_label` = "DIAGNOSIS-GRID /
unconstrained-contaminant-search" and `config.penning_mode`/`penning_rp` —
a downstream consumer should check that field, not just this README or the
filename, before trusting a number from this grid.
