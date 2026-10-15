# Three beams, one office floor, one standard fire

Three simply supported 5.0 m floor beams — **aluminium (EN AW-6082 T6), steel (S355) and glulam (GL24h)** — are each
sized so that their ambient design bending resistance carries the same office floor under their own Eurocode
(slab + finishes 6.0 kN/m plus self-weight, imposed 6.0 kN/m; M_Ed = (1.35G + 1.5Q)L²/8). They are then heated from
below and both sides by the ISO 834 standard fire, unprotected. The result is the first time each beam's modelled
**sectional bending resistance** drops below the fire-situation demand M_fi = (G + 0.5Q)L²/8.

| beam | section | resistance < demand after |
|---|---|---|
| aluminium | IPE-proportioned, h ≈ 0.205 m | ≈ 4.6 min |
| steel | IPE-proportioned, h ≈ 0.179 m | ≈ 13.5 min |
| glulam | 175 × 349 mm | ≈ 64 min |

```
python reproduce.py
```

## Method (all from the cited standards; no fitted constant)
* ISO 834 gas temperature; steel EN 1993-1-2 eq. 4.25 lumped heating (c_a(θ), ε 0.7, α_c 25, shadow factor
  k_sh computed from the section), resistance k_y,θ W_pl f_y / κ₁ with κ₁ = 0.7 (three-sided, slab on top).
* Aluminium EN 1999-1-2: c = 903 + 0.41θ, ε 0.3 (clean surface), Table 1 reduction factors for 6082-T6.
* Glulam EN 1995-1-2 reduced cross-section: β_n 0.7 mm/min, d₀ 7 mm (k₀ ramp below 20 min), k_fi 1.15, three sides.
* Input-sensitivity index u (not a statistical uncertainty or a bounding interval): one-at-a-time sweep (emissivities, convection, shadow factor, charring rate, d₀, timber density,
  imposed load, time step), half-ranges combined in quadrature; every adjacent gap must exceed max(10 %, 3u).
* Independent checks: integrator vs exact solution; steel failure temperature vs EN eq. 4.22; reduced section by
  hand at 60 min; time-step halving; section classes (fire ε for steel).

## Worked-example benchmarks
`python benchmarks.py` reproduces published worked examples with this solver's own functions: steel (SCI Fire Design
Example 1, the JRC 2012 fire workshop report, SBE MS6, UPT course tables), glulam (JRC report, Frangi example; EN 1995-1-2
eq. 4.1) and aluminium (TALAT 2713 resistance; Maljaars et al. HERON 2010 heating, approximate; two independent integrators;
EN 1999-1-2 Table 1a). No published EN worked example of unprotected-aluminium heating with stated inputs and result was found;
the aluminium heating rests on the independent integrators and the HERON comparison. 43/43 checks pass.

## Scope
Unprotected members; the slab restrains the top face continuously (no lateral-torsional buckling); no composite
action; bending only (deflection, shear and connections are not checked); continuous (non-catalogue) sizing.
This compares three idealised beams under one disclosed rule. It is **not a fire rating** of real beams or materials.
