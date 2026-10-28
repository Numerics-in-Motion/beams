# One glulam beam, two or three, in a standard fire

The office floor of [`fire-three-beams`](../fire-three-beams/) (span 5.0 m, slab + finishes 6.0 kN/m plus self-weight,
imposed 6.0 kN/m, simply supported) sits on the same spruce glulam (GL24h), cut three ways:

| contender | members | total width × depth |
|---|---|---|
| one | 1 × 240 × 320 mm | 240 × 320 mm |
| two | 2 × 120 × 320 mm | 240 × 320 mm |
| three | 3 × 80 × 320 mm | 240 × 320 mm |

The same total width, depth and section modulus, so the same cold design bending resistance (utilisation 0.871 for all
three, **when the beams are held sideways**: see "Restraint"). The model assigns each member M_fi/n. Every member is
heated on three sides (bottom and both sides; the top is shielded) by the ISO 834 standard fire. The endpoint is the
first time the **summed calculated sectional bending resistance** of a contender falls below the fire-situation load
M_fi. It is a code calculation of sectional resistance only.

```
python reproduce.py            # everything (~6 min): 54 scenes restrained and unrestrained, sensitivities, checks
python reproduce.py --quick    # the restrained race, the controls and the outside checks (~20 s)
```

## Result (frozen in `reference/result.json`)

| | one | two | three |
|---|---|---|---|
| RCS, code inputs (min) | 89.4 | 46.9 | 29.5 |
| RPM, code inputs (min) | 90.7 | 50.1 | 32.5 |

Over all **54 scenes** (method RCS / RPM; load by the EN 1990 fire combination with ψ 0.5 or 0.3, or by the simplified
η_fi = 0.6; γ_M,fi 1.0 / 1.1 / 1.2; imposed load 5 / 6 / 7 kN/m) and at the **joint adverse corner** of the uncertain
inputs (β_n 0.65–0.75 mm/min, β_0 0.60–0.70 mm/min, d_0 7–10 mm, density 380–460 kg/m³; the single beam at its
least time, the other contender at its greatest):

* one keeps its calculated resistance **at least 1.5×** as long as two (54/54 scenes clear max(10 %, 3u));
* one keeps it **at least 2.3×** as long as three (54/54);
* the smallest joint-corner ratios over all scenes are in `reference/result.json` (`restrained`, `min_ratio`);
* two against three is **unresolved** (27/54).

u is the one-at-a-time RSS half-range of the uncertain inputs, maximised over their grid (a sensitivity index, not a
statistical bound).

## Method — EN 1995-1-2:2004 (historical edition)

* **RCS** (4.2.2): effective section, d_ef = β_n t + k_0 d_0, k_0 = t/20 below 20 min (Table 4.1), d_0 7 mm;
  strength and stiffness with k_mod,fi = 1.0 (4.2.2(5)); k_fi 1.15 (glulam, Table 2.1); γ_M,fi 1.0.
* **RPM** (4.2.3): residual section (β_n t); k_mod,fi = 1 − p/(200 A_r) for bending strength, linear from t = 0 to
  20 min; 1 − p/(330 A_r) for the modulus of elasticity (used only in the unrestrained check).
* **Charring** (3.4.2): β_n 0.7 mm/min (glulam, Table 3.1). The one-dimensional β_0 = 0.65 mm/min with corner
  rounding r = d_char,0 is allowed only while the ORIGINAL width b ≥ b_min (8.15 d_char,0 below 13 mm, else
  2 d_char,0 + 80 mm); each member's resistance at each time is the smaller of the permitted representations
  (RPM; in RCS only as a sensitivity, `rcs_beta0_sensitivity`, with the same conclusions).
* The adhesive is assumed to keep bond integrity for the period (3.5); all three use the same product.
* Shear and bearing: disregarded for unnotched rectangular sections (4.3.1 simplified method). Deflection,
  connections and the slab are not checked.

## Restraint

The registered race assumes each beam's compression edge is held sideways continuously and its supports
torsionally (EN 1995-1-1 6.3.3(5): k_crit = 1). **This is an assumption, not a verified detail.** The study also
measures the **unrestrained** case (EN 1995-1-2:2004 4.3.2: when bracing is lost in the fire, lateral-torsional stability is
checked without it):

* k_crit of EN 1995-1-1 6.3.3 on the charred effective section at every time step: σ_m,crit = π √(E_0.05 I_z G_0.05
  I_tor) / (l_ef W_y) (6.31); λ_rel,m = √(f_m,k / σ_m,crit) (6.30); k_crit = 1 (λ ≤ 0.75), 1.56 − 0.75 λ
  (0.75 < λ ≤ 1.4), 1/λ² (λ > 1.4) (6.34);
* l_ef = 0.9 L (Table 6.1, simply supported, uniform load) + 2h (load on the compression edge) = 5.14 m, torsion held
  at the supports only;
* GL24h E_0.05 = 9.4 GPa (EN 1194 class table) and G_0.05 = 0.54 GPa (5/6 of G_mean 650 MPa, EN 14080:2013 5.1.3):
  the lower of the two printed tables each;
* RCS: the fire factors cancel in λ (strength and stiffness both k_mod,fi = 1, k_fi and γ_M,fi on both); RPM:
  λ² = k_mod,fi,m f_m,k / (k_mod,fi,E σ_m,crit), with E and G both reduced by 1 − p/(330 A_r);
* the rounded β_0 section uses the inscribed rectangle for I_z and I_tor (a lower bound) and its own W_y.

Unrestrained, the order holds and both claimed pairs clear in 54/54 scenes (one/two ≥ ×1.89, one/three ≥ ×3.89;
nominal prescribed RCS 89 / 33 / 14 min; with the rounded-β_0 RCS sensitivity the single beam's nominal time is 86 min). The single beam's λ reaches at most 0.775 at its endpoint. A bound that uses only the
code's first range (k_crit = 1 for λ ≤ 0.75) — the single beam credited only until its λ first exceeds 0.75, the thin
beams given their restrained times — gives the same minima as the restrained race (`branch_free_bound`), so the conclusion does not
depend on the 1.56 − 0.75 λ and 1/λ² branches.

**Equal cold strength assumes the beams are held sideways; unrestrained, the three thin beams would not pass the
cold check** (λ_rel,m 0.97, k_crit 0.83, utilisation 1.05).

## Outside checks (`benchmarks()`)

| check | printed | here |
|---|---|---|
| JRC 2012 fire workshop report, ch. 6 (Frangi) 6.1.4.3: GL24h 160 × 735, 3 sides, 30 min, W_fi | 8664 × 10³ mm³ | 8664.0 × 10³ |
| Kodur, CE808 lecture Example 10.1 (RPM): 130 × 720, 60 min, 0.64 mm/min, k_mod,fi | 0.805 | 0.80469 |
| same, residual W | 4.12 × 10⁶ mm³ | 4.119 × 10⁶ |
| Timoshenko & Goodier torsion factor β of a rectangle, h/b 1 / 2 / 4 / 10 | 0.141 / 0.229 / 0.281 / 0.312 | 0.1407 / 0.2289 / 0.2809 / 0.3123 |

`fire-three-beams/benchmarks.py` holds the wider set of worked examples for the reduced cross-section method.
All checks above read the same formulas the race uses; none of it is validated against fire tests.

## Sources

* EN 1995-1-2:2004 text: HTML transcription at dnaop.com (pages 60570_4–60570_8; a national-adoption draft copy,
  secondary): k_fi Table 2.1, β Table 3.1, 3.4.2, 4.2.2 (incl. (5): k_mod,fi = 1.0 for strength and stiffness),
  4.2.3, 4.3.2 (bracing lost in the fire). The RPM factor for the modulus of elasticity (1/330) was not legible
  in that transcription; it is the standard's expression (4.5) as commonly reproduced, used only in the
  unrestrained check, where it affects only the thin beams.
* EN 1995-1-1 6.3.3 (secondary, HTML): structuralbasics.com, "Lateral Torsional Buckling of Timber Beams" (6.31,
  6.30, l_ef = 0.9 s + 2h, I_tor of a rectangle); FIN EC help, "Lateral torsional buckling" (6.3.3, k_crit = 1 for
  λ ≤ 0.75, 6.32). **The 1.56 − 0.75 λ and 1/λ² branches were not found in a clean HTML source**; the conclusion does
  not depend on them (branch-free bound above).
* GL24h stiffness (secondary, HTML): studiomadera.it, "Classi di resistenza legno lamellare" (E_0,g,05 9.4,
  G_g,mean 0.72 kN/mm²); Edinburgh Napier CWST, "5th percentile shear modulus" (EN 14080:2013 5.1.3: G_g,05 =
  5/6 G_g,mean).
* Member sizes: Stora Enso glulam product page (spruce GL24h, widths 80–260 mm in 20 mm steps, heights in 40 mm
  steps, mean density 420 kg/m³, β_n 0.7 mm/min) — a manufactured size range, not a stock claim.
* JRC report: Vassart et al., "Eurocodes: Background & Applications — Structural Fire Design — Worked examples" (2012).
* Kodur, CE808 lecture notes, Example 10.1 (HTML transcript, secondary).

## Scope

Idealised members under one floor, one code edition (2004; later editions are not claimed), three-sided exposure,
bending only. No statement about real buildings or design choices.
