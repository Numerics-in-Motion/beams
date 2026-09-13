# Less live load next door, more sag here

A steel floor girder runs continuously over four columns: three **3.0 m** bays,
the rig's published I-section, the same **6 kN/m** live load on every bay. Then
the live load comes off the **middle** bay only. The end bay keeps exactly the
same load.

| the end bay's largest sag | |
|---|---|
| live load on every bay | **8.72 mm** |
| live load taken off the middle bay | **12.57 mm** (× 1.44) |

```
python reproduce.py            # the whole study and every gate, seconds
python reproduce.py --quick    # the displayed girder and the control
```

## What is claimed

In the registered rigid, level, equal-span linear-elastic continuous-girder
models, the governing live-load pattern loaded the target end bay and every
second bay from it. The largest end-bay sag was 1.39 to 1.76 times the
all-bays-loaded sag; permanent load reduces the ratio, and flexible supports
can shrink or reverse the effect.

**Why.** In a continuous girder the load in the next bay bends the girder
upward beside the column they share, and that hogging reduces this bay's sag.
Take that load away and the bay sags further, although nothing was added to
it. For two equal spans on pin/roller supports under a uniform load the
coefficients are exact: every span loaded 0.005416 wL⁴/EI, one span loaded
0.009151 wL⁴/EI.

**The registered cases** — rigid, level supports, equal spans:

| | |
|---|---|
| number of spans | 2, 3, 4, every target span |
| live-load arrangements | all 2ⁿ subsets of loaded spans |
| outer ends | pinned, fixed |
| load | uniform (primary), equal-resultant midspan point |
| mesh | 12, 24, 48, 96 elements per span |

Over 36 cases and 384 arrangements, the largest sag of every span comes from
loading **fewer** spans than all of them — in every case the span and its
alternates — and the end-span ratio runs from **× 1.39** (two spans, fixed ends,
point loads) to **× 1.76** (four spans, fixed ends, uniform load); × 1.44–1.76
under uniform load alone. Interior spans are not quoted: with every span
loaded they hardly sag, so their ratio is not a comparable number.

## What changes it

| choice | effect on the end-span ratio |
|---|---|
| permanent load on every span, D/(D+L) = 0 / 0.25 / 0.50 / 0.75 / 0.90 | smallest × 1.44 / 1.33 / 1.22 / 1.11 / 1.04 under uniform load: it shrinks, it does not reverse |
| interior support stiffness κ_v = k_v L³/EI from 0 to rigid | **reverses**: at κ_v ≤ 100 unloading the neighbour no longer increases the sag (ratio ≤ 1.0) |
| outer-end rotational restraint κ_θ = k_θ L/EI from 0 to rigid | × 1.44–1.76, no reversal |
| unequal spans, interior/end ratio 0.5–1.5 | no reversal; where a short interior span lifts under full load the ratio is not reported |
| interior support settlement, ± 0.5 and ± 1 of the full-load sag | the load-induced ratio is unchanged (linear); the total-sag ratio runs × 1.25–2.47 |

## What is not

**Deflection, not strength.** With every bay loaded the interior support moment
governs first yield; the displayed girder stays at 59 % of yield in every
arrangement. Linear elastic, small displacement, **quasi-static load cases**:
removing the load is not a dynamic event. A 2-D beam: no
slab, no composite action, no connection or column flexibility in the core.
Real supports are not rigid, and flexible supports can shrink or reverse the
effect.

## How it is solved

The rig's own frame solver (`fem.assemble_global`) assembles the girder; the
supports, rotational and vertical springs, prescribed settlements and
consistent nodal loads are built in `solver/pattern_beam.py`, and `fem.py` is
not modified. Deflection inside each element is recovered **exactly** for an
Euler–Bernoulli beam — the cubic Hermite interpolation of the nodal solution
plus the element's own fixed-end deflection under its uniform load — and the
bending moment the elastic check reads is EI times the curvature of that same
field. Every arrangement is a linear superposition of single-span solutions.

## Verification, not validation

| check | result |
|---|---|
| two spans, pin/roller, uniform load: 0.005416 and 0.009151 wL⁴/EI from beam theory | within 1.2e-5 |
| the same girder's support moment wL²/8 and sagging moment 9wL²/128 | within 5e-11 |
| every reported ratio, 48 → 96 elements per span | u = 2.7e-6 |
| every rigid equal-span case: the largest sag comes from fewer loaded spans, by more than 5u | 36 / 36 |
| the governing arrangement is the span and its alternates | 36 / 36 |
| the displayed girder below yield in every arrangement | 59 % |

All of it reads the same idealisation the solver reads. No measured girder is
reproduced here.

## Files

```
solver/pattern_beam.py         the study: the girder, the supports, every arrangement
solver/pattern_beam_gates.py   the gates and every shipped string
data/pattern_beam.json         the whole run
data/gates.json                the gates as they were passed
reference/                     the frozen answer and every module's SHA-256
```

The worker in the video is a generated illustration drawn 1.70 m tall on the
girder's own scale; it carries no number.
