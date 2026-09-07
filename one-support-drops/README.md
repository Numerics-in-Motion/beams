# One support drops

**Nothing is on the beam. The ground under one end goes down. Two of the four
beams do not care, and two of them reach yield — at 83.3 mm and 41.7 mm.**

```
python reproduce.py           # every gate        (a few seconds)
python reproduce.py --quick   # the four supports only
```

## The four

The same 3.0 m span, the same section, the same steel. Only the restraint at
the ends changes, and the only input is the left support moving down.

| support | settles at first yield | as a fraction of the span |
|---|---|---|
| cantilever | **never** | the beam goes down with it |
| pinned + roller | **never** | it tilts; nothing is stretched |
| fixed + roller | **83.3 mm** | L/36 |
| fixed + fixed | **41.7 mm** | L/72 |

The two that can move freely develop no internal force at all — measured at
4.6 × 10⁻¹² of the fixed case, not asserted from a reaction count. The two
held at both ends cannot move freely, so **compatibility turns the movement
into stress**.

### "Nothing is applied" is not the claim

A prescribed displacement does work and produces reactions. It is a boundary
condition of the other kind, not the absence of one. What is true is that
**no load is placed on the beam**.

## The part that is not obvious

```
d_yield = sigma_y * L^2 / (k * E * c)
```

**The second moment of area is not in it.**

| section | stiffness | depth | yields at |
|---|---|---|---|
| solid circle | ×1.0 | 71 mm | **52.5 mm** |
| I-section | ×1.5 | 90 mm | 41.7 mm |
| hollow tube | ×10.3 | 170 mm | **22.1 mm** |
| box | ×14.1 | 170 mm | **22.1 mm** |

Four sections of the **same steel area**. The stiffness spreads by 194 %, and
`d_yield × depth` is one number to 1 part in 10¹². The tube and the box are
**36 % apart in stiffness and give the identical answer**, because they are
the same depth. **The least stiff section is the most tolerant.**

A stiffer beam buys nothing against the ground moving.

## So which one is fragile?

Span squared over depth — and it runs against intuition.

| | span / depth | yields at |
|---|---|---|
| 3 m span, 300 mm deep | 10 | **12.5 mm** |
| 12 m span, 90 mm deep | 133 | **667 mm** |

**53×, and the long floppy one wins.** The stubby deep beam is the fragile
one.

## The claim is narrow on purpose

**More redundancy does not order the answer.** `scope_gate` holds the
counterexample: a 1-degree-redundant 2 m beam yields at **37 mm** while a
3-degree-redundant 12 m beam takes **667 mm**. Everything above is stated for
**one span and one section with the end restraint as the only variable**, and
a test fails if the wording ever widens.

## Yield is not the first thing that goes wrong

41.7 mm is **L/72**. A serviceability reference is nearer **L/500**, so
cracking, tilt and doors that stop closing arrive about **seven times
sooner**. `serviceability_gate` exists so that a number on screen is never
left to read as a safe allowance.

## Gates

| gate | what it forbids |
|---|---|
| `replication_gate` | **NO-GO.** a solve that does not reproduce `M = k EI d / L²` |
| `determinacy_gate` | **NO-GO.** calling a case free by counting reactions instead of solving it |
| `section_gate` | **NO-GO.** letting `d_yield` depend on I, which is the claim |
| `scope_gate` | **NO-GO.** stating "more redundancy is worse" without its counterexample |
| `serviceability_gate` | **NO-GO.** a yield settlement reported with no serviceability reference |
| `mesh_gate` | an answer that moves when the element size is halved |
| `slenderness_gate` | reporting a case with span/depth below 10 |

Worst replication error **6.15 × 10⁻¹²**; mesh drift **1.36 × 10⁻¹²** over
12 → 24 → 48 elements.

## The solver

`fem.py` and `beams_published.py` at the repository root are the published
beam study's own rig — the beam builder, the four support conditions and the
four sections — checked by SHA-256 before `reproduce.py` runs, and checked
against the study's own module by test so the extraction cannot drift. New here is `settlement.py`, which adds what the
frame solver did not have: a **prescribed support displacement**, solved in
the partitioned form with the `−K_fc u_c` term rather than by zeroing the
constrained degrees of freedom.

## Verification, not validation

Every check reads the same E, section and geometry the solver reads, and the
closed form it reproduces is the textbook one for the same idealisation — so
agreement is verification. **Validating this would need a real beam and a
real settlement. There is neither.**

## What this is not

Linear elastic, small displacement, quasi-static, ideal supports. A real
fixed end rotates a little and a real bearing has friction and gaps, and
**every one of those relieves the secondary moment** — so these are upper
bounds, not predictions. First yield, not collapse. One span, one settling
support, no creep, no relaxation, no cracking, no soil.
