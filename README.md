# Beams

Beams, solved. One repository for the rig, and one folder per study inside it.

```
fem.py                 the frame solver, ONCE
beams_published.py     the four supports and the beam builder
column_geometry.py     the four sections, at one fixed steel area
one-support-drops/     what the ground moving does to a beam nobody loaded
less-live-load-next-door/  what the next bay's live load does to this one
```

## Why one repository and not one per study

This channel published twenty-seven repositories, one per video, and the same
solver was copied into repository after repository. **Every copy is a place
that rots on its own**, and nothing notices when one drifts from another.

Here the rig's modules sit at the root, once. Each study folder ships only
what is new to it, and its `reproduce.py` checks the root modules and its own
against a frozen SHA-256 before running anything. A change to a shared module
that moves an old study's numbers makes **that study's own check fail**,
which is exactly what twenty-seven copies could not do.

The earlier per-video repositories are archived, not deleted: everything ever
linked from a caption still resolves.

## Running one

```
cd one-support-drops
python reproduce.py           # the whole study and every gate
python reproduce.py --quick   # the headline only
```

`requirements.txt` is numpy and scipy. Nothing else is needed.

## The studies

| folder | asks | answers with |
|---|---|---|
| [`one-support-drops`](one-support-drops/) | a support settles and nothing is loaded — which beam is hurt? | the settlement at first yield |
| [`less-live-load-next-door`](less-live-load-next-door/) | a continuous girder loses the live load on the next bay — does this bay sag less, or more? | the largest sag of the loaded bay |

## What the rig is

A 2D frame solver: three degrees of freedom a node, Euler–Bernoulli elements,
linear elastic. The four support conditions (cantilever, pinned + roller,
fixed + roller, fixed + fixed) and the four sections at one fixed steel area
come from the published beam study and are not re-typed here.

**Verification, not validation.** Every check inside a study reads the same
material and geometry its solver reads. Validating any of it would need real
beams and real measurements, and there are none.
