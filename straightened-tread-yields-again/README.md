# Straightened tread, earlier renewed plastic flow

An **idealized** aluminium strip is bent down beyond first yield, unloaded,
pressed upward, unloaded and then retested under the same downward central
patch. Its support fixture is a captive pin/roller laboratory model, not a
loose installed stair bearing. The strip is 0.8 m long, 0.28 m wide and
12 mm thick. The central full-width patch spans 0.1 m. Self-weight is retained
at every stage.

The registered 51-value sweep varies only the reverse press force. The
closest return to the virgin **gravity-loaded** shape is `q=1.175`, selected
from bow alone. Its full-span residual deviation is **0.051433 mm**, against
**0.671483 mm** for the bent, unstraightened control. First **new plastic
flow** on downward reload begins at **3.297 kN**, or **0.84445×** the virgin
first-yield patch load of **3.904 kN**. This is not a drop in material yield
stress, elastic stiffness, collapse load, fatigue life or safe load.

The six fixed-q sensitivity legs (isotropic/kinematic hardening assumptions,
central-patch widths 0.09/0.11 m, and yield stresses 99/121 MPa) preserve
earlier renewed flow in this model. They do not validate a particular alloy.
The [Kosel et al. 2011 primary paper](https://www.lanem.si/papers/brojan/2011_JMEP.pdf)
supports the general springback/residual-stress concept, but does **not**
validate this patch/gravity/re-yield threshold.

## Reproduce

Install the dependencies in `requirements.txt`, then run:

```text
python reproduce.py --quick
python reproduce.py
```

The full command recomputes the 51-cell baseline and three independent
51-cell refinements, six fixed-q sensitivity legs with their own controls and
refinements, independent pure-bending/reverse-stress-block checks, signed
history symmetry, equilibrium, scope and event brackets—**252 histories**.
It runs in a temporary isolated directory and compares the result with
`reference/reference_canonical.json`, leaving the published data unchanged.
`solver/` contains the full source; `protocol/` preserves both frozen
pre-result documents; `data/` preserves the original records.

The summed refinement shifts are **operational estimates, not certified
bounds**. Independent dense analytic integration gives the unstraightened
control bow **0.672006 mm**, just outside its saved operational upper estimate
**0.671974 mm**. This small discrepancy does not change the registered
selection or editorial materiality result, but prohibits treating the
interval as guaranteed error coverage. The video uses calculated physical
motion with an explicitly labelled vertical magnification; its generated
press image is a setup illustration, not evidence.

No installed-stair rating, repair instruction or safety recommendation follows
from this model.
