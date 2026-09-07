# -*- coding: utf-8 -*-
"""One support drops. Nothing is loaded. Which beam is hurt?

THE QUESTION

The published beam study put one beam on four different supports and loaded
it. This asks
the opposite: **put no load on at all, and let one support settle.** The four
beams are the same four, the span is the same 3.0 m, the section is the same,
and the only input is a downward movement of the left-hand support.

WHAT COMES OUT

    cantilever        0     the beam just goes down with it
    pinned + roller   0     it tilts; nothing is stretched
    fixed + roller    3EI d / L^2
    fixed + fixed     6EI d / L^2

The two that are held at both ends cannot move freely, so **compatibility
turns the movement into stress.** The two that can move develop none.

**AND THE CLAIM IS NOT "REDUNDANT MEANS WORSE".** The design review was
explicit: with an arbitrary layout, more redundancy does NOT guarantee a
larger secondary action or an earlier yield. What is defensible is the
narrower statement this study actually measures -- **one span, one section,
one settling support, and only the restraint at the ends changed.**

THE SINGLE MEASURED QUANTITY

**The settlement at first yield**: the drop `d` at which the extreme fibre
first reaches the yield stress. `M = k EI d / L^2` and `sigma = M c / I`, so

    d_yield = sigma_y L^2 / (k E c)

which does not contain I at all -- **the second moment of area cancels.**
A stiffer section does not buy you any tolerance to the ground moving. What
sets it is the span squared over the depth.

"MOVEMENT" IS AN INPUT, NOT AN ABSENCE OF INPUT

The video may not say "nothing is applied". A prescribed displacement does
work and produces reactions; it is a boundary condition of the other kind,
not the lack of one. What is true is that **no load is placed on the beam.**

WHAT IS NOT CLAIMED

Linear elastic, small displacement, quasi-static. Ideal supports: a real
fixed end rotates a little and a real bearing has friction and gaps, and
every one of those RELIEVES the secondary moment -- so these numbers are the
upper bound, not a prediction. First yield, not collapse. No creep,
relaxation, cracking or soil behaviour. One span, one settling support.

WHY THE SOLVER IS NEW

`fem.solve` zeroes the constrained degrees of freedom and solves the free
partition only, so it cannot express "this support is somewhere else now".
The partitioned form with the `-K_fc u_c` term is what this needs; the
continuum solver `bracket_fem.solve_constrained` already had it and the frame
solver did not. **`fem.py` is not modified** -- its assembly is imported and
the constrained solve is done here.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import column_geometry as CG                              # noqa: E402
import beams_published as G19                             # noqa: E402
from fem import assemble_global                           # noqa: E402
from utils import load_config                             # noqa: E402

# ---- the freeze block ----------------------------------------------------
SUPPORTS = ("cantilever", "pinned", "propped", "fixed")
SECTION = "i_section"          # 019's own primary section
N_ELEMENTS = 24                # 019 uses 12; refinement is checked, not assumed
SETTLING_NODE = 0              # the LEFT support, in every case
D_PROBE = 0.010                # 10 mm, a probe only -- the answer is d_yield
# the closed form for `M = k E I d / L^2`, one end settling
K_CLOSED = {"cantilever": 0.0, "pinned": 0.0, "propped": 3.0, "fixed": 6.0}


def config():
    return load_config(G19.config_path())


def section(cfg, name=SECTION):
    for s in CG.build_all_columns(cfg):
        if s.name == name:
            return s
    raise KeyError(name)


def solve_settlement(support, sec, delta, cfg, n=N_ELEMENTS, span=None):
    """K u = 0 with one support DOF prescribed at `-delta`.

    Returns (max |M| over the beam, the full displacement vector, the span).
    The load vector is genuinely zero: every internal force here comes from
    the prescribed movement, through the `-K_fc u_c` coupling term.
    """
    E = float(cfg["material"]["youngs_modulus_pa"])
    nodes, members, span = G19.build_beam(n, span=span, cfg=cfg)
    K, cache = assemble_global(nodes, members, E, sec.area, sec.inertia)
    ndof = 3 * len(nodes)

    fixed = sorted(G19.constrained(support, len(nodes)))
    free = [d for d in range(ndof) if d not in fixed]

    dof = 3 * SETTLING_NODE + 1                       # vertical at that node
    if dof not in fixed:
        raise ValueError("node %d is not vertically supported under %r"
                         % (SETTLING_NODE, support))

    u = np.zeros(ndof)
    u[dof] = -float(delta)                            # downward
    rhs = -K[np.ix_(free, fixed)] @ u[fixed]          # no applied load at all
    u[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)

    # The moment is read exactly the way the published study reads it:
    # local end forces from the cached (T, k_local, L), entries 2 and 5.
    m_max = 0.0
    for k, (a, b) in enumerate(members):
        T, kl, _L = cache[k]
        idx = [3 * a, 3 * a + 1, 3 * a + 2, 3 * b, 3 * b + 1, 3 * b + 2]
        fl = kl @ (T @ u[idx])
        m_max = max(m_max, abs(fl[2]), abs(fl[5]))
    return m_max, u, span


def moment_closed_form(support, sec, delta, cfg, span=None):
    """k E I d / L^2 -- the textbook value this study must reproduce."""
    E = float(cfg["material"]["youngs_modulus_pa"])
    span = float(span or cfg["geometry"]["span_m"])
    return K_CLOSED[support] * E * sec.inertia * float(delta) / span ** 2


def settlement_at_first_yield(support, sec, cfg, span=None, n=N_ELEMENTS):
    """The drop at which the extreme fibre reaches yield. inf if never.

    Solved from the model, not from the formula: one probe solve fixes the
    linear constant M/d, and yield is then a division. That keeps the answer
    tied to what the solver produced.
    """
    k = K_CLOSED[support]
    if k == 0.0:
        return float("inf")
    sig_y = float(cfg["material"]["failure_stress_pa"])
    m_probe, _u, _span = solve_settlement(support, sec, D_PROBE, cfg,
                                          n=n, span=span)
    m_per_delta = m_probe / D_PROBE
    m_yield = sig_y * sec.inertia / sec.fibre_c
    return m_yield / m_per_delta
