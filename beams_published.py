# -*- coding: utf-8 -*-
"""The beam rig, as published.

This is the publishable face of the beam study "You Asked for Beams. The
Loser Won." and "Not Same Support" -- the beam builder, the four support
conditions and the four sections at one fixed steel area.

WHY IT IS A SEPARATE FILE

The study module those functions live in carries internal video numbering and
a viewer's name in its comments, neither of which may appear in a published
package. Extracting them here is the same arrangement the roof and floor rigs
use, and the extraction is not taken on trust: `tests/test_beams_published_is
_the_published_rig.py` asserts that this module and the study module return
**identical** node sets and identical constrained-degree-of-freedom sets for
every support, so a drift between the two fails a test rather than shipping.

The bodies below are the study's own, unchanged.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import column_geometry as CG                              # noqa: E402
from utils import load_config                             # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(ROOT, "config", "beams_config.json")

N_ELEMENTS = 12
SUPPORTS = ["cantilever", "pinned", "propped", "fixed"]
SUPPORT_LABEL = dict(cantilever="CANTILEVER", pinned="PINNED BOTH ENDS",
                     propped="FIXED + PINNED", fixed="FIXED BOTH ENDS")


def config_path():
    """The beam config, in the development tree or in the rig package."""
    if os.path.exists(CFG):
        return CFG
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(4):
        p = os.path.join(here, "config", "beams_config.json")
        if os.path.exists(p):
            return p
        here = os.path.dirname(here)
    raise SystemExit("beams_config.json not found from %s"
                     % os.path.abspath(__file__))


def config():
    return load_config(config_path())


def sections(cfg=None):
    """The four sections, from the column builder, unchanged."""
    return CG.build_all_columns(cfg or config())


def build_beam(n=N_ELEMENTS, span=None, cfg=None):
    """A straight beam of `n` elements, with an exact node at midspan."""
    span = span or float(cfg["geometry"]["span_m"])
    if n % 2:
        n += 1
    xs = np.linspace(0.0, span, n + 1)
    nodes = np.column_stack([xs, np.zeros(n + 1)])
    return nodes, [(i, i + 1) for i in range(n)], span


def constrained(support, n_nodes):
    """Which degrees of freedom the support condition holds.

    dof order per node is (u_x, u_y, theta). A "fixed" end holds all three; a
    pin holds both translations; a roller holds only the vertical.

    One horizontal restraint is enough for equilibrium in every case, and
    adding a second at the other end would put the beam in axial restraint --
    which is a different structure (it would develop membrane action) and not
    what anybody means by these names.
    """
    last = n_nodes - 1
    if support == "cantilever":
        return {0, 1, 2}                                  # built in, far end free
    if support == "pinned":
        return {0, 1, 3 * last + 1}                       # pin + roller
    if support == "propped":
        return {0, 1, 2, 3 * last + 1}                    # built in + roller
    if support == "fixed":
        return {0, 1, 2, 3 * last + 1, 3 * last + 2}      # built in + built in
    raise ValueError(support)
