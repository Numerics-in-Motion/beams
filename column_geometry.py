"""Geometry generation for the four column cross-section candidates
(Solid Circular Rod / Circular Hollow Tube / Box Tube / I-Section).

All four columns share the SAME length, target cross-sectional area
(hence same material volume and mass), envelope, end conditions, and
initial-imperfection rule (spec section 4). They differ only in how the
shared area is distributed across the cross-section -- solved by
`section_properties.py` from a fixed proportioning rule chosen before any
simulation is run (spec section 5).

Each column is meshed into `n_elements` two-node frame elements along its
length (y axis, matching the tower convention: y = height/length, x =
lateral). The initial geometric imperfection `w0(y) = delta0*sin(pi*y/L)`
is baked directly into the node x-coordinates -- this needs zero new code
in `fem.py`, since `fem.element_global_stiffness` already computes each
element's own angle from its node coordinates and never assumed a member
is axis-aligned.

Reused from geometry.py: `Structure` (node/member storage, single
area/inertia/fibre_c per structure -- already exactly the right shape for
a column, which has one uniform section along its whole length).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from geometry import Structure
import section_properties as sp


@dataclass
class Column(Structure):
    """A `Structure` plus the column-specific metadata the nonlinear
    solver and renderer need: which nodes are the pinned base / laterally
    -restrained top, the section's own dimensioned properties, and the
    imperfection amplitude actually used (for reporting)."""
    bottom_node: int = 0
    top_node: int = -1
    section: dict = field(default_factory=dict)
    imperfection_amplitude: float = 0.0


def build_column(name: str, section: dict, length: float, n_elements: int,
                 imperfection_ratio: float) -> Column:
    """Build one meshed, imperfect column for the given (already-solved)
    `section` dict from `section_properties.py`."""
    delta0 = imperfection_ratio * length
    n_nodes = n_elements + 1
    ys = np.linspace(0.0, length, n_nodes)
    xs = delta0 * np.sin(np.pi * ys / length)
    nodes = np.column_stack([xs, ys])
    members = [(i, i + 1) for i in range(n_elements)]

    col = Column(
        name=name, nodes=nodes, members=members,
        width=float(np.max(xs) - np.min(xs)), height=length,
        area=section["area"], inertia=section["I_min"],
        # the fibre distance that goes WITH I_min, not the other one
        fibre_c=sp.fibre_for(section, section["I_min"]),
        bottom_node=0, top_node=n_elements,
        section=section, imperfection_amplitude=delta0,
    )
    # Cosmetic only (visualization.draw_structure icon selection) -- the
    # real boundary conditions used by the solver come from
    # `column_fem.column_constrained_dofs`, entirely independent of this.
    # A ground/"fixed" icon at the base and a downward load arrow at the
    # top (both mechanisms already built into draw_structure) read far
    # more clearly for a column than the pin/roller icon pair meant for
    # horizontal shapes/bridges.
    col.fixed_nodes = [0]
    col.load_nodes = [n_elements]
    col.load_weights = [1.0]
    return col


def build_all_columns(config: dict) -> list[Column]:
    """Build the four column candidates from `config["geometry"]`, using
    the fixed proportioning rules documented in report.md -- nothing here
    is tuned after seeing a simulation result."""
    g = config["geometry"]
    L = g["length_m"]
    A = g["target_area_m2"]
    B_max = g["envelope_width_m"]
    H_max = g["envelope_height_m"]
    n_elem = g["n_elements"]
    imp_ratio = g["imperfection_ratio"]

    sections = dict(
        solid_circle=sp.solid_circle_properties(A),
        hollow_tube=sp.hollow_tube_properties(A, g["tube_diameter_ratio"] * B_max),
        box=sp.box_properties(A, g["box_width_ratio"] * B_max),
        i_section=sp.i_section_properties(
            A, g["i_section_height_ratio"] * H_max,
            g["i_section_flange_width_ratio"] * (g["i_section_height_ratio"] * H_max),
            g["i_section_tf_tw_ratio"]),
    )

    for name, sec in sections.items():
        if sec["outer_dim"] > max(B_max, H_max) + 1e-9:
            raise ValueError(f"{name}: outer dimension {sec['outer_dim']} "
                             f"exceeds envelope {max(B_max, H_max)}")

    return [build_column(name, sec, L, n_elem, imp_ratio)
           for name, sec in sections.items()]
