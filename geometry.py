"""Geometry generation for the three competing structures.

All three structures occupy the SAME bounding rectangle (width_m x height_m),
share the SAME support layout (pin at bottom-left, roller at bottom-right) and
the SAME load point (top centre).  They differ only in cell topology:

    triangle : a rectangular node grid, fully triangulated (each quad cell gets
               one diagonal) -> load carried mostly as member axial force.
    square   : the same rectangular node grid with horizontal + vertical members
               only (quadrilateral cells, no diagonals) -> load carried by frame
               bending, the classic "weak" un-braced panel.
    hexagon  : a honeycomb (hexagonal-cell) lattice filling the same rectangle.

Member cross-sectional areas are chosen per structure so that every structure
uses exactly the same total material volume  (sum(L_i * A_i) = target volume),
with a single uniform area inside each structure (spec 4.2).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from utils import circular_section_properties


@dataclass
class Structure:
    name: str
    nodes: np.ndarray                 # (Nn, 2) coordinates [m]
    members: list                     # list of (i, j) node index pairs
    width: float
    height: float
    # assigned after volume matching:
    area: float = 0.0                 # uniform member area [m^2]
    inertia: float = 0.0              # section second moment [m^4]
    fibre_c: float = 0.0              # outer-fibre distance [m]
    # WHICH SECTION SHAPE the inertia above describes. `assign_matched_area`
    # fills in a SOLID round bar, because that is the rule the static
    # strength studies use. For a dynamic problem that is the wrong shape and it is
    # the difference between a real building and a floppy one, so the choice
    # has to be made out loud: `sections.apply_sections` stamps it here and
    # `dynamics.simulate_earthquake` refuses a structure that never made it.
    # A shipped study once had towers that were solid bars because nobody
    # chose; leaving this None is what makes that impossible now.
    section_model: str = None         # "solid" | "tube", set explicitly
    d_over_t: float = None            # tube proportion, when section_model is
                                      # "tube"
    # boundary conditions (filled by set_boundary_conditions):
    pin_node: int = -1
    roller_node: int = -1
    # if True, the "roller" node is actually a second pin (both translations
    # fixed) -- used for a two-hinged arch, which must resist horizontal
    # thrust at both supports. Default False reproduces the original
    # pin+roller scheme exactly for every existing structure.
    roller_is_pin: bool = False
    # nodes with ALL 3 dof fixed (e.g. every base-column node of a tower on
    # a rigid foundation). If non-empty, this REPLACES the pin/roller scheme
    # entirely for BC purposes; empty (default) reproduces every existing
    # structure's behaviour exactly.
    fixed_nodes: list = field(default_factory=list)
    load_nodes: list = field(default_factory=list)     # node indices
    load_weights: list = field(default_factory=list)   # fraction of total load

    # ---- derived quantities -------------------------------------------------
    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_members(self) -> int:
        return len(self.members)

    def member_lengths(self) -> np.ndarray:
        p = self.nodes
        return np.array([np.linalg.norm(p[j] - p[i]) for (i, j) in self.members])

    def total_length(self) -> float:
        return float(self.member_lengths().sum())

    def total_volume(self) -> float:
        return float((self.member_lengths() * self.area).sum())

    def total_mass(self, density: float) -> float:
        return self.total_volume() * density


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _dedupe_nodes(coords: list, tol: float = 1e-9):
    """Merge coincident points; return (unique_coords array, index_map dict)."""
    unique = []
    keys = {}

    def key(p):
        return (round(p[0] / tol), round(p[1] / tol))

    index_of = []
    for p in coords:
        k = key(p)
        if k not in keys:
            keys[k] = len(unique)
            unique.append(p)
        index_of.append(keys[k])
    return np.array(unique, dtype=float), index_of


def _unique_members(pairs) -> list:
    """Remove duplicate and zero-length (i==i) member definitions."""
    seen = set()
    out = []
    for i, j in pairs:
        if i == j:
            continue
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        out.append((key[0], key[1]))
    return out


# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------
def generate_square(nx: int, ny: int, W: float, H: float) -> Structure:
    """Rectangular grid, horizontal + vertical members only (quad cells)."""
    xs = np.linspace(0.0, W, nx + 1)
    ys = np.linspace(0.0, H, ny + 1)
    coords = [(x, y) for y in ys for x in xs]
    nodes = np.array(coords, dtype=float)

    def idx(ix, iy):
        return iy * (nx + 1) + ix

    members = []
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            if ix < nx:
                members.append((idx(ix, iy), idx(ix + 1, iy)))   # horizontal
            if iy < ny:
                members.append((idx(ix, iy), idx(ix, iy + 1)))   # vertical
    return Structure("square", nodes, _unique_members(members), W, H)


def generate_triangle(nx: int, ny: int, W: float, H: float) -> Structure:
    """Rectangular grid, fully triangulated (each quad cell + one diagonal)."""
    base = generate_square(nx, ny, W, H)

    def idx(ix, iy):
        return iy * (nx + 1) + ix

    members = list(base.members)
    for iy in range(ny):
        for ix in range(nx):
            # alternate diagonal direction so the truss is symmetric-ish
            if (ix + iy) % 2 == 0:
                members.append((idx(ix, iy), idx(ix + 1, iy + 1)))
            else:
                members.append((idx(ix + 1, iy), idx(ix, iy + 1)))
    return Structure("triangle", base.nodes, _unique_members(members), W, H)


def generate_hexagon(cols: int, rows: int, W: float, H: float) -> Structure:
    """Honeycomb (flat-top hexagonal cells) rescaled to fill W x H exactly.

    Boundary cells may be incomplete; the generation rule is: tile flat-top
    hexagons on an offset column grid, keep every hexagon edge, then rescale the
    raw vertex cloud so its bounding box maps onto [0,W] x [0,H].
    """
    R = 1.0  # circumradius cancels out after rescaling
    ang = np.deg2rad([0, 60, 120, 180, 240, 300])
    raw = []
    for c in range(cols):
        for r in range(rows):
            cx = 1.5 * R * c
            cy = np.sqrt(3.0) * R * r + (np.sqrt(3.0) * R / 2.0 if c % 2 else 0.0)
            verts = [(cx + R * np.cos(a), cy + R * np.sin(a)) for a in ang]
            for k in range(6):
                raw.append((verts[k], verts[(k + 1) % 6]))

    # collect vertices + edges with de-duplication
    all_pts = []
    for a, b in raw:
        all_pts.append(a)
        all_pts.append(b)
    nodes, imap = _dedupe_nodes(all_pts, tol=1e-6)

    members = []
    for e in range(len(raw)):
        members.append((imap[2 * e], imap[2 * e + 1]))
    members = _unique_members(members)

    # rescale bounding box -> [0,W] x [0,H]
    xmin, ymin = nodes.min(axis=0)
    xmax, ymax = nodes.max(axis=0)
    nodes = nodes.copy()
    nodes[:, 0] = (nodes[:, 0] - xmin) / (xmax - xmin) * W
    nodes[:, 1] = (nodes[:, 1] - ymin) / (ymax - ymin) * H

    return Structure("hexagon", nodes, members, W, H)


# ---------------------------------------------------------------------------
# boundary conditions & material volume
# ---------------------------------------------------------------------------
def _nearest_node(nodes: np.ndarray, target) -> int:
    d = np.linalg.norm(nodes - np.asarray(target), axis=1)
    return int(np.argmin(d))


def _split_load_at_x(nodes: np.ndarray, candidate_idx: np.ndarray, x_target: float):
    """Given a set of candidate node indices, return (load_nodes, load_weights)
    such that the resultant of a unit load acts exactly at x = x_target: an
    exact-match node if one exists, otherwise a linear split between the
    nearest candidate left and right of x_target (or the single nearest
    candidate if all lie on one side). Shared by every structure's boundary-
    condition setup (shapes and bridges alike) so the "resultant acts at
    centre" guarantee is implemented in exactly one place."""
    xs = nodes[candidate_idx, 0]

    exact = candidate_idx[np.isclose(xs, x_target, atol=1e-6)]
    if len(exact) > 0:
        return [int(exact[0])], [1.0]

    left = candidate_idx[xs < x_target]
    right = candidate_idx[xs >= x_target]
    if len(left) > 0 and len(right) > 0:
        li = int(left[np.argmax(nodes[left, 0])])   # closest-from-left
        ri = int(right[np.argmin(nodes[right, 0])])  # closest-from-right
        xl, xr = nodes[li, 0], nodes[ri, 0]
        if abs(xr - xl) < 1e-9:
            return [li], [1.0]
        # linear interpolation so resultant acts at x_target
        wr = (x_target - xl) / (xr - xl)
        wl = 1.0 - wr
        return [li, ri], [wl, wr]

    ni = int(candidate_idx[np.argmin(np.abs(xs - x_target))])
    return [ni], [1.0]


def set_boundary_conditions(s: Structure) -> None:
    """Assign pin, roller and load nodes identically for every structure:
    pin  -> nearest node to bottom-left  (0, 0)
    roller-> nearest node to bottom-right (W, 0)
    load -> top-region node(s) nearest to top centre (W/2, H); if two top nodes
            straddle the centre they share the load weighted so the resultant
            acts at x = W/2 (spec 4.4)."""
    W, H = s.width, s.height
    s.pin_node = _nearest_node(s.nodes, (0.0, 0.0))
    s.roller_node = _nearest_node(s.nodes, (W, 0.0))

    # candidate top nodes: within 2% of H of the maximum node height
    ymax = s.nodes[:, 1].max()
    top_mask = s.nodes[:, 1] >= ymax - 0.02 * H
    top_idx = np.where(top_mask)[0]
    s.load_nodes, s.load_weights = _split_load_at_x(s.nodes, top_idx, W / 2.0)


def assign_matched_area(s: Structure, target_volume: float) -> None:
    """Choose the uniform member area so sum(L_i * A) == target_volume."""
    L = s.total_length()
    area = target_volume / L
    r, inertia, c = circular_section_properties(area)
    s.area = area
    s.inertia = inertia
    s.fibre_c = c


def build_all_structures(config: dict) -> list:
    """Generate the three structures, assign BCs and matched areas."""
    g = config["geometry"]
    W = float(g["width_m"])
    H = float(g["height_m"])
    nx = int(g["cell_count_x"])
    ny = int(g["cell_count_y"])
    hc = int(g.get("hex_cols", 4))
    hr = int(g.get("hex_rows", 2))
    target_v = float(g["target_material_volume_m3"])

    structures = [
        generate_triangle(nx, ny, W, H),
        generate_square(nx, ny, W, H),
        generate_hexagon(hc, hr, W, H),
    ]
    for s in structures:
        set_boundary_conditions(s)
        assign_matched_area(s, target_v)
    return structures
