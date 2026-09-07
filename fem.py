"""2D frame (Euler-Bernoulli beam-column) finite element analysis.

Each node has 3 dof: (u_x, u_y, theta_z).  A frame model is used (rather than a
pure truss) because the un-braced square panel is a mechanism under a truss
model; frame elements carry the un-braced panel through bending stiffness so all
three topologies can be compared with the SAME element formulation (spec 5.1,
recommended case B).  Justification is recorded in repot.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MemberResult:
    index: int
    node_i: int
    node_j: int
    length: float
    axial_force: float      # N, tension positive [N]
    axial_stress: float     # N / A [Pa]
    max_moment: float       # max|M| at the two ends [N m]
    bending_stress: float   # |M| c / I [Pa]
    combined_stress: float  # |axial_stress| + bending_stress [Pa]


def local_stiffness(E: float, A: float, I: float, L: float) -> np.ndarray:
    """6x6 local frame stiffness for dof order [ui,vi,ti, uj,vj,tj]."""
    ea = E * A / L
    ei = E * I
    L2, L3 = L * L, L * L * L
    k = np.array([
        [ ea,        0,          0,        -ea,       0,          0        ],
        [ 0,   12*ei/L3,    6*ei/L2,        0, -12*ei/L3,    6*ei/L2       ],
        [ 0,    6*ei/L2,     4*ei/L,        0,  -6*ei/L2,     2*ei/L       ],
        [-ea,        0,          0,         ea,       0,          0        ],
        [ 0,  -12*ei/L3,   -6*ei/L2,        0,  12*ei/L3,   -6*ei/L2       ],
        [ 0,    6*ei/L2,     2*ei/L,        0,  -6*ei/L2,     4*ei/L       ],
    ], dtype=float)
    return k


def transform_matrix(c: float, s: float) -> np.ndarray:
    """6x6 transformation from global to local (c=cos, s=sin of member angle)."""
    T = np.zeros((6, 6))
    R = np.array([[ c, s, 0],
                  [-s, c, 0],
                  [ 0, 0, 1]])
    T[0:3, 0:3] = R
    T[3:6, 3:6] = R
    return T


def element_global_stiffness(E, A, I, xi, yi, xj, yj):
    """Return (Kglobal 6x6, T, L, c, s) for one frame element."""
    dx, dy = xj - xi, yj - yi
    L = np.hypot(dx, dy)
    c, s = dx / L, dy / L
    kl = local_stiffness(E, A, I, L)
    T = transform_matrix(c, s)
    Kg = T.T @ kl @ T
    return Kg, T, L, c, s


def assemble_global(nodes, members, E, A, I):
    """Assemble the global stiffness matrix (3 dof/node) for the given member
    list.  Returns (K, elem_cache) where elem_cache[m] = (T, kl, L)."""
    n = len(nodes)
    ndof = 3 * n
    K = np.zeros((ndof, ndof))
    cache = {}
    for m, (i, j) in enumerate(members):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        Kg, T, L, c, s = element_global_stiffness(E, A, I, xi, yi, xj, yj)
        dofs = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
        for a in range(6):
            for b in range(6):
                K[dofs[a], dofs[b]] += Kg[a, b]
        cache[m] = (T, local_stiffness(E, A, I, L), L)
    return K, cache


def constrained_dofs(structure) -> list:
    """Return the list of fixed global dof indices for the BCs.

    If `structure.fixed_nodes` is non-empty (towers: every base-column node
    on a rigid foundation), ALL 3 dof are fixed at each of those nodes and
    the pin/roller scheme below is not used at all.

    Otherwise (shapes and bridges):
    pin   : u_x, u_y fixed (rotation free)
    roller: u_y fixed (u_x free) -- unless `structure.roller_is_pin` is True,
            in which case this support is also a full pin (u_x, u_y fixed).
            The latter is used for a two-hinged arch, which needs horizontal-
            thrust restraint at both ends; every other structure leaves this
            attribute unset/False and gets the original pin+roller scheme.
    """
    fixed_nodes = getattr(structure, "fixed_nodes", None)
    if fixed_nodes:
        fixed = []
        for n in fixed_nodes:
            fixed.extend([3 * n, 3 * n + 1, 3 * n + 2])
        return sorted(set(fixed))

    fixed = []
    p = structure.pin_node
    r = structure.roller_node
    fixed.extend([3 * p, 3 * p + 1])
    if getattr(structure, "roller_is_pin", False):
        fixed.extend([3 * r, 3 * r + 1])
    else:
        fixed.append(3 * r + 1)
    return sorted(set(fixed))


def solve(structure, members, E, A, I, load_vector, singular_tol=1e-12):
    """Solve K u = F with the structure's BCs applied to the given member set.

    Returns dict with keys:
        u            : full displacement vector (ndof,)
        member_results : list[MemberResult]
        reactions    : reaction vector at constrained dofs (ndof,)
        singular     : bool (True if the free system is (near) singular)
        min_pivot_ratio : smallest/largest eigenvalue ratio of K_ff
    """
    n = len(structure.nodes)
    ndof = 3 * n
    K, cache = assemble_global(structure.nodes, members, E, A, I)

    fixed = set(constrained_dofs(structure))
    free = [d for d in range(ndof) if d not in fixed]

    Kff = K[np.ix_(free, free)]
    Ff = load_vector[free]

    # singularity / mechanism check via symmetric eigenvalues
    if len(free) == 0:
        singular = True
        u = np.zeros(ndof)
        ratio = 0.0
        return dict(u=u, member_results=[], reactions=np.zeros(ndof),
                    singular=True, min_pivot_ratio=0.0)

    evals = np.linalg.eigvalsh(Kff)
    lam_max = evals[-1]
    lam_min = evals[0]
    ratio = lam_min / lam_max if lam_max > 0 else 0.0
    singular = (lam_min <= singular_tol * lam_max) or (lam_min <= 0.0)

    u = np.zeros(ndof)
    if not singular:
        uf = np.linalg.solve(Kff, Ff)
        u[free] = uf

    reactions = K @ u - load_vector  # nonzero only at fixed dofs

    member_results = []
    if not singular:
        for m, (i, j) in enumerate(members):
            T, kl, L = cache[m]
            ue = u[[3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]]
            ul = T @ ue
            fl = kl @ ul            # local end forces
            axial = fl[3]           # tension positive (force on node j, local x)
            Mi, Mj = fl[2], fl[5]
            member_results.append((m, i, j, L, axial, Mi, Mj))

    return dict(u=u, member_results=member_results, reactions=reactions,
                singular=singular, min_pivot_ratio=float(ratio))
