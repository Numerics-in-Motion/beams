"""Spinoff 048: where one ideal rotational restraint works.

One uniform simply-supported Euler--Bernoulli beam carries a symmetric UDL.
A point rotational spring connects the beam rotation to a fixed reference at
one of four pre-registered nodes.  It is not an internal splice.

The finite-element result is checked against an independent force-method
solution.  Maximum sag is found from every exact polynomial stationary point,
not from a plotting grid.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy.linalg import solve
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output_spinoff" / "Spinoff_048" / "data"

E = 200.0e9
I = 4.0e-4
L = 6.0
Q = 10.0e3
KAPPA = 1.0
POSITIONS = (1 / 8, 1 / 4, 3 / 8, 1 / 2)
MIRRORS = (7 / 8, 3 / 4, 5 / 8)
MESHES = (64, 128, 256)
PRODUCTION_MESH = 256
DEFLECTION_MAGNIFICATION = 100
U_FLOOR = 1.0e-10
PARENT_FROZEN = {"all": 0.005416121605827781,
                 "one": 0.009150556074890303}


def element_stiffness(EI: float, h: float) -> np.ndarray:
    return EI / h ** 3 * np.array([
        [12, 6 * h, -12, 6 * h],
        [6 * h, 4 * h * h, -6 * h, 2 * h * h],
        [-12, -6 * h, 12, -6 * h],
        [6 * h, 2 * h * h, -6 * h, 4 * h * h],
    ], dtype=float)


def element_load(q: float, h: float) -> np.ndarray:
    """Consistent downward UDL; v is positive upward."""
    return np.array([-q * h / 2, -q * h * h / 12,
                     -q * h / 2, q * h * h / 12], dtype=float)


def _candidate_node(a_frac: float, n: int) -> int:
    raw = float(a_frac) * n
    node = int(round(raw))
    if abs(raw - node) > 1.0e-15 * n:
        raise ValueError("registered restraint coordinate is not an exact node")
    return node


def validate_registered_inputs(a_frac: float, kappa: float, n: int) -> None:
    if a_frac not in POSITIONS + MIRRORS:
        raise ValueError("coordinate is outside the registered candidates and mirrors")
    if not (kappa >= 0.0 or math.isinf(kappa)):
        raise ValueError("rotational stiffness must be non-negative")
    _candidate_node(a_frac, n)


def _banded_ldlt_high_precision(A: np.ndarray, b: np.ndarray,
                                bandwidth: int = 3) -> np.ndarray:
    """Solve one symmetric positive band system at 60 decimal digits.

    The fourth-order beam matrix has a condition number that makes 1e-10
    equilibrium and mirror gates inappropriate for a dense float64 solve at
    256 elements.  A narrow-band LDLT is cheap enough to keep the registered
    mesh and the pre-frozen tolerances unchanged.
    """
    n = len(b)
    with mp.workdps(60):
        d = [mp.mpf("0") for _ in range(n)]
        lower = [dict() for _ in range(n)]
        for i in range(n):
            lo = max(0, i - bandwidth)
            for j in range(lo, i):
                s = mp.mpf(str(A[i, j]))
                k0 = max(0, i - bandwidth, j - bandwidth)
                for k in range(k0, j):
                    if k in lower[i] and k in lower[j]:
                        s -= lower[i][k] * d[k] * lower[j][k]
                lower[i][j] = s / d[j]
            diag = mp.mpf(str(A[i, i]))
            for k, lik in lower[i].items():
                diag -= lik * lik * d[k]
            d[i] = diag
        y = [mp.mpf("0") for _ in range(n)]
        for i in range(n):
            y[i] = mp.mpf(str(b[i])) - sum((v * y[j] for j, v in lower[i].items()),
                                            mp.mpf("0"))
        z = [y[i] / d[i] for i in range(n)]
        x = [mp.mpf("0") for _ in range(n)]
        for i in range(n - 1, -1, -1):
            s = z[i]
            for j in range(i + 1, min(n, i + bandwidth + 1)):
                if i in lower[j]:
                    s -= lower[j][i] * x[j]
            x[i] = s
    return np.array([float(v) for v in x])


def assemble_single(a_frac: float, kappa: float, n: int,
                    *, validate: bool = True, q_by_element=None):
    if validate:
        validate_registered_inputs(a_frac, kappa, n)
    node = _candidate_node(a_frac, n)
    # Solve the registered dimensionless system on x/L in [0, 1].  Carrying
    # E=200 GPa and metre DOFs in one matrix needlessly loses symmetry and
    # equilibrium digits; dimensional scales are restored in solve_case().
    h = 1.0 / n
    ndof = 2 * (n + 1)
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    ke = element_stiffness(1.0, h)
    if q_by_element is None:
        q_by_element = np.ones(n)
    q_by_element = np.asarray(q_by_element, dtype=float)
    if q_by_element.shape != (n,):
        raise ValueError("one UDL value is required per element")
    if validate:
        validate_symmetric_load(q_by_element)
    for e, qe in enumerate(q_by_element):
        d = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(d, d)] += ke
        F[d] += element_load(float(qe), h)
    K_beam = K.copy()
    fixed = {0, 2 * n}
    k_theta = math.inf
    if math.isinf(kappa):
        fixed.add(2 * node + 1)
    else:
        k_theta = kappa
        K[2 * node + 1, 2 * node + 1] += k_theta
    free = np.array([j for j in range(ndof) if j not in fixed], dtype=int)
    u = np.zeros(ndof)
    u[free] = _banded_ldlt_high_precision(K[np.ix_(free, free)], F[free])
    R = K @ u - F
    rot_dof = 2 * node + 1
    with mp.workdps(60):
        beam_moment_bar = (sum((mp.mpf(str(K_beam[rot_dof, j])) *
                                mp.mpf(str(u[j])) for j in range(ndof)), mp.mpf("0"))
                           - mp.mpf(str(F[rot_dof])))
        spring_balance_bar = (mp.mpf("0") if math.isinf(kappa) else
                              beam_moment_bar + mp.mpf(str(k_theta)) *
                              mp.mpf(str(u[rot_dof])))
    return dict(K=K, F=F, u=u, R=R, n=n, h=h, node=node,
                a_frac=float(a_frac), kappa=float(kappa), k_theta=k_theta,
                q_by_element=q_by_element, fixed=sorted(fixed), EI_field=1.0,
                beam_moment_bar=float(beam_moment_bar),
                spring_balance_bar=float(spring_balance_bar))


def _element_coefficients(v1, t1, v2, t2, q, h, EI=E * I):
    """Ascending powers of xi for the exact element deflection."""
    b = -q * h ** 4 / (24.0 * EI)
    return np.array([
        v1,
        h * t1,
        -3 * v1 - 2 * h * t1 + 3 * v2 - h * t2 + b,
        2 * v1 + h * t1 - 2 * v2 + h * t2 - 2 * b,
        b,
    ], dtype=float)


def _polyval(c, x):
    return sum(float(ci) * x ** i for i, ci in enumerate(c))


def exact_extrema(solution) -> dict:
    u = solution["u"]
    h = solution["h"]
    qv = solution["q_by_element"]
    points = []
    for e, q in enumerate(qv):
        c = _element_coefficients(u[2 * e], u[2 * e + 1],
                                  u[2 * e + 2], u[2 * e + 3], q, h,
                                  EI=solution.get("EI_field", E * I))
        dc = np.array([c[1], 2 * c[2], 3 * c[3], 4 * c[4]])
        roots = np.roots(dc[::-1])
        xis = [0.0, 1.0]
        xis += [float(z.real) for z in roots
                if abs(z.imag) < 1.0e-10 and -1.0e-12 <= z.real <= 1 + 1.0e-12]
        for xi in xis:
            xi = min(1.0, max(0.0, xi))
            points.append((e * h + xi * h, _polyval(c, xi)))
    x_min, v_min = min(points, key=lambda p: p[1])
    return {"sag_m": float(-v_min), "x_at_max_m": float(x_min),
            "min_v_m": float(v_min)}


def sample_field(solution, samples=1001):
    xs_bar = np.linspace(0.0, 1.0, samples)
    u, h, n, qv = (solution[k] for k in ("u", "h", "n", "q_by_element"))
    vs = np.empty_like(xs_bar)
    for j, x in enumerate(xs_bar):
        e = min(int(x / h), n - 1)
        xi = (x - e * h) / h
        c = _element_coefficients(u[2 * e], u[2 * e + 1],
                                  u[2 * e + 2], u[2 * e + 3], qv[e], h,
                                  EI=solution.get("EI_field", E * I))
        vs[j] = _polyval(c, xi)
    return xs_bar * L, vs * (Q * L**4 / (E * I))


def force_method(a_frac: float, kappa: float) -> dict:
    """Independent simply-supported compatibility solution, downward positive."""
    a = a_frac * L
    EI = E * I

    def v0(x):
        return Q * x * (L ** 3 - 2 * L * x * x + x ** 3) / (24 * EI)

    theta0 = Q * (L ** 3 - 6 * L * a * a + 4 * a ** 3) / (24 * EI)
    alpha = (a ** 3 + (L - a) ** 3) / (3 * EI * L ** 2)
    if math.isinf(kappa):
        couple = -theta0 / alpha
        theta = 0.0
    else:
        k_theta = kappa * EI / L
        couple = -k_theta * theta0 / (1.0 + k_theta * alpha)
        theta = theta0 + couple * alpha

    c1 = L / 3 - a + a * a / (2 * L)

    def g(x):
        if x <= a:
            return (x ** 3 / (6 * L) + c1 * x) / EI
        return (x ** 3 / (6 * L) - x * x / 2
                + (c1 + a) * x - a * a / 2) / EI

    def v(x):
        return v0(x) + couple * g(x)

    candidates = [(0.0, v(0.0)), (L, v(L)), (a, v(a))]
    for lo, hi in ((0.0, a), (a, L)):
        if hi - lo > 0:
            opt = minimize_scalar(lambda z: -v(z), bounds=(lo, hi),
                                  method="bounded", options={"xatol": 1e-13})
            candidates.append((float(opt.x), float(v(opt.x))))
    x_max, sag = max(candidates, key=lambda p: p[1])
    return {"sag_m": float(sag), "x_at_max_m": float(x_max),
            "theta_rad": float(theta), "couple_nm": float(couple),
            "alpha_rad_per_nm": float(alpha)}


def solve_case(a_frac: float, kappa: float = KAPPA, n: int = PRODUCTION_MESH,
               *, validate: bool = True) -> dict:
    s = assemble_single(a_frac, kappa, n, validate=validate)
    ex = exact_extrema(s)
    node = s["node"]
    v_scale = Q * L**4 / (E * I)
    theta_scale = Q * L**3 / (E * I)
    moment_scale = Q * L**2
    theta_bar = float(s["u"][2 * node + 1])
    theta = theta_bar * theta_scale
    spring_moment = (0.0 if math.isinf(kappa) else
                     -s["k_theta"] * theta_bar * moment_scale)
    vertical_reaction = float(s["R"][0] + s["R"][2 * n]) * Q * L
    balance_scale = max(abs(spring_moment / moment_scale), 1.0e-12)
    return {
        "a_frac": float(a_frac), "kappa": "inf" if math.isinf(kappa) else float(kappa),
        "mesh": int(n), "sag_m": ex["sag_m"] * v_scale,
        "x_at_max_m": ex["x_at_max_m"] * L,
        "theta_rad": theta, "spring_moment_nm": float(spring_moment),
        "vertical_reaction_n": vertical_reaction,
        "vertical_equilibrium_rel": abs(vertical_reaction - Q * L) / (Q * L),
        "spring_compatibility_rel": (0.0 if math.isinf(kappa) or spring_moment == 0 else
                                      abs(spring_moment + (kappa * E * I / L) * theta)
                                      / abs(spring_moment)),
        "spring_equilibrium_rel": (0.0 if math.isinf(kappa) else
                                    abs(s["spring_balance_bar"]) / balance_scale),
    }


def _assemble_two_span(n_per_span: int, loaded: str) -> dict:
    n = 2 * n_per_span
    h = 1.0 / n_per_span
    EI = 1.0
    K = np.zeros((2 * (n + 1), 2 * (n + 1)))
    F = np.zeros(2 * (n + 1))
    ke = element_stiffness(EI, h)
    qv = np.zeros(n)
    for e in range(n):
        q = 1.0 if loaded == "all" or e < n_per_span else 0.0
        qv[e] = q
        d = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(d, d)] += ke
        F[d] += element_load(q, h)
    fixed = {0, 2 * n_per_span, 2 * n}
    free = [i for i in range(len(F)) if i not in fixed]
    u = np.zeros_like(F)
    u[free] = solve(K[np.ix_(free, free)], F[free], assume_a="sym")
    # Reuse polynomial extremum search with unit EI and only the first span.
    pts = []
    for e in range(n_per_span):
        c = _element_coefficients(u[2 * e], u[2 * e + 1], u[2 * e + 2],
                                  u[2 * e + 3], qv[e], h, EI=1.0)
        roots = np.roots(np.array([c[1], 2*c[2], 3*c[3], 4*c[4]])[::-1])
        xis = [0.0, 1.0] + [float(z.real) for z in roots
                            if abs(z.imag) < 1e-10 and 0 <= z.real <= 1]
        pts.extend((e*h + xi*h, _polyval(c, xi)) for xi in xis)
    return {"coefficient": float(-min(v for _, v in pts))}


def parent_control() -> dict:
    got = {k: _assemble_two_span(96, k)["coefficient"] for k in ("all", "one")}
    rel = {k: abs(got[k] / PARENT_FROZEN[k] - 1.0) for k in got}
    return {"frozen_038": PARENT_FROZEN, "reproduced": got,
            "worst_relative": max(rel.values())}


def validate_symmetric_load(q_by_element) -> None:
    qv = np.asarray(q_by_element, dtype=float)
    if (not np.allclose(qv, qv[::-1], rtol=0.0, atol=1.0e-12)
            or not np.allclose(qv, 1.0, rtol=0.0, atol=1.0e-12)):
        raise ValueError("the registered load must be uniform and mirror-symmetric")


def run_study(write: bool = True) -> dict:
    meshes = {str(n): [solve_case(a, n=n) for a in POSITIONS] for n in MESHES}
    cases = meshes[str(PRODUCTION_MESH)]
    baseline = force_method(0.5, 0.0)["sag_m"]
    for row in cases:
        row["sag_mm"] = row["sag_m"] * 1000
        row["reduction_percent"] = 100 * (1 - row["sag_m"] / baseline)
        row["theta_mrad"] = row["theta_rad"] * 1000
        row["spring_moment_knm"] = row["spring_moment_nm"] / 1000
    movements = []
    for i in range(len(POSITIONS)):
        s128 = meshes["128"][i]["sag_m"]
        s256 = meshes["256"][i]["sag_m"]
        movements.append(abs(s128 / s256 - 1.0))
    u = max(U_FLOOR, max(movements))
    ordered = sorted(cases, key=lambda r: r["sag_m"])
    gap = abs(ordered[1]["sag_m"] / ordered[0]["sag_m"] - 1.0)
    mirrors = [{"left": solve_case(a), "right": solve_case(1-a)}
               for a in POSITIONS[:-1]]
    analytic = [{"a_frac": a, "fe": solve_case(a), "force": force_method(a, KAPPA)}
                for a in POSITIONS]
    limits = {
        "pinned": [{"a_frac": a, "fe": solve_case(a, 0.0),
                    "force": force_method(a, 0.0)} for a in POSITIONS],
        "rigid": [{"a_frac": a, "fe": solve_case(a, math.inf),
                   "force": force_method(a, math.inf)} for a in POSITIONS],
    }
    out = {
        "study": "ideal rotational restraint to a fixed reference",
        "model": {"E_pa": E, "I_m4": I, "L_m": L, "q_n_per_m": Q,
                  "kappa_theta": KAPPA, "positions_x_over_L": list(POSITIONS),
                  "production_mesh": PRODUCTION_MESH,
                  "deflection_magnification": DEFLECTION_MAGNIFICATION},
        "parent": parent_control(), "baseline_sag_m": baseline,
        "cases": cases, "meshes": meshes, "mesh_movements": movements,
        "u": u, "best_runner_gap": gap,
        "winner_x_over_L": ordered[0]["a_frac"],
        "mirrors": mirrors, "analytic": analytic, "limits": limits,
        "claim": ("Of the four registered points, x/L = 1/8 gives the smallest maximum sag; "
                  "the midpoint restraint does nothing under symmetric loading because the "
                  "beam rotation there is already zero."),
        "scope": ("An ideal dimensionless rotational restraint to a fixed reference, not a "
                  "splice, real connection, strength check or installation recommendation."),
    }
    if write:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "rotational_restraint.json").write_text(
            json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    return out


def main():
    r = run_study(write=True)
    print("IDEAL ROTATIONAL RESTRAINT POSITION")
    for row in r["cases"]:
        print("x/L=%5.3f  sag=%7.4f mm  reduction=%6.2f%%  theta=%+.4f mrad" %
              (row["a_frac"], row["sag_mm"], row["reduction_percent"],
               row["theta_mrad"]))
    print("winner: x/L = %.3f; gap %.3f%%; u %.3e" %
          (r["winner_x_over_L"], 100*r["best_runner_gap"], r["u"]))


if __name__ == "__main__":
    main()
