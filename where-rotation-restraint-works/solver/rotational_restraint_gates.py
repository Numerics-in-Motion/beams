"""Registered gates for the ideal-rotational-restraint study."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import rotational_restraint as R


OUT = R.OUT
TOL = {
    "parent": 2.0e-5,
    "pinned": 1.0e-9,
    "finite": 1.0e-8,
    "rigid": 1.0e-8,
    "mirror": 1.0e-9,
    "equilibrium": 1.0e-10,
    "compatibility": 1.0e-10,
}


def rel(a, b):
    if a == b == 0:
        return 0.0
    if b == 0:
        return abs(a)
    return abs(a / b - 1.0)


def _force_errors(rows):
    e = []
    for row in rows:
        e.append(rel(row["fe"]["sag_m"], row["force"]["sag_m"]))
        if not math.isinf(float(row["fe"]["kappa"]) if row["fe"]["kappa"] != "inf" else math.inf):
            # FE uses upward-positive displacement; the independent force method
            # uses downward-positive displacement.  Preserve and test that sign.
            fe_down = -row["fe"]["theta_rad"]
            fm_down = row["force"]["theta_rad"]
            e.append(abs(fe_down - fm_down) / max(abs(fm_down), 1.0e-12))
            # The force method's applied couple is the negative of the FE
            # spring reaction moment under the two sign conventions.
            fe_m = row["fe"]["spring_moment_nm"]
            fm_m = row["force"]["couple_nm"]
            e.append(abs(fe_m + fm_m) / max(abs(fm_m), 1.0e-9))
    return e


def run(data=None, write=True):
    d = data or R.run_study(write=write)
    gates = []

    def add(name, passes, **detail):
        gates.append({"name": name, "passes": bool(passes), **detail})

    add("parent", d["parent"]["worst_relative"] <= TOL["parent"],
        worst_relative=d["parent"]["worst_relative"], tolerance=TOL["parent"])

    pinned_errors = _force_errors(d["limits"]["pinned"])
    pinned_sags = [x["fe"]["sag_m"] for x in d["limits"]["pinned"]]
    pinned_closed = 5.0 / 384.0 * R.Q * R.L**4 / (R.E * R.I)
    pinned_errors += [rel(v, pinned_closed) for v in pinned_sags]
    add("pinned_limit", max(pinned_errors) <= TOL["pinned"],
        worst_relative=max(pinned_errors), tolerance=TOL["pinned"])

    finite_errors = _force_errors(d["analytic"])
    add("finite_force_method", max(finite_errors) <= TOL["finite"],
        worst_relative=max(finite_errors), tolerance=TOL["finite"])

    rigid_errors = [rel(x["fe"]["sag_m"], x["force"]["sag_m"])
                    for x in d["limits"]["rigid"]]
    add("rigid_limit", max(rigid_errors) <= TOL["rigid"],
        worst_relative=max(rigid_errors), tolerance=TOL["rigid"])

    mirror_errors = [rel(x["left"]["sag_m"], x["right"]["sag_m"])
                     for x in d["mirrors"]]
    add("mirror", max(mirror_errors) <= TOL["mirror"],
        worst_relative=max(mirror_errors), tolerance=TOL["mirror"])

    eq = max(x["vertical_equilibrium_rel"] for x in d["cases"])
    cp = max(x["spring_compatibility_rel"] for x in d["cases"])
    se = max(x["spring_equilibrium_rel"] for x in d["cases"])
    add("equilibrium", eq <= TOL["equilibrium"],
        worst_relative=eq, tolerance=TOL["equilibrium"])
    add("spring_compatibility", cp <= TOL["compatibility"] and
        se <= TOL["compatibility"], worst_relative=cp,
        independent_equilibrium_relative=se, tolerance=TOL["compatibility"])

    add("resolution", d["best_runner_gap"] >= 5*d["u"], u=d["u"],
        gap=d["best_runner_gap"], headroom=d["best_runner_gap"]/(5*d["u"]))

    midpoint = next(x for x in d["cases"] if x["a_frac"] == 0.5)
    add("midpoint_null",
        abs(midpoint["theta_rad"]) <= 1e-12 and
        rel(midpoint["sag_m"], d["baseline_sag_m"]) <= 1e-9,
        theta_rad=midpoint["theta_rad"],
        sag_relative_to_unrestrained=rel(midpoint["sag_m"], d["baseline_sag_m"]))

    known_bad = {}
    try:
        R.solve_case(1/8, -1.0)
        known_bad["negative_stiffness"] = False
    except ValueError:
        known_bad["negative_stiffness"] = True
    # Faults injected through the numerical result path must miss the frozen
    # force-method result by more than its gate tolerance.
    fm = R.force_method(1/8, R.KAPPA)
    bad_coord = R.solve_case(33/256, R.KAPPA, validate=False)
    known_bad["coordinate_mapping"] = rel(bad_coord["sag_m"], fm["sag_m"]) > TOL["finite"]
    bad_operator = R.solve_case(1/8, 1.02, validate=False)
    known_bad["spring_operator"] = rel(bad_operator["sag_m"], fm["sag_m"]) > TOL["finite"]
    try:
        qv = np.ones(R.PRODUCTION_MESH)
        qv[3] *= 1.01
        R.assemble_single(1/8, R.KAPPA, R.PRODUCTION_MESH,
                          q_by_element=qv, validate=True)
        known_bad["load_symmetry"] = False
    except ValueError:
        known_bad["load_symmetry"] = True
    add("known_bad", all(known_bad.values()), checks=known_bad)

    winner_ok = d["winner_x_over_L"] == 1/8
    reductions = [x["reduction_percent"] for x in d["cases"]]
    payoff_ok = winner_ok and abs(reductions[-1]) < 1e-7 and reductions[0] > 10
    add("registered_payoff", payoff_ok, winner=d["winner_x_over_L"],
        reductions_percent=reductions)

    rows = d["cases"]
    content = {
        "title": ["WHERE SHOULD ONE", "ROTATION RESTRAINT GO?"],
        "setup": ["SAME 6.0 m STEEL BEAM", "SAME 10.0 kN/m FLOOR LOAD",
                  "SAME IDEAL RESTRAINT · κθ = 1", "ONLY x/L CHANGES"],
        "labels": ["x/L = 1/8", "x/L = 1/4", "x/L = 3/8", "x/L = 1/2"],
        "values_mm": ["%.3f mm" % x["sag_mm"] for x in rows],
        "twist": ["MIDPOINT: MAXIMUM BENDING MOMENT",
                  "BUT ROTATION = 0", "SPRING MOMENT = 0"],
        "answer": ["BEST OF THESE FOUR: x/L = 1/8",
                   "MAXIMUM SAG: %.3f mm" % rows[0]["sag_mm"],
                   "%.1f%% BELOW THE UNRESTRAINED SAG" % rows[0]["reduction_percent"]],
        "disclosure": ["IDEAL RESTRAINT TO A FIXED REFERENCE",
                       "NOT A SPLICE OR REAL CONNECTION",
                       "DEFLECTION ×100 · LINEAR ELASTIC MODEL",
                       "NOT AN INSTALLATION RECOMMENDATION"],
        "repo": "REPRODUCTION: github.com/Numerics-in-Motion/beams",
        "illustration": "OPENING CLOSE-UP IS AN ILLUSTRATION",
    }
    narration = {
        "n00": "A six-metre steel beam, carrying the same floor load.",
        "n01": "One ideal rotational restraint can go at four registered positions.",
        "n02": "It connects beam rotation to a fixed reference. It is not a splice.",
        "n03": "Every case uses the same stiffness: kappa theta equals one.",
        "n04": "We measure the greatest downward deflection anywhere on the beam.",
        "n05": "At one eighth of the span, it falls to %.3f millimetres." % rows[0]["sag_mm"],
        "n06": "At one quarter: %.3f. At three eighths: %.3f." %
               (rows[1]["sag_mm"], rows[2]["sag_mm"]),
        "n07": "But at the midpoint, the spring moment is zero.",
        "n08": "The beam bends most there, yet its rotation is already zero.",
        "n09": "So the rotational restraint has nothing to resist.",
        "n10": "Of these four registered positions, one eighth gives the least maximum downward deflection.",
        "n11": "An ideal model, not a real connection or installation recommendation.",
    }
    all_pass = all(g["passes"] for g in gates)
    result = {"passed": sum(g["passes"] for g in gates), "total": len(gates),
              "gates": gates, "tolerances": TOL,
              "all_pass": all_pass,
              "no_go": {g["name"]: g["passes"] for g in gates},
              "kill_line_survives": all_pass,
              "claim": d["claim"], "scope": d["scope"],
              "on_screen": content, "narration": narration}
    if write:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "gates.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main():
    g = run(write=True)
    for row in g["gates"]:
        print("%-24s %s" % (row["name"], "PASS" if row["passes"] else "FAIL"))
    print("%d/%d gates pass" % (g["passed"], g["total"]))
    raise SystemExit(0 if g["kill_line_survives"] else 1)


if __name__ == "__main__":
    main()
