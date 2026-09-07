# -*- coding: utf-8 -*-
"""Gates for "one support drops". Seven, five of them NO-GO.

    python src/settlement_gates.py --out <dir>

A NO-GO gate is one whose failure means the video cannot be made at all, not
one whose wording would change.

  replication_gate   NO-GO  the prescribed-displacement solve must reproduce
                            `M = k EI d / L^2` for all four supports.
  determinacy_gate   NO-GO  the two free cases must come out zero IN THE MODEL
                            ACTUALLY SOLVED, not by counting reactions. The
                            design review was explicit: the 027 truss counts as
                            pin-jointed determinate while its elements carry
                            moment, so determinacy has to be measured.
  section_gate       NO-GO  `d_yield` must be independent of I. This is the
                            claim: a stiffer section buys no tolerance. Four
                            sections of the SAME area, and `d_yield * c` must
                            be the same number for all of them.
  scope_gate         NO-GO  a counterexample must exist to "more redundancy is
                            worse", or the claim would be stated too widely.
                            The study may only say it for a fixed span and
                            section with the end restraint as the variable.
  mesh_gate                 halving the element size must not move `d_yield`.
  slenderness_gate          every reported case must have L/depth >= 10, where
                            a Euler-Bernoulli beam is a defensible model.
  serviceability_gate NO-GO the yield settlement must be reported against a
                            serviceability reference, so the video cannot leave
                            "42 mm is fine" on screen.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import column_geometry as CG                              # noqa: E402
import settlement as S                                    # noqa: E402

# Where the result is written is the CALLER's business. A module that ships
# inside a public package must not carry this repository's own layout.

REPLICATION_TOL = 1e-9      # the closed form is exact; this is round-off only
MESH_TOL = 1e-6
SLENDERNESS_MIN = 10.0      # below this, shear deformation is not negligible
SERVICE_REF = 500.0         # a common serviceability limit is around L/500


def echo(s=""):
    print(s)


def gate(label, ok, detail=""):
    echo("   [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                           ("   %s" % detail) if detail else ""))
    return bool(ok)


def run(write=False, out_dir=None):
    t0 = time.time()
    cfg = S.config()
    sec = S.section(cfg)
    span = float(cfg["geometry"]["span_m"])
    E = float(cfg["material"]["youngs_modulus_pa"])
    sig = float(cfg["material"]["failure_stress_pa"])
    out = {"gates": {}}

    echo("ONE SUPPORT DROPS -- nothing is on the beam")
    echo("   %s, span %.2f m, section %s (A %.4f m2, I %.3e m4, c %.4f m)"
         % (cfg["material"]["name"], span, sec.name, sec.area, sec.inertia,
            sec.fibre_c))
    echo()

    # ---- gate 1 (NO-GO): the solve reproduces the closed form -------------
    rep = {}
    worst = 0.0
    m_fixed = S.moment_closed_form("fixed", sec, S.D_PROBE, cfg)
    for sup in S.SUPPORTS:
        m, _u, _s = S.solve_settlement(sup, sec, S.D_PROBE, cfg)
        c = S.moment_closed_form(sup, sec, S.D_PROBE, cfg)
        err = abs(m - c) / c if c else m / m_fixed
        rep[sup] = {"solver_nm": float(m), "closed_nm": float(c),
                    "rel_err": float(err)}
        worst = max(worst, err)
    out["gates"]["replication_gate"] = {
        "demands": "NO-GO. The prescribed-displacement solve must reproduce "
                   "M = k EI d / L^2 for every support, including the two "
                   "that must be exactly zero",
        "per_support": rep, "worst_rel_err": float(worst),
        "tol": REPLICATION_TOL, "pass": bool(worst < REPLICATION_TOL)}
    gate("replication gate", out["gates"]["replication_gate"]["pass"],
         "worst %.2e" % worst)

    # ---- gate 2 (NO-GO): determinacy measured, not counted ----------------
    free_cases = [s for s in S.SUPPORTS if S.K_CLOSED[s] == 0.0]
    held_cases = [s for s in S.SUPPORTS if S.K_CLOSED[s] != 0.0]
    ratios = {s: rep[s]["solver_nm"] / m_fixed for s in free_cases}
    out["gates"]["determinacy_gate"] = {
        "demands": "NO-GO. The free cases must develop no moment in the model "
                   "SOLVED, not in a reaction count. Counting can call a frame "
                   "determinate while its elements carry moment",
        "free": free_cases, "held": held_cases,
        "moment_over_fixed": {k: float(v) for k, v in ratios.items()},
        "pass": bool(all(v < REPLICATION_TOL for v in ratios.values()))}
    gate("determinacy gate", out["gates"]["determinacy_gate"]["pass"],
         "free cases at %.1e of the fixed case" % max(ratios.values()))

    # ---- gate 3 (NO-GO): the second moment of area cancels ----------------
    per_sec = {}
    for s in CG.build_all_columns(cfg):
        d = S.settlement_at_first_yield("fixed", s, cfg)
        per_sec[s.name] = {"area_m2": float(s.area), "I_m4": float(s.inertia),
                           "c_m": float(s.fibre_c), "d_yield_m": float(d),
                           "d_times_c": float(d * s.fibre_c)}
    prod = [v["d_times_c"] for v in per_sec.values()]
    spread_dc = (max(prod) - min(prod)) / np.mean(prod)
    ivals = [v["I_m4"] for v in per_sec.values()]
    spread_I = (max(ivals) - min(ivals)) / np.mean(ivals)
    out["gates"]["section_gate"] = {
        "demands": "NO-GO. d_yield must not depend on I. Four sections of the "
                   "same area: d_yield * c must be one number, while I itself "
                   "varies. A stiffer section buys no tolerance to settlement",
        "per_section": per_sec,
        "spread_of_d_times_c": float(spread_dc),
        "spread_of_I": float(spread_I),
        "pass": bool(spread_dc < 1e-9 and spread_I > 0.10)}
    gate("section gate", out["gates"]["section_gate"]["pass"],
         "I spreads %.0f%%, d*c spreads %.1e" % (100 * spread_I, spread_dc))

    # ---- gate 4 (NO-GO): the claim has a boundary, and it is shown --------
    # A LESS redundant short beam against a MORE redundant long one. If the
    # less redundant one yields first, "more redundancy is worse" is false as
    # a general statement and the video may only say it for a fixed span.
    short_span, long_span = 2.0, 12.0
    d_less = S.settlement_at_first_yield("propped", sec, cfg, span=short_span)
    d_more = S.settlement_at_first_yield("fixed", sec, cfg, span=long_span)
    out["gates"]["scope_gate"] = {
        "demands": "NO-GO. A counterexample to 'more redundancy is worse' must "
                   "exist, so the claim is stated for a fixed span and section "
                   "with the end restraint as the only variable",
        "less_redundant": {"support": "propped", "degree": 1,
                           "span_m": short_span, "d_yield_m": float(d_less)},
        "more_redundant": {"support": "fixed", "degree": 3,
                           "span_m": long_span, "d_yield_m": float(d_more)},
        "counterexample": bool(d_less < d_more),
        "pass": bool(d_less < d_more)}
    gate("scope gate", out["gates"]["scope_gate"]["pass"],
         "1 deg at %.1f m yields at %.0f mm, 3 deg at %.1f m at %.0f mm"
         % (short_span, 1000 * d_less, long_span, 1000 * d_more))

    # ---- gate 5: the mesh ------------------------------------------------
    mesh = {}
    for n in (12, 24, 48):
        mesh["n=%d" % n] = float(S.settlement_at_first_yield(
            "fixed", sec, cfg, n=n))
    vals = list(mesh.values())
    drift = (max(vals) - min(vals)) / np.mean(vals)
    out["gates"]["mesh_gate"] = {
        "demands": "d_yield must not move when the element size is halved",
        "per_mesh": mesh, "drift": float(drift), "tol": MESH_TOL,
        "pass": bool(drift < MESH_TOL)}
    gate("mesh gate", out["gates"]["mesh_gate"]["pass"], "drift %.2e" % drift)

    # ---- gate 6: the beam theory is applicable ----------------------------
    depth = 2.0 * sec.fibre_c
    cases = {"study": span / depth,
             "scope_short": short_span / depth,
             "scope_long": long_span / depth,
             # the fragility sweep is reported too, and it changes the DEPTH.
             # Checking only the study section left its deepest case unchecked.
             "fragility_deep": 3.0 / 0.30,
             "fragility_slender": 12.0 / 0.09}
    out["gates"]["slenderness_gate"] = {
        "demands": "every reported case must have span/depth >= %.0f, where "
                   "shear deformation is small enough for a Euler-Bernoulli "
                   "beam to be a defensible model" % SLENDERNESS_MIN,
        "depth_m": float(depth), "span_over_depth": cases,
        "min": float(min(cases.values())),
        "pass": bool(min(cases.values()) >= SLENDERNESS_MIN)}
    gate("slenderness gate", out["gates"]["slenderness_gate"]["pass"],
         "min L/d = %.1f" % min(cases.values()))

    # ---- gate 7 (NO-GO): yield is not the first thing that goes wrong -----
    d_fixed = S.settlement_at_first_yield("fixed", sec, cfg)
    as_ratio = span / d_fixed
    out["gates"]["serviceability_gate"] = {
        "demands": "NO-GO. d_yield must be reported against a serviceability "
                   "reference. Yield is far from the first consequence of a "
                   "support moving, and a number on screen without that "
                   "context reads as a safe allowance",
        "d_yield_m": float(d_fixed), "as_span_over": float(as_ratio),
        "service_reference_span_over": SERVICE_REF,
        "times_beyond_service": float(SERVICE_REF / as_ratio),
        "pass": bool(as_ratio < SERVICE_REF)}
    gate("serviceability gate", out["gates"]["serviceability_gate"]["pass"],
         "yields at L/%.0f, service reference L/%.0f -- %.0fx beyond"
         % (as_ratio, SERVICE_REF, SERVICE_REF / as_ratio))

    # ---- the numbers the video is allowed to quote ------------------------
    yields = {s: S.settlement_at_first_yield(s, sec, cfg) for s in S.SUPPORTS}
    fragile = {}
    for L, dep in ((3.0, 0.30), (12.0, 0.09)):
        # a section of that depth, same area, so only L^2/c moves
        cc = dep / 2.0
        fragile["L%.0f_d%.2f" % (L, dep)] = float(
            sig * L ** 2 / (S.K_CLOSED["fixed"] * E * cc))
    out["meta"] = {
        "span_m": span, "section": sec.name, "area_m2": float(sec.area),
        "I_m4": float(sec.inertia), "fibre_c_m": float(sec.fibre_c),
        "E_pa": E, "yield_pa": sig,
        "d_yield_m": {k: (None if np.isinf(v) else float(v))
                      for k, v in yields.items()},
        "d_yield_as_span_over": {
            k: (None if np.isinf(v) else float(span / v))
            for k, v in yields.items()},
        "fragility_sweep_m": fragile,
        "fragility_ratio": float(max(fragile.values()) / min(fragile.values())),
        "quantity": "The settlement at first yield: the downward movement of "
                    "ONE support at which the extreme fibre first reaches the "
                    "yield stress, with no load on the beam. Linear elastic, "
                    "small displacement, quasi-static, ideal supports -- every "
                    "real imperfection RELIEVES the secondary moment, so these "
                    "are upper bounds, not predictions.",
        "seconds": float(time.time() - t0)}
    out["all_pass"] = all(v["pass"] for v in out["gates"].values())
    out["no_go"] = {k: out["gates"][k]["pass"]
                    for k in ("replication_gate", "determinacy_gate",
                              "section_gate", "scope_gate",
                              "serviceability_gate")}

    echo()
    echo("   NO-GO gates: %s"
         % ", ".join("%s %s" % (k, "ok" if v else "FAIL")
                     for k, v in out["no_go"].items()))
    echo()
    echo("   %-14s %18s %14s" % ("support", "settles at yield", "as span /"))
    for s in S.SUPPORTS:
        v = yields[s]
        echo("   %-14s %18s %14s"
             % (s, "never" if np.isinf(v) else "%.1f mm" % (1000 * v),
                "-" if np.isinf(v) else "L/%.0f" % (span / v)))
    echo()
    echo("   the section does not matter: I spreads %.0f%% and d*c is one "
         "number to %.0e" % (100 * spread_I, spread_dc))
    echo("   what matters is L^2 / depth: %.1f mm at 3 m deep 0.30 m, "
         "%.0f mm at 12 m deep 0.09 m -- %.0fx"
         % (1000 * min(fragile.values()), 1000 * max(fragile.values()),
            out["meta"]["fragility_ratio"]))
    echo("   %.0f s" % out["meta"]["seconds"])

    if write:
        # `out` is the result dict -- the destination needs its own name
        dest = out_dir or os.path.join(os.getcwd(), "data")
        os.makedirs(dest, exist_ok=True)
        p = os.path.join(dest, "gates.json")
        with open(p, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        echo()
        echo("   wrote %s" % p)
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None,
                    help="directory to write gates.json into")
    a = ap.parse_args(argv)
    res = run(write=True, out_dir=a.out)
    return 0 if res["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
