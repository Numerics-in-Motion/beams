"""One support drops, nothing is on the beam, and two of the four are hurt.

    python reproduce.py           # every gate            (a few seconds)
    python reproduce.py --quick   # the four supports only

WHAT THIS ASKS

The published beam study put one 3 m beam on four different supports and
loaded it. This puts NO load on at all and moves the left-hand support down instead.

    cantilever        never yields   the beam goes down with it
    pinned + roller   never yields   it tilts; nothing is stretched
    fixed + roller     83.3 mm
    fixed + fixed      41.7 mm

The two that can move freely develop nothing. The two held at both ends
cannot, so compatibility turns the movement into stress.

"NOTHING IS APPLIED" IS NOT THE CLAIM

A prescribed displacement does work and produces reactions; it is a boundary
condition of the other kind, not the absence of one. What is true is that no
LOAD is placed on the beam.

THE PART THAT IS NOT OBVIOUS

    d_yield = sigma_y L^2 / (k E c)

**The second moment of area is not in it.** Four sections of the same steel
area, whose I spreads 194 %, all give the same d_yield times depth to 1 part
in 1e12 -- and the hollow tube and the box, 36 % apart in stiffness, give the
identical 22.1 mm because they are the same depth. The LEAST stiff section is
the most tolerant. A stiffer beam buys nothing against the ground moving.

What buys tolerance is span squared over depth, and it runs against
intuition: a 3 m beam 300 mm deep yields at 12.5 mm while a 12 m beam 90 mm
deep takes 667 mm. **53x, and the floppy one wins.**

THE CLAIM IS NARROW ON PURPOSE

More redundancy does NOT order the answer in general. `scope_gate` holds the
counterexample: a 1-degree-redundant 2 m beam yields at 37 mm while a
3-degree-redundant 12 m beam takes 667. The result here is stated for ONE
span and ONE section with the end restraint as the only variable.

YIELD IS NOT THE FIRST THING THAT GOES WRONG

41.7 mm is L/72. A serviceability reference is nearer L/500, so cracking,
tilt and doors that stop closing arrive about seven times sooner. A number on
screen without that context reads as a safe allowance, and it is not one.

VERIFICATION, NOT VALIDATION

Every check reads the same E, section and geometry the solver reads, and the
closed form it reproduces is the textbook one for the same idealisation.
**Validating this would need a real beam and a real settlement. There is
neither.**

WHAT IS NOT CLAIMED

Linear elastic, small displacement, quasi-static, ideal supports. A real
fixed end rotates a little and a real bearing has friction and gaps, and
every one of those RELIEVES the secondary moment -- so these are upper
bounds, not predictions. First yield, not collapse. One span, one settling
support, no creep, no relaxation, no cracking, no soil.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "solver"))
sys.path.insert(0, RIG)

import numpy as np                      # noqa: E402
import settlement as S                  # noqa: E402
import settlement_gates as G            # noqa: E402

REF = os.path.join(HERE, "reference", "reference_canonical.json")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _where(k):
    """This study's modules sit beside it; the rig's own sit at the root."""
    p = os.path.join(HERE, "solver", k + ".py")
    return p if os.path.exists(p) else os.path.join(RIG, k + ".py")


def check_sources(ref):
    bad = []
    for k, want in ref["source_sha256"].items():
        got = sha(_where(k))
        if got != want:
            bad.append("%s differs: %s != %s" % (k, got[:12], want[:12]))
    if bad:
        raise SystemExit("the modules do not match the reference:\n  "
                         + "\n  ".join(bad))
    print("   %d modules match the reference by SHA-256 -- the beam rig's own "
          "solver, not a retyped one" % len(ref["source_sha256"]))


def quick():
    cfg = S.config()
    sec = S.section(cfg)
    span = float(cfg["geometry"]["span_m"])
    print("   %-16s %18s %12s" % ("support", "settles at yield", "as span /"))
    for sup in S.SUPPORTS:
        d = S.settlement_at_first_yield(sup, sec, cfg)
        print("   %-16s %18s %12s"
              % (sup, "never" if np.isinf(d) else "%.1f mm" % (1000 * d),
                 "-" if np.isinf(d) else "L/%.0f" % (span / d)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="the four supports only, no gates")
    a = ap.parse_args(argv)

    print("ONE SUPPORT DROPS -- nothing is on the beam")
    ref = json.load(open(REF, encoding="utf-8"))
    check_sources(ref)
    print()

    if a.quick:
        quick()
        return 0

    res = G.run(write=False)
    m = res["meta"]

    print()
    print("   NO-GO gates: %s"
          % ", ".join("%s %s" % (k, "ok" if v else "FAIL")
                      for k, v in res["no_go"].items()))

    # the frozen reference is the check, not the printout
    want = ref["meta"]["d_yield_m"]
    worst = 0.0
    for k, v in want.items():
        got = m["d_yield_m"][k]
        if v is None or got is None:
            assert v is None and got is None, k
            continue
        worst = max(worst, abs(got - v) / v)
    print("   reproduces the frozen reference to %.2e" % worst)
    if worst > 1e-12:
        raise SystemExit("the answer moved against the frozen reference")

    print()
    quick()
    print()
    print("   the section does not decide it: I spreads %.0f %% and "
          "d_yield x depth is one number to %.0e"
          % (100 * res["gates"]["section_gate"]["spread_of_I"],
             res["gates"]["section_gate"]["spread_of_d_times_c"]))
    print("   span^2 / depth does: %.1f mm at 3 m deep 0.30 m against "
          "%.0f mm at 12 m deep 0.09 m -- %.0fx"
          % (1000 * min(m["fragility_sweep_m"].values()),
             1000 * max(m["fragility_sweep_m"].values()),
             m["fragility_ratio"]))
    print()
    print("   NOT a prediction: ideal supports, first yield, one span. "
          "Yield is L/%.0f, about %.0fx beyond an L/%.0f serviceability "
          "reference."
          % (res["gates"]["serviceability_gate"]["as_span_over"],
             res["gates"]["serviceability_gate"]["times_beyond_service"],
             res["gates"]["serviceability_gate"]
             ["service_reference_span_over"]))
    return 0 if res["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
