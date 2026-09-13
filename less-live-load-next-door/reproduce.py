"""Less live load next door, more sag here.

    python reproduce.py            # the whole study and every gate (seconds)
    python reproduce.py --quick    # the displayed girder and the control only

WHAT IS CLAIMED

In the registered rigid, level, equal-span linear-elastic continuous-girder
models, the governing live-load pattern loaded the target end bay and every
second bay from it. The largest end-bay sag was 1.39 to 1.76 times the
all-bays-loaded sag; permanent load reduces the ratio, and flexible supports
can shrink or reverse the effect.

WHAT IS NOT

Deflection, not strength. Quasi-static load cases, not a dynamic event. No slab,
composite action, connection or column flexibility. Interior spans are not
quoted.

WHAT THIS DOES

Checks every module and the config against the SHA-256 frozen with the
reference, re-runs the study (every span count, end condition, load type,
target span and all 2^n live-load subsets, at four meshes, plus every
sensitivity), re-applies every gate, and fails if a gated number has moved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "solver"))
sys.path.insert(0, RIG)
RUN = tempfile.mkdtemp(prefix="pattern-")
os.environ["PATTERN_DATA"] = RUN

import pattern_beam as PB                 # noqa: E402
import pattern_beam_gates as G            # noqa: E402

REF = os.path.join(HERE, "reference", "reference_canonical.json")
TOL = 1e-9


def sha(path):
    """SHA-256 of the file with line endings normalised to LF, so the check
    means the same on a Linux clone, a Windows checkout that converts to
    CRLF, and the tree it was frozen from."""
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


def _where(k):
    p = os.path.join(HERE, "solver", k + ".py")
    return p if os.path.exists(p) else os.path.join(RIG, k + ".py")


def check_sources(ref):
    bad = [k for k, want in ref["source_sha256"].items()
           if sha(_where(k)) != want]
    if sha(os.path.join(RIG, "config", "beams_config.json")) != ref["config_sha256"]:
        bad.append("config/beams_config.json")
    if bad:
        raise SystemExit("these do not match the reference: %s" % ", ".join(bad))
    print("   %d modules and the config match the reference by SHA-256"
          % len(ref["source_sha256"]))


def moved(a, b):
    return abs(a / b - 1.0)


def quick(ref):
    cfg = PB.config()
    c = PB.control(cfg)
    d = PB.display_case(cfg)
    print("   two spans, one loaded / all loaded: %.6f / %.6f wL^4/EI"
          % (c["fe"][PB.PRODUCTION_MESH]["one"], c["fe"][PB.PRODUCTION_MESH]["all"]))
    print("   displayed girder: end bay %.2f mm with every bay loaded, %.2f mm "
          "with the middle bay unloaded (x%.3f)"
          % (d["sag_all_mm"], d["sag_alternate_mm"],
             d["sag_alternate_mm"] / d["sag_all_mm"]))
    worst = max(moved(d["sag_all_mm"], ref["display"]["sag_all_mm"]),
                moved(d["sag_alternate_mm"], ref["display"]["sag_alternate_mm"]))
    print("   reproduce the reference to %.1e" % worst)
    return 0 if worst <= TOL else 1


def full(ref):
    PB.main()
    res = G.run(write=True)
    print()
    for g in res["gates"]:
        print("   %-11s %s" % (g["name"], "PASS" if g["passes"] else "FAIL"))
    gates = {g["name"]: g for g in res["gates"]}
    pairs = [(gates["core"]["end_span_range"][0], ref["end_span_range_exact"][0]),
             (gates["core"]["end_span_range"][1], ref["end_span_range_exact"][1]),
             (gates["alternate"]["display_ratio"], ref["display_ratio"])]
    worst = max(moved(a, b) for a, b in pairs)
    print()
    print("   end-span sag ratio %.3f to %.3f over %d cases and %d subsets"
          % (gates["core"]["end_span_range"][0], gates["core"]["end_span_range"][1],
             gates["core"]["cases"], gates["core"]["subsets"]))
    print("   the gated numbers reproduce the reference to %.1e" % worst)
    if worst > TOL:
        raise SystemExit("the answer moved against the frozen reference")
    print()
    print("   claim: %s" % res["claim"])
    return 0 if res["passed"] == res["total"] else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    print("LESS LIVE LOAD NEXT DOOR, MORE SAG HERE")
    ref = json.load(open(REF, encoding="utf-8"))
    check_sources(ref)
    return quick(ref) if a.quick else full(ref)


if __name__ == "__main__":
    raise SystemExit(main())
