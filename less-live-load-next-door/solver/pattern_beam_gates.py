# -*- coding: utf-8 -*-
"""Less live load next door, more sag here -- the gates.

    python src/pattern_beam_gates.py      # reads data/pattern_beam.json, writes gates.json

Every line is the design review's, registered before the claim was written.

  control_gate      the two-span pin/roller UDL coefficients from beam theory
                    (0.005416 and 0.009151 wL^4/EI) reproduced within 0.5 %
  resolution_gate   u = the largest change of any reported ratio from 48 to
                    96 elements per span; u < 1 %, and the smallest effect
                    (ratio - 1) exceeds 5u
  core_gate         rigid, level, equal spans; 2/3/4 spans x pinned/fixed ends
                    x UDL/point x every target span x all 2^n live-load
                    subsets: the largest target-span sag comes from a subset
                    loading FEWER spans than all, and exceeds the all-spans
                    sag by more than the resolution margin
  alternate_gate    wherever "alternating spans" is named, exhaustive
                    enumeration confirms it is the governing subset
  elastic_gate      the displayed girder stays below yield under every
                    live-load subset
  flexible_gate     the on-screen caveat matches the stiffness sweep: if a
                    flexible interior support reverses the effect, the words
                    say "shrink or reverse"
  wording_gate      nothing shipped says less load on this span, always, never,
                    all continuous beams, safer, stronger, failure, collapse or
                    code-required, and no interior-span ratio is quoted
"""
from __future__ import annotations

import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import pattern_beam as PB             # noqa: E402

CONTROL_LINE = 0.005
RESOLUTION_LINE = 0.01
EFFECT_MULTIPLE = 5.0

TITLE = "LESS LIVE LOAD NEXT DOOR, MORE SAG HERE"
CLAIM = None                      # written by run() from the measured range
ON_SCREEN = dict(
    setup="A STEEL FLOOR GIRDER, CONTINUOUS OVER FOUR COLUMNS",
    same="SAME LIVE LOAD ON THIS BAY",
    removed="NEIGHBOURING LIVE LOAD REMOVED",
    model="RIGID, LEVEL, EQUAL-SPAN MODEL",
    headline_template="GOVERNING END-SPAN SAG: {lo:.2f}–{hi:.2f}× — LIVE LOAD ONLY",
    why="THE LOAD NEXT DOOR BENDS THIS BAY UP BESIDE THE SHARED COLUMN",
    why_sub="NEXT-BAY LOAD, ISOLATED",
    rise_template="PEAK RISE HERE: {mm:.1f} mm",
    permanent="Permanent load reduces the ratio.",
    flexible="Flexible supports can shrink or reverse the effect.",
    linear="LINEAR ELASTIC 2-D BEAM MODEL",
    magnified_template="VERTICAL DEFLECTION ×{factor:d}",
)
SILENT = dict(
    result_template="end bay {lo:.2f}–{hi:.2f}× across the cases",
    limits=["live load only; permanent load reduces it",
            "rigid, level, equal spans",
            "flexible supports can reverse it"],
)
NARRATION = dict(
    a1="A steel floor girder, continuous over four columns.",
    a2="The same live load on every bay.",
    b1="Now take the live load off the middle bay only.",
    c1="The end bay sags more, with nothing added to it.",
    d1="The load next door was lifting this bay.",
    e1="Alternating bays governed every rigid, level, equal-span end-bay case.",
    e2="One point three nine to one point seven six times.",
    z1="Permanent load shrinks it.",
    z2="Flexible supports can reverse it.",
)
IDENTITY = ["WORKER DRAWN TO SCALE, 1.70 m - AN ILLUSTRATION"]
MAGNIFICATION = 50
FORBIDDEN = [r"less load on this (span|bay)", r"\balways\b", r"\bnever\b",
             r"all continuous beams", r"\bsafer\b", r"\bstronger\b",
             r"\bfail", r"\bcollaps", r"code[- ]required", r"\bbounce",
             r"\b13\b|13\.\d|\b13×", r"times more"]


def load(path=None):
    return json.load(io.open(path or os.path.join(PB.OUT, "pattern_beam.json"),
                             encoding="utf-8"))


def _key(r):
    return (r["n_spans"], r["ends"], r["load"], r["target"])


def control_gate(doc):
    cf = doc["control"]["closed_form"]
    worst = 0.0
    for per, v in doc["control"]["fe"].items():
        worst = max(worst, abs(v["all"] / cf["all"] - 1), abs(v["one"] / cf["one"] - 1))
    for per, v in doc["control"]["moment"].items():
        worst = max(worst, abs(v["support_over_wL2_8"] - 1),
                    abs(v["sagging_over_9wL2_128"] - 1))
    return dict(name="control", closed_form=cf,
                fe_production=doc["control"]["fe"][str(PB.PRODUCTION_MESH)],
                worst=worst, passes=bool(worst <= CONTROL_LINE))


def _by(doc):
    out = {}
    for r in doc["core"]:
        out.setdefault(_key(r), {})[r["mesh"]] = r
    return out


def resolution_gate(doc):
    by = _by(doc)
    u = max(abs(v[PB.PRODUCTION_MESH]["ratio"] / v[PB.REFINED_MESH]["ratio"] - 1)
            for v in by.values())
    smallest = min(v[PB.REFINED_MESH]["ratio"] - 1 for v in by.values())
    return dict(name="resolution", u=u, smallest_effect=smallest,
                passes=bool(u < RESOLUTION_LINE and smallest > EFFECT_MULTIPLE * u))


def core_gate(doc, u):
    by = _by(doc)
    bad = []
    for k, v in sorted(by.items()):
        r = v[PB.REFINED_MESH]
        if not (len(r["governing"]) < k[0] and r["ratio"] - 1 > EFFECT_MULTIPLE * u):
            bad.append(dict(case=k, governing=r["governing"], ratio=r["ratio"]))
    ends = [v[PB.REFINED_MESH]["ratio"] for v in by.values()
            if v[PB.REFINED_MESH]["end_span"]]
    ends_udl = [v[PB.REFINED_MESH]["ratio"] for k, v in by.items()
                if v[PB.REFINED_MESH]["end_span"] and k[2] == "udl"]
    return dict(name="core", cases=len(by), subsets=sum(2 ** k[0] for k in by),
                bad=bad, end_span_range=[min(ends), max(ends)],
                end_span_range_udl=[min(ends_udl), max(ends_udl)],
                passes=not bad)


def alternate_gate(doc):
    by = _by(doc)
    bad = [dict(case=k, governing=v[PB.REFINED_MESH]["governing"],
                alternate=v[PB.REFINED_MESH]["alternate"])
           for k, v in by.items()
           if v[PB.REFINED_MESH]["governing"] != v[PB.REFINED_MESH]["alternate"]]
    d = PB.DISPLAY
    shown = by[(d["n_spans"], d["ends"], d["load"], d["target"])][PB.REFINED_MESH]
    return dict(name="alternate", not_alternate=bad,
                display_governing=shown["governing"],
                display_ratio=shown["ratio"], passes=not bad)


def elastic_gate(doc):
    s = doc["display"]["worst_stress_over_yield"]
    return dict(name="elastic", worst_stress_over_yield=s, passes=bool(s < 0.8))


def flexible_gate(doc):
    rows = doc["support_stiffness"]
    reverses = any(r["alone_ratio"] < 1.0 or r["ratio"] <= 1.0 + 1e-9
                   for r in rows if r["kappa_v"] != "inf")
    says = "reverse" in ON_SCREEN["flexible"].lower()
    return dict(name="flexible", reversal_found=reverses,
                caveat=ON_SCREEN["flexible"], passes=bool(reverses == says))


def strings(rng):
    lo, hi = rng
    s = [TITLE] + IDENTITY + list(NARRATION.values())
    for k, v in ON_SCREEN.items():
        if k == "headline_template":
            s.append(v.format(lo=lo, hi=hi))
        elif k == "magnified_template":
            s.append(v.format(factor=MAGNIFICATION))
        elif k == "rise_template":
            s.append(v.format(mm=4.1))
        else:
            s.append(v)
    s.append(SILENT["result_template"].format(lo=lo, hi=hi))
    s += SILENT["limits"]
    return s


def hits(text):
    return [(m.group(0), text[max(0, m.start() - 30):m.end() + 30])
            for pat in FORBIDDEN for m in re.finditer(pat, text, re.I)]


def wording_gate(texts):
    found = [h for t in texts for h in hits(t)]
    return dict(name="wording", hits=found, passes=not found)


def run(write=True):
    doc = load()
    res = resolution_gate(doc)
    core = core_gate(doc, res["u"])
    rng = [round(core["end_span_range"][0], 2), round(core["end_span_range"][1], 2)]
    texts = strings(rng)
    gates = [control_gate(doc), res, core, alternate_gate(doc),
             elastic_gate(doc), flexible_gate(doc), wording_gate(texts)]
    on = dict(ON_SCREEN)
    on["headline"] = ON_SCREEN["headline_template"].format(lo=rng[0], hi=rng[1])
    on["magnified"] = ON_SCREEN["magnified_template"].format(factor=MAGNIFICATION)
    silent = dict(SILENT)
    silent["result"] = SILENT["result_template"].format(lo=rng[0], hi=rng[1])
    claim = ("In the registered rigid, level, equal-span linear-elastic "
             "continuous-girder models, the governing live-load pattern loaded "
             "the target end bay and every second bay from it. The largest "
             "end-bay sag was %.2f to %.2f times the all-bays-loaded sag; "
             "permanent load reduces the ratio, and flexible supports can "
             "shrink or reverse the effect." % (rng[0], rng[1]))
    out = dict(title=TITLE, claim=claim, on_screen=on, silent=silent,
               narration=NARRATION, identity=IDENTITY,
               magnification=MAGNIFICATION, display=doc["display"],
               end_span_range=rng, gates=gates,
               passed=sum(bool(g["passes"]) for g in gates), total=len(gates))
    if write:
        with io.open(os.path.join(PB.OUT, "gates.json"), "w", encoding="utf-8",
                     newline="\n") as fh:
            json.dump(out, fh, indent=1, ensure_ascii=False, default=str)
    return out


def main():
    out = run()
    for g in out["gates"]:
        extra = {k: v for k, v in g.items() if k not in ("name", "passes")}
        print("%-11s %s  %s" % (g["name"], "PASS" if g["passes"] else "FAIL",
                                json.dumps(extra, default=str)[:230]))
    print("%d / %d" % (out["passed"], out["total"]))
    print("claim:", out["claim"])
    return 0 if out["passed"] == out["total"] else 1


if __name__ == "__main__":
    sys.exit(main())
