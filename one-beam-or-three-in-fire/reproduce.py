"""Reproduce the one / two / three glulam beam fire comparison from frozen sources.

    python reproduce.py            # everything: 54 scenes restrained and unrestrained, sensitivities, checks (~6 min)
    python reproduce.py --quick    # the restrained race, the controls and the outside checks (~20 s)

Compares every number with reference/result.json and exits nonzero, naming the quantity, if anything disagrees or a
source file changed. Hashes are taken over LF-normalised bytes, so a CRLF checkout reproduces too.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solver"))
SOURCES = ("solver/split_glulam.py",)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def flat(d, pre=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{pre}{k}."))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(flat(v, f"{pre}{i}."))
    else:
        out[pre[:-1]] = d
    return out


def same(a, b):
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) or isinstance(b, str) or a is None or b is None:
        return a == b
    if isinstance(a, float) and math.isinf(a) or isinstance(b, float) and math.isinf(b):
        return a == b
    return abs(a - b) <= 1e-7 * max(1.0, abs(b))


def main():
    quick = "--quick" in sys.argv
    ref = json.loads((HERE / "reference" / "result.json").read_text(encoding="utf-8"))
    for rel, want in ref["source_sha256"].items():
        if sha(HERE / rel) != want:
            raise SystemExit(f"source changed: {rel}")
    import split_glulam as sg
    got = json.loads(json.dumps(sg.study(quick=quick)))
    want = {k: v for k, v in ref.items() if k != "source_sha256"}
    if quick:
        want = {k: want[k] for k in got}
        want["claims"] = {k: v for k, v in want["claims"].items() if k in got["claims"]}
    g, w = flat(got), flat(want)
    bad = sorted(k for k in set(g) | set(w) if k not in g or k not in w or not same(g[k], w[k]))
    if bad:
        raise SystemExit("MISMATCH (%d): %s" % (len(bad), ", ".join(bad[:12])))
    r = got["restrained"]
    print("restrained (k_crit = 1), %s:" % got["edition"])
    for m in ("RCS", "RPM"):
        print("  nominal %s: one %.1f, two %.1f, three %.1f min" % ((m,) + tuple(r["nominal_" + m][n] for n in "123")))
    print("  one vs two  : %d/54 scenes, joint adverse corner >= x%.3f" % (r["1v2"]["passes"], r["1v2"]["min_ratio"]))
    print("  one vs three: %d/54 scenes, joint adverse corner >= x%.3f" % (r["1v3"]["passes"], r["1v3"]["min_ratio"]))
    print("  two vs three: %d/54 scenes -> %s" % (r["2v3"]["passes"], got["claims"]["two_vs_three"]))
    if not quick:
        u, bf = got["unrestrained"], got["branch_free_bound"]
        print("unrestrained (EN 1995-1-1 6.3.3 k_crit on the charred section, l_ef %.2f m):" % got["controls"]["l_ef"])
        print("  one vs two %d/54 (>= x%.3f), one vs three %d/54 (>= x%.3f); branch-free bound x%.3f / x%.3f"
              % (u["1v2"]["passes"], u["1v2"]["min_ratio"], u["1v3"]["passes"], u["1v3"]["min_ratio"],
                 bf["one_over_two_min"], bf["one_over_three_min"]))
    print("outside checks: %d/%d" % (sum(b["ok"] for b in got["benchmarks"]), len(got["benchmarks"])))
    if not all(b["ok"] for b in got["benchmarks"]):
        raise SystemExit("an outside check disagrees")
    print("Every number matches the frozen reference.")


if __name__ == "__main__":
    if "--freeze" in sys.argv:
        import split_glulam as sg
        res = json.loads(json.dumps(sg.study()))
        res["source_sha256"] = {rel: sha(HERE / rel) for rel in SOURCES}
        with open(HERE / "reference" / "result.json", "w", encoding="utf-8", newline="\n") as f:
            json.dump(res, f, indent=1)
            f.write("\n")
        print("frozen")
    else:
        main()
