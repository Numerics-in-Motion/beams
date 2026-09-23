"""Reproduce the ideal rotational-restraint position study."""
from __future__ import annotations
import argparse, hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "solver"))
import rotational_restraint as R
import rotational_restraint_gates as G

REF = os.path.join(HERE, "reference", "reference_canonical.json")

def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()

def rel(a, b):
    return abs(a / b - 1.0) if b else abs(a)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    ref = json.load(open(REF, encoding="utf-8"))
    bad = [m for m, h in ref["source_sha256"].items()
           if sha(os.path.join(HERE, "solver", m + ".py")) != h]
    if bad:
        raise SystemExit("source hash mismatch: " + ", ".join(bad))
    d = R.run_study(write=False)
    if args.quick:
        rows = d["cases"]
        print("x/L and maximum sag (mm):")
        for row in rows:
            print("  %.3f  %.4f" % (row["a_frac"], row["sag_m"] * 1000))
        print("midpoint rotation: %.3e rad" % rows[-1]["theta_rad"])
    else:
        g = G.run(d, write=False)
        for row in g["gates"]:
            print("  %-24s %s" % (row["name"], "PASS" if row["passes"] else "FAIL"))
        if not g["kill_line_survives"]:
            raise SystemExit("a registered gate failed")
    worst = max(rel(a["sag_m"], b["sag_m"])
                for a, b in zip(d["cases"], ref["cases"]))
    if worst > 1e-12 or d["winner_x_over_L"] != ref["winner_x_over_L"]:
        raise SystemExit("the frozen answer moved")
    print("reproduce the frozen reference to %.1e" % worst)
    print(ref["claim"])
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
