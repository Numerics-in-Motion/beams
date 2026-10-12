"""Reproduce Spinoff 065: one ideal lateral brace at three heights on an IPE 300, load on top or hanging.

    python reproduce.py

Re-runs both elastic methods (Ritz series and thin-walled beam FE) on the registered grid (spans 4/6/8 m, two
load heights, four brace positions), the validation checks and the EN 1993-1-1 design resistances, and compares
every published number with reference/result.json. Exits nonzero, naming the quantity, on any mismatch or if a
source changed.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solver"))


def sha(p):
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    ref = json.loads((HERE / "reference" / "result.json").read_text(encoding="utf-8"))
    for rel, want in ref["source_sha256"].items():
        if sha(HERE / rel) != want:
            raise SystemExit(f"source changed: {rel}")
    import ltb065 as M
    g = M.run_grid()
    v = M.run()
    bad = []
    for key, row in ref["grid"].items():
        for b, vals in row.items():
            for q in ("P_cr_fe", "P_Rd"):
                got = g["grid"][key][b][q]
                if abs(got - vals[q]) > 1e-6 * vals[q]:
                    bad.append(f"{key} {b} {q}: {got} vs {vals[q]}")
    if g["claims"] != ref["claims"] or not all(g["claims"].values()):
        bad.append("claims")
    if not v["validation_pass"]:
        bad.append("validation")
    if bad:
        raise SystemExit("MISMATCH: " + "; ".join(bad[:5]))
    for load in ("on_top", "hanging"):
        r = g["grid"][f"6|{load}"]
        print(f"6 m, load {load}:", ", ".join(f"{b} Pcr {r[b]['P_cr_fe'] / 1e3:.1f} kN / PRd {r[b]['P_Rd'] / 1e3:.1f} kN"
                                             for b in ("none", "bottom", "centre", "top")))
    print("claims:", ", ".join(k for k, ok in g["claims"].items() if ok), "-> all hold")
    print(f"Ritz vs FE max relative difference {g['max_ritz_fe_rel']:.2e}; validation pass")


if __name__ == "__main__":
    main()
