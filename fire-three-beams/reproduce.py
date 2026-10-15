"""Reproduce the three-beam standard-fire comparison from frozen sources.

    python reproduce.py

Re-runs the sizing of all three beams, their heating / charring histories, the input-sensitivity sweep that
defines u, the gap test against max(10 %, 3u) and the five independent checks, then compares every published
number with reference/result.json. Exits nonzero, naming the quantity, if anything disagrees or a source changed.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solver"))


def sha(p):
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()   # line-ending independent


def main():
    ref = json.loads((HERE / "reference" / "result.json").read_text(encoding="utf-8"))
    for rel, want in ref["source_sha256"].items():
        if sha(HERE / rel) != want:
            raise SystemExit(f"source changed: {rel}")
    import fire_beams as fb
    with tempfile.TemporaryDirectory() as tmp:
        o = fb.run(tmp)
    bad = []
    for k, v in ref["t_fail_s"].items():
        if abs(o["t_fail_s"][k] - v) > 1e-6 * max(1.0, v):
            bad.append(f"failure time {k}: {o['t_fail_s'][k]} vs {v}")
    for k, v in ref["u_s"].items():
        if abs(o["u_s"][k] - v) > 1e-6 * max(1.0, v):
            bad.append(f"u {k}")
    if o["claims"] != ref["claims"] or not all(o["claims"].values()):
        bad.append("claims")
    if not o["validation_pass"]:
        bad.append("validation")
    if bad:
        raise SystemExit("MISMATCH: " + "; ".join(bad))
    for k in ("aluminium", "steel", "glulam"):
        print(f"{k:9s} bending resistance < floor load after {o['t_fail_s'][k] / 60:5.1f} min  (u {o['u_s'][k] / 60:.1f} min)")
    for g in o["gaps"]:
        print(f"{g['pair'][0]} -> {g['pair'][1]}: gap {g['gap']:.0%} > max(10%, 3u) = {g['threshold']:.0%} -> holds")
    print("all validation checks pass; reproduced", ", ".join(ref["editions"].values()))


if __name__ == "__main__":
    main()
