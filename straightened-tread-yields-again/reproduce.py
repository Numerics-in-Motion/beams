"""Standalone, isolated reproduction of the 051 idealized tread history."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REF=json.loads((ROOT/"reference"/"reference_canonical.json").read_text(encoding="utf-8"))


def verify_sources():
    expected=json.loads((ROOT/"reference"/"source_sha256.json").read_text(encoding="utf-8"))
    for rel,want in expected.items():
        source=ROOT/rel
        if not source.exists() or hashlib.sha256(source.read_bytes()).hexdigest()!=want:
            raise RuntimeError(f"source changed: {rel}")


def close(actual,expected,tol,label):
    if abs(actual-expected)>tol:
        raise RuntimeError(f"{label}: {actual} vs frozen {expected}, tolerance {tol}")


def quick():
    sys.path.insert(0,str(ROOT/"solver"))
    from tread_history051 import Spec,one_case
    spec=Spec()
    q0=one_case(spec,0.,diagnostic=True)
    qstar=one_case(spec,REF["q_star"],diagnostic=True)
    close(q0["bow_mm"],REF["control_bow_mm"],1e-6,"control bow")
    close(q0["event"]["ratio"],1.25,1e-6,"control reload event")
    close(qstar["bow_mm"],REF["straightened_bow_mm"],1e-6,"straightened bow")
    close(qstar["event"]["ratio"],REF["re_yield_over_virgin"],1e-6,"re-yield ratio")
    print("QUICK PASS: q=0 and q*=1.175 fresh-specimen histories")


def full():
    # Keep the published reference and copied source immutable. All scratch
    # calculations happen inside a disposable isolated checkout.
    with tempfile.TemporaryDirectory(prefix="tread051_repro_") as temp:
        work=Path(temp)
        shutil.copytree(ROOT/"solver",work/"solver")
        shutil.copytree(ROOT/"protocol",work/"protocol")
        (work/"data").mkdir()
        shutil.copy2(ROOT/"data"/"pilot.json",work/"data"/"pilot.json")
        commands=[
          ["solver/tread_history051.py","sweep","--out","data/sweep_base.json"],
          ["solver/tread_history051.py","sweep","--nx","161","--out","data/sweep_nx161.json"],
          ["solver/tread_history051.py","sweep","--nz","201","--out","data/sweep_nz201.json"],
          ["solver/tread_history051.py","sweep","--steps","2","--out","data/sweep_steps2.json"],
          ["solver/tread_history051_sensitivity.py"],
          ["solver/tread_history051_validate.py"],
        ]
        for args in commands:
            p=subprocess.run([sys.executable,*args],cwd=work,capture_output=True,text=True)
            if p.returncode:
                raise RuntimeError(f"{' '.join(args)} failed:\n{p.stdout[-2500:]}\n{p.stderr[-2500:]}")
            print("PASS",' '.join(args))
        v=json.loads((work/"data"/"validation.json").read_text(encoding="utf-8"))
        d=json.loads((work/"data"/"sensitivity.json").read_text(encoding="utf-8"))
        if not v["all_pass"]:raise RuntimeError("registered validation gate failed")
        if v["q_star"]!=[REF["q_star"]]:raise RuntimeError("co-minimum changed")
        zero=d["cells"][0]
        cell=d["cells"][d["q_star_indices"][0]]
        close(zero["bow_mm"],REF["control_bow_mm"],1e-6,"control bow")
        close(cell["bow_mm"],REF["straightened_bow_mm"],1e-6,"straightened bow")
        close(cell["r"],REF["re_yield_over_virgin"],1e-6,"re-yield ratio")
        if len(d["sensitivity"])!=6:raise RuntimeError("missing sensitivity leg")
        print("FULL PASS: 51×4 primary cells, six sensitivity legs, 252 histories and all gates")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    args=parser.parse_args()
    verify_sources()
    quick() if args.quick else full()
