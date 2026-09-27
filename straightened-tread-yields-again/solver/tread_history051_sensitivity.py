"""Registered 051 selection, per-cell error envelopes, and fixed-q sensitivities."""
from __future__ import annotations

import json
import time
from pathlib import Path

from tread_history051 import A,E,FY,QGRID,Spec,one_case

DATA=Path(__file__).resolve().parents[1]/"data"


def read(name):
    return json.loads((DATA/name).read_text(encoding="utf-8"))["result"]


def uncertainty(base,alts):
    ub=sum(abs(alt["bow_mm"]-base["bow_mm"]) for alt in alts)
    ur=sum(abs(alt["event"]["ratio"]-base["event"]["ratio"]) for alt in alts)+.002
    return ub,ur


def choose(base,refinements):
    cells=[]
    for i,b in enumerate(base):
        ub,ur=uncertainty(b,[r[i] for r in refinements])
        cells.append({"q":b["q"],"bow_mm":b["bow_mm"],"u_b_mm":ub,
                      "r":b["event"]["ratio"],"u_r":ur,
                      "bow_interval_mm":[b["bow_mm"]-ub,b["bow_mm"]+ub]})
    minimum=min(range(len(cells)),key=lambda i:cells[i]["bow_mm"])
    included={minimum}
    while True:
        before=len(included)
        for i,c in enumerate(cells):
            if any(c["bow_interval_mm"][0]<=cells[j]["bow_interval_mm"][1]
                   and c["bow_interval_mm"][1]>=cells[j]["bow_interval_mm"][0]
                   for j in included):
                included.add(i)
        if len(included)==before: break
    return cells,sorted(included)


def run_case(spec,q):
    variants=[one_case(spec,q,81,101,1,True),
              one_case(spec,q,161,101,1,True),
              one_case(spec,q,81,201,1,True),
              one_case(spec,q,81,101,2,True)]
    ub,ur=uncertainty(variants[0],variants[1:])
    return {"base":variants[0],"refinements":variants[1:],"u_b_mm":ub,"u_r":ur}


def main():
    start=time.perf_counter()
    base=read("sweep_base.json")
    refs=[read("sweep_nx161.json"),read("sweep_nz201.json"),read("sweep_steps2.json")]
    cells,selected=choose(base,refs)
    legs={"isotropic_0p01E":Spec(hiso=.01*E),
          "kinematic_0p01E":Spec(hkin=.01*E),
          "patch_0p090":Spec(a=.09),"patch_0p110":Spec(a=.11),
          "fy_99MPa":Spec(fy=99e6),"fy_121MPa":Spec(fy=121e6)}
    outcomes={}
    for name,spec in legs.items():
        control=run_case(spec,0.)
        chosen={}
        for i in selected:
            selected_case=run_case(spec,float(QGRID[i]))
            b=selected_case["base"]["bow_mm"]
            b0=control["base"]["bow_mm"]
            upper_b=(b+selected_case["u_b_mm"])/(b0-control["u_b_mm"])
            upper_r=selected_case["base"]["event"]["ratio"]+selected_case["u_r"]
            selected_case["upper_bow_ratio"]=upper_b
            selected_case["upper_reyield_ratio"]=upper_r
            selected_case["pass"]=(selected_case["base"]["event"]["pre"]>0
                                   and not selected_case["base"]["event"]["immediate"]
                                   and upper_b<=.20 and upper_r<1.)
            chosen[str(QGRID[i])]=selected_case
        outcomes[name]={"spec":vars(spec),"py_N":spec.py,"control":control,
                        "selected":chosen,"pass":all(v["pass"] for v in chosen.values())}
        print(name,"PASS" if outcomes[name]["pass"] else "FAIL",
              [(q,round(v["upper_bow_ratio"],5),round(v["upper_reyield_ratio"],5)) for q,v in chosen.items()],flush=True)
    target=(DATA/"sensitivity.json")
    data={"protocol":"051 v2+v3","q_star_indices":selected,
          "q_star_values":[float(QGRID[i]) for i in selected],"cells":cells,
          "sensitivity":outcomes,"elapsed_s":time.perf_counter()-start}
    target.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print("selected",data["q_star_values"],"all_sensitivity_pass",all(v["pass"] for v in outcomes.values()),
          "elapsed_s",round(data["elapsed_s"],2))


if __name__=="__main__":main()
