"""Independent physical/numerical gates for registered 051 results."""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from tread_history051 import B,E,FY,L,T,QGRID,Beam,Spec,moment,shape_from_curvature,hermite_max_abs
from tread_history051_sensitivity import choose,uncertainty

DATA=Path(__file__).resolve().parents[1]/"data"


def load(name):
    return json.loads((DATA/name).read_text(encoding="utf-8"))


def block_moment(k1:float,k2:float,fy:float) -> float:
    """Closed-piecewise integral of reverse stress, independent of Simpson/return map.

    Initial virgin loading makes a rectangular elastic core of half-depth z0.
    The reversed trial stress is affine on each of three initial-strain zones;
    introduce its ±fy crossing points and integrate each affine/clipped zone.
    """
    h=T/2
    z0=fy/(E*k1)
    segments=[(-h,-z0,E*(k2-k1),-fy),
              (-z0,z0,E*k2,0.),
              (z0,h,E*(k2-k1),fy)]
    total=0.
    for low,high,slope,intercept in segments:
        roots=[low,high]
        if abs(slope)>1e-16:
            for boundary in (-fy,fy):
                z=(boundary-intercept)/slope
                if low<z<high: roots.append(z)
        roots=sorted(set(roots))
        for lo,hi in zip(roots[:-1],roots[1:]):
            mid=(lo+hi)/2
            trial=slope*mid+intercept
            if trial>fy: aa,bb=0.,fy
            elif trial< -fy: aa,bb=0.,-fy
            else: aa,bb=slope,intercept
            total += aa*(hi**3-lo**3)/3+bb*(hi**2-lo**2)/2
    return B*total


def block_stress(z:np.ndarray,k1:float,k2:float,fy:float):
    ep=np.sign(z)*np.maximum(k1*np.abs(z)-fy/E,0.)
    return np.clip(E*(k2*z-ep),-fy,fy)


def mirrored_histories():
    """Gravity-off, sign-reversed production return-map histories for EPP/iso/kin."""
    results={}
    for name,s in {"epp":Spec(),"isotropic":Spec(hiso=.01*E),
                   "kinematic":Spec(hkin=.01*E)}.items():
        plus,minus=Beam(s),Beam(s)
        plus.k[:]=0.
        minus.k[:]=0.
        ky=s.fy/(E*T/2)
        path=np.r_[np.linspace(0,1.5*ky,61),np.linspace(1.5*ky,0,61)[1:],
                   np.linspace(0,-1.2*ky,61)[1:],np.linspace(-1.2*ky,0,61)[1:]]
        worst=0.
        for kp in path:
            for beam,k in ((plus,kp),(minus,-kp)):
                kval=np.full(len(beam.x),k)
                _,_,sigma,gamma,sign,_=beam._trial(kval)
                beam.ep+=gamma*sign
                beam.acc+=gamma
                beam.alpha+=s.hkin*gamma*sign
                beam.k=kval
            _,_,sp,_,_,_=plus._trial(plus.k)
            _,_,sm,_,_,_=minus._trial(minus.k)
            worst=max(worst,float(np.max(np.abs(sp+sm)))/s.fy,
                      float(np.max(np.abs(plus.ep+minus.ep))),
                      float(np.max(np.abs(plus.acc-minus.acc))),
                      float(np.max(np.abs(plus.alpha+minus.alpha)))/s.fy)
        results[name]=worst
    return results


def dense_analytic_control_bow():
    """An independent dense x-integration for the q=0 control only."""
    s=Spec()
    x=np.unique(np.r_[np.linspace(0,L,16001),(L-s.a)/2,(L+s.a)/2])
    m=moment(x,1.25*s.py,s)
    h=T/2
    core=np.sqrt(np.maximum(3-2*m/s.my,1e-15))*h
    kplastic=s.fy/(E*core)
    kfirst=np.where(m<=s.my,m/(E*s.iy),kplastic)
    delta_k=kfirst-m/(E*s.iy)
    dy,dth=shape_from_curvature(x,delta_k)
    return hermite_max_abs(x,dy,dth)*1000


def event_evidence(record,spec):
    ev=record["event"]
    pre=ev.get("pre",float("nan"))
    ratio=ev.get("ratio",float("nan"))
    if ev.get("immediate"):
        return (pre==0 and ratio==0 and ev.get("one_sided_confirmed") is True)
    if not (math.isfinite(pre) and 0<pre<=1.35*spec.py
            and math.isfinite(ratio) and abs(ratio-pre/spec.py)<=1e-10):
        return False
    brackets=ev.get("brackets",[])
    if len(brackets)!=2:return False
    for bracket,half in zip(brackets,(.002,.001)):
        keys=("lower_ratio","upper_ratio","below_acc_increment","above_acc_increment")
        if not all(key in bracket and isinstance(bracket[key],(int,float))
                   and math.isfinite(bracket[key]) for key in keys):
            return False
        if bracket["below_acc_increment"]<0 or bracket["above_acc_increment"]<=0:
            return False
        lower=max(0.,ratio-half)
        upper=min(1.35,ratio+half)
        if (abs(bracket.get("lower_ratio",float("nan"))-lower)>1e-10
            or abs(bracket.get("upper_ratio",float("nan"))-upper)>1e-10
            or bracket.get("below_acc_increment",float("nan"))>1e-12
            or bracket.get("above_acc_increment",float("nan"))<=1e-12):
            return False
    return True


def anchors(q:float):
    s=Spec()
    my_analytic=s.fy*B*T*T/6
    py_analytic=(my_analytic-s.weight*L*L/8)/(L/4-s.a/8)
    virgin=Beam(s)
    virgin_numeric_event=virgin.event()
    idx=int(np.argmin(np.abs(virgin.x-L/2)))
    virgin.step(py_analytic)
    _,_,sig,_,_,_=virgin._trial(virgin.k)
    yield_gap=abs(float(np.max(np.abs(sig)))-s.fy)/s.fy
    linear_gap=abs(virgin.k[idx]-s.my/(E*s.iy))/(s.my/(E*s.iy))
    pure=[]
    z=virgin.z
    w=virgin.wz
    ky=s.fy/(E*T/2)
    for factor in (1.2,1.5,2.0):
        sig=np.clip(E*(factor*ky)*z,-s.fy,s.fy)
        numeric=B*float(np.sum(sig*z*w))/s.my
        exact=1.5-.5/factor**2
        pure.append({"kappa_over_ky":factor,"numeric":numeric,"analytic":exact,
                     "relative_gap":abs(numeric-exact)/exact})
    # Centre reverse case: first max is analytically solvable for ideal-plastic rectangle.
    m1=float(moment(np.array([L/2]),1.25*s.py,s)[0])
    core=(T/2)*math.sqrt(3-2*m1/s.my)
    k1=s.fy/(E*core)
    m2=float(moment(np.array([L/2]),-q*s.py,s)[0])
    lo,hi=-1.,k1
    if not block_moment(k1,lo,s.fy)<m2<block_moment(k1,hi,s.fy):
        raise RuntimeError("independent reverse stress-block root not bracketed")
    for _ in range(90):
        mid=(lo+hi)/2
        if block_moment(k1,mid,s.fy)<m2: lo=mid
        else: hi=mid
    k2=(lo+hi)/2
    history=Beam(s)
    history.ramp(1.25*s.py,100)
    k1_numeric=float(history.k[idx])
    acc_at_max=history.acc.copy()
    history.ramp(0,100)
    k_unload=float(history.k[idx])
    history.step(-.001*s.py)
    unload_no_flow=float(np.max(history.acc-acc_at_max))
    dm=float(moment(np.array([L/2]),-.001*s.py,s)[0]-moment(np.array([L/2]),0,s)[0])
    elastic_slope=abs((history.k[idx]-k_unload)/dm-1/(E*s.iy))/(1/(E*s.iy))
    history.ramp(-q*s.py,int(round((q-.001)/.025))*4)
    _,_,sigma_numeric,_,_,_=history._trial(history.k)
    sigma_analytic=block_stress(z,k1,k2,s.fy)
    block={"q":q,"analytic_k1":k1,"numeric_k1":k1_numeric,
           "k1_relative_gap":abs(k1_numeric-k1)/abs(k1),
           "analytic_reverse_kappa":k2,"numeric_reverse_kappa":float(history.k[idx]),
           "reverse_kappa_relative_gap":abs(float(history.k[idx])-k2)/max(abs(k2),1e-12),
           "max_stress_gap_over_fy":float(np.max(np.abs(sigma_numeric[idx]-sigma_analytic)))/s.fy,
           "analytic_block_moment_gap_over_my":abs(block_moment(k1,k2,s.fy)-m2)/s.my,
           "elastic_unload_slope_relative_gap":elastic_slope,
           "elastic_unload_acc_increment":unload_no_flow}
    sign_symmetry=mirrored_histories()
    # Total post-history gravity moment and excess history resultant.
    history.ramp(0,100)
    _,_,post_sigma,_,_,_=history._trial(history.k)
    virgin_g_sigma=E*moment(history.x,0,s)[:,None]/(E*s.iy)*z[None,:]
    excess=post_sigma-virgin_g_sigma
    excess_n=float(np.max(np.abs(B*np.sum(excess*w,axis=1))))/(s.fy*B*T)
    excess_m=float(np.max(np.abs(B*np.sum(excess*z*w,axis=1))))/s.my
    return {"analytic_py_N":py_analytic,"solver_py_N":virgin_numeric_event["pre"],"yield_gap":yield_gap,
            "virgin_linear_gap":linear_gap,"pure_bending":pure,"reverse_block":block,
            "gravity_off_mirrored_histories":sign_symmetry,
            "dense_analytic_control_bow_mm":dense_analytic_control_bow(),
            "history_excess_axial_rel":excess_n,"history_excess_moment_rel":excess_m}


def main():
    start=time.perf_counter()
    sweep_names=("sweep_base.json","sweep_nx161.json","sweep_nz201.json","sweep_steps2.json")
    sweeps=[load(name) for name in sweep_names]
    expected_resolutions=((81,101,1),(161,101,1),(81,201,1),(81,101,2))
    sens=load("sensitivity.json")
    rows=sweeps[0]["result"]
    completeness=True
    for run,(nx,nz,steps) in zip(sweeps,expected_resolutions):
        completeness &= (run.get("mode")=="sweep" and run.get("protocol")=="051 v2+v3"
                         and (run.get("nx"),run.get("nz"),run.get("step_multiplier"))==(nx,nz,steps)
                         and run.get("spec")==vars(Spec()) and len(run.get("result",[]))==51
                         and abs(run.get("py_N",0)-Spec().py)<1e-9)
        completeness &= all(abs(row.get("q",-99)-float(QGRID[i]))<1e-12
                            for i,row in enumerate(run.get("result",[])))
    if not completeness:raise RuntimeError("exact ordered q grid / sweep metadata incomplete")
    cells,selected=choose(rows,[run["result"] for run in sweeps[1:]])
    completeness &= (sens.get("q_star_indices")==selected
                     and sens.get("q_star_values")==[float(QGRID[i]) for i in selected]
                     and len(sens.get("cells",[]))==51)
    if not selected:raise RuntimeError("no selected cell")
    anchor=anchors(cells[selected[0]]["q"])
    u_py=max(abs(anchor["analytic_py_N"]-anchor["solver_py_N"])/anchor["analytic_py_N"],
             anchor["yield_gap"])
    control=cells[0]
    primary=[]
    for i in selected:
        c=cells[i]
        upper_b=(c["bow_mm"]+c["u_b_mm"])/(control["bow_mm"]-control["u_b_mm"])
        upper_r=(c["r"]+c["u_r"])/(1-u_py)
        primary.append({"q":c["q"],"upper_bow_ratio":upper_b,
                        "upper_reyield_ratio":upper_r,"pass":c["r"]>0
                        and not rows[i]["event"]["immediate"]
                        and upper_b<=.10 and upper_r<=.90
                        and c["u_b_mm"]<=.005 and c["u_r"]<=.005})
    anchors_pass=(anchor["yield_gap"]<1e-6 and anchor["virgin_linear_gap"]<1e-6
                  and max(p["relative_gap"] for p in anchor["pure_bending"])<.005
                  and anchor["reverse_block"]["k1_relative_gap"]<.005
                  and anchor["reverse_block"]["reverse_kappa_relative_gap"]<.005
                  and anchor["reverse_block"]["max_stress_gap_over_fy"]<.005
                  and anchor["reverse_block"]["analytic_block_moment_gap_over_my"]<1e-8
                  and anchor["reverse_block"]["elastic_unload_slope_relative_gap"]<1e-6
                  and anchor["reverse_block"]["elastic_unload_acc_increment"]<1e-12
                  and max(anchor["gravity_off_mirrored_histories"].values())<1e-12
                  and anchor["history_excess_axial_rel"]<1e-7
                  and anchor["history_excess_moment_rel"]<1e-7)
    expected_specs={"isotropic_0p01E":Spec(hiso=.01*E),
                    "kinematic_0p01E":Spec(hkin=.01*E),
                    "patch_0p090":Spec(a=.09),"patch_0p110":Spec(a=.11),
                    "fy_99MPa":Spec(fy=99e6),"fy_121MPa":Spec(fy=121e6)}
    sensitivity_pass=set(sens.get("sensitivity",{}))==set(expected_specs)
    sensitivity_checks={}
    all_records=[(r,Spec()) for run in sweeps for r in run["result"]]
    q0_pass=all(abs(run["result"][0]["event"]["ratio"]-1.25)<=.005 for run in sweeps)
    for name,spec in expected_specs.items():
        leg=sens.get("sensitivity",{}).get(name,{})
        correct_spec=leg.get("spec")==vars(spec) and abs(leg.get("py_N",0)-spec.py)<1e-9
        control_case=leg.get("control",{})
        selected_cases=leg.get("selected",{})
        expected_keys={str(float(QGRID[i])) for i in selected}
        correct_keys=set(selected_cases)==expected_keys
        sensitivity_pass &= correct_spec and correct_keys
        if not correct_spec or not correct_keys:
            sensitivity_checks[name]={"pass":False,"reason":"missing or wrong spec/cell"}
            continue
        batches=[(0.,control_case)]+[(float(QGRID[i]),selected_cases[str(float(QGRID[i]))])
                                     for i in selected]
        complete_batches=True
        for q,case in batches:
            variants=[case.get("base")]+case.get("refinements",[])
            complete=(len(variants)==4 and all(
                isinstance(v,dict) and abs(v.get("q",-99)-q)<1e-12
                and (v.get("nx"),v.get("nz"),v.get("step_multiplier"))==expected_resolutions[j]
                and v.get("spec")==vars(spec)
                for j,v in enumerate(variants)))
            complete_batches &= complete
            if not complete:continue
            all_records.extend((v,spec) for v in variants)
        sensitivity_pass &= complete_batches
        if not complete_batches:
            sensitivity_checks[name]={"pass":False,"reason":"missing refinements"}
            continue
        cb=control_case["base"]
        cub,cur=uncertainty(cb,control_case["refinements"])
        q0_pass &= abs(cb["event"]["ratio"]-1.25)<=.005
        numeric_py=Beam(spec).event()["pre"]
        leg_u_py=abs(numeric_py-spec.py)/spec.py
        leg_checks=[]
        for i in selected:
            case=selected_cases[str(float(QGRID[i]))]
            b=case["base"]
            ub,ur=uncertainty(b,case["refinements"])
            b_ratio=(b["bow_mm"]+ub)/(cb["bow_mm"]-cub)
            r_ratio=(b["event"]["ratio"]+ur)/(1-leg_u_py)
            passed=(b["event"]["pre"]>0 and not b["event"]["immediate"]
                    and b_ratio<=.20 and r_ratio<1 and ub<=.005 and ur<=.005)
            leg_checks.append({"q":b["q"],"upper_bow_ratio":b_ratio,
                               "upper_reyield_ratio":r_ratio,"pass":passed})
        leg_pass=all(c["pass"] for c in leg_checks)
        sensitivity_pass &= leg_pass
        sensitivity_checks[name]={"pass":leg_pass,"cells":leg_checks,
                                  "numeric_virgin_py_N":numeric_py,
                                  "virgin_py_relative_gap":leg_u_py}
    completeness &= len(all_records)==252
    scope_pass=all(r["max_strain"]<=.005 and r["max_slope_rad"]<=.12
                   and r["max_y_over_L"]<=.04 and r["max_n_rel"]<=1e-7
                   and r["max_m_rel"]<=1e-7 for r,_ in all_records)
    event_pass=all(event_evidence(r,spec) for r,spec in all_records)
    precision_pass=all(c["u_b_mm"]<=.005 for c in cells)
    pilot=load("pilot.json")
    elapsed_primary=sweeps[0]["elapsed_s"]
    elapsed_other=sum(run["elapsed_s"] for run in sweeps[1:])+sens["elapsed_s"]+time.perf_counter()-start
    elapsed_total=elapsed_primary+elapsed_other
    compute_pass=(pilot["elapsed_s"]<=60 and elapsed_primary<=1800
                  and elapsed_other<=1800 and elapsed_total<=3600)
    dense=anchor["dense_analytic_control_bow_mm"]
    control_operational_interval=[control["bow_mm"]-control["u_b_mm"],
                                  control["bow_mm"]+control["u_b_mm"]]
    gates={"complete_registered_matrix":completeness,"anchors":anchors_pass,
           "scope_all_252":scope_pass,"event_brackets_all_252":event_pass,
           "q0_exact":q0_pass,"all_cell_bow_precision":precision_pass,
           "primary":all(p["pass"] for p in primary),"six_sensitivities":sensitivity_pass,
           "compute_budget":compute_pass}
    sources=[Path("solver/tread_history051.py"),Path("solver/tread_history051_sensitivity.py"),
             Path("solver/tread_history051_validate.py"),
             Path("protocol/051_registration_tread_history_v2.md"),
             Path("protocol/051_registration_tread_history_v3.md")]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    result={"protocol":"051 v2+v3","anchors":anchor,"primary":primary,
            "sensitivity_recomputed":sensitivity_checks,"virgin_py_relative_gap":u_py,
            "control_operational_bow_interval_mm":control_operational_interval,
            "dense_analytic_control_bow_outside_operational_interval":
                not control_operational_interval[0]<=dense<=control_operational_interval[1],
            "uncertainty_note":"Summed refinement shifts are operational estimates, not certified bounds; independent dense analytic control bow may lie outside the saved operational interval. Do not present interval as guaranteed coverage.",
            "gates":gates,"all_pass":all(gates.values()),
            "q_star":sens["q_star_values"],"compute_elapsed_s":elapsed_total,
            "source_sha256":hashes}
    (DATA/"validation.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({"gates":gates,"primary":primary,"anchors":anchor,
                      "compute_elapsed_s":elapsed_total},indent=2))


if __name__=="__main__":main()
