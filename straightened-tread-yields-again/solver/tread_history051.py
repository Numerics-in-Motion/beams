"""Registered Spinoff 051 elastic-plastic tread-strip experiment.

Protocol: analyze/260919/051_registration_tread_history_v2.md plus v3.md.
The section return map is the production model; analytic anchors live separately.
No media code may reinterpret failed scientific gates as a released result.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


L, B, T, E, FY, RHO, G, A = .8, .28, .012, 70e9, 110e6, 2700., 9.81, .1
QGRID = np.arange(51, dtype=float) * .025


@dataclass(frozen=True)
class Spec:
    a: float = A
    fy: float = FY
    hiso: float = 0.
    hkin: float = 0.

    @property
    def iy(self) -> float:
        return B * T**3 / 12

    @property
    def my(self) -> float:
        return self.fy * B * T**2 / 6

    @property
    def weight(self) -> float:
        return RHO * B * T * G

    @property
    def mg(self) -> float:
        return self.weight * L**2 / 8

    @property
    def c(self) -> float:
        return L / 4 - self.a / 8

    @property
    def py(self) -> float:
        return (self.my - self.mg) / self.c


def simpson_weights(n: int, depth: float = T) -> np.ndarray:
    if n % 2 != 1:
        raise ValueError("Simpson nodes must be odd")
    h = depth / (n - 1)
    w = np.ones(n)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    return w * h / 3


def grid(nx: int, a: float) -> np.ndarray:
    return np.unique(np.r_[np.linspace(0, L, nx), (L-a)/2, (L+a)/2])


def moment(x: np.ndarray, p: float, spec: Spec) -> np.ndarray:
    left, right = (L-spec.a)/2, (L+spec.a)/2
    reaction = (spec.weight * L + p) / 2
    return reaction*x - spec.weight*x*x/2 - p/(2*spec.a)*(
        np.maximum(x-left, 0)**2 - np.maximum(x-right, 0)**2
    )


def shape_from_curvature(x: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Integrate linearly interpolated curvature exactly on each segment."""
    h = np.diff(x)
    theta = np.zeros_like(x)
    y = np.zeros_like(x)
    theta[1:] = np.cumsum(h * (k[:-1]+k[1:])/2)
    y[1:] = np.cumsum(h * theta[:-1] + h*h*(k[:-1]/3+k[1:]/6))
    rotation = -y[-1]/L
    return y+rotation*x, theta+rotation


def hermite_max_abs(x: np.ndarray, y: np.ndarray, slopes: np.ndarray) -> float:
    """Global |y| extrema of cubic Hermite interpolant, including interior roots."""
    top = float(np.max(np.abs(y)))
    for i, h in enumerate(np.diff(x)):
        a = 2*y[i]-2*y[i+1]+h*(slopes[i]+slopes[i+1])
        b = -3*y[i]+3*y[i+1]-h*(2*slopes[i]+slopes[i+1])
        c = h*slopes[i]
        roots = np.roots([3*a, 2*b, c]) if abs(a)>1e-16 else np.roots([2*b,c]) if abs(b)>1e-16 else []
        for r in roots:
            if abs(np.imag(r)) < 1e-12 and 0 < np.real(r) < 1:
                u = float(np.real(r))
                top = max(top, abs(((a*u+b)*u+c)*u+y[i]))
    return top


class Beam:
    def __init__(self, spec: Spec, nx: int = 81, nz: int = 101):
        self.spec, self.x, self.z = spec, grid(nx,spec.a), np.linspace(-T/2,T/2,nz)
        self.wz = simpson_weights(nz)
        self.ep = np.zeros((len(self.x),nz))
        self.acc = np.zeros_like(self.ep)
        self.alpha = np.zeros_like(self.ep)
        self.k = moment(self.x,0,spec)/(E*spec.iy)
        self.p = 0.
        self.max_strain = float(np.max(np.abs(self.k[:,None]*self.z)))
        self.max_slope = 0.
        self.max_y = 0.
        self.max_n_rel = 0.
        self.max_m_rel = 0.
        self.stages = 0
        self._record_scope()

    def clone(self) -> Beam:
        return copy.deepcopy(self)

    def _trial(self, k: np.ndarray):
        s = self.spec
        strain = k[:,None]*self.z[None,:]
        trial = E*(strain-self.ep)
        eta = trial-self.alpha
        signed = np.where(eta>=0,1.,-1.)
        f = np.abs(eta)-(s.fy+s.hiso*self.acc)
        active = f > 0
        gamma = np.maximum(f,0)/(E+s.hiso+s.hkin)
        sigma = trial-E*gamma*signed
        tangent = np.where(active,E*(s.hiso+s.hkin)/(E+s.hiso+s.hkin),E)
        m = B*np.sum(sigma*self.z*self.wz,axis=1)
        dm = B*np.sum(tangent*self.z*self.z*self.wz,axis=1)
        return m,dm,sigma,gamma,signed,strain

    def _record_scope(self):
        y,theta = shape_from_curvature(self.x,self.k)
        self.max_slope=max(self.max_slope,float(np.max(np.abs(theta))))
        self.max_y=max(self.max_y,float(np.max(np.abs(y))))
        if self.max_strain > .005+1e-9 or self.max_slope>.12+1e-9 or self.max_y/L>.04+1e-9:
            raise RuntimeError(f"scope: strain={self.max_strain:.6g}, slope={self.max_slope:.6g}, y/L={self.max_y/L:.6g}")

    def step(self,p:float):
        s=self.spec
        target=moment(self.x,p,s)
        k=self.k+(target-moment(self.x,self.p,s))/(E*s.iy)
        for _ in range(50):
            m,dm,sigma,gamma,sign,strain=self._trial(k)
            err=m-target
            if np.max(np.abs(err)) <= 1e-9*s.my:
                break
            if np.min(dm) <= 0:
                raise RuntimeError("nonpositive section tangent")
            update=np.clip(err/dm,-.02,.02)
            k=k-update
        else:
            raise RuntimeError(f"section Newton failed at P/Py={p/s.py:.6f}; moment residual {np.max(np.abs(err))}")
        self.k=k
        self.ep += gamma*sign
        self.acc += gamma
        self.alpha += s.hkin*gamma*sign
        self.p=p
        self.stages+=1
        n=B*np.sum(sigma*self.wz,axis=1)
        self.max_n_rel=max(self.max_n_rel,float(np.max(np.abs(n))/(s.fy*B*T)))
        self.max_m_rel=max(self.max_m_rel,float(np.max(np.abs(m-target))/s.my))
        self.max_strain=max(self.max_strain,float(np.max(np.abs(strain))))
        self._record_scope()
        if self.max_n_rel>1e-7 or self.max_m_rel>1e-7:
            raise RuntimeError(f"equilibrium: N={self.max_n_rel}, M={self.max_m_rel}")

    def ramp(self,to:float,n:int):
        for p in np.linspace(self.p,to,n+1)[1:]:
            self.step(float(p))

    def bow(self, virgin_y:np.ndarray,virgin_th:np.ndarray) -> float:
        y,theta=shape_from_curvature(self.x,self.k)
        return hermite_max_abs(self.x,y-virgin_y,theta-virgin_th)

    def event(self) -> dict:
        """Frozen-state outward elastic search, in units of applied patch load."""
        s=self.spec
        _,_,sigma,_,_,_=self._trial(self.k)
        eta=sigma-self.alpha
        f=np.abs(eta)-(s.fy+s.hiso*self.acc)
        dmdp=moment(self.x,1,s)-moment(self.x,0,s)
        dsdp=E*self.z[None,:]*dmdp[:,None]/(E*s.iy)
        # At an initially inward-moving surface the trial stress can cross
        # the *opposite* yield surface later; current-sign filtering misses it.
        direction=np.where(dsdp>=0,1.,-1.)
        gap=(s.fy+s.hiso*self.acc)-direction*eta
        eligible=np.abs(dsdp)>0
        immediate=eligible & (gap<=1e-8*s.fy)
        if np.any(immediate):
            probe=self.clone()
            before=probe.acc.copy()
            probe.step(1e-5*s.py)
            if float(np.max(probe.acc-before)) <= 1e-12:
                raise RuntimeError("immediate event lacks one-sided plastic-flow confirmation")
            return {"pre":0.,"ratio":0.,"immediate":True,"one_sided_confirmed":True,
                    "location_m":float(self.x[np.where(immediate)[0][0]])}
        roots=np.where(eligible,gap/np.maximum(np.abs(dsdp),1e-30),np.inf)
        idx=np.unravel_index(np.argmin(roots),roots.shape)
        pre=float(roots[idx])
        if not math.isfinite(pre) or pre>1.35*s.py:
            raise RuntimeError("re-yield event not bracketed below 1.35 Py")
        probes=[]
        for half in (.002,.001):
            lo=max(0.,pre-half*s.py)
            hi=min(1.35*s.py,pre+half*s.py)
            outside=[]
            for p_probe in (lo,hi):
                probe=self.clone()
                before=probe.acc.copy()
                probe.step(p_probe)
                outside.append(float(np.max(probe.acc-before)))
            if outside[0]>1e-12 or outside[1]<=1e-12:
                raise RuntimeError(f"event bracket failed at {pre/s.py:.6f}: {outside}")
            probes.append({"lower_ratio":lo/s.py,"upper_ratio":hi/s.py,
                           "below_acc_increment":outside[0],"above_acc_increment":outside[1]})
        return {"pre":pre,"ratio":pre/s.py,"immediate":False,"location_m":float(self.x[idx[0]]),
                "z_m":float(self.z[idx[1]]),"brackets":probes}

    def diagnostic(self,steps:int):
        self.ramp(1.35*self.spec.py,steps)


def one_case(spec:Spec,q:float,nx=81,nz=101,step_mult=1,diagnostic=False) -> dict:
    beam=Beam(spec,nx,nz)
    virgin_y,virgin_th=shape_from_curvature(beam.x,beam.k)
    beam.ramp(1.25*spec.py,100*step_mult)
    beam.ramp(0.,100*step_mult)
    if q:
        beam.ramp(-q*spec.py,max(1,int(round(q/.025))*4*step_mult))
        beam.ramp(0.,100*step_mult)
    bow=beam.bow(virgin_y,virgin_th)
    ev=beam.event()
    if diagnostic:
        beam.diagnostic(100*step_mult)
    return {"q":q,"nx":nx,"nz":nz,"step_multiplier":step_mult,
            "spec":vars(spec),"bow_mm":bow*1000,"event":ev,"max_strain":beam.max_strain,
            "max_slope_rad":beam.max_slope,"max_y_over_L":beam.max_y/L,
            "max_n_rel":beam.max_n_rel,"max_m_rel":beam.max_m_rel,"stages":beam.stages}


def sweep(spec:Spec,nx=81,nz=101,step_mult=1,qs=QGRID) -> list[dict]:
    beam=Beam(spec,nx,nz)
    virgin_y,virgin_th=shape_from_curvature(beam.x,beam.k)
    beam.ramp(1.25*spec.py,100*step_mult)
    beam.ramp(0.,100*step_mult)
    records=[]
    qlast=0.
    for q in qs:
        q=float(q)
        if q < qlast-1e-12:
            raise ValueError("q must be sorted")
        if q>qlast:
            n=max(1,int(round((q-qlast)/.025))*4*step_mult)
            beam.ramp(-q*spec.py,n)
        unloaded=beam.clone()
        unloaded.ramp(0.,100*step_mult)
        bow_mm=unloaded.bow(virgin_y,virgin_th)*1000
        event=unloaded.event()
        unloaded.diagnostic(100*step_mult)
        records.append({"q":q,"bow_mm":bow_mm,
                        "event":event,"max_strain":unloaded.max_strain,
                        "max_slope_rad":unloaded.max_slope,"max_y_over_L":unloaded.max_y/L,
                        "max_n_rel":unloaded.max_n_rel,"max_m_rel":unloaded.max_m_rel,
                        "stages":unloaded.stages})
        qlast=q
    return records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=["pilot","sweep","single"])
    parser.add_argument("--out",type=Path)
    parser.add_argument("--nx",type=int,default=81)
    parser.add_argument("--nz",type=int,default=101)
    parser.add_argument("--steps",type=int,default=1)
    parser.add_argument("--q",type=float,default=0.)
    parser.add_argument("--a",type=float,default=A)
    parser.add_argument("--fy",type=float,default=FY)
    parser.add_argument("--hiso",type=float,default=0.)
    parser.add_argument("--hkin",type=float,default=0.)
    args=parser.parse_args()
    spec=Spec(args.a,args.fy,args.hiso,args.hkin)
    start=time.perf_counter()
    if args.mode=="pilot":
        result=one_case(spec,0.,args.nx,args.nz,args.steps,True)
    elif args.mode=="single":
        result=one_case(spec,args.q,args.nx,args.nz,args.steps,True)
    else:
        result=sweep(spec,args.nx,args.nz,args.steps)
    data={"protocol":"051 v2+v3","mode":args.mode,"spec":vars(spec),
          "nx":args.nx,"nz":args.nz,"step_multiplier":args.steps,
          "py_N":spec.py,"my_Nm":spec.my,"mg_Nm":spec.mg,
          "elapsed_s":time.perf_counter()-start,"result":result}
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps({"mode":args.mode,"elapsed_s":data["elapsed_s"],"py_N":spec.py,
                      "result":result if args.mode!="sweep" else f"{len(result)} cells"}))


if __name__=="__main__":
    main()
