"""One office floor, three beams, one standard fire: aluminium, steel or glulam?

Each beam is sized (continuously, no catalogue rounding) so that its AMBIENT design bending resistance equals the
same office-floor demand under its own Eurocode (rule R1, the only rule the claim uses):
    G = 6.0 kN/m slab + finishes + the beam's own weight (iterated), Q = 6.0 kN/m (office, psi_1 = 0.5), L = 5.0 m
    M_Ed = (1.35 G + 1.5 Q) L^2 / 8          M_fi = (G + 0.5 Q) L^2 / 8
Simply supported, top face shielded and continuously restrained by the slab (no lateral-torsional buckling), no
composite action, free to expand. Heated from below and both sides (three sides) by the ISO 834 standard fire.
Endpoint: the first time the modelled SECTIONAL bending resistance falls below M_fi (not floor collapse).

  steel      S355, IPE proportions (IPE 240 scaled), EN 1993-1-2: eq. 4.25 lumped heating, c_a(theta) 3.4.1.2,
             eps_m 0.7, alpha_c 25, k_sh = 0.9 [A_m/V]_b/[A_m/V]; resistance k_y,theta W_pl f_y / (kappa_1 kappa_2),
             kappa_1 = 0.7 (unprotected, three sides, concrete slab), kappa_2 = 1.
  aluminium  EN AW-6082 T6, f_o 260 MPa, gamma_M1 1.1, same proportions, EN 1999-1-2: c = 903 + 0.41 theta,
             rho 2700, eps 0.3 (clean surface), Table 1 k_o,theta; resistance k_o,theta W_pl f_o.
  glulam     GL24h, b = h/2, k_mod 0.8, gamma_M 1.25 (ambient); EN 1995-1-2 reduced cross-section,
             beta_n 0.7 mm/min, d_0 7 mm with k_0 = t/20 for t < 20 min, k_fi 1.15, three-sided charring.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

EDITIONS = {"steel": "EN 1993-1-2:2005", "aluminium": "EN 1999-1-2:2007", "timber": "EN 1995-1-2:2004"}
L = 5.0
G_SLAB, Q_IMP, PSI1 = 6.0e3, 6.0e3, 0.5
SIGMA = 5.67e-8
T_END, DT = 7200.0, 1.0
ORDER = ("aluminium", "steel", "glulam")

# EN 1993-1-2 Table 3.1
ST_T = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
ST_KY = [1, 1, 1, 1, 1, .78, .47, .23, .11, .06, .04, .02, 0]
# EN 1999-1-2 Table 1, EN AW-6082 T6
AL_T = [20, 100, 150, 200, 250, 300, 350, 550]
AL_K = [1.00, 0.90, 0.79, 0.65, 0.38, 0.20, 0.11, 0.0]
# IPE 240 (m, m^2, m^3)
IPE = dict(h=0.240, b=0.120, tw=0.0062, tf=0.0098, r=0.015, A=3.912e-3, Wpl=366.6e-6)


def iso834(t_s):
    return 20.0 + 345.0*math.log10(8.0*t_s/60.0 + 1.0)


def c_steel(th):
    if th < 600:
        return 425 + 7.73e-1*th - 1.69e-3*th**2 + 2.22e-6*th**3
    if th < 735:
        return 666 + 13002/(738 - th)
    if th < 900:
        return 545 + 17820/(th - 731)
    return 650.0


def c_al(th):
    return 903.0 + 0.41*th


def ipe(scale):
    g = {k: IPE[k]*scale for k in ("h", "b", "tw", "tf", "r")}
    g["A"] = IPE["A"]*scale**2
    g["Wpl"] = IPE["Wpl"]*scale**3
    perim3 = (2*g["b"] - g["tw"]) + 2*(g["h"] - 2*g["tf"]) + g["b"] + 2*g["tf"]
    g["AmV"] = perim3/g["A"]
    g["AmV_box"] = (2*g["h"] + g["b"])/g["A"]
    g["ksh"] = 0.9*g["AmV_box"]/g["AmV"]
    return g


def demands(self_weight, Q=Q_IMP):
    G = G_SLAB + self_weight
    return (1.35*G + 1.5*Q)*L*L/8, (G + PSI1*Q)*L*L/8


def size_steel(Q=Q_IMP, fy=355e6, rho=7850.0):
    sw = 300.0
    for _ in range(40):
        MEd, Mfi = demands(sw, Q)
        g = ipe((MEd/fy/IPE["Wpl"])**(1/3)); sw = g["A"]*rho*9.81
    g.update(MEd=MEd, Mfi=Mfi, fy=fy, mu0=Mfi/(g["Wpl"]*fy), self_weight=sw)
    return g


def size_aluminium(Q=Q_IMP, fo=260e6, gm1=1.1, rho=2700.0):
    sw = 100.0
    for _ in range(40):
        MEd, Mfi = demands(sw, Q)
        g = ipe((MEd*gm1/fo/IPE["Wpl"])**(1/3)); sw = g["A"]*rho*9.81
    g.update(MEd=MEd, Mfi=Mfi, fo=fo, mu0=Mfi/(g["Wpl"]*fo), self_weight=sw)
    return g


def size_glulam(Q=Q_IMP, fmk=24e6, kmod=0.8, gm=1.25, rho=420.0):
    sw = 200.0
    for _ in range(40):
        MEd, Mfi = demands(sw, Q)
        W = MEd/(kmod*fmk/gm); h = (12*W)**(1/3); b = h/2; sw = b*h*rho*9.81
    return dict(b=b, h=h, MEd=MEd, Mfi=Mfi, fmk=fmk, self_weight=sw, mu0=Mfi/(b*h*h/6*1.15*fmk))


def class_checks(st, al):
    eps_fi = 0.85*math.sqrt(235e6/st["fy"])
    cf = (st["b"] - st["tw"] - 2*st["r"])/2/st["tf"]; cw = (st["h"] - 2*st["tf"] - 2*st["r"])/st["tw"]
    eps_al = math.sqrt(250e6/al["fo"])
    bf = (al["b"] - al["tw"] - 2*al["r"])/2/al["tf"]; bw = 0.4*(al["h"] - 2*al["tf"] - 2*al["r"])/al["tw"]
    return dict(steel_flange=(cf, 9*eps_fi), steel_web=(cw, 72*eps_fi), steel_ok=cf <= 9*eps_fi and cw <= 72*eps_fi,
                al_flange=(bf, 4.5*eps_al), al_web=(bw, 16*eps_al), al_ok=bf <= 4.5*eps_al and bw <= 16*eps_al)


def heat(amv, ksh, rho, c_of, eps, alpha_c=25.0, dt=DT, t_end=T_END):
    n = int(round(t_end/dt))
    t = np.arange(n + 1)*dt
    th = np.empty(n + 1); th[0] = 20.0
    for i in range(n):
        tg = iso834(t[i])
        h = alpha_c*(tg - th[i]) + eps*SIGMA*((tg + 273)**4 - (th[i] + 273)**4)
        th[i + 1] = th[i] + ksh*amv/(c_of(th[i])*rho)*h*dt
    return t, th


def crossing(t, ratio):
    """First time ratio = resistance/demand falls below 1 (linear interpolation)."""
    idx = np.flatnonzero(ratio < 1.0)
    if not len(idx):
        return math.inf
    i = idx[0]
    if i == 0:
        return 0.0
    f = (ratio[i - 1] - 1.0)/(ratio[i - 1] - ratio[i])
    return float(t[i - 1] + f*(t[i] - t[i - 1]))


def steel_history(eps=0.7, alpha_c=25.0, ksh_f=1.0, dt=DT, Q=Q_IMP):
    g = size_steel(Q)
    t, th = heat(g["AmV"], g["ksh"]*ksh_f, 7850.0, c_steel, eps, alpha_c, dt)
    ky = np.interp(th, ST_T, ST_KY)
    ratio = ky*g["Wpl"]*g["fy"]/0.7/g["Mfi"]
    return g, t, th, ratio


def aluminium_history(eps=0.3, alpha_c=25.0, ksh_f=1.0, dt=DT, Q=Q_IMP):
    g = size_aluminium(Q)
    t, th = heat(g["AmV"], g["ksh"]*ksh_f, 2700.0, c_al, eps, alpha_c, dt)
    k = np.interp(th, AL_T, AL_K)
    ratio = k*g["Wpl"]*g["fo"]/g["Mfi"]
    return g, t, th, ratio


def char_depth(t_s, beta=0.7e-3, d0=0.007):
    tm = t_s/60.0
    return beta*tm + np.minimum(tm/20.0, 1.0)*d0


def glulam_history(beta=0.7e-3, d0=0.007, rho=420.0, dt=DT, Q=Q_IMP, t_end=T_END):
    g = size_glulam(Q, rho=rho)
    t = np.arange(int(round(t_end/dt)) + 1)*dt
    d = char_depth(t, beta, d0)
    bw = np.maximum(g["b"] - 2*d, 0.0); hw = np.maximum(g["h"] - d, 0.0)
    ratio = bw*hw*hw/6*1.15*g["fmk"]/g["Mfi"]
    return g, t, d, ratio


def failure_times(**kw):
    s = {"aluminium": aluminium_history(**{k: v for k, v in kw.items() if k in ("eps_al", "alpha_c", "ksh_f", "dt", "Q")
                                            and k != "eps_al"}, **({"eps": kw["eps_al"]} if "eps_al" in kw else {})),
         "steel": steel_history(**{k: v for k, v in kw.items() if k in ("alpha_c", "ksh_f", "dt", "Q")},
                                **({"eps": kw["eps_st"]} if "eps_st" in kw else {})),
         "glulam": glulam_history(**{k: v for k, v in kw.items() if k in ("beta", "d0", "rho", "dt", "Q")})}
    return {k: crossing(v[1], v[3]) for k, v in s.items()}


SWEEP = {"eps_st": (0.6, 0.8), "eps_al": (0.3, 0.7), "alpha_c": (25.0, 35.0), "ksh_f": (0.9, 1.1),
         "beta": (0.65e-3, 0.75e-3), "d0": (0.007, 0.010), "rho": (380.0, 460.0), "Q": (5.0e3, 7.0e3), "dt": (1.0, 0.5)}


def uncertainty(base):
    contrib = {k: [] for k in ORDER}
    rows = {}
    for name, (lo, hi) in SWEEP.items():
        a, b = failure_times(**{name: lo}), failure_times(**{name: hi})
        rows[name] = {k: (a[k], b[k]) for k in ORDER}
        for k in ORDER:
            contrib[k].append(abs(b[k] - a[k])/2)
    u = {k: math.sqrt(sum(x*x for x in contrib[k])) for k in ORDER}
    return u, rows


# ---------------------------------------------------------------- independent checks

def v_integrator_exact():
    """Explicit Euler (production step) vs the exact solution for constant gas temperature, linear h, constant c."""
    amv, h, c, rho, tg = 250.0, 40.0, 600.0, 7850.0, 800.0
    a = amv*h/(c*rho); th = 20.0
    for _ in range(900):
        th += a*(tg - th)*DT
    exact = tg - (tg - 20.0)*math.exp(-a*900*DT)
    return abs(th - exact)/(exact - 20.0)


def v_steel_eq422(st):
    """The steel failure temperature (k_y = mu0 * kappa1) must match EN eq. 4.22 evaluated at mu0*kappa1."""
    mu = st["mu0"]*0.7
    th_tab = float(np.interp(-mu, [-k for k in ST_KY], ST_T))
    th_eq = 39.19*math.log(1/(0.9674*mu**3.833) - 1) + 482
    return th_tab, th_eq


def v_char_closed_form(gl):
    """Reduced section at 60 min by hand: d_ef = 0.7*60 + 7 = 49 mm."""
    d = 0.049
    W = (gl["b"] - 2*d)*(gl["h"] - d)**2/6
    hand = W*1.15*gl["fmk"]/gl["Mfi"]
    t = np.array([3600.0])
    code = float(((gl["b"] - 2*char_depth(t))*(gl["h"] - char_depth(t))**2/6*1.15*gl["fmk"]/gl["Mfi"])[0])
    return hand, code


def v_dt_halving():
    a, b = failure_times(dt=1.0), failure_times(dt=0.5)
    return {k: abs(a[k] - b[k])/b[k] for k in ORDER}


def run(out_dir):
    out_dir = Path(out_dir)
    (out_dir).mkdir(parents=True, exist_ok=True)
    al, ta, tha, ra = aluminium_history()
    st, ts, ths, rs = steel_history()
    gl, tg_, dch, rg = glulam_history()
    tfail = {"aluminium": crossing(ta, ra), "steel": crossing(ts, rs), "glulam": crossing(tg_, rg)}
    u, rows = uncertainty(tfail)
    order = sorted(ORDER, key=lambda k: tfail[k])
    gaps = []
    for a, b in zip(order, order[1:]):
        ur = max(u[a]/tfail[a], u[b]/tfail[b])
        gap = (tfail[b] - tfail[a])/tfail[a]
        gaps.append(dict(pair=[a, b], gap=gap, threshold=max(0.10, 3*ur), passes=gap > max(0.10, 3*ur)))
    cls = class_checks(st, al)
    integ = v_integrator_exact()
    th_tab, th_eq = v_steel_eq422(st)
    hand, code = v_char_closed_form(gl)
    dth = v_dt_halving()
    th_at_fail_steel = float(np.interp(tfail["steel"], ts, ths))
    validation = {
        "integrator_vs_exact_rel": integ, "integrator_ok": integ < 0.01,
        "steel_theta_table": th_tab, "steel_theta_eq422": th_eq, "steel_eq422_ok": abs(th_tab - th_eq) < 15.0,
        "steel_theta_at_failure": th_at_fail_steel, "steel_failure_temp_consistent": abs(th_at_fail_steel - th_tab) < 2.0,
        "char_hand": hand, "char_code": code, "char_ok": abs(hand - code) < 1e-9,
        "dt_halving_rel": dth, "dt_ok": max(dth.values()) < 0.01,
        "classes": cls, "classes_ok": cls["steel_ok"] and cls["al_ok"],
    }
    validation["all"] = all(v for k, v in validation.items() if k.endswith("_ok"))
    claims = {
        "order_al_steel_glulam": order == ["aluminium", "steel", "glulam"],
        "gaps_exceed_max_10pct_3u": all(g["passes"] for g in gaps),
    }
    res = {
        "editions": EDITIONS, "span_m": L, "loads": dict(G_slab=G_SLAB, Q=Q_IMP, psi1=PSI1),
        "sections": {
            "aluminium": {k: al[k] for k in ("h", "b", "A", "Wpl", "AmV", "ksh", "MEd", "Mfi", "mu0", "self_weight")},
            "steel": {k: st[k] for k in ("h", "b", "A", "Wpl", "AmV", "ksh", "MEd", "Mfi", "mu0", "self_weight")},
            "glulam": {k: gl[k] for k in ("b", "h", "MEd", "Mfi", "mu0", "self_weight")},
        },
        "t_fail_s": tfail, "u_s": u, "sweep": rows, "gaps": gaps, "order": order,
        "validation": validation, "validation_pass": validation["all"], "claims": claims,
    }
    (out_dir/"results.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    step = 5
    np.savez_compressed(out_dir/"frames.npz", time_s=ta[::step], gas_C=np.array([iso834(x) for x in ta[::step]]),
                        aluminium_C=tha[::step], steel_C=ths[::step], char_m=dch[::step],
                        aluminium_ratio=ra[::step], steel_ratio=rs[::step], glulam_ratio=rg[::step])
    return res


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]/"out"
    r = run(target)
    print(json.dumps({"t_fail_min": {k: round(v/60, 2) for k, v in r["t_fail_s"].items()},
                      "u_min": {k: round(v/60, 2) for k, v in r["u_s"].items()},
                      "gaps": [(g["pair"], round(g["gap"], 3), round(g["threshold"], 3), g["passes"]) for g in r["gaps"]],
                      "validation_pass": r["validation_pass"], "claims": r["claims"]}, indent=1))
    print({k: v for k, v in r["validation"].items() if k != "classes"})
