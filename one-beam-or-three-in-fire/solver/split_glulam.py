"""One glulam beam, two or three: the same wood under the same floor in the ISO 834 standard fire.

The floor is fire-three-beams' office floor (L 5.0 m, slab + finishes 6.0 kN/m plus self-weight, imposed 6.0 kN/m,
simply supported). It is carried by ONE 240 x 320, TWO 120 x 320 or THREE 80 x 320 mm GL24h members -- the same total
width, depth and section modulus, so the same cold design resistance -- heated on three sides (top shielded).

EN 1995-1-2:2004 (historical edition), each method a separate race, never mixed across contenders:
  RCS  4.2.2  effective section, d_ef = beta_n t + k_0 d_0 (k_0 = t/20 below 20 min), k_mod,fi = 1, k_fi 1.15
  RPM  4.2.3  residual section (beta_n t), k_mod,fi = 1 - p / (200 A_r) for bending, linear from t = 0 to 20 min
  3.4.2       the one-dimensional rate beta_0 with corner rounding r = d_char,0 is permitted only while the ORIGINAL
              width b >= b_min (8.15 d_char,0 below 13 mm, else 2 d_char,0 + 80 mm); each member's resistance at each
              time is the ADVERSE (smaller) of the permitted representations (RPM; RCS as a sensitivity only).
Endpoint: the first time the SUMMED sectional bending resistance of the n members falls below M_fi (minutes). This is
the calculated sectional resistance, not floor collapse.

Restraint. The registered race assumes the compression edge is held continuously and the supports torsionally
(EN 1995-1-1 6.3.3(5), k_crit = 1). The UNRESTRAINED case is measured too: k_crit of EN 1995-1-1 6.3.3 on the charred
section at every time step, l_ef = 0.9 L + 2 h (Table 6.1, uniform load, load on the compression edge), torsion held at
the supports only; sigma_m,crit by (6.31) with E_0.05 9.4 GPa (GL24h, EN 1194 class table) and G_0.05 0.54 GPa
(5/6 of G_mean 650 MPa, EN 14080:2013 5.1.3) -- the adverse pair of the two printed tables. In RCS the fire factors
cancel in lambda (4.2.2(5): strength and stiffness with k_mod,fi = 1); in RPM, E and G take 1 - p / (330 A_r).
"""
from __future__ import annotations

import itertools
import math

import numpy as np

L, G_SLAB, Q_NOM = 5.0, 6.0e3, 6.0e3       # fire-three-beams' office floor
B, H = 0.240, 0.320                         # total width, depth (m)
FMK, KFI = 24e6, 1.15                       # GL24h f_m,k; k_fi glulam (EN 1995-1-2 Table 2.1)
KMOD, GAMMA_M = 0.8, 1.25                   # ambient check only (service class 1, medium-term; glulam gamma_M)
E05, G05 = 9.4e9, 0.54e9
LEF = 0.9 * L + 2 * H
EDITION = "EN 1995-1-2:2004"
NS = (1, 2, 3)


# ---------------------------------------------------------------------------------------------------- sections
def b_min(dc):
    return 8.15 * dc if dc < 0.013 else 2 * dc + 0.080


def W_rounded(bw, hw, r, ny=200):
    """section modulus, area of a bw x hw residual section whose two lower corners are rounded with radius r."""
    y = (np.arange(ny) + 0.5) / ny * hw; dy = hw / ny
    w = np.full(ny, bw); m = y < r
    w[m] = bw - 2 * (r - np.sqrt(np.maximum(r * r - (r - y[m]) ** 2, 0)))
    w = np.maximum(w, 0)
    A = (w * dy).sum(); yc = (w * y * dy).sum() / A
    I = (w * (y - yc) ** 2 * dy).sum() + (w * dy ** 3 / 12).sum()
    return I / max(yc, hw - yc), A, w, y


def i_tor(b, h):
    """St Venant torsion constant of a rectangle: (1/3) l s^3 (1 - 0.63 s/l + 0.052 (s/l)^5), s the short side."""
    s, l = min(b, h), max(b, h)
    q = s / l
    return l * s ** 3 / 3 * (1 - 0.63 * q + 0.052 * q ** 5)


def k_crit_of(lam):
    if lam <= 0.75:
        return 1.0
    if lam <= 1.4:
        return 1.56 - 0.75 * lam
    return 1.0 / lam ** 2


def sigma_crit(bz, hz, Wy, lef=None):
    lef = LEF if lef is None else lef
    Iz = hz * bz ** 3 / 12
    return math.pi * math.sqrt(E05 * Iz * G05 * i_tor(bz, hz)) / (lef * Wy)


def ltb(bz, hz, Wy, f_over_E=1.0, lef=None):
    """(k_crit, lambda_rel,m) of a section whose lateral / torsional stiffness is that of the rectangle bz x hz."""
    if bz <= 0 or hz <= 0 or Wy <= 0:
        return 0.0, float("inf")
    lam = math.sqrt(FMK * f_over_E / sigma_crit(bz, hz, Wy, lef))
    return k_crit_of(lam), lam


# ---------------------------------------------------------------------------------------------------- resistance
def resistance(n, tm, method, p, restrained=True, detail=False, b_tot=B, h=H):
    """summed sectional bending resistance (N m) of n members at time tm (min); adverse permitted representation.
    restrained=False multiplies each representation by its own k_crit (EN 1995-1-1 6.3.3)."""
    bn = b_tot / n
    k0 = min(tm / 20, 1)
    out, info = [], []
    dn = p["beta_n"] * tm
    if method == "RCS":
        de = dn + k0 * p["d0"]; bw, hw = bn - 2 * de, h - de
        if bw <= 0 or hw <= 0:
            out.append(0.0); info.append(dict(kc=0.0, lam=float("inf"), rep="notional", bw=0.0, hw=0.0, kmod=1.0))
        else:
            Wy = bw * hw * hw / 6
            kc, lam = (1.0, 0.0) if restrained else ltb(bw, hw, Wy)
            out.append(kc * Wy * KFI * FMK / p["gMfi"])
            info.append(dict(kc=kc, lam=lam, rep="notional", bw=bw, hw=hw, kmod=1.0))
    else:
        bw, hw = bn - 2 * dn, h - dn
        if bw <= 0 or hw <= 0:
            out.append(0.0); info.append(dict(kc=0.0, lam=float("inf"), rep="notional", bw=0.0, hw=0.0, kmod=0.0))
        else:
            d20 = p["beta_n"] * 20; b20, h20 = bn - 2 * d20, h - d20
            pa20 = (b20 + 2 * h20) / (b20 * h20)
            pa = (bw + 2 * hw) / (bw * hw)
            if tm >= 20:
                k, kE = 1 - pa / 200, 1 - pa / 330
            else:
                k, kE = 1 + tm / 20 * (-pa20 / 200), 1 + tm / 20 * (-pa20 / 330)
            Wy = bw * hw * hw / 6
            kc, lam = (1.0, 0.0) if restrained else ltb(bw, hw, Wy, max(k, 1e-9) / max(kE, 1e-9))
            out.append(kc * max(k, 0) * Wy * KFI * FMK / p["gMfi"])
            info.append(dict(kc=kc, lam=lam, rep="notional", bw=bw, hw=hw, kmod=max(k, 0)))
    d0c = p["beta_0"] * tm
    if bn >= b_min(d0c) and (method == "RPM" or p.get("rcs_beta0", False)):
        if method == "RCS":
            de = d0c + k0 * p["d0"]; bw, hw = bn - 2 * de, h - de
            W = W_rounded(bw, hw, d0c)[0] if bw > 0 and hw > 0 else 0.0
            kc, lam = (1.0, 0.0) if restrained else ltb(bw, hw - d0c, W)
            out.append(kc * W * KFI * FMK / p["gMfi"])
            info.append(dict(kc=kc, lam=lam, rep="rounded", bw=bw, hw=hw, kmod=1.0))
        else:
            bw, hw = bn - 2 * d0c, h - d0c
            W, A, w, y = W_rounded(bw, hw, d0c)
            per = bw + 2 * hw - (4 - math.pi) * d0c
            if tm >= 20:
                k, kE = 1 - per / (200 * A), 1 - per / (330 * A)
            else:
                d20 = p["beta_0"] * 20; W20, A20, _, _ = W_rounded(bn - 2 * d20, h - d20, d20)
                per20 = (bn - 2 * d20) + 2 * (h - d20) - (4 - math.pi) * d20
                k, kE = 1 + tm / 20 * (-(per20 / (200 * A20))), 1 + tm / 20 * (-(per20 / (330 * A20)))
            kc, lam = (1.0, 0.0) if restrained else ltb(bw, hw - d0c, W, max(k, 1e-9) / max(kE, 1e-9))
            out.append(kc * max(k, 0) * W * KFI * FMK / p["gMfi"])
            info.append(dict(kc=kc, lam=lam, rep="rounded", bw=bw, hw=hw, kmod=max(k, 0)))
    i = int(np.argmin(out))
    return (n * out[i], info[i]) if detail else n * out[i]


def m_fi(p):
    sw = B * H * p["rho"] * 9.81; G = G_SLAB + sw
    if p["load"] == "comb":
        return (G + p["psi"] * p["Q"]) * L ** 2 / 8
    return 0.6 * (1.35 * G + 1.5 * p["Q"]) * L ** 2 / 8


def m_ed(p):
    sw = B * H * p["rho"] * 9.81
    return (1.35 * (G_SLAB + sw) + 1.5 * p["Q"]) * L ** 2 / 8


def t_fail(n, method, p, restrained=True):
    M = m_fi(p)
    if resistance(n, 0.0, method, p, restrained) < M:
        return 0.0
    lo, hi = 0.0, 200.0
    if resistance(n, hi, method, p, restrained) >= M:
        return hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if resistance(n, mid, method, p, restrained) >= M:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lam_exceeds(n, method, p, lim=0.75, dt=0.05, t_end=200.0):
    """first time the governing representation's lambda_rel,m exceeds lim (k_crit < 1 from then on), by scan."""
    for t in np.arange(0, t_end, dt):
        _, inf = resistance(n, t, method, p, False, detail=True)
        if inf["lam"] > lim:
            return float(t)
    return t_end


# ---------------------------------------------------------------------------------------------------- the sweep
UNC = dict(beta_n=(0.65e-3, 0.7e-3, 0.75e-3), beta_0=(0.60e-3, 0.65e-3, 0.70e-3), d0=(0.007, 0.0085, 0.010),
           rho=(380.0, 420.0, 460.0))
NOM = dict(beta_n=0.7e-3, beta_0=0.65e-3, d0=0.007, rho=420.0)
SCENES = [dict(method=m, load=l, gMfi=g, psi=s, Q=q)
          for m in ("RCS", "RPM") for l in ("comb", "eta") for g in (1.0, 1.1, 1.2)
          for s in ((0.3, 0.5) if l == "comb" else (0.5,)) for q in (5e3, 6e3, 7e3)]
NOMINAL_SCENE = dict(load="comb", gMfi=1.0, psi=0.5, Q=6e3)


def scene_grid(S):
    names = [k for k in UNC if not (S["method"] == "RPM" and k == "d0")]
    return names, [dict(NOM, **dict(zip(names, v))) for v in itertools.product(*[UNC[k] for k in names])]


def sweep(restrained=True, rcs_beta0=False):
    """every scene: nominal times, u (one-at-a-time RSS half-range, max over the uncertain grid), joint adverse corner
    ratio min(t_a) / max(t_b) per pair, need = max(10 %, 3u), pass if ratio - 1 > need."""
    res = []
    for S in SCENES:
        names, grid = scene_grid(S)
        P = lambda u: dict(u, load=S["load"], gMfi=S["gMfi"], psi=S["psi"], Q=S["Q"], rcs_beta0=rcs_beta0)  # noqa
        tf = lambda n, g: t_fail(n, S["method"], P(g), restrained)  # noqa
        T = {n: [tf(n, g) for g in grid] for n in NS}

        def u_at(n, g):
            c = []
            for k in names:
                lo, hi = dict(g), dict(g); lo[k], hi[k] = UNC[k][0], UNC[k][-1]
                c.append(abs(tf(n, hi) - tf(n, lo)) / 2)
            base = tf(n, g)
            return math.sqrt(sum(x * x for x in c)) / base if base > 0 else float("inf")
        um = {n: max(u_at(n, g) for g in grid + [NOM]) for n in NS}
        row = dict(S, nominal={n: tf(n, NOM) for n in NS}, u=um, tmin={n: min(T[n]) for n in NS},
                   tmax={n: max(T[n]) for n in NS}, pairs={},
                   arg_min={n: grid[int(np.argmin(T[n]))] for n in NS},
                   arg_max={n: grid[int(np.argmax(T[n]))] for n in NS})
        for a, b in ((1, 2), (1, 3), (2, 3)):
            r = min(T[a]) / max(T[b]) if max(T[b]) > 0 else float("inf")
            need = max(0.10, 3 * max(um[a], um[b]))
            row["pairs"]["%dv%d" % (a, b)] = [r, need, bool(r - 1 > need)]
        res.append(row)
    return res


def branch_free_bound():
    """Unrestrained claim without the 1.56 - 0.75 lambda and 1 / lambda^2 branches: the single beam is credited only
    until the earlier of its restrained endpoint and the first time its lambda exceeds 0.75 (k_crit = 1 there by the
    code's first range); the thin beams get their RESTRAINED times (an upper bound for any k_crit <= 1)."""
    worst = {"1v2": float("inf"), "1v3": float("inf")}
    lam_max = 0.0
    for S in SCENES:
        _, grid = scene_grid(S)
        P = lambda u: dict(u, load=S["load"], gMfi=S["gMfi"], psi=S["psi"], Q=S["Q"])  # noqa
        t1, t2, t3 = [], [], []
        for g in grid:
            p = P(g)
            tr = t_fail(1, S["method"], p, True)
            tu = t_fail(1, S["method"], p, False)
            lam = resistance(1, tu, S["method"], p, False, detail=True)[1]["lam"]
            lam_max = max(lam_max, lam)
            t1.append(min(tr, lam_exceeds(1, S["method"], p)) if lam > 0.75 else tr)
            t2.append(t_fail(2, S["method"], p, True)); t3.append(t_fail(3, S["method"], p, True))
        worst["1v2"] = min(worst["1v2"], min(t1) / max(t2))
        worst["1v3"] = min(worst["1v3"], min(t1) / max(t3))
    return dict(one_over_two_min=worst["1v2"], one_over_three_min=worst["1v3"], single_lambda_max_at_endpoint=lam_max)


def summary(res):
    out = {}
    for k in ("1v2", "1v3", "2v3"):
        out[k] = dict(passes=sum(r["pairs"][k][2] for r in res), scenes=len(res),
                      min_ratio=min(r["pairs"][k][0] for r in res),
                      min_ratio_gamma1=min(r["pairs"][k][0] for r in res if r["gMfi"] == 1.0))
    out["scenes"] = [dict({k: r[k] for k in ("method", "load", "gMfi", "psi", "Q")},
                          t_min={str(n): r["tmin"][n] for n in NS}, t_max={str(n): r["tmax"][n] for n in NS},
                          arg_min={str(n): r["arg_min"][n] for n in NS}, arg_max={str(n): r["arg_max"][n] for n in NS},
                          ratio={k: r["pairs"][k][0] for k in ("1v2", "1v3", "2v3")},
                          need={k: r["pairs"][k][1] for k in ("1v2", "1v3", "2v3")}) for r in res]
    for m in ("RCS", "RPM"):
        r = [x for x in res if x["method"] == m and all(x[k] == v for k, v in NOMINAL_SCENE.items())][0]
        out["nominal_" + m] = {str(n): r["nominal"][n] for n in NS}
        out["u_" + m] = {str(n): r["u"][n] for n in NS}
    return out


# ---------------------------------------------------------------------------------------------------- outside checks
def benchmarks():
    """published worked examples, recomputed with this module's functions."""
    rows = []
    # JRC 2012 fire workshop report ch. 6 (Frangi) 6.1.4.3: GL24h 160 x 735, 3 sides, 30 min -> W_fi 8664e3 mm3
    p = dict(NOM, gMfi=1.0)
    R = resistance(1, 30.0, "RCS", p, True, b_tot=0.160, h=0.735)
    rows.append(dict(name="JRC 6.1.4.3 W_fi [1e3 mm3], RCS 30 min, 160 x 735", printed=8664.0,
                     got=R / (KFI * FMK) * 1e9 / 1e3, tol_rel=0.001))
    # Kodur CE808 lecture, Example 10.1 (RPM): 130 x 720, 60 min, rate 0.64 -> k_mod,fi 0.805, W 4.12e6 mm3
    pk = dict(beta_n=0.64e-3, beta_0=1.0, d0=0.0, rho=420.0, gMfi=1.0)
    Rk, inf = resistance(1, 60.0, "RPM", pk, True, detail=True, b_tot=0.130, h=0.720)
    rows.append(dict(name="Kodur Ex. 10.1 k_mod,fi (RPM, 60 min)", printed=0.805, got=inf["kmod"], tol_abs=0.0005))
    rows.append(dict(name="Kodur Ex. 10.1 W_r [1e6 mm3]", printed=4.12, got=inf["bw"] * inf["hw"] ** 2 / 6 * 1e9 / 1e6,
                     tol_abs=0.005))
    # torsion constant of a rectangle vs the tabulated beta of Timoshenko & Goodier (h/b 1, 2, 4, 10)
    for hb, beta in ((1.0, 0.141), (2.0, 0.229), (4.0, 0.281), (10.0, 0.312)):
        rows.append(dict(name="I_tor/(h b^3) at h/b %g (Timoshenko & Goodier)" % hb, printed=beta,
                         got=i_tor(1.0, hb) / hb, tol_abs=0.002))
    for r in rows:
        r["ok"] = bool(abs(r["got"] - r["printed"]) <= (r["tol_abs"] if "tol_abs" in r else r["tol_rel"] * r["printed"]))
    return rows


def controls():
    p = dict(NOM, **NOMINAL_SCENE)
    util = {str(n): m_ed(p) / (n * (B / n) * H * H / 6 * KMOD * FMK / GAMMA_M) for n in NS}
    cold = {}
    for n in NS:
        bn = B / n
        kc, lam = ltb(bn, H, bn * H * H / 6)
        cold[str(n)] = dict(lam=lam, k_crit=kc, utilisation_unrestrained=util[str(n)] / kc)
    return dict(cold_utilisation=util, cold_unrestrained=cold, l_ef=LEF)


# ---------------------------------------------------------------------------------------------------- histories
def history(method="RCS", t_end=100.0, dt=0.5):
    """nominal scene, nominal inputs: per contender, the char depth, the zero-strength layer, the residual section,
    the RPM factor and the summed resistance / M_fi, restrained and unrestrained."""
    p = dict(NOM, **NOMINAL_SCENE)
    M = m_fi(p)
    ts = np.arange(0, t_end + 1e-9, dt)
    out = dict(method=method, t=[float(t) for t in ts], m_fi=M, members={})
    for n in NS:
        rows = []
        for t in ts:
            R, inf = resistance(n, t, method, p, True, detail=True)
            Ru, infu = resistance(n, t, method, p, False, detail=True)
            rows.append(dict(char=p["beta_n"] * t, zero=min(t / 20, 1) * p["d0"] if method == "RCS" else 0.0,
                             bw=inf["bw"], hw=inf["hw"], rep=inf["rep"], kmod=inf["kmod"], r=R / M,
                             r_unres=Ru / M, k_crit=infu["kc"], lam=infu["lam"]))
        out["members"][str(n)] = dict(b=B / n, rows=rows, t_end=t_fail(n, method, p, True),
                                     t_end_unres=t_fail(n, method, p, False))
    return out


def study(quick=False):
    R = sweep(True)
    out = dict(edition=EDITION, floor=dict(L=L, G_slab=G_SLAB, Q=Q_NOM), sections_mm={"1": [240, 320], "2": [120, 320],
               "3": [80, 320]}, restrained=summary(R), controls=controls(), benchmarks=benchmarks(),
               histories={m: history(m) for m in ("RCS", "RPM")})
    out["claims"] = dict(one_over_two=out["restrained"]["1v2"]["passes"] == 54,
                         one_over_three=out["restrained"]["1v3"]["passes"] == 54,
                         two_vs_three="unresolved" if out["restrained"]["2v3"]["passes"] < 54 else "resolved")
    if quick:
        return out
    out["rcs_beta0_sensitivity"] = summary(sweep(True, rcs_beta0=True))
    out["unrestrained"] = summary(sweep(False))
    out["unrestrained_rcs_beta0"] = summary(sweep(False, rcs_beta0=True))
    out["branch_free_bound"] = branch_free_bound()
    u = out["unrestrained"]
    out["claims"]["order_holds_unrestrained"] = (u["1v2"]["passes"] == 54 and u["1v3"]["passes"] == 54
                                                 and out["branch_free_bound"]["one_over_two_min"] > 1
                                                 and out["branch_free_bound"]["one_over_three_min"] > 1)
    return out
