"""Spinoff 065 -- one lateral brace, three heights, one steel beam: lateral-torsional buckling.

A simply supported IPE 300 (fork supports), span 6.0 m, one midspan point load applied ON THE TOP FLANGE
(gravity load from something resting on the beam). One ideal lateral brace at midspan, attached to the bottom
flange, the shear centre or the top flange (plus no brace). How much load before the beam twists sideways?

Two independent elastic methods (they share only the section constants and the loading):
  * RITZ  -- sine series in lateral displacement u and twist phi, N harmonics, the moment-coupling integral
             evaluated numerically, the load-height term and the brace spring added directly;
  * FE    -- thin-walled beam elements, Hermite cubics for (u, u') and (phi, phi') per node: lateral bending
             E I_z, St Venant torsion G I_t and warping E I_w stiffness; geometric stiffness from M(x) u'' phi
             by Gauss quadrature; load height and brace at the midspan node.
Design resistance: EN 1993-1-1:2005 6.3.2.2 general case (rolled I, h/b <= 2 -> curve a, alpha_LT = 0.21),
M_b,Rd = chi_LT W_pl f_y (gamma_M1 = 1.0), applied to each M_cr.

Section IPE 300 (Euronorm 19-57 tables): h 300, b 150, t_w 7.1, t_f 10.7 mm; I_z 603.8 cm4, I_t 20.12 cm4,
I_w 125.9e3 cm6, W_pl,y 628.4 cm3. Steel S355, E 210 GPa, G 81 GPa.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import scipy.linalg as la

ROOT = Path(__file__).resolve().parents[1]
E, G = 210e9, 81e9
FY = 355e6
H, B, TW, TF = 0.300, 0.150, 0.0071, 0.0107
IZ, IT, IW, WPL = 603.8e-8, 20.12e-8, 125.9e-9, 628.4e-6
HS = H - TF                                   # distance between flange centroids
SPAN = 6.0
TOP, CENTRE, BOTTOM = HS / 2, 0.0, -HS / 2
BRACES = {"none": None, "bottom": BOTTOM, "centre": CENTRE, "top": TOP}
RIGID = 1e12                                  # N/m, "ideal" brace


def mcr_uniform(L):
    """Classical closed form, uniform moment, fork supports."""
    return math.pi / L * math.sqrt(E * IZ * G * IT * (1 + math.pi ** 2 * E * IW / (G * IT * L ** 2)))


# ------------------------------------------------------------------ Ritz
def p_ritz(L=SPAN, e=None, a=TOP, kb=RIGID, N=25, nx=4001, uniform=False):
    x = np.linspace(0, L, nx)
    Mx = np.ones_like(x) if uniform else np.where(x <= L / 2, x / 2, (L - x) / 2)
    n = np.arange(1, N + 1)
    k = n * math.pi / L
    S = np.sin(np.outer(k, x))
    K = np.zeros((2 * N, 2 * N)); Gm = np.zeros((2 * N, 2 * N))
    K[:N, :N] = np.diag(E * IZ * k ** 4 * L / 2)
    K[N:, N:] = np.diag((E * IW * k ** 4 + G * IT * k ** 2) * L / 2)
    C = np.trapezoid((k ** 2)[:, None, None] * S[:, None, :] * S[None, :, :] * Mx, x, axis=2)
    Gm[:N, N:] = C; Gm[N:, :N] = C.T
    s = np.sin(n * math.pi / 2)
    if not uniform:
        Gm[N:, N:] += a * np.outer(s, s)          # load point falls a*phi^2/2: destabilising for a > 0
    if e is not None:
        v = np.concatenate([s, e * s]); K = K + kb * np.outer(v, v)
    w = la.eigvals(Gm, K)
    w = np.real(w[np.abs(np.imag(w)) < 1e-9]); w = w[w > 1e-30]
    return 1.0 / w.max()


# ------------------------------------------------------------------ FE
def _hermite(xi, l):
    N = np.array([1 - 3 * xi ** 2 + 2 * xi ** 3, l * (xi - 2 * xi ** 2 + xi ** 3),
                  3 * xi ** 2 - 2 * xi ** 3, l * (-xi ** 2 + xi ** 3)])
    d1 = np.array([-6 * xi + 6 * xi ** 2, l * (1 - 4 * xi + 3 * xi ** 2), 6 * xi - 6 * xi ** 2,
                   l * (-2 * xi + 3 * xi ** 2)]) / l
    d2 = np.array([-6 + 12 * xi, l * (-4 + 6 * xi), 6 - 12 * xi, l * (-2 + 6 * xi)]) / l ** 2
    return N, d1, d2


def p_fe(L=SPAN, e=None, a=TOP, kb=RIGID, ne=40, uniform=False):
    """Thin-walled beam FE; DOFs per node (u, u', phi, phi')."""
    assert ne % 2 == 0
    l = L / ne
    nd = 4 * (ne + 1)
    K = np.zeros((nd, nd)); Gm = np.zeros((nd, nd))
    gx, gw = np.polynomial.legendre.leggauss(6)
    for el in range(ne):
        x0 = el * l
        iu = [4 * el, 4 * el + 1, 4 * el + 4, 4 * el + 5]
        ip = [4 * el + 2, 4 * el + 3, 4 * el + 6, 4 * el + 7]
        for g_, w_ in zip(gx, gw):
            xi = (g_ + 1) / 2
            x = x0 + xi * l
            N, d1, d2 = _hermite(xi, l)
            jw = w_ * l / 2
            M = 1.0 if uniform else (x / 2 if x <= L / 2 else (L - x) / 2)
            K[np.ix_(iu, iu)] += E * IZ * np.outer(d2, d2) * jw
            K[np.ix_(ip, ip)] += (E * IW * np.outer(d2, d2) + G * IT * np.outer(d1, d1)) * jw
            c = -M * np.outer(d2, N) * jw              # -M u'' phi  (sign matches the Ritz -u'' basis)
            Gm[np.ix_(iu, ip)] += c; Gm[np.ix_(ip, iu)] += c.T
    mid = ne // 2
    if not uniform:
        Gm[4 * mid + 2, 4 * mid + 2] += a
    if e is not None:
        v = np.zeros(nd); v[4 * mid] = 1.0; v[4 * mid + 2] = e
        K += kb * np.outer(v, v)
    fixed = [0, 2, 4 * ne, 4 * ne + 2]               # fork supports: u = phi = 0 at both ends
    free = np.setdiff1d(np.arange(nd), fixed)
    w = la.eigvals(Gm[np.ix_(free, free)], K[np.ix_(free, free)])
    w = np.real(w[np.abs(np.imag(w)) < 1e-9]); w = w[w > 1e-30]
    return 1.0 / w.max()


# ------------------------------------------------------------------ design resistance
def mb_rd(mcr):
    lam = math.sqrt(WPL * FY / mcr)
    al = 0.21
    ph = 0.5 * (1 + al * (lam - 0.2) + lam ** 2)
    chi = min(1.0, 1 / (ph + math.sqrt(max(ph ** 2 - lam ** 2, 0.0))))
    return chi * WPL * FY, chi, lam


def run():
    out = {"section": "IPE 300", "span_m": SPAN, "fy": FY, "E": E, "G": G,
           "IZ": IZ, "IT": IT, "IW": IW, "WPL": WPL, "load_height_m": TOP, "cases": {}}
    for name, e in BRACES.items():
        pr = p_ritz(e=e)
        pf = p_fe(e=e)
        mcr = pf * SPAN / 4
        mb, chi, lam = mb_rd(mcr)
        out["cases"][name] = dict(P_cr_ritz=pr, P_cr_fe=pf, rel=abs(pr - pf) / pf, M_cr=mcr,
                                  M_bRd=mb, chi_LT=chi, lam_LT=lam, P_Rd=4 * mb / SPAN)
    c = out["cases"]
    base = c["none"]
    for name in c:
        c[name]["gain_elastic"] = c[name]["P_cr_fe"] / base["P_cr_fe"]
        c[name]["gain_design"] = c[name]["P_Rd"] / base["P_Rd"]
    # validation
    v = {}
    v["uniform_ritz_vs_closed"] = abs(p_ritz(uniform=True) - mcr_uniform(SPAN)) / mcr_uniform(SPAN)
    v["uniform_fe_vs_closed"] = abs(p_fe(uniform=True) - mcr_uniform(SPAN)) / mcr_uniform(SPAN)
    v["C1_point_load_shear_centre"] = p_fe(a=0.0) * SPAN / 4 / mcr_uniform(SPAN)
    v["fe_mesh"] = {str(n): p_fe(e=TOP, ne=n) for n in (20, 40, 80)}
    v["ritz_N"] = {str(n): p_ritz(e=TOP, N=n) for n in (9, 25, 41)}
    v["ritz_fe_max_rel"] = max(x["rel"] for x in c.values())
    out["validation"] = v
    out["validation_pass"] = (v["uniform_ritz_vs_closed"] < 1e-6 and v["uniform_fe_vs_closed"] < 1e-3
                              and 1.30 < v["C1_point_load_shear_centre"] < 1.40 and v["ritz_fe_max_rel"] < 0.01)
    return out


if __name__ == "__main__":
    o = run()
    for k, v in o["cases"].items():
        print(f"{k:7s} Ritz {v['P_cr_ritz']/1e3:7.1f}  FE {v['P_cr_fe']/1e3:7.1f} kN (rel {v['rel']:.2e})  "
              f"elastic x{v['gain_elastic']:.2f}   P_Rd {v['P_Rd']/1e3:6.1f} kN  design x{v['gain_design']:.2f}  chi {v['chi_LT']:.2f}")
    print(json.dumps(o["validation"], indent=1, default=float))
    print("validation_pass", o["validation_pass"])


# ------------------------------------------------------------------ v2: the load-height grid and mode shapes
SPANS = (4.0, 6.0, 8.0)
LOADS = {"on_top": TOP, "hanging": BOTTOM}


def mode_fe(L=SPAN, e=None, a=TOP, kb=RIGID, ne=40):
    """Critical load and the normalised mode (x, u, phi) from the FE eigenvector (max |u| + |hs/2 phi| = 1)."""
    l = L / ne
    nd = 4 * (ne + 1)
    K = np.zeros((nd, nd)); Gm = np.zeros((nd, nd))
    gx, gw = np.polynomial.legendre.leggauss(6)
    for el in range(ne):
        x0 = el * l
        iu = [4 * el, 4 * el + 1, 4 * el + 4, 4 * el + 5]
        ip = [4 * el + 2, 4 * el + 3, 4 * el + 6, 4 * el + 7]
        for g_, w_ in zip(gx, gw):
            xi = (g_ + 1) / 2; x = x0 + xi * l
            N, d1, d2 = _hermite(xi, l); jw = w_ * l / 2
            M = x / 2 if x <= L / 2 else (L - x) / 2
            K[np.ix_(iu, iu)] += E * IZ * np.outer(d2, d2) * jw
            K[np.ix_(ip, ip)] += (E * IW * np.outer(d2, d2) + G * IT * np.outer(d1, d1)) * jw
            c = -M * np.outer(d2, N) * jw
            Gm[np.ix_(iu, ip)] += c; Gm[np.ix_(ip, iu)] += c.T
    mid = ne // 2
    Gm[4 * mid + 2, 4 * mid + 2] += a
    if e is not None:
        v = np.zeros(nd); v[4 * mid] = 1.0; v[4 * mid + 2] = e; K += kb * np.outer(v, v)
    fixed = [0, 2, 4 * ne, 4 * ne + 2]
    free = np.setdiff1d(np.arange(nd), fixed)
    w, V = la.eig(Gm[np.ix_(free, free)], K[np.ix_(free, free)])
    w = np.real(w); ok = np.flatnonzero((np.abs(np.imag(la.eigvals(Gm[np.ix_(free, free)], K[np.ix_(free, free)]))) < 1e-9) & (w > 1e-30))
    i = ok[np.argmax(w[ok])]
    full = np.zeros(nd); full[free] = np.real(V[:, i])
    x = np.linspace(0, L, ne + 1); u = full[0::4]; ph = full[2::4]
    s = np.max(np.abs(u) + np.abs(HS / 2 * ph))
    if u[np.argmax(np.abs(u))] < 0:
        s = -s
    return 1.0 / w[i], x, u / s, ph / s


def run_grid():
    rec = {"section": "IPE 300", "fy": FY, "E": E, "G": G, "hs": HS, "grid": {}}
    for L in SPANS:
        for lname, a in LOADS.items():
            row = {}
            for bname, e in BRACES.items():
                pr = p_ritz(L=L, e=e, a=a); pf = p_fe(L=L, e=e, a=a)
                mb, chi, lam = mb_rd(pf * L / 4)
                row[bname] = dict(P_cr_ritz=pr, P_cr_fe=pf, rel=abs(pr - pf) / pf, P_Rd=4 * mb / L, chi=chi, lam=lam)
            for b in row:
                row[b]["gain_el"] = row[b]["P_cr_fe"] / row["none"]["P_cr_fe"]
                row[b]["gain_rd"] = row[b]["P_Rd"] / row["none"]["P_Rd"]
            rec["grid"][f"{L:.0f}|{lname}"] = row
    g = rec["grid"]
    c1 = all(min(r["top"][m] / r["centre"][m], r["centre"][m] / r["none"][m]) > 1.10
             for k, r in g.items() if k.endswith("on_top") for m in ("gain_el", "gain_rd"))
    c2 = all(abs(r["top"]["P_cr_fe"] - r["centre"]["P_cr_fe"]) / r["top"]["P_cr_fe"] < 1e-6
             for k, r in g.items() if k.endswith("hanging"))
    c3 = all(max(r["bottom"]["gain_el"], r["bottom"]["gain_rd"]) <= 1.10 for k, r in g.items() if k.endswith("on_top"))
    rec["claims"] = {"1_on_top_top_beats_centre_beats_none": c1, "2_hanging_centre_equals_top": c2,
                     "3_on_top_bottom_adds_le_10pct": c3}
    rec["max_ritz_fe_rel"] = max(r[b]["rel"] for r in g.values() for b in r)
    return rec
