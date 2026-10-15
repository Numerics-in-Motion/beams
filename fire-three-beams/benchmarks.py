"""Benchmarks of src/solver/fire_beams.py against published, independent worked examples (agentF, 2026-10-13).

Every check calls the module's OWN code path (heat(), c_steel, c_al, char_depth, steel_history,
aluminium_history, glulam_history, crossing, the ST_/AL_ tables). Where an example's geometry/load differs from
the module's sized beams, the module's size_*() function (and, to isolate the resistance formula at a published
temperature, heat()) is temporarily monkey-patched IN THIS PROCESS ONLY -- solver/fire_beams.py is not edited.

Sources (transcribed by hand; page/sheet numbers as printed in the source):
  [SCI]  SCI, "Fire Resistance design - Example 1: Two-storey office building", job BCF 196, sheets 9-12
         https://www.SteelConstruction.info/images/9/9f/Fire_Design_Example-1.pdf
  [JRC]  Vassart, Zhao, Cajot, Robert, Meyer, Frangi, "Eurocodes: Background & Applications - Structural Fire
         Design - Worked examples", JRC report (workshop 27-28 Nov 2012), ch. 2 (B. Zhao) 2.8.3 pp.103-108 and
         ch. 6 (A. Frangi) 6.1.4 pp.220-223.
         https://eurocodes.jrc.ec.europa.eu/doc/2012_11_WS_fire/report/2012_11_WS_fire.pdf
  [SBE]  Steel Buildings in Europe, "Multi-Storey Steel Buildings Part 6: Fire Engineering", Worked Example
         "Fire safety strategies and design approach of steel floor beam", sheets 5-7 (pp. 6-53 .. 6-55)
         https://www.SteelConstruction.info/images/c/c1/SBE_MS6.pdf
  [UPT]  Politehnica University Timisoara fire-design course, "Examples.pdf" (tube in tension, beam with
         concrete slab on top, and Table 1 'Temperatures on unprotected cross sections under ISO fire')
         https://dan.ct.upt.ro/fire/files/Examples.pdf
  [TALAT2713] EAA TALAT Lecture 2713 "Fire Design Example based on ENV 1999-2 (Eurocode 9)", 1999, sec. 6
         http://aluminium-guide.com/wp-content/uploads/2020/01/2713.pdf
  [HERON] Maljaars, Twilt, Fellinger, Snijder, Soetens, "Aluminium structures exposed to fire conditions - an
         overview", HERON 55(2) 2010, sec. 5.2 design example 2 (pp. 113-114)  https://heronjournal.nl/55-2/1.pdf
  [EN1999] EN 1999-1-2:2007 (online copy), Table 1a and 3.3.1.2
         https://www.phd.eng.br/wp-content/uploads/2014/11/en.1999.1.2.2007.pdf
  [EN1995] EN 1995-1-2:2004 (online copy), 4.2.2 eq. (4.1), Table 4.1 / Fig. 4.2a
         https://www.phd.eng.br/wp-content/uploads/2015/12/en.1995.1.2.2004.pdf
"""
from __future__ import annotations

import contextlib
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).absolute().parent / "solver"))
import fire_beams as fb  # noqa: E402

RESULTS = []


def check(method, ex_id, quantity, published, module, tol, tol_kind="abs", note="", ref="published"):
    diff = module - published
    lim = tol if tol_kind == "abs" else tol*abs(published)
    ok = abs(diff) <= lim
    rel_s = f" ({100*diff/published:+.2f}%)" if published else ""
    RESULTS.append((method, ex_id, ok))
    tol_s = f"+/-{tol:g}" if tol_kind == "abs" else f"+/-{100*tol:g}%"
    print(f"  [{'PASS' if ok else 'FAIL'}] {quantity}: {ref} {published:.4g} | module {module:.4g} | "
          f"diff {diff:+.4g}{rel_s} | tol {tol_s}" + (f"\n         note: {note}" if note else ""))
    return ok


@contextlib.contextmanager
def patched(**attrs):
    old = {k: getattr(fb, k) for k in attrs}
    try:
        for k, v in attrs.items():
            setattr(fb, k, v)
        yield
    finally:
        for k, v in old.items():
            setattr(fb, k, v)


def const_heat(theta):
    """Stand-in for fb.heat that returns a fixed member temperature, to evaluate the module's resistance path
    (table interpolation + resistance formula inside *_history) at a temperature quoted by the publication."""
    return lambda *a, **k: (np.array([0.0]), np.array([float(theta)]))


def t_reach(t, th, target):
    """Time [min] at which the module's temperature history first reaches target (linear interpolation)."""
    i = int(np.flatnonzero(th >= target)[0])
    return float(t[i - 1] + (target - th[i - 1])/(th[i] - th[i - 1])*(t[i] - t[i - 1]))/60.0


def at_time(t, th, t_min):
    return float(np.interp(t_min*60.0, t, th))


# =================================================================================================== STEEL
def steel():
    print("\n=== METHOD 1: STEEL, EN 1993-1-2:2005 =====================================================")

    print("\n[S1] SCI Example 1 (BCF 196) sheets 11-12: 533x210x122 UKB S275, 3-sided, Am/V = 108 m-1, eps_m 0.7,"
          "\n     alpha_c 25, c_a(theta), ISO 834; published time to theta_cr = 669 C")
    for ksh, pub in ((1.0, 18.8), (0.698, 22.7)):
        t, th = fb.heat(108.0, ksh, 7850.0, fb.c_steel, 0.7, 25.0, dt=1.0, t_end=3600.0)
        check("steel", "S1", f"time to 669 C, k_sh={ksh} [min]", pub, t_reach(t, th, 669.0), 0.3,
              note="publication dt = 5 s (EN max), module dt = 1 s")

    print("\n[S2] SCI Example 1 end-to-end through fb.steel_history (W_pl 3196 cm3, f_y 265 MPa -> M_Rd 847 kNm,"
          "\n     M_fi,d 352 kNm, kappa_1 0.7, k_sh 0.698, Am/V 108)")
    g = dict(AmV=108.0, ksh=0.698, Wpl=3196e-6, fy=265e6, Mfi=352e3)
    with patched(size_steel=lambda Q=None, **k: dict(g)):
        _, t, th, ratio = fb.steel_history()
    check("steel", "S2", "M_fi,t,Rd at t=0 = M_Rd/kappa_1 [kNm]", 1210.0, ratio[0]*g["Mfi"]/1e3, 0.005, "rel")
    tf = fb.crossing(t, ratio)/60
    check("steel", "S2", "failure time (resistance < M_fi,d) [min]", 22.7, tf, 0.5,
          note=f"module theta at failure {np.interp(tf*60, t, th):.1f} C (Table 3.1 linear interpolation of "
               f"k_y = mu0 = 0.291) vs 669 C used by SCI (EN eq. 4.22); the ~6 C gap is the eq.4.22-vs-table "
               f"difference, not a heating error (S1 matches)")

    print("\n[S3] JRC worked examples, ch.2 Example 1 (pp.106-108): IPE (H 300) S275 secondary beam, 4-sided,"
          "\n     Am/V 186, k_sh 0.706, eps 0.7, alpha_c 25, mu0 0.308, kappa = 1")
    th_cr = float(np.interp(-0.308, [-k for k in fb.ST_KY], fb.ST_T))
    check("steel", "S3", "theta_cr from Table 3.1 at mu0=0.308 [C]", 667.0, th_cr, 1.0)
    t, th = fb.heat(186.0, 0.706, 7850.0, fb.c_steel, 0.7, 25.0, dt=1.0, t_end=3600.0)
    check("steel", "S3", "time to 667 C, c_a(theta) [min]", 17.0, t_reach(t, th, 667.0), 0.5,
          note="publication reports to 0.5 min resolution ('17 minutes')")
    t, th = fb.heat(186.0, 0.706, 7850.0, lambda x: 600.0, 0.7, 25.0, dt=3.0, t_end=3600.0)
    check("steel", "S3", "time to 667 C, c_a = 600 const, dt 3 s [min]", 16.5, t_reach(t, th, 667.0), 0.5,
          note="publication averages gas temperature over each 3 s step; module uses start-of-step value")

    print("\n[S4] SBE MS6 worked example sheets 5-7: IPE300 S275, 3-sided, Am/V 200, k_sh 0.652, c_a = 600 const,"
          "\n     eps 0.7, dt 5 s; W_pl 628 cm3, kappa_1 kappa_2 = 0.7")
    t, th = fb.heat(200.0, 0.652, 7850.0, lambda x: 600.0, 0.7, 25.0, dt=5.0, t_end=3600.0)
    check("steel", "S4", "theta_a at 60 min [C]", 935.0, float(th[-1]), 15.0,
          note="published value is read off the sheet-6 graph")
    g = dict(AmV=200.0, ksh=0.652, Wpl=628e-6, fy=275e6, Mfi=80.9e3)
    with patched(size_steel=lambda Q=None, **k: dict(g), heat=const_heat(935.0)):
        _, _, _, ratio = fb.steel_history()
    check("steel", "S4", "M_fi,Rd at 935 C [kNm]", 12.83, ratio[0]*g["Mfi"]/1e3, 0.025, "rel",
          note="DISAGREEMENT IN SOURCE: linear interpolation of Table 3.1 at 935 C gives k_y = 0.053; the sheet "
               "prints 0.052 (rounding/slip). Module k_y = %.4f." % np.interp(935.0, fb.ST_T, fb.ST_KY))

    print("\n[S5] UPT Table 1 'Temperatures on unprotected cross sections under ISO fire' ([Am/V]sh, i.e. k_sh"
          "\n     included; table assumptions not printed -- taken as eps 0.7, alpha_c 25, c_a(theta))")
    amv = [10, 15, 20, 25, 30, 40, 60, 100, 200, 400]
    rows = {10: [73, 97, 120, 142, 164, 204, 277, 392, 553, 640],
            20: [159, 219, 273, 323, 369, 448, 561, 676, 734, 754],
            30: [257, 351, 431, 499, 554, 636, 721, 767, 828, 837],
            40: [359, 479, 570, 639, 689, 734, 798, 865, 878, 882],
            48: [439, 570, 660, 717, 737, 792, 875, 901, 907, 910]}
    worst = 0.0
    for a, col in zip(amv, range(len(amv))):
        t, th = fb.heat(float(a), 1.0, 7850.0, fb.c_steel, 0.7, 25.0, dt=1.0, t_end=48*60.0)
        for tm, vals in rows.items():
            d = at_time(t, th, tm) - vals[col]
            worst = max(worst, abs(d))
            if abs(d) > 5.0:
                print(f"     cell Am/V={a}, t={tm} min: table {vals[col]} module {at_time(t, th, tm):.1f}")
    check("steel", "S5", f"max |module - table| over {len(amv)*len(rows)} cells [C]", 0.0, worst, 5.0, ref="target")

    print("\n[S6] UPT 'Fire design of a beam with a concrete slab on top': HE160A S355, W_pl 245.1 cm3,"
          "\n     theta_a = 772 C at R30 -> k_y 0.144 -> mid-span M_fi,t,Rd = k_y f_y W_pl / 0.7")
    g = dict(AmV=108.0, ksh=1.0, Wpl=245.1e-6, fy=355e6, Mfi=6.33e3)
    with patched(size_steel=lambda Q=None, **k: dict(g), heat=const_heat(772.0)):
        _, _, _, ratio = fb.steel_history()
    check("steel", "S6", "k_y at 772 C", 0.144, float(np.interp(772.0, fb.ST_T, fb.ST_KY)), 0.0005)
    check("steel", "S6", "M_fi,t,Rd mid-span [kNm]", 17.9, ratio[0]*g["Mfi"]/1e3, 0.01, "rel")
    check("steel", "S6", "k_y at 828 C (tube example)", 0.096, float(np.interp(828.0, fb.ST_T, fb.ST_KY)), 0.0005)


# =============================================================================================== ALUMINIUM
def rk4_al(amv_ksh, eps, alpha_c=25.0, dt=0.25, t_end=1800.0):
    """INDEPENDENT integrator for EN 1999-1-2 eq. (4.10): classical RK4, own gas curve and own c_al."""
    sig = 5.67e-8
    gas = lambda s: 20.0 + 345.0*math.log10(8.0*s/60.0 + 1.0)
    cal = lambda th: 0.41*th + 903.0
    f = lambda s, th: amv_ksh/(cal(th)*2700.0)*(alpha_c*(gas(s) - th)
                                                 + eps*sig*((gas(s) + 273.0)**4 - (th + 273.0)**4))
    n = int(round(t_end/dt))
    ts = np.arange(n + 1)*dt
    th = np.empty(n + 1); th[0] = 20.0
    for i in range(n):
        s, y = ts[i], th[i]
        k1 = f(s, y); k2 = f(s + dt/2, y + dt/2*k1); k3 = f(s + dt/2, y + dt/2*k2); k4 = f(s + dt, y + dt*k3)
        th[i + 1] = y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    return ts, th


def aluminium():
    print("\n=== METHOD 2: ALUMINIUM, EN 1999-1-2:2007 =================================================")
    print("  NOTE: no published EN 1999-1-2 worked example with fully stated inputs AND numbers for the heating of an"
          "\n  UNPROTECTED aluminium member was found online (searched JRC, TALAT 2503/2704/2713, HERON 2005/2010,"
          "\n  EN 1999-1-2 copies, course notes). A1 is a published resistance example; A2 is a published but"
          "\n  approximate unprotected-heating example with partly inferred inputs; A3-A5 are the fallback checks.")

    print("\n[A1] TALAT 2713 sec. 6 (ENV 1999-1-2 based, 6082-T6 values at 200/250 C identical to EN Table 1a):"
          "\n     Beam F theta 232 C -> k_0.2 = 0.48, M_fi,t,Rd = k M_Rd gamma_M1/gamma_M,fi = 0.48*341*1.10 = 180.0 kNm")
    MRd, gm1 = 341e3, 1.10
    g = dict(AmV=200.0, ksh=1.0, Wpl=MRd*gm1/260e6, fo=260e6, Mfi=139.8e3)
    with patched(size_aluminium=lambda Q=None, **k: dict(g), heat=const_heat(232.0)):
        _, _, _, ratio = fb.aluminium_history()
    check("aluminium", "A1", "k_o,theta at 232 C", 0.48, float(np.interp(232.0, fb.AL_T, fb.AL_K)), 0.005,
          note="publication rounds 0.4772 to 0.48")
    check("aluminium", "A1", "M_fi,t,Rd beam F [kNm]", 180.0, ratio[0]*g["Mfi"]/1e3, 0.01, "rel",
          note="module keeps k unrounded (0.4772); 0.6% gap is the publication's rounding of k")
    check("aluminium", "A1", "k_o,theta at 225 C (column B)", 0.515, float(np.interp(225.0, fb.AL_T, fb.AL_K)),
          0.0005)

    print("\n[A2] HERON 55(2) 2010 sec. 5.2 example 2: unprotected RHS 75x37.5x4 (6060-T66), 3 sides exposed, ISO 834,"
          "\n     EN 1999-1-2 eq.(1): 'reaches the melt temperature (appr. 600 C) after 11 minutes'."
          "\n     INFERRED inputs (not printed in source): Am/V = (2*75+37.5)/(2*4*(75+37.5-8)) = 224 m-1, k_sh 1,"
          "\n     eps 0.3 (EN 2.2(2) clean surface), alpha_c 25")
    amv = (2*75 + 37.5)*1e-3/(2*4*(75 + 37.5 - 2*4)*1e-6)   # 187.5 mm / 836 mm2
    t, th = fb.heat(amv, 1.0, 2700.0, fb.c_al, 0.3, 25.0, dt=1.0, t_end=1800.0)
    check("aluminium", "A2", f"time to 600 C, Am/V={amv:.0f}, eps 0.3 [min]", 11.0, t_reach(t, th, 600.0), 1.0,
          note="approximate published value ('appr.'), inputs partly inferred -> semi-quantitative only")
    t7, th7 = fb.heat(amv, 1.0, 2700.0, fb.c_al, 0.7, 25.0, dt=1.0, t_end=1800.0)
    print(f"     (discrimination: with eps 0.7 the module would give {t_reach(t7, th7, 600.0):.1f} min, i.e. the "
          f"published 11 min is consistent with eps 0.3)")

    print("\n[A3] Second, independent integrator (own RK4, dt 0.25 s, own gas curve and c_al) and scipy solve_ivp,"
          "\n     same EN inputs, on the module's OWN sized aluminium beam (fb.size_aluminium())")
    al = fb.size_aluminium()
    amv_ksh = al["AmV"]*al["ksh"]
    t, th = fb.heat(al["AmV"], al["ksh"], 2700.0, fb.c_al, 0.3, 25.0, dt=1.0, t_end=1800.0)
    tr, thr = rk4_al(amv_ksh, 0.3)
    print(f"     section: Am/V {al['AmV']:.1f} m-1, k_sh {al['ksh']:.3f}, k_sh*Am/V {amv_ksh:.1f} m-1")
    for tm in (5, 10, 15, 20, 30):
        check("aluminium", "A3", f"theta at {tm} min, module heat() vs RK4 [C]", at_time(tr, thr, tm),
              at_time(t, th, tm), 0.5, ref="RK4")
    for target in (250.0, 300.0):
        check("aluminium", "A3", f"time to {target:.0f} C, module vs RK4 [min]", t_reach(tr, thr, target),
              t_reach(t, th, target), 0.002, "rel", ref="RK4")
    try:
        from scipy.integrate import solve_ivp
        sig = 5.67e-8
        gas = lambda s: 20.0 + 345.0*math.log10(8.0*s/60.0 + 1.0)
        rhs = lambda s, y: [amv_ksh/((0.41*y[0] + 903.0)*2700.0)*(25.0*(gas(s) - y[0]) + 0.3*sig*(
            (gas(s) + 273.0)**4 - (y[0] + 273.0)**4))]
        sol = solve_ivp(rhs, (0.0, 1800.0), [20.0], method="LSODA", rtol=1e-10, atol=1e-8,
                        t_eval=[300.0, 600.0, 900.0, 1200.0, 1800.0])
        for s, y in zip(sol.t, sol.y[0]):
            check("aluminium", "A3", f"theta at {s/60:.0f} min, module heat() vs scipy LSODA [C]", float(y),
                  at_time(t, th, s/60), 0.5, ref="LSODA")
    except ImportError:
        print("     scipy not available - LSODA cross-check skipped")

    print("\n[A4] Closed form: module heat() with gas temperature patched to a constant 800 C, eps = 0 (convection"
          "\n     only), c = 903 + 0.41 theta  =>  t(theta) = rho/(alpha k_sh Am/V) [(a + b Tg) ln((Tg-20)/(Tg-theta))"
          "\n     - b (theta - 20)], a = 903, b = 0.41")
    Tg, a, b, alpha = 800.0, 903.0, 0.41, 25.0
    exact = lambda thv: 2700.0/(alpha*amv_ksh)*((a + b*Tg)*math.log((Tg - 20.0)/(Tg - thv)) - b*(thv - 20.0))
    with patched(iso834=lambda s: Tg):
        t, th = fb.heat(al["AmV"], al["ksh"], 2700.0, fb.c_al, 0.0, alpha, dt=1.0, t_end=7200.0)
    for target in (250.0, 300.0, 500.0):
        check("aluminium", "A4", f"time to {target:.0f} C, module vs exact [min]", exact(target)/60,
              t_reach(t, th, target), 0.005, "rel", ref="exact")

    print("\n[A5] Module data vs online copy of EN 1999-1-2:2007 Table 1a, EN AW-6082 T6 (and clause 3.3.1.2)")
    en_T = [20, 100, 150, 200, 250, 300, 350, 550]
    en_k = [1.00, 0.90, 0.79, 0.65, 0.38, 0.20, 0.11, 0.0]
    print("     transcribed row: 1,00 0,90 0,79 0,65 0,38 0,20 0,11 0 at 20/100/150/200/250/300/350/550 C"
          "\n     (the scanned copy's OCR fuses the 250/300 C cells as '0,381;20'; 0,65/0,38 at 200/250 C are"
          "\n     corroborated by TALAT 2713 sec. 6.1)")
    check("aluminium", "A5", "max |AL_K - EN Table 1a|", 0.0,
          max(abs(x - y) for x, y in zip(fb.AL_K, en_k)) + (0.0 if list(fb.AL_T) == en_T else 1.0), 1e-9,
          ref="target")
    check("aluminium", "A5", "c_al(300 C) vs 0.41*300+903 [J/kgK]", 1026.0, fb.c_al(300.0), 1e-9, ref="EN")
    print("     emissivity 0.3 (clean) / 0.7 (painted, sooted): EN 1999-1-2 2.2(2); the online scan renders the"
          "\n     numbers as an image, so this is taken from Maljaars HERON 50(4) 2005 sec. 2.3, which quotes it"
          "\n     (https://heronjournal.nl/50-4/6.pdf) -- not a check of the code text itself.")


# ================================================================================================== GLULAM
def glulam():
    print("\n=== METHOD 3: GLULAM, EN 1995-1-2:2004 reduced cross-section ===============================")

    print("\n[G1] JRC worked examples ch.6 sec. 6.1.4.3 (p.222): main beam GL24h 160x735, 3 sides, t = 30 min,"
          "\n     beta_n 0.7, d0 7, k_fi 1.15: b_fi 104, h_fi 707, W_fi 8664e3 mm3, f_m,d,fi 27.6, M_d,fi 118.1 kNm,"
          "\n     sigma_d,fi 13.6 N/mm2")
    g = dict(b=0.160, h=0.735, fmk=24e6, Mfi=118.1e3)
    with patched(size_glulam=lambda Q=None, fmk=None, kmod=None, gm=None, rho=None: dict(g)):
        _, t, d, ratio = fb.glulam_history(t_end=1800.0)
    MR = ratio[-1]*g["Mfi"]
    check("glulam", "G1", "d_ef = char_depth(30 min) [mm]", 28.0, 1e3*float(fb.char_depth(np.array([1800.0]))[0]),
          1e-9)
    check("glulam", "G1", "M_fi,Rd = f_m,d,fi W_fi = 27.6*8664e3 [kNm]", 27.6*8664e3/1e6, MR/1e3, 0.001, "rel")
    W = MR/(1.15*g["fmk"])
    check("glulam", "G1", "sigma_d,fi = M_d,fi/W_fi [N/mm2]", 13.6, g["Mfi"]/W/1e6, 0.05)
    check("glulam", "G1", "utilisation sigma/f = 13.6/27.6", 13.6/27.6, 1/ratio[-1], 0.005)

    print("\n[G2] JRC sec. 6.1.4.2 (pp.221-222): secondary beam C24 solid 120x260, 3 sides, 30 min, beta_n 0.8:"
          "\n     b_fi 58, h_fi 229, W_fi 506.9e3 mm3  (geometry only: the module's k_fi = 1.15 is glulam-specific,"
          "\n     so W is recovered from the module ratio by dividing out 1.15 f_m,k)")
    g = dict(b=0.120, h=0.260, fmk=24e6, Mfi=7.0e3)
    with patched(size_glulam=lambda Q=None, fmk=None, kmod=None, gm=None, rho=None: dict(g)):
        _, t, d, ratio = fb.glulam_history(beta=0.8e-3, t_end=1800.0)
    check("glulam", "G2", "W_fi [1e3 mm3]", 506.9, ratio[-1]*g["Mfi"]/(1.15*g["fmk"])*1e6, 0.001, "rel")

    print("\n[G3] k0 ramp (t < 20 min) -- CLAUSE check only, no published numeric worked example with t < 20 min"
          "\n     was found: EN 1995-1-2 eq.(4.1) d_ef = beta_n t + k0 d0, d0 = 7 mm, Table 4.1 k0 = t/20 (t<20),"
          "\n     1.0 (t>=20) (also stated in JRC report Fig. 6.1.3, p.220)")
    for tm, pub in ((10, 0.7*10 + 7*10/20), (15, 0.7*15 + 7*15/20), (20, 0.7*20 + 7), (60, 0.7*60 + 7)):
        check("glulam", "G3", f"d_ef at {tm} min [mm]", pub, 1e3*float(fb.char_depth(np.array([tm*60.0]))[0]),
              1e-9, ref="EN clause")


def main():
    print("Benchmarks of", Path(fb.__file__).absolute())
    print("Editions in module:", fb.EDITIONS)
    steel(); aluminium(); glulam()
    print("\n=== SUMMARY ================================================================================")
    for m in ("steel", "aluminium", "glulam"):
        rs = [r for r in RESULTS if r[0] == m]
        exs = sorted({r[1] for r in rs})
        print(f"  {m:10s} checks {len(rs):3d}  pass {sum(r[2] for r in rs):3d}  fail {sum(not r[2] for r in rs):3d}"
              f"   examples: {', '.join(exs)}")
    nfail = sum(not r[2] for r in RESULTS)
    print(f"  OVERALL: {'PASS' if nfail == 0 else 'FAIL'} ({len(RESULTS) - nfail}/{len(RESULTS)})")
    return nfail


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
