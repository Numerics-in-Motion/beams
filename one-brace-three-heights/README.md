# One brace, three heights: where should it hold an I-beam?

A simply supported IPE 300 (S355, fork supports), one midspan point load, and **one ideal (rigid) lateral brace**
at midspan — on the bottom flange, at the shear centre ("middle") or on the top flange — or none. The beam fails by
twisting sideways: lateral-torsional buckling.

| 6 m span | none | bottom | middle | top |
|---|---|---|---|---|
| load **on the top flange**: elastic P_cr | 60.1 kN | 61.3 | 119.0 | 308.9 |
| design P_Rd (EN 1993-1-1, curve a) | 51.1 kN | 51.9 | 86.8 | 126.5 |
| load **hanging** from the bottom flange: P_cr | 111.4 kN | 196.5 | 308.9 | 308.9 |
| design P_Rd | 83.1 kN | 112.0 | 126.5 | 126.5 |

Claims, checked on spans 4, 6 and 8 m (a 10 % bar; the checks were run after the grid had been probed):

1. With the load on top, the top brace beats the middle brace (design +30–48 %), which beats no brace.
2. Hang the load below and, **with an ideal brace**, the middle brace equals the top one — the governing buckled
   shape becomes antisymmetric and keeps the midspan still. This is conditional, not general: with a 100 kN/m brace
   the tie disappears (8 m, hanging: top +58 % vs middle +33 % over no brace).
3. With the load on top, the bottom-flange brace adds at most 10 % (0–8 %).

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

**Method.** Two discretisations of the same elastic energy (lateral bending, St Venant torsion, warping, the moment
coupling, the load-height term and the brace spring on u + e·φ): a Ritz sine series and a thin-walled beam finite
element model. They agree within 3.2 × 10⁻⁵ — a numerical consistency check, not physical validation. Checks: the
uniform-moment closed form (exact / 3 × 10⁻⁸), the C1 factor for a midspan point load at the shear centre (1.361), and
the NCCI SN003 approximate formula for the unbraced beam (0.8–1.1 % at the shear centre, 2.1–5.5 % for top and
hanging loads; SN003 is itself approximate). Design resistance: EN 1993-1-1:2005 6.3.2.2 general case, curve a.

**Limits.** Ideal brace, ideal fork supports, one load case, elastic buckling plus the design curve; not a design for
any real floor.
