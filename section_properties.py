"""Closed-form cross-section property formulas for the four column
candidates, all solved to hit the same target area `A_target` under a
proportioning rule that is fixed *before* any simulation is run (spec
section 5) -- nothing here is tuned after seeing a result.

Every function returns a dict with at least:
    area, I_x, I_y, I_min, c (outer-fibre distance for the governing axis),
    plus the section's own dimensions for rendering/reporting.
"""
from __future__ import annotations

import numpy as np

from utils import circular_section_properties


def fibre_for(section, inertia):
    """The extreme-fibre distance that BELONGS to the inertia being used.

    Bending stress is `M c / I` with both referenced to the same axis. Passing
    an inertia and letting this pick the matching `c` is the only way to use
    these dictionaries that cannot be got wrong, which is why the ambiguous
    `c` key no longer exists.
    """
    if abs(inertia - section["I_x"]) <= 1e-18:
        return section["c_x"]
    if abs(inertia - section["I_y"]) <= 1e-18:
        return section["c_y"]
    raise ValueError("%s: inertia %.6e is neither I_x (%.6e) nor I_y (%.6e), "
                     "so there is no fibre distance that goes with it"
                     % (section.get("kind", "section"), inertia,
                        section["I_x"], section["I_y"]))


def solid_circle_properties(area: float) -> dict:
    """Solid circular rod: A = pi r^2, I = pi r^4 / 4 = A^2/(4 pi)."""
    r, inertia, c = circular_section_properties(area)
    return dict(kind="solid_circle", area=area, r=r, I_x=inertia, I_y=inertia,
               I_min=inertia, c_x=c, c_y=c, outer_dim=2 * r)


def hollow_tube_properties(area: float, D: float) -> dict:
    """Circular hollow tube, fixed outer diameter `D` (a design choice, not
    tuned after the fact -- see report.md).  Solve inner diameter `d` from
    A = pi/4 * (D^2 - d^2) = area, then I = pi/64 * (D^4 - d^4)."""
    d2 = D * D - 4.0 * area / np.pi
    if d2 <= 0:
        raise ValueError(
            f"hollow_tube_properties: area {area} too large for outer "
            f"diameter D={D} (would need d^2={d2} <= 0)")
    d = np.sqrt(d2)
    t = (D - d) / 2.0
    inertia = np.pi / 64.0 * (D ** 4 - d ** 4)
    c = D / 2.0
    return dict(kind="hollow_tube", area=area, D=D, d=d, t=t,
               I_x=inertia, I_y=inertia, I_min=inertia, c_x=c, c_y=c,
               outer_dim=D)


def box_properties(area: float, B: float) -> dict:
    """Square hollow (box) tube, fixed outer width `B`.  Solve wall
    thickness `t` from A = B^2 - (B-2t)^2 = area, then
    I = (B^4 - (B-2t)^4) / 12 (equal about both axes by symmetry)."""
    inner = B * B - area
    if inner <= 0:
        raise ValueError(
            f"box_properties: area {area} too large for outer width B={B}")
    b_inner = np.sqrt(inner)
    t = (B - b_inner) / 2.0
    if t <= 0:
        raise ValueError(f"box_properties: solved wall thickness t={t} <= 0")
    inertia = (B ** 4 - b_inner ** 4) / 12.0
    c = B / 2.0
    return dict(kind="box", area=area, B=B, b_inner=b_inner, t=t,
               I_x=inertia, I_y=inertia, I_min=inertia, c_x=c, c_y=c,
               outer_dim=B)


def i_section_properties(area: float, h: float, b_f: float,
                         tf_tw_ratio: float = 2.0) -> dict:
    """Symmetric I-section.  Fixed proportioning rule (chosen once, before
    any simulation): overall height `h`, flange width `b_f`, and a fixed
    flange/web thickness ratio `tf_tw_ratio` (t_f = tf_tw_ratio * t_w).  The
    one free scale parameter is `t_w`, solved so the resulting area equals
    `area`.

    Area = 2 * b_f * t_f + (h - 2*t_f) * t_w
         = 2 * b_f * tf_tw_ratio * t_w + (h - 2*tf_tw_ratio*t_w) * t_w
    -> quadratic in t_w:  -2*tf_tw_ratio*t_w^2 + (2*b_f*tf_tw_ratio + h)*t_w - area = 0
    """
    a_coef = -2.0 * tf_tw_ratio
    b_coef = 2.0 * b_f * tf_tw_ratio + h
    c_coef = -area
    disc = b_coef ** 2 - 4 * a_coef * c_coef
    if disc < 0:
        raise ValueError("i_section_properties: no real solution for t_w "
                         f"(area={area}, h={h}, b_f={b_f})")
    t_w_candidates = [(-b_coef + np.sqrt(disc)) / (2 * a_coef),
                      (-b_coef - np.sqrt(disc)) / (2 * a_coef)]
    t_w_options = [t for t in t_w_candidates if 0 < t < h / (2 * tf_tw_ratio)]
    if not t_w_options:
        raise ValueError(f"i_section_properties: no valid t_w root in "
                         f"{t_w_candidates}")
    t_w = min(t_w_options)
    t_f = tf_tw_ratio * t_w
    h_web = h - 2 * t_f

    check_area = 2 * b_f * t_f + h_web * t_w
    if abs(check_area - area) > 1e-9 * max(area, 1.0):
        raise ValueError(f"i_section_properties: area check failed "
                         f"({check_area} vs {area})")

    # Strong axis (bending about z, in-plane of the web height h):
    i_web = t_w * h_web ** 3 / 12.0
    i_flange_own = b_f * t_f ** 3 / 12.0
    d_flange = (h_web + t_f) / 2.0  # centroid offset of each flange
    i_flange_total = 2 * (i_flange_own + b_f * t_f * d_flange ** 2)
    I_x = i_web + i_flange_total

    # Weak axis (bending about y, perpendicular to the web):
    i_web_y = h_web * t_w ** 3 / 12.0
    i_flange_y = 2 * (t_f * b_f ** 3 / 12.0)
    I_y = i_web_y + i_flange_y

    I_min = min(I_x, I_y)
    # The extreme fibre is a DIFFERENT distance on each axis, and the
    # bending stress M c / I only means anything when the two come from
    # the same one. Returning a single `c` here paired the weak-axis I_y
    # with the strong-axis h/2 for every caller -- exactly twice the right
    # bending term. Both are returned now and `c` is gone, so the pairing
    # has to be made on purpose.
    c_x = h / 2.0            # bending about the strong axis, I_x
    c_y = b_f / 2.0          # bending about the weak axis, I_y
    return dict(kind="i_section", area=check_area, h=h, b_f=b_f, t_f=t_f,
               t_w=t_w, h_web=h_web, I_x=I_x, I_y=I_y, I_min=I_min,
               c_x=c_x, c_y=c_y, outer_dim=max(h, b_f))
