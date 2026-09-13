# -*- coding: utf-8 -*-
"""Less live load next door, more sag here.

THE QUESTION

A steel floor girder runs continuously over several columns. Put the same
live load on every bay, then take it off the neighbouring bay and leave it on
this one. Does this bay sag less, because there is less load on the girder,
or more?

WHAT COMES OUT

More. In a continuous beam the load in the next bay bends the girder UP over
the shared support, and that upward bending is holding this bay's sag down.
Remove it and the bay sags further, although nothing was added to it.

For two equal spans on pin/roller supports under a uniform load the textbook
values are exact: every span loaded, 0.00541 wL^4/EI; one span loaded,
0.00915 wL^4/EI. The study reproduces both from beam theory, independently
of the finite elements (`closed_form_two_span`).

THE SINGLE MEASURED QUANTITY

The largest downward deflection anywhere in the target bay. Deflection
inside each element is recovered EXACTLY for an Euler-Bernoulli beam: the
cubic Hermite interpolation of the nodal solution plus the element's own
fixed-end deflection under its uniform load. So the answer does not depend on
where the nodes happen to fall.

WHAT IS SWEPT (registered with the design review before any claim)

  core (rigid, level supports, equal spans)
      spans                2, 3, 4, every target span
      live-load subsets    all 2^n
      outer ends           pinned, fixed
      load                 uniform (primary), equal-resultant midspan point
      mesh                 12, 24, 48, 96 elements per span
  sensitivities (reported, and each changes the wording it bears on)
      permanent load       D/(D+L) = 0, 0.25, 0.50, 0.75, 0.90 on every span
      unequal spans        interior/end span ratio 0.5, 0.75, 1, 1.25, 1.5
      interior support     vertical stiffness kappa_v = k_v L^3/EI
                           0, 10, 30, 100, 300, 1000, 3000, rigid
      outer end rotation   kappa_theta = k_theta L/EI = 0, 0.1, 1, 10, rigid
      settlement           interior supports, +/-0.5 and +/-1 of the
                           all-bays-loaded sag, uniform and alternating

WHAT IS NOT CLAIMED

Linear elastic, small displacement, quasi-static load cases (removing the
load is not a dynamic event). A 2-D beam: no slab, no composite action, no
connection or column flexibility in the core; flexible supports can shrink or
reverse the effect, which is why the core is stated for rigid, level,
equal-span supports. Deflection, not strength: under every bay loaded the
interior support moment governs first yield.

`fem.py` IS NOT MODIFIED. Its element assembly is imported; the supports,
springs, prescribed settlements and load vectors are built here.
"""
from __future__ import annotations

import io
import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import beams_published as BP                              # noqa: E402
import column_geometry as CG                              # noqa: E402
from fem import assemble_global                           # noqa: E402
from utils import load_config                             # noqa: E402

# ---- the freeze block ----------------------------------------------------
SECTION = "i_section"
SPANS = (2, 3, 4)
ENDS = ("pinned", "fixed")
LOADS = ("udl", "point")
MESHES = (12, 24, 48, 96)
PRODUCTION_MESH = 48
REFINED_MESH = 96
SAMPLES = 17                       # per element, on the exact deflection
DEAD_FRACTIONS = (0.0, 0.25, 0.50, 0.75, 0.90)
SPAN_RATIOS = (0.5, 0.75, 1.0, 1.25, 1.5)
KAPPA_V = (0.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, np.inf)
KAPPA_THETA = (0.0, 0.1, 1.0, 10.0, np.inf)
SETTLEMENTS = (-1.0, -0.5, 0.0, 0.5, 1.0)
DISPLAY = dict(n_spans=3, ends="pinned", load="udl", target=0)
LIVE_W = 6.0e3                     # N/m on the displayed girder (elastic with margin)
OUT = os.environ.get("PATTERN_DATA") or os.path.join(ROOT, "data")


def config():
    return load_config(BP.config_path())


def section(cfg, name=SECTION):
    for s in CG.build_all_columns(cfg):
        if s.name == name:
            return s
    raise KeyError(name)


# ---- the model -----------------------------------------------------------
class Beam:
    """One continuous beam: geometry, supports, and its unit load cases."""

    def __init__(self, lengths, per_span, cfg=None, end_kappa=(0.0, 0.0),
                 kappa_v=np.inf, sec=None):
        cfg = cfg or config()
        self.cfg = cfg
        self.E = float(cfg["material"]["youngs_modulus_pa"])
        self.sec = sec or section(cfg)
        self.A, self.I = self.sec.area, self.sec.inertia
        self.EI = self.E * self.I
        self.lengths = [float(x) for x in lengths]
        self.n_spans = len(self.lengths)
        self.per = int(per_span)
        self.L_ref = self.lengths[0]
        xs = [0.0]
        self.support_nodes = [0]
        self.span_of_element = []
        for s, Ls in enumerate(self.lengths):
            for k in range(1, self.per + 1):
                xs.append(xs[-1] + Ls / self.per)
                self.span_of_element.append(s)
            self.support_nodes.append(len(xs) - 1)
        self.x = np.asarray(xs)
        self.nodes = np.column_stack([self.x, np.zeros_like(self.x)])
        self.members = [(i, i + 1) for i in range(len(xs) - 1)]
        K, self.cache = assemble_global(self.nodes, self.members, self.E,
                                        self.A, self.I)
        ndof = 3 * len(xs)
        fixed = {0}                                # one horizontal restraint
        first, last = self.support_nodes[0], self.support_nodes[-1]
        # outer supports: vertically rigid; rotation spring kappa_theta
        for nd, kap in ((first, end_kappa[0]), (last, end_kappa[1])):
            fixed.add(3 * nd + 1)
            if np.isinf(kap):
                fixed.add(3 * nd + 2)
            elif kap > 0:
                K[3 * nd + 2, 3 * nd + 2] += kap * self.EI / self.L_ref
        # interior supports: vertical spring kappa_v (rigid = fixed)
        self.interior = self.support_nodes[1:-1]
        for nd in self.interior:
            if np.isinf(kappa_v):
                fixed.add(3 * nd + 1)
            elif kappa_v > 0:
                K[3 * nd + 1, 3 * nd + 1] += kappa_v * self.EI / self.L_ref ** 3
        self.K = K
        self.fixed = sorted(fixed)
        self.free = [d for d in range(ndof) if d not in fixed]
        self.ndof = ndof
        self._chol = cho_factor(K[np.ix_(self.free, self.free)])
        ij = np.asarray(self.members)
        self._i, self._j = ij[:, 0], ij[:, 1]
        self._le = self.x[self._j] - self.x[self._i]

    # --- solves ---
    def _solve(self, F, prescribed=None):
        u = np.zeros(self.ndof)
        rhs = F[self.free].copy()
        if prescribed:
            for dof, val in prescribed.items():
                u[dof] = val
            fx = self.fixed
            rhs -= self.K[np.ix_(self.free, fx)] @ u[fx]
        u[self.free] = cho_solve(self._chol, rhs)
        return u

    def load_vector(self, span, kind, w=1.0):
        F = np.zeros(self.ndof)
        q = np.zeros(len(self.members))
        if kind == "udl":
            for e, s in enumerate(self.span_of_element):
                if s != span:
                    continue
                i, j = self.members[e]
                le = self.x[j] - self.x[i]
                F[3 * i + 1] -= w * le / 2
                F[3 * i + 2] -= w * le * le / 12
                F[3 * j + 1] -= w * le / 2
                F[3 * j + 2] += w * le * le / 12
                q[e] = w
        elif kind == "point":
            mid = self.support_nodes[span] + self.per // 2
            if self.per % 2:
                raise ValueError("point load needs a node at midspan")
            F[3 * mid + 1] -= w * self.lengths[span]      # equal resultant
        else:
            raise KeyError(kind)
        return F, q

    def field(self, u, q):
        """Exact deflection (positive = up) at SAMPLES points per element."""
        xi = np.linspace(0.0, 1.0, SAMPLES)
        N1 = 1 - 3 * xi ** 2 + 2 * xi ** 3
        N2 = xi - 2 * xi ** 2 + xi ** 3
        N3 = 3 * xi ** 2 - 2 * xi ** 3
        N4 = -xi ** 2 + xi ** 3
        bub = xi ** 2 * (1 - xi) ** 2
        i, j, le = self._i, self._j, self._le[:, None]
        return (N1 * u[3 * i + 1][:, None] + le * N2 * u[3 * i + 2][:, None]
                + N3 * u[3 * j + 1][:, None] + le * N4 * u[3 * j + 2][:, None]
                - np.asarray(q)[:, None] * le ** 4 * bub / (24.0 * self.EI))

    def unit_fields(self, kind):
        """The deflection field of each span loaded alone, per unit w."""
        out = []
        for s in range(self.n_spans):
            F, q = self.load_vector(s, kind)
            out.append(self.field(self._solve(F), q))
        return np.asarray(out)                     # (spans, elements, samples)

    def settlement_field(self, deltas):
        """Deflection from prescribed interior-support movements (rigid only)."""
        pres = {3 * nd + 1: float(d) for nd, d in zip(self.interior, deltas)}
        u = self._solve(np.zeros(self.ndof), prescribed=pres)
        return self.field(u, np.zeros(len(self.members)))

    def span_mask(self, span):
        return np.asarray(self.span_of_element) == span

    def sag(self, field, span):
        """Largest downward deflection in `span` (positive number)."""
        return float(-field[self.span_mask(span)].min())

    def moment_field(self, u, q):
        """Bending moment M = EI v'' from the same exact deflection, at the
        same samples: no end-force sign conventions to carry by hand."""
        xi = np.linspace(0.0, 1.0, SAMPLES)
        d1 = -6 + 12 * xi
        d2 = -4 + 6 * xi
        d3 = 6 - 12 * xi
        d4 = -2 + 6 * xi
        fpp = 2 - 12 * xi + 12 * xi ** 2
        i, j, le = self._i, self._j, self._le[:, None]
        vpp = ((d1 * u[3 * i + 1][:, None] + le * d2 * u[3 * i + 2][:, None]
                + d3 * u[3 * j + 1][:, None] + le * d4 * u[3 * j + 2][:, None]) / le ** 2
               - np.asarray(q)[:, None] * le ** 2 * fpp / (24.0 * self.EI))
        return self.EI * vpp

    def moments(self, u, q):
        """The largest |M| anywhere on the girder."""
        return float(np.abs(self.moment_field(u, q)).max())


def subsets(n):
    return [tuple(s for s in range(n) if mask >> s & 1) for mask in range(2 ** n)]


def alternate(n, target):
    return tuple(s for s in range(n) if (s - target) % 2 == 0)


def enumerate_case(beam, kind, fields=None):
    """Every target span x every live-load subset: sag of the target span."""
    fields = beam.unit_fields(kind) if fields is None else fields
    n = beam.n_spans
    rows = []
    all_s = tuple(range(n))
    for t in range(n):
        sags = {}
        for S in subsets(n):
            f = fields[list(S)].sum(axis=0) if S else np.zeros(fields.shape[1:])
            sags[S] = beam.sag(f, t)
        gov = max(sags, key=sags.get)
        # the all-bays sag of a short bay between long ones can cancel to
        # nothing (the bay lifts); a ratio over it is not a number to report
        scale = max(abs(v) for v in sags.values())
        cancels = sags[all_s] <= 1e-9 * scale

        def over(v):
            return float("inf") if cancels else v / sags[all_s]
        rows.append(dict(target=t, all=sags[all_s], governing=list(gov),
                         governing_sag=sags[gov], ratio=over(sags[gov]),
                         denominator_cancels=bool(cancels),
                         alternate=list(alternate(n, t)),
                         alternate_ratio=over(sags[alternate(n, t)]),
                         alone_ratio=over(sags[(t,)]),
                         end_span=t in (0, n - 1),
                         sags={",".join(map(str, S)) or "-": v
                               for S, v in sags.items()}))
    return rows


def closed_form_two_span(L=1.0, EI=1.0, w=1.0):
    """Two equal spans, pin/roller, UDL: beam theory, no finite elements.

    Left span, with the hogging moment M at the interior support:
    v(x) = w x (L^3 - 2 L x^2 + x^3) / 24EI - M x (L^2 - x^2) / 6EIL,
    M = wL^2/8 with both spans loaded, wL^2/16 with the left span alone.
    """
    x = np.linspace(0.0, L, 200001)
    base = w * x * (L ** 3 - 2 * L * x ** 2 + x ** 3) / (24 * EI)
    out = {}
    for key, M in (("all", w * L * L / 8.0), ("one", w * L * L / 16.0)):
        v = base - M * x * (L * L - x * x) / (6 * EI * L)
        out[key] = float(v.max()) * EI / (w * L ** 4)
    return out


# ---- the study -----------------------------------------------------------
def core(cfg):
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n, ends, kind, per in itertools.product(SPANS, ENDS, LOADS, MESHES):
        kap = (0.0, 0.0) if ends == "pinned" else (np.inf, np.inf)
        b = Beam([L] * n, per, cfg, end_kappa=kap)
        for r in enumerate_case(b, kind):
            r.update(n_spans=n, ends=ends, load=kind, mesh=per,
                     all_coefficient=r["all"] * b.EI / L ** 4)
            rows.append(r)
    return rows


def control(cfg):
    L = float(cfg["geometry"]["span_m"])
    cf = closed_form_two_span()
    out = dict(closed_form=cf, fe={}, moment={})
    for per in MESHES:
        b = Beam([L, L], per, cfg)
        f = b.unit_fields("udl")
        out["fe"][per] = dict(all=b.sag(f.sum(axis=0), 0) * b.EI / L ** 4,
                              one=b.sag(f[0], 0) * b.EI / L ** 4)
        # the moment the elastic check reads, against wL^2/8 at the support
        F = np.zeros(b.ndof)
        q = np.zeros(len(b.members))
        for s in (0, 1):
            Fs, qs = b.load_vector(s, "udl")
            F += Fs
            q += qs
        M = b.moment_field(b._solve(F), q)
        out["moment"][per] = dict(
            support_over_wL2_8=float(-M[b.span_mask(0)][-1, -1] / (L * L / 8.0)),
            sagging_over_9wL2_128=float(M.max() / (9.0 * L * L / 128.0)))
    return out


def dead_fraction(cfg):
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n, ends, kind in itertools.product(SPANS, ENDS, LOADS):
        kap = (0.0, 0.0) if ends == "pinned" else (np.inf, np.inf)
        b = Beam([L] * n, PRODUCTION_MESH, cfg, end_kappa=kap)
        f = b.unit_fields(kind)
        full = f.sum(axis=0)
        for t in range(n):
            base = enumerate_case(b, kind, f)[t]
            gov = f[base["governing"]].sum(axis=0)
            for d in DEAD_FRACTIONS:
                pat = d * full + (1 - d) * gov
                rows.append(dict(n_spans=n, ends=ends, load=kind, target=t,
                                 end_span=t in (0, n - 1), dead_fraction=d,
                                 ratio=b.sag(pat, t) / b.sag(full, t)))
    return rows


def unequal(cfg):
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n, ends, r in itertools.product(SPANS, ENDS, SPAN_RATIOS):
        lengths = [L] + [r * L] * (n - 2) + [L] if n > 2 else [L, r * L]
        kap = (0.0, 0.0) if ends == "pinned" else (np.inf, np.inf)
        b = Beam(lengths, PRODUCTION_MESH, cfg, end_kappa=kap)
        f = b.unit_fields("udl")
        for row in enumerate_case(b, "udl", f):
            alone_single = b.sag(f[row["target"]], row["target"])
            rows.append(dict(n_spans=n, ends=ends, span_ratio=r,
                             target=row["target"], end_span=row["end_span"],
                             ratio=row["ratio"], governing=row["governing"],
                             reverses=row["ratio"] <= 1.0,
                             all_over_alone=row["all"] / alone_single))
    return rows


def support_stiffness(cfg):
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n, kv in itertools.product(SPANS, KAPPA_V):
        b = Beam([L] * n, PRODUCTION_MESH, cfg, kappa_v=kv)
        f = b.unit_fields("udl")
        full = f.sum(axis=0)
        t = 0
        alone = b.sag(f[0], 0)
        alt = b.sag(f[list(alternate(n, 0))].sum(axis=0), 0)
        best = enumerate_case(b, "udl", f)[0]
        rows.append(dict(n_spans=n, kappa_v=float(kv), target=t,
                         alone_ratio=alone / b.sag(full, 0),
                         alternate_ratio=alt / b.sag(full, 0),
                         governing=best["governing"], ratio=best["ratio"]))
    return rows


def end_rotation(cfg):
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n, kt in itertools.product(SPANS, KAPPA_THETA):
        b = Beam([L] * n, PRODUCTION_MESH, cfg, end_kappa=(kt, kt))
        best = enumerate_case(b, "udl")[0]
        rows.append(dict(n_spans=n, kappa_theta=float(kt), target=0,
                         ratio=best["ratio"], governing=best["governing"],
                         alternate_ratio=best["alternate_ratio"]))
    return rows


def settlement(cfg):
    """Total sag and load-induced sag, with interior supports displaced."""
    L = float(cfg["geometry"]["span_m"])
    rows = []
    for n in SPANS:
        b = Beam([L] * n, PRODUCTION_MESH, cfg)
        f = b.unit_fields("udl")
        full = f.sum(axis=0)
        ref = b.sag(full, 0)
        gov = f[list(alternate(n, 0))].sum(axis=0)
        k = len(b.interior)
        for mode, s in itertools.product(("uniform", "alternating"), SETTLEMENTS):
            signs = [1.0] * k if mode == "uniform" else [(-1.0) ** i for i in range(k)]
            deltas = [-s * ref * g for g in signs]          # positive s = down
            fs = b.settlement_field(deltas)
            total_all = b.sag(full + fs, 0)
            total_pat = b.sag(gov + fs, 0)
            rows.append(dict(n_spans=n, mode=mode, settlement_over_sag=s,
                             total_ratio=(total_pat / total_all
                                          if total_all > 0 else float("nan")),
                             incremental_ratio=b.sag(gov, 0) / b.sag(full, 0),
                             total_all_sag=total_all, total_pattern_sag=total_pat))
    return rows


def display_case(cfg):
    """The girder on screen, at its real load: elastic check and sags in mm."""
    L = float(cfg["geometry"]["span_m"])
    d = DISPLAY
    b = Beam([L] * d["n_spans"], PRODUCTION_MESH, cfg)
    fy = float(cfg["material"]["failure_stress_pa"])
    c = b.sec.fibre_c
    worst_stress = 0.0
    for S in subsets(d["n_spans"]):
        F = np.zeros(b.ndof)
        q = np.zeros(len(b.members))
        for s in S:
            Fs, qs = b.load_vector(s, d["load"], LIVE_W)
            F += Fs
            q += qs
        u = b._solve(F)
        worst_stress = max(worst_stress, b.moments(u, q) * c / b.I)
    f = b.unit_fields(d["load"]) * LIVE_W
    full = f.sum(axis=0)
    alt = f[list(alternate(d["n_spans"], d["target"]))].sum(axis=0)
    return dict(n_spans=d["n_spans"], ends=d["ends"], load=d["load"],
                target=d["target"], live_w_n_per_m=LIVE_W, span_m=L,
                section=b.sec.name, sag_all_mm=b.sag(full, 0) * 1e3,
                sag_alternate_mm=b.sag(alt, 0) * 1e3,
                worst_stress_over_yield=worst_stress / fy)


def main():
    cfg = config()
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    doc = dict(frozen=dict(section=SECTION, span_m=float(cfg["geometry"]["span_m"]),
                           spans=SPANS, ends=ENDS, loads=LOADS, meshes=MESHES,
                           production_mesh=PRODUCTION_MESH,
                           refined_mesh=REFINED_MESH, samples=SAMPLES,
                           dead_fractions=DEAD_FRACTIONS, span_ratios=SPAN_RATIOS,
                           kappa_v=[float(k) for k in KAPPA_V],
                           kappa_theta=[float(k) for k in KAPPA_THETA],
                           settlements=SETTLEMENTS, display=DISPLAY,
                           live_w=LIVE_W),
               control=control(cfg), core=core(cfg),
               dead_fraction=dead_fraction(cfg), unequal=unequal(cfg),
               support_stiffness=support_stiffness(cfg),
               end_rotation=end_rotation(cfg), settlement=settlement(cfg),
               display=display_case(cfg))
    doc["runtime_s"] = time.time() - t0

    def clean(o):
        if isinstance(o, float) and np.isinf(o):
            return "inf"
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, np.floating):
            return clean(float(o))
        return o
    with io.open(os.path.join(OUT, "pattern_beam.json"), "w", encoding="utf-8",
                 newline="\n") as fh:
        json.dump(clean(doc), fh, indent=1)
    print("wrote %s  %.1f s" % (os.path.join(OUT, "pattern_beam.json"),
                                doc["runtime_s"]))


if __name__ == "__main__":
    main()
