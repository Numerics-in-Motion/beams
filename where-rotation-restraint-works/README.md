# Where one rotational restraint works

A uniform **6.000 m** simply supported steel beam carries the same symmetric
**10.000 kN/m** distributed load in every case. One ideal rotational spring,
with `kappa_theta = k_theta L / EI = 1`, restrains the beam's **absolute local
rotation to a fixed reference** at one of four registered points.

| x/L | maximum downward deflection (mm) | reduction from unrestrained (%) |
|---|---:|---:|
| 1/8 | 1.8148 | 13.96 |
| 1/4 | 1.9196 | 9.00 |
| 3/8 | 2.0478 | 2.92 |
| 1/2 | 2.1094 | 0.00 |

**Of these four registered positions, x/L = 1/8 gives the least maximum
downward deflection. At x/L = 1/2, the ideal rotational restraint is inactive
under symmetric UDL because the unrestrained rotation is zero there, despite
the unrestrained beam's maximum bending moment there.**

This is not a continuous-position optimum. It is an ideal dimensionless
restraint, not a splice, connection-strength calculation, real connection or
installation recommendation. Deflection, not strength.

## Why the midpoint does nothing

A rotational spring reacts to **rotation**, not bending moment. For this
symmetric beam and load, the unrestrained slope is exactly zero at midspan.
Adding a finite rotational spring there therefore creates zero spring moment
and leaves the maximum sag unchanged at 2.1094 mm.

## Reproduce

```
python reproduce.py
python reproduce.py --quick
```

The full run recomputes the four cases, mirror controls, pinned and rigid
limits, 64/128/256 meshes, force-method solution, vertical equilibrium,
spring compatibility and four known-bad injections. It requires all 11 gates.

The finite-element solve uses a 60-digit banded LDLT so the pre-registered
symmetry and equilibrium tolerances are tested without relaxing them. Maximum
deflection is found from every exact within-element polynomial stationary
point, not a plotting grid. `reference/reference_canonical.json` is the frozen
release reference and records LF-normalised source hashes.

## Scope

Linear-elastic Euler--Bernoulli beam, small displacement, quasi-static load.
Verification against beam theory and an independent force-method solution, not
validation against a physical beam.
