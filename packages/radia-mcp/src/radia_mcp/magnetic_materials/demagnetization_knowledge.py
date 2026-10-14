"""Demagnetization factor (反磁場係数) for ellipsoids + numerical methods.

Source papers (public-safe curated corpus):
- 反磁場係数.pdf
- 円柱状磁心の反磁界係数と磁気力.pdf
- 磁気回路の反磁場係数.pdf (17 MB)
"""

OVERVIEW = r"""
# 反磁場係数 N (demagnetization factor)

## Why it matters

A finite magnetic body in an external field H_ext develops a
demagnetizing field opposing the external field:
```
H_internal = H_ext - N * M
```
where N is the demagnetization tensor.  For isotropic bodies, N is
diagonal in the body's principal axes:
- N_x + N_y + N_z = 1   (sum rule, SI units)
- N_x = N_y = N_z = 1/3   (sphere, isotropic limit)

For a long cylinder: N_axial → 0, N_radial → 1/2.

## Why FE analysts care

When using a measured BH-curve for FE material data, the measurement
was usually done with a (close to N=0) ring-core or strip-tester
geometry.  Using this BH on a finite body in FE without correcting for
the geometric demag gives wrong results.

However, in FE simulation, the demag effect emerges NATURALLY from the
geometry — the geometric N is implicit in the mesh.  So you should:
- USE BH measured on N=0 geometry (ring core, Epstein frame)
- LET FE solve for the geometric demag

PITFALL: if BH was measured on a finite-N specimen (rare but possible),
the data already includes that N → using it in FE double-counts the
demag, giving low-permeability artifact.

## Reference

The original closed-form solutions for ellipsoid demag are:
- Osborn, "Demagnetizing factors of the general ellipsoid", Phys. Rev.
  67(11-12):351, 1945
- Stoner, Phil. Mag. 36:803, 1945

For non-ellipsoid (cube, cylinder, prism), N is NOT well-defined as a
single tensor — it varies point-to-point.  Numerical computation
required.
"""


OSBORN_FORMULAS = r"""
# Osborn 1945 ellipsoid closed-form demag factors

## Sphere (radius R)
N_x = N_y = N_z = 1/3

## Infinite cylinder (axis = z, radius R)
N_x = N_y = 1/2,  N_z = 0

## Thin disk / sheet (axis = z, perpendicular to plane)
N_x = N_y = 0,  N_z = 1

## General prolate ellipsoid (semi-axes a > b = c)
Eccentricity e = sqrt(1 - (b/a)²) ∈ (0, 1)

```
N_a = (1 - e²)/e³ * (1/(2e) * ln((1+e)/(1-e)) - 1)
N_b = N_c = (1 - N_a) / 2
```

In the limit a/b → 1 (sphere): N_a → 1/3.
In the limit a/b → ∞ (needle): N_a → 0, N_b → 1/2.

## General oblate ellipsoid (semi-axes a = b > c)
Eccentricity e = sqrt(1 - (c/a)²) ∈ (0, 1)

```
N_c = (1/e²) * (1 - sqrt(1-e²)/e * arcsin(e))
N_a = N_b = (1 - N_c) / 2
```

In the limit a/c → 1 (sphere): N_c → 1/3.
In the limit a/c → ∞ (thin disk): N_c → 1.

## Triaxial ellipsoid (a > b > c)
Elliptic integrals (Carlson form):
```
N_a, N_b, N_c — see Osborn 1945 Table I-III for tabulated values
```
Pre-computed tables for common ratios are in the lab folder.

## Python helper (numerical, for tabulation)

```python
import scipy.special as sp
import math
def osborn_prolate(a_over_b):
    if abs(a_over_b - 1) < 1e-9:
        return (1/3, 1/3, 1/3)
    e = math.sqrt(1 - 1/a_over_b**2)
    N_a = (1 - e**2) / e**3 * (0.5/e * math.log((1+e)/(1-e)) - 1)
    N_b = (1 - N_a) / 2
    return (N_a, N_b, N_b)

# Sphere check: a_over_b → 1
# Needle:       a_over_b = 10 → N_a ≈ 0.020, N_b ≈ 0.490
# Long-rod:     a_over_b = 100 → N_a ≈ 0.0004, N_b ≈ 0.4998
```

## Reference table (常用値)

| Shape | Aspect ratio | N_axial | N_radial |
|-------|--------------|---------|----------|
| Sphere | 1 | 0.333 | 0.333 |
| Prolate ellipsoid (a/b=2) | 2 | 0.174 | 0.413 |
| Prolate ellipsoid (a/b=5) | 5 | 0.056 | 0.472 |
| Prolate ellipsoid (a/b=10) | 10 | 0.020 | 0.490 |
| Prolate ellipsoid (a/b=100) | 100 | 0.0004 | 0.4998 |
| Long cylinder | ∞ | 0 | 0.500 |
| Oblate (a/c=2) | 2 | 0.527 | 0.236 |
| Thin disk | 0 (c/a → 0) | 1.000 | 0 |

These match Osborn 1945 Table I (within rounding).
"""


NON_ELLIPSOID = r"""
# Non-ellipsoid demag: cylinders, prisms, cubes

For non-ellipsoid bodies, N is NOT a constant tensor — it varies with
position inside the body (because the internal field is not uniform).

Two common approximations:

## "Magnetometric" demag factor (volume-averaged)

```
<N_axial> = (1/V) ∫∫∫ N_axial(r) dV
```
This is the value Most useful for "average" magnetization in
single-domain analysis.

## "Fluxmetric" demag factor (cross-section averaged)

```
N_axial_flux = average of N_axial over the midplane cross-section
```
This is what corresponds to flux-coil measurement (e.g., Epstein frame).

The two differ by ~5-15% for typical aspect ratios.

## Closed-form for short cylinders (Brown 1962, Sato-Ishii 1989)

For a SOLID cylinder of length L, diameter D:
- Magnetometric N at center axis: tabulated in `磁気回路の反磁場係数.pdf`
- Common engineering approximation:
```
N_axial ≈ 1 / (1 + 2 * sqrt(L/D + L²/(2*D²)))   (very rough)
```
Better: use the Sato-Ishii tabulated values from the lab paper.

## Lab papers

- `public-safe curated corpus`
  — cylindrical core demag + corresponding magnetic force
- `public-safe curated corpus` (17 MB)
  — magnetic circuit perspective, includes air gap effect

## When to compute N analytically vs let FE do it

| Scenario | Do this |
|----------|---------|
| Sizing study, no FE yet | Use Osborn ellipsoid approximation |
| Material BH measurement | Use ring core or Epstein → N≈0 → no correction |
| Final design with FE | Let FE compute the field; N is implicit |
| Cross-validate FE | Compare FE-computed center H to N_analytical predicts |

## Radia practical recommendation

For finite PM body design:
1. Use Osborn for INITIAL sizing (Br/(1+N(μ_r-1)) on operating line)
2. Move to FE for FINAL design — let calc_em_force.py / Radia compute
   the actual demag field
3. Check against Osborn limit values as sanity check
"""


def get_demagnetization_knowledge(topic: str = "overview") -> str:
    """Dispatch by topic.

    Topics:
        overview         - Why N matters, sum rule, FE caveats
        osborn           - Osborn 1945 ellipsoid closed-form formulas
        non_ellipsoid    - Cylinder, prism, cube (numerical / approx)
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "why", "intro"):
        return OVERVIEW
    if topic in ("osborn", "ellipsoid", "closed_form", "formulas"):
        return OSBORN_FORMULAS
    if topic in ("non_ellipsoid", "cylinder", "prism", "cube"):
        return NON_ELLIPSOID
    if topic == "all":
        return "\n\n".join([OVERVIEW, OSBORN_FORMULAS, NON_ELLIPSOID])
    return (f"Unknown topic '{topic}'. Available: overview, osborn, "
            "non_ellipsoid, all.")
