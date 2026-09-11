"""Kelvin inversion pullback for every form degree and orientation.

The sphere Kelvin inversion is the lab's open-boundary map: the exterior of a
ball is solved on a second ball glued to the first by a periodic
identification.  Every quantity crossing that identification has to be
transformed, and the transform depends on TWO independent properties of the
quantity, its form degree and its orientation type.  Getting either wrong
produces a plausible-looking field that is wrong by orders of magnitude.

Distilled from:
  - docs/kelvin/KELVIN_TRANSFORMATION.md 2.2-2.3 (the lab derivation)
  - Bossavit 1998 Ch.2 (straight vs twisted forms in premetric EM)
  - radia.kelvin_source (the implementation these notes describe)
  - measured incidents, 2026-09-08 (see the traps topic)
"""

KELVIN_MAP = """
# The map, and why it is its own inverse

## Two spheres, not one
The production geometry has the PHYSICAL ball of radius R at `c_P` and its
Kelvin image, a second ball of the same radius, at `c_K`.  The mesh glues
their surfaces by a periodic identification.  The map is

    x' = c_K + R^2 (x - c_P) / |x - c_P|^2 ,     rho  = |x  - c_P|
    x  = c_P + R^2 (x' - c_K) / |x' - c_K|^2 ,    rho' = |x' - c_K|

so `rho rho' = R^2`.  A single-centre helper CANNOT express this: with
`c_K - c_P` of 0.66 m on the C-type mesh, using one centre lands the query
in the wrong element, or outside the mesh entirely.

## Jacobian
    J = dx'/dx = (R/rho)^2 H(n) ,      H(n) v = v - 2 (v.n) n ,  n = (x-c_P)/rho

`H` is the Householder reflection about the tangent plane: it flips the radial
component and keeps the tangential ones, `H^2 = I`, `det H = -1`.  Hence

    |det J| = (R/rho)^6 ,     sgn(det J) = -1   (orientation reversing)

## Involution
Applying the map twice is the identity, and the metric factors compose to one
because `R/rho' = rho/R`.  So ONE table serves both directions, provided the
factor `s = R / rho_target` is evaluated in the frame being mapped INTO:
Kelvin -> real uses the physical radius, real -> Kelvin uses the
computational radius.
"""

KELVIN_FORM_TABLE = """
# The complete pullback table

Form degree and orientation type are INDEPENDENT axes.  Write `s = R/rho` with
`rho` the radius in the target frame, and `H` the Householder reflection.

| degree | proxy | straight | twisted |
|--------|-------|----------|---------|
| 0 | scalar f | `f` | `-f` |
| 1 | covector a | `+s^2 H a` | `-s^2 H a` |
| 2 | pseudovector b | `-s^4 H b` | `+s^4 H b` |
| 3 | density u | `-s^6 u` | `+s^6 u` |

## Where each minus sign comes from
The two minus signs have DIFFERENT origins and must not be conflated.

1. Degrees 2 and 3 carry a `sgn(det J)` from converting the form to a vector or
   scalar PROXY: a 2-form proxy transforms with `det(J) J^{-1}`, a 3-form proxy
   with `det(J)`.  This is the Levi-Civita/Hodge identification, not physics.
2. TWISTED (outer-oriented, density-like) forms carry one additional
   `sgn(det J) = -1` because the map reverses orientation.  This is physics:
   it is what distinguishes `H` from `E` and `B` from `D`.

A straight 2-form and a twisted 1-form therefore look superficially similar
(both end up with a bare minus in some places) for entirely unrelated reasons.

## Sanity check that catches a wrong exponent
Energy density is a twisted 3-form, so `U' = +s'^6 U`.  Combined with
`nu' = nu_0 (rho'/R)^2` and the twisted 1-form rule for `H`, the compactified
energy integral reproduces the physical one exactly:

    1/2 int mu' |H'|^2 dV'  =  1/2 int mu_0 |H|^2 dV

because `mu' |H'|^2 = mu_0 s'^2 * s'^4 |H|^2 = mu_0 s'^6 |H|^2` and
`dV = s'^6 dV'`.  Any exponent error breaks this identity immediately, so it
is the cheapest verification available.
"""

KELVIN_PHYSICS_ASSIGNMENT = """
# Which electromagnetic quantity is which form

| quantity | symbol | degree | orientation | transform |
|----------|--------|--------|-------------|-----------|
| electric scalar potential | V | 0 | straight | `V` |
| MAGNETIC scalar potential | phi_m | 0 | **twisted** | `-phi_m` |
| magnetic vector potential | A | 1 | straight | `+s^2 H A` |
| electric field | E | 1 | straight | `+s^2 H E` |
| magnetic field strength | H | 1 | **twisted** | `-s^2 H H` |
| magnetic flux density | B | 2 | straight | `-s^4 H B` |
| electric displacement | D | 2 | **twisted** | `+s^4 H D` |
| current density | J | 2 | **twisted** | `+s^4 H J` |
| charge density | rho_e | 3 | **twisted** | `+s^6 rho_e` |
| magnetic energy density | U_m | 3 | **twisted** | `+s^6 U_m` |

## The one that catches people
The magnetic scalar potential is a 0-form AND twisted.  Its transform carries a
minus but NO metric factor.  Treating it as a plain composition (the natural
guess, since 0-forms usually just compose) silently flips the sign of the whole
exterior field.

## Consequence at the identified spheres
At `rho = rho' = R` the factor is 1 and only `H` acts, so:

  - the NORMAL component of a 1-form flips, the tangential survives;
  - the NORMAL component of a 2-form proxy survives, the tangential flips.

That is why an HCurl unknown (`A`, straight 1-form, tangential trace) glues
straight across a plain periodic identification, while the magnetic scalar
potential does NOT: it is anti-continuous there.  Solve for `-phi_m'` in the
exterior and the unknown becomes continuous again, which is what lets one
periodic H1 space hold the source enclosure and the exterior together.
"""

KELVIN_API = """
# The implementation

```python
from radia.kelvin_source import (
    kelvin_physical_to_computational,   # coordinates, offset aware
    kelvin_computational_to_physical,
    kelvin_transform_form,              # the table, either direction
    kelvin_solution_to_physical,        # Kelvin -> real, guards the domain
    kelvin_solution_to_computational,   # real -> Kelvin, guards the domain
    evaluate_kelvin_exterior,           # evaluate a solved field at far points
    KELVIN_FORMS,                       # name -> (degree, twisted)
)
```

`form` accepts a physics name from `KELVIN_FORMS` (`"scalar_potential"`,
`"vector_potential"`, `"field_strength"`, `"flux_density"`,
`"current_density"`, `"charge_density"`, ...) or an explicit
`(degree, twisted)` pair, so a quantity without a name is still expressible.

Far-field probe of a Kelvin-transformed solve:

```python
B_far = evaluate_kelvin_exterior(
    result["B_cf"], mesh, points_outside_the_ball,
    kelvin_center=kelvin_center, physical_center=physical_center,
    radius=R, form="flux_density")
```

Forward direction, for a source the solver needs INSIDE the Kelvin ball, the
native coefficients do it in C++ with the same conventions:
`radia.KelvinRadiaScalarPotential`, `KelvinRadiaVectorPotential`,
`KelvinRadiaFieldStrength`, `KelvinRadiaFluxDensity`.
"""

KELVIN_TRAPS = """
# Traps, all of them measured rather than imagined (2026-09-08)

## 1. Evaluating a physical coefficient at a Kelvin coordinate
A Radia field coefficient is a PHYSICAL-coordinate object.  Calling it at a
point inside the Kelvin ball returns the field of a completely different
location.  Measured on the C-type coarse mesh at half the Kelvin radius: the
plain coefficient gives 4e-6 T where the pulled-back source is 3.3e-4 T, two
orders of magnitude, with nothing to indicate a problem.  Either supply the
pulled-back source or make that region NaN; never let it return a number.

## 2. A single-centre map on a two-sphere geometry
`kelvin_map_3d(points, center, R)` inverts about ONE centre.  The production
mesh has the physical ball and its image at different centres, so the
single-centre helpers are only valid on concentric research meshes.

## 3. Forgetting that phi_m is twisted
See the physics table.  The magnetic scalar potential's minus sign is the
difference between a correct exterior and a sign-flipped one.  A test that only
checks magnitudes will not catch it; compare against the analytic source field
with its direction.

## 4. Assuming the exterior is coupled because the mesh has the material
A Kelvin block can be present, meshed, and completely inert.  The decisive test
is a PERTURBATION: scale the exterior permeability and see whether the physical
solution moves.  A live exterior responds at the 1e-3 level; a dead one at
1e-15.  Integrating the exterior energy is NOT sufficient evidence, and neither
is a nonzero potential there.

## 5. Trusting an energy that contains the source self energy
The vacuum self energy of a coil has an integrand with a kink on the conductor
surface, which cuts arbitrarily through mesh elements.  On the C-type coarse
mesh that integral swung between 1.90 and 2.75 J as the quadrature order went
8, 10, 12, 14, 18.  Subtract the same constant from both formulations
POINTWISE, so the large common part cancels inside the quadrature rather than
between two large numbers.
"""

_TOPICS = {
    "map": KELVIN_MAP,
    "table": KELVIN_FORM_TABLE,
    "physics": KELVIN_PHYSICS_ASSIGNMENT,
    "api": KELVIN_API,
    "traps": KELVIN_TRAPS,
}


def get_kelvin_pullback_documentation(topic: str = "all") -> str:
    """Return the Kelvin pullback knowledge for ``topic``."""
    key = (topic or "all").strip().lower()
    if key in ("all", ""):
        return "\n\n".join(_TOPICS[name] for name in
                           ("map", "table", "physics", "api", "traps"))
    if key not in _TOPICS:
        raise ValueError(
            "unknown topic %r; choose from %s or 'all'"
            % (topic, sorted(_TOPICS)))
    return _TOPICS[key]
