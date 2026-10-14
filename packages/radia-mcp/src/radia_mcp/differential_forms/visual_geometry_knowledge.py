"""Visual differential geometry as a reasoning layer for Radia.

The geometric progression is distilled from Tristan Needham's *Visual
Differential Geometry and Forms* (2021; Japanese translation, 2026).  The
finite-element and computational-electromagnetics synthesis remains Radia's
own, grounded in the Bossavit and FEEC references listed by this package.
"""

FIVE_ACTS = r"""
# The five-act map: from space to executable calculus

Needham's five acts form a useful engineering decision tree.  Each act adds
one structure; later structures must not be smuggled into earlier ones.

| Act | Added structure | Invariant question | Radia consequence |
|---|---|---|---|
| I. Space | manifold, chart, orientation | What survives a change of coordinates? | Geometry and labels are not component arrays. |
| II. Metric | lengths, angles, areas, volume | What requires a ruler? | Metric and material enter mass/Hodge operators. |
| III. Curvature | local failure of flatness | Is the geometry intrinsically curved? | Curved mesh geometry must be evaluated through mapped elements. |
| IV. Parallel transport | connection and holonomy | How are directions compared at different points? | Rotating frames and periodic sectors need an explicit transport convention. |
| V. Forms | k-forms, wedge, d, integration | What is the natural integrand and boundary law? | H1/HCurl/HDiv/L2 realize one de Rham chain. |

The progression is

    chart -> metric -> curvature -> connection/holonomy -> forms.

Forms then compress the preceding geometry into a coordinate-independent
calculus.  A formula is not considered understood merely because its component
terms cancel.  Radia should be able to name the geometric object, the structure
it depends on, and the invariant that tests it.
"""


INTRINSIC_METRIC = r"""
# Intrinsic geometry, metric, and the Hodge boundary

## Intrinsic versus extrinsic

An intrinsic measurement can be made by an observer living on the manifold:
distance, angle, geodesic curvature, Gaussian curvature, and topology.  An
extrinsic measurement depends on an embedding into a larger space: a chosen
normal, bending in ambient space, or the second fundamental form.

This distinction is operational in CAE.  Two differently embedded meshes may
represent the same intrinsic surface problem, while a shell-bending or shape
optimization problem may deliberately depend on the embedding.  A solver and
its validation must say which problem it is solving.

## The metric is an operator

The metric g converts tangent vectors into covectors and defines length,
angle, area, and volume.  With an orientation it induces the Hodge star

    star_g : Omega^k(M) -> Omega^(n-k)(M).

The exterior derivative d does not require g.  Therefore keep the split

    topology / incidence / d          metric + material / Hodge star
    ------------------------          ------------------------------
    exact combinatorial identities    approximate constitutive operator

under every coordinate map.  For a smooth map phi,

    d(phi^* omega) = phi^*(d omega),

but the Hodge star changes with the pulled-back metric and material.  This is
the geometric reason that Kelvin maps, transformation optics, curved elements,
and anisotropic media modify weights while leaving the de Rham incidence chain
intact.

## Radia rule

Let NGSolve own Piola maps, orientation transforms, quadrature, and mapped
GridFunction evaluation.  Do not reconstruct high-order HCurl/HDiv basis
values from raw component formulas in Python.  A coordinate component is not
the field itself, and a vector proxy for a k-form already contains a metric and
orientation choice.
"""


CURVATURE_HOLONOMY = r"""
# Curvature, holonomy, and topology are three different loop effects

## Curvature from transport

Parallel transport compares tangent vectors at different points.  Transport
around a closed loop generally returns a rotated vector.  On an oriented
surface, with a consistent sign convention,

    holonomy angle = integral_over_region K dA          (mod 2*pi).

For a small loop this is the local curvature test.  Globally, Gauss-Bonnet
adds boundary geodesic curvature and corner turning:

    integral K dA + integral_boundary k_g ds + sum corner_turning
        = 2*pi*chi(M).

For a closed triangulated surface, the sum of vertex angle defects is the
curvature integral.  This gives a mesh QA invariant independent of element
size: it must equal 2*pi times the Euler characteristic.

## Do not conflate the three loop phenomena

1. Levi-Civita holonomy is geometric.  It changes with curvature and metric.
2. A de Rham period integral of a closed non-exact form is topological.  It can
   remain nonzero on a flat domain with a hole.
3. The circulation of the electromagnetic potential A is gauge-theoretic;
   Stokes relates it to magnetic flux where a spanning surface is available.

All three involve loops, but they live in different complexes and answer
different questions.  In particular, a nonzero period is not evidence of
Riemann curvature, and a curved coordinate chart is not evidence that the
underlying space is intrinsically curved.

## Radia uses

- Check surface-mesh orientation and Euler characteristic before using a
  stream-function or surface-current space.
- Compare transported frames, not raw Cartesian components, across rotating
  sectors or moving work coordinates.
- Use cohomology generators for holes in conductors; do not try to repair them
  by changing the material/Hodge matrix.
"""


FORMS_STOKES = r"""
# Forms are native integrands

A k-form is best understood by what it integrates over:

| Degree | Geometric carrier | Typical EM meaning | FE space in 3-D |
|---:|---|---|---|
| 0 | point | scalar potential | H1 |
| 1 | oriented line | E, H, vector potential circulation | HCurl |
| 2 | oriented surface | B, D, current flux | HDiv |
| 3 | oriented volume | charge or source density | L2 |

The wedge product builds oriented area and volume.  The exterior derivative
raises degree by one and obeys

    d(d omega) = 0,
    d(alpha wedge beta)
      = d alpha wedge beta + (-1)^degree(alpha) alpha wedge d beta.

The fundamental theorem of exterior calculus is

    integral_chain d omega = integral_boundary(chain) omega.

It contains the fundamental theorem of calculus, classical Stokes, and the
divergence theorem in one statement.  Its discrete dual is the incidence
identity boundary(boundary(chain)) = 0, hence

    C G = 0,       D C = 0

for discrete grad, curl, and div maps.  These identities are exact algebraic
QA gates, not convergence observations.

Locally, the Poincare lemma says closed forms are exact on contractible
domains.  Globally, de Rham cohomology records the obstruction.  A solver that
assumes every curl-free field is a gradient must first prove that the relevant
cohomology group vanishes or explicitly add cuts/harmonic representatives.
"""


CARTAN = r"""
# Cartan's moving-frame calculus

Choose an orthonormal coframe theta^i and connection 1-forms omega^i_j.
In one common convention Cartan's structure equations are

    T^i       = d theta^i + omega^i_j wedge theta^j,
    Omega^i_j = d omega^i_j + omega^i_k wedge omega^k_j.

For the torsion-free Levi-Civita connection, T^i = 0.  Metric compatibility
makes the connection matrix skew-symmetric in an orthonormal frame.  The
curvature 2-form Omega measures the failure of parallel transport to close;
on a surface its single independent component is K theta^1 wedge theta^2.

## Convention guard

Signs change if the frame matrix, connection, orientation, or index placement
is defined differently.  Needham's moving-frame convention writes the flat
Euclidean Maurer-Cartan relations in the form d[theta]=[omega] wedge [theta]
and d[omega]=[omega] wedge [omega].  Do not copy either sign pattern into code
without recording the convention and verifying a known sphere and a flat
frame.  The invariant content is:

- torsion is the first-structure residual;
- curvature is the second-structure residual;
- a rotating frame may have nonzero connection but zero curvature;
- Bianchi identities are covariant closure statements.

## Radia boundary

Use a connection when comparing vector or tensor components at different
points or times.  Do not add connection coefficients to the exterior
derivative of ordinary scalar-valued forms: d is already coordinate invariant.
For finite elements, mapped form evaluation and Piola transforms carry this
geometry; an extra hand-written Christoffel correction usually double-counts
it.
"""


MAXWELL_GEOMETRY = r"""
# Maxwell as geometry rather than three unrelated vector identities

On a spatial slice, use

    e, h, a in Omega^1,       b, d_field, j in Omega^2,
    rho in Omega^3.

On spacetime, combine electric and magnetic fields into an electromagnetic
2-form F.  With one common sign convention,

    F = b + e wedge dt.

Then the single equation

    dF = 0

contains both magnetic-flux closure d b = 0 and Faraday induction
partial_t b + d e = 0.  The sourced half is written dG = J after choosing the
appropriate twisted excitation form G and current form J.  Constitutive laws,
including material anisotropy and temperature dependence, connect F and G and
carry the metric/Hodge information.

## Vector-proxy guard

In three Euclidean dimensions a metric and orientation identify a 1-form with
a vector and a 2-form with another vector.  This makes dot/cross notation
convenient, but the identification is not free.  Under curvilinear maps or
material transformations, keep the Hodge operation explicit before comparing
components.

## Discrete translation

    Omega^0 --d--> Omega^1 --d--> Omega^2 --d--> Omega^3
       H1           HCurl          HDiv           L2

This is why scalar nodal interpolation is not a valid substitute for an HCurl
electric field or an HDiv magnetic flux.  Polynomial order does not repair a
wrong form degree.  Radia's MQS/Darwin field propagation still uses the
Laplace kernel; this geometric classification does not imply a full-wave
Helmholtz formulation.
"""


RADIA_WORKFLOW = r"""
# Geometry-first workflow for Radia development

1. Name the physical quantity and its integration carrier; assign its form
   degree before choosing an array layout.
2. Write conservation and kinematic laws with d.  Check d squared equals zero
   at the continuous and incidence-matrix levels.
3. Put metric, material, and coordinate-map effects in Hodge/constitutive
   operators.  Check symmetry/positivity when physics requires it.
4. Audit topology.  If the domain has holes, compute periods, Betti numbers,
   cuts, or harmonic representatives explicitly.
5. Declare orientation, pullback/Piola, and Cartan sign conventions at every
   interface between meshes, frames, or rotating sectors.
6. Cross-check a mapped result with d(phi^*omega)=phi^*(domega), not by raw
   component equality.
7. Validate geometry with the executable `differential_forms_geometry_gate`.

## Symptom-to-structure table

| Symptom | Inspect first | Typical mistake |
|---|---|---|
| div(curl) is nonzero | de Rham incidence chain | incompatible spaces or orientation |
| result changes under a coordinate map | pullback and Hodge | comparing components instead of forms |
| a loop mode remains after gauge fixing | cohomology | assuming closed implies exact globally |
| rotating-sector components jump | connection/transport | no frame convention |
| a curved surface has wrong total curvature | orientation and Euler characteristic | flipped/duplicated faces |
| constitutive matrix is nonsymmetric/indefinite | metric/Hodge assembly | bad Jacobian, label, or material tensor |
| interpolated B develops sources | representation | nodal smoothing of an HDiv 2-form |
| A changes after gauge fixing | gauge policy | treating a potential as a physical invariant |
"""


QA_CONTRACT = r"""
# Executable geometry QA contract

`differential_forms_geometry_gate(summary_json)` accepts schema
`radia-visual-geometry-validation/v1` and one of four profiles:

- `de_rham`: nilpotency plus Hodge symmetry/positivity;
- `mapped_em`: de Rham, pullback commutation, Hodge, and dF closure;
- `surface`: Hodge, Cartan connection, Gauss-Bonnet, and holonomy;
- `full`: every check above.

All operator residuals are dimensionless relative residuals.  Surface angles
are in radians.  The gate checks:

    ||C G||, ||D C||,
    ||d phi^* - phi^* d||,
    Hodge symmetry and minimum eigenvalue,
    Cartan metric compatibility, torsion, curvature, and Bianchi residuals,
    integral K + boundary terms - 2*pi*chi,
    wrapped(holonomy - sign*integral K),
    ||dF||.

Missing required sections, non-finite values, undeclared normalization, and
undeclared orientation/Cartan conventions fail loudly.  Numerical violations
return `needs_attention` with individual checks and metrics.

The `connection` object reports
`metric_compatibility_relative_residual`, `torsion_relative_residual`,
`curvature_form_relative_residual`, and `bianchi_relative_residual`.  The
`surface` object reports the curvature integral, boundary geodesic curvature,
corner turning, and integer Euler characteristic.  The `holonomy` object adds
the measured angle, the signed curvature integral, and an orientation sign.
"""


SOURCE_SCOPE = r"""
# Source scope and attribution boundary

Primary geometric source:

Tristan Needham, *Visual Differential Geometry and Forms: A Mathematical
Drama in Five Acts*, Princeton University Press, 2021.  Japanese translation:
Sumio Yamada (supervising translator) and Haruyuki Kawabe, Maruzen Publishing,
2026, ISBN 978-4-621-31240-7.

The five-act progression, visual interpretations of metric/curvature,
holonomy as integrated curvature, forms as integrands, generalized Stokes,
and the moving-frame route to Cartan's structure equations are distilled from
that book.  The FE spaces, Piola rules, discrete gates, NGSolve ownership
boundary, and Radia workflow are a computational-electromagnetics synthesis
from this package's Bossavit, Whitney, and FEEC sources.  Do not attribute
those implementation details to Needham.

No textbook images or extended passages are embedded in this public knowledge
module.  It records derived concepts, equations, and engineering checks.
"""


_TOPICS = {
    "five_acts": FIVE_ACTS,
    "intrinsic_metric": INTRINSIC_METRIC,
    "curvature_holonomy": CURVATURE_HOLONOMY,
    "forms_stokes": FORMS_STOKES,
    "cartan": CARTAN,
    "maxwell": MAXWELL_GEOMETRY,
    "radia_workflow": RADIA_WORKFLOW,
    "qa": QA_CONTRACT,
    "source_scope": SOURCE_SCOPE,
}


def get_visual_geometry_documentation(topic: str = "all") -> str:
    """Return a visual differential-geometry topic for Radia."""

    key = topic.lower().strip().replace("-", "_")
    aliases = {
        "acts": "five_acts",
        "metric": "intrinsic_metric",
        "intrinsic": "intrinsic_metric",
        "curvature": "curvature_holonomy",
        "holonomy": "curvature_holonomy",
        "forms": "forms_stokes",
        "stokes": "forms_stokes",
        "workflow": "radia_workflow",
        "gate": "qa",
        "source": "source_scope",
    }
    key = aliases.get(key, key)
    if key == "all":
        return "\n\n".join(_TOPICS.values())
    if key in _TOPICS:
        return _TOPICS[key]
    return (
        f"Unknown topic '{topic}'. Available: all, "
        + ", ".join(_TOPICS)
        + "."
    )
