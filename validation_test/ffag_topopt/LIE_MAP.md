# Canonical Lie-map topology contract

## Scope and coordinates

The canonical coordinate order is

```text
z = (x, px/p0, y, py/p0, ell, delta).
```

The Poisson signs are `(+1,+1,-1)` for the `(x,px)`, `(y,py)`, and
`(ell,delta)` pairs.  For a piecewise-constant body segment, Radia expands

```text
H = -(1+h*x)*sqrt((1+delta)^2-px^2-py^2) - a_s + H_ref(delta)
```

through homogeneous degree five.  The covariant longitudinal vector
potential `a_s` contains source-free normal and skew multipoles through
decapole.  `H_ref` supplies the finite-reference-speed slip term.  This gives
symmetric `H2/H3/H4/H5` tensors and the Hamiltonian vector-field jet
`A=J H2`, `F2=J H3`, `F3=J H4`, `F4=J H5`.

The runtime value kernel is C++ and is exposed by
`radia.beam.canonical_body_hamiltonian_jet` and the standalone-MEX MATLAB
entry `radia.beam.canonicalBodyHamiltonianJet`.  The topology module carries
analytic parameter tangents for all nine multipole components and checks its
values against that native kernel.  The tracked Wolfram Language script
`lie_map_symbolic_oracle.wls` independently expands the exact parent
Hamiltonian and writes the committed symbolic coefficients in
`lie_map_symbolic_reference.json`; Wolfram Language is not needed at runtime.

## RK reference orbit and Bishop/RMF frame

The reference path is acquired before the Lie expansion.  For a realized
three-dimensional field,
`radia.ffag_topopt.recover_periodic_planar_closed_orbit` integrates position
and unit tangent with DOP853 and solves the one-cell periodic closure.  Its
station values form a `PlanarDesignOrbit`.  Since the present FFAG contract is
planar, torsion is zero and the unit tangent together with the fixed bend axis
defines the straightened transverse frame.  This is geometrically equivalent
to the planar Frenet frame away from zero curvature, without using a
curvature-derived normal as the numerical transport rule.

`fourth_order_lie_map_from_hcurl_transverse` is the production bridge from the
tracked orbit to optics. Its field argument is a real vector-potential
GridFunction in `HCurl(order=p)`. At every segment it constructs the moving
frame and evaluates A itself throughout the transverse volume; `curl(A)` is
not substituted into the canonical Hamiltonian. No HDiv projection or regular
field grid is introduced, and `p` is read from the retained FESpace. Thus the
first RK pass obtains the reference trajectory and frame; the later A-map RK
propagates deviations using the same unexpanded A Hamiltonian. A caller-
supplied fixed `PlanarDesignOrbit` remains supported, but is not mislabelled as
an automatically recovered closed orbit.

`PlanarDesignOrbit` also exposes the RK stations as global
`(X(s),Y(s),Z(s))`, a cubic-Hermite `position_at(s)`, `tangent_at(s)`, the
local `(x,y,s)` frame, and the segment turning curvature `h(s)`.  The public
`track_canonical_hamiltonian_s` path integrates canonical coordinates with
this same `h(s)` using fixed classical RK4 or adaptive DOP853/RK45/RK23.  A
matched uniform-field circle is an exact zero-deviation trajectory in both
fixed and adaptive paths.

The RK verification surface accepts two explicit physical routes.  Canonical
A-RK (`track_hcurl_vector_potential_canonical_s`) evaluates `A_s,A_y` and their
NGSolve element-interior derivatives from the HCurl GridFunction in the exact
unexpanded Hamiltonian.  Cartesian B-RK
(`track_b_coefficient_cartesian_s`) evaluates HDiv-MMM's direct B
CoefficientFunction.  `compare_hcurl_a_map_to_b_coefficient_rk` runs both and
compares their exit coordinates.  HDiv-MMM may generate both A and B as
CoefficientFunctions, but `project_earlytimes_grid_function_maps` establishes
the explicit EarlyTimes boundary by copying A into an HCurl p=5 GridFunction.
The Lie and canonical A-RK routes reject an unprojected A CoefficientFunction.
The independent Cartesian B-RK route may instead evaluate HDiv-MMM's original
B CoefficientFunction directly; an HDiv p=4 GridFunction is an optional
projection check.  This dual input does not weaken the Lie contract: Lie
construction remains HCurl-A-only.

## Measured-field topology calibration

Measured median-plane B is an inverse-design target, not a field-map input.
`MeasuredMedianPlaneFieldTarget` records physical `(s,x,y=0)` samples, selected
local components `(B_x,B_y,B_s)`, and their measurement bands.
`build_measured_median_plane_field_response_matrix` constructs native HDiv-MMM
observation rows at exactly those points, and
`optimize_hdiv_mmm_magnet_to_measured_median_plane` changes whole pole elements
until complete three-dimensional re-solves match the measurements.  The
accepted pole topology is then re-solved to generate the HCurl A-map and the
independent direct B-map used below.

There is deliberately no B-spline or polynomial continuation from the measured
plane to an off-plane B field or vector potential.  Measurement samples that
are not used by the topology objective remain validation observations only.

## Direct HCurl A-map Hamiltonian

`sample_transverse_vector_potential` evaluates a live HCurl GridFunction
throughout the upper/lower moving-frame volume and records `(A_x,A_y,A_s)`
without first discarding A through `curl(A)`.
`canonical_vector_potential_hamiltonian_rhs` implements the exact static parent

```text
H = -(1+h*x) sqrt((1+delta)^2-(px-a_x)^2-(py-a_y)^2)
    - a_s + H_ref(delta).
```

The axial-gauge production contract is `a_x=0`; the RHS still accepts nonzero
`a_x` so the reduction can be checked.  On an open cell, a local scalar gauge
can impose `A_x=0` and `A_s=A_y=0` on the design orbit.  On a complete closed
ring, a single-valued periodic gauge can set `A_s=0` everywhere on the orbit
only when its longitudinal circulation vanishes.  Otherwise each cell uses a
local gauge and the seam gauge transition is a canonical generating function.

`fourth_order_lie_map_from_hcurl_transverse` is the public map boundary.  It
accepts only the constrained HCurl GridFunction, evaluates A itself, certifies
`A_x=0`, the design-orbit gauge, the declared median-plane symmetry, and HCurl
order at least five, then recovers the degree-five Hamiltonian jet internally.
The coefficient-level Hamiltonian and propagation routines are private
implementation kernels; arbitrary `A_y,A_s` polynomial arrays are not a field
map input.  The internal forward-mode calculation returns `R/T/U/V`,
`f3/f4/f5`, and their HCurl-sample response.  It also rejects a homogeneous map
when the Hamiltonian linear term `H1` exceeds the declared reference-orbit
tolerance; gauge conditions alone do not prove that the supplied curve is a
physical design orbit.

The deleted median-plane A reduction is not a diagnostic or fallback path.
General `x-y` coupling always retains the HCurl off-plane information instead
of fabricating it from a two-dimensional slice.

The general production bridge samples the complete HCurl field on a local
`(x,y)` tensor patch.  Its private jet recovery retains every symmetry-allowed
triangular monomial `x**i*y**j`, `1 <= i+j <= 5`, for both `A_y` and `A_s`.
Independent fits on `x<=0` and `x>=0` report
the aperture-scaled jet mismatch across the design-orbit face; they are not
averaged into a fabricated HCurl normal trace.

The scaled least-squares recovery stores its exact linear sample-to-coefficient
Jacobian.  `differentiate_hcurl_transverse_lie_map` chains an arbitrary set of
sampled-A design responses through that Jacobian and the forward Hamiltonian
Lie AD, returning responses of `R/T/U/V` and `f3/f4/f5`.  Its
`objective_gradient` reverse-contracts tensor cotangents into a scalar-objective
gradient over the supplied design modes.  The returned `H1` response is the
reference-orbit equality-constraint Jacobian: fixed-orbit map gradients are
valid on its null space unless the design orbit is recomputed.

## p=5 sampled aperture certificate

`certify_p5_lie_aperture_against_b_map` finds the largest consecutively passing
set of transverse entrance rings for a fourth-order map built from
`HCurl(order=5)`.  Its validation input is the independent magnetic-flux-density
GridFunction obtained from the HDiv-MMM field route.  That B-map is passed
directly to Cartesian Runge--Kutta with
`field_representation="magnetic_flux_density"`; it is never reconstructed from
the A-map used by LIE.

For every requested spatial ring, normalized transverse-momentum ring, and
momentum deviation, the certificate separately gates LIE truncation against
the unexpanded canonical A-map RK, the canonical A-map RK against the direct
B-map Cartesian RK, and the total LIE-versus-B-map discrepancy.  It also gates
exit-plane residual and static-field momentum conservation.  The result is a
deterministic sampled certificate over nested rings, not an interval enclosure
between angular samples.  HDiv-MMM integration convergence and p=5/p=6 field
projection convergence remain independent upstream evidence.

For a supplied HCurl result and one or more entrance states,
`compare_hcurl_lie_map_to_direct_rk` materializes this triangle directly.  Its
three reported arrays are `lie_minus_a_rk`, `a_rk_minus_b_rk`, and
`lie_minus_b_rk`; the first retains all six canonical coordinates, while the
two Cartesian-B comparisons exclude the unavailable arrival coordinate
`ell`.  The B trajectory uses A only for canonical/mechanical momentum
conversion at the entrance and exit, never to generate the Lorentz force.

The numerical frame convention is the Bishop rotation-minimizing frame (RMF),
not a differentiated Frenet normal.  The native field adapter uses the
fourth-order double-reflection construction of Wang, Juttler, Zheng, and Liu,
*ACM Transactions on Graphics* 27(1), 2008,
[doi:10.1145/1330511.1330513](https://doi.org/10.1145/1330511.1330513).
It discretizes the parallel frame introduced by Bishop,
*American Mathematical Monthly* 82(3), 1975,
[doi:10.1080/00029890.1975.11993807](https://doi.org/10.1080/00029890.1975.11993807).
This keeps the transverse axes defined at zero-curvature stations and prevents
frame roll from appearing as a physical normal/skew multipole mix.  For the
current planar orbit, fixed-binormal and Bishop/RMF framing agree up to the
entrance-frame rotation.  For stations covering a complete spatial closed
loop, `periodic_frame=True` computes the one-turn holonomy and distributes its
compensating roll uniformly in chord arc length.  This is the constant-twist
periodic minimal-twist boundary-value construction; it is deliberately
separate from the zero-twist open-path RMF.  A symmetry-reduced cell must first
apply its cell symmetry to the closure edge and must not set this flag as if
the cell endpoints were an ordinary closed loop.  See Farouki and Moon,
*Advances in Computational Mathematics* 44, 2018,
[doi:10.1007/s10444-018-9599-3](https://doi.org/10.1007/s10444-018-9599-3),
and Farouki, Kim, and Moon, *Computer Aided Geometric Design* 76, 2020,
[doi:10.1016/j.cagd.2019.101802](https://doi.org/10.1016/j.cagd.2019.101802).

## Dragt--Finn representation

The high-order map uses the factorial Taylor convention

```text
z_out = R z + T[z,z]/2 + U[z,z,z]/6 + V[z,z,z,z]/24 + O(z^5)
```

and is factorized as

```text
M = R o exp(:f3:) o exp(:f4:) o exp(:f5:) + O(z^5).
```

`dragt_finn_factorize_third_order` removes the cubic-generator self-cascade
before extracting `f4`.  `dragt_finn_factorize_fourth_order` additionally
removes the third repeated `f3` action and the ordered `f3`--`f4` cross term
before extracting `f5`.  It symmetrizes the homogeneous rank-five generator,
reconstructs `V`, and reports formal symplectic residual coefficients through
cubic Jacobian degree.  A target that fails the declared formal symplectic
tolerance is rejected rather than projected silently.

`apply_dragt_finn_map` evaluates the factors at finite amplitude with an
implicit-midpoint Hamiltonian flow.  This application path is symplectic as a
finite map; the stored `R/T/U/V` remains its formal fourth-order truncation.

## Topology-optimization flow

The complete batch path is

```text
HCurl(order=p) vector potential A
  -> NGSolve-native curl(A)
  -> normal/skew multipoles through decapole
  -> native H2/H3/H4/H5 and A/F2/F3/F4
  -> forward-AD R/T/U/V and f3/f4/f5
  -> target-minus-realized map residual
  -> TSVD reachable/unreachable split
  -> ACA + thin QR + TSVD material inverse
  -> whole-HEX binary proposal
  -> complete active-system solve and Lie-map acceptance gate.
```

The material inverse below the field correction may still use the independent
HDiv-MMM response basis. The HCurl GridFunction is the canonical realized-field
boundary presented to EarlyTimes; direct-B GridFunctions exist only as an
explicit compatibility and cross-check mode.

Skew quadrupoles control linear `x-y` coupling, normal/skew sextupoles control
second-order cross terms, normal/skew octupoles control third-order terms, and
normal/skew decapoles plus lower-order cascades control fourth-order cross
terms.  The local TSVD residual is a certificate for
the current field basis and linearization; it is not a claim that an arbitrary
numerical transfer map is globally realizable by a binary Maxwell design.

## Deliberate boundary

This implementation completes the fourth-order map for source-free,
piecewise-constant body multipoles through decapole (`f3/f4/f5`).  Direct
fifth-degree kinematic terms are expanded from the exact square-root parent,
so the stored fourth-order map has fifth-order trajectory error against the
untruncated parent Hamiltonian.  Longitudinal fringe/edge vector potentials
and arbitrary-order normal forms remain separate work.  Planar closed-orbit recovery itself is available
through the RK bridge above; nonplanar torsion needs a richer frame contract.

The Hamiltonian and factorization conventions follow the published
Lie-algebraic accelerator-map formulations in Iselin,
[CERN SL-Note-2000-001](https://cds.cern.ch/record/702566/files/sl-note-2000-001.pdf),
and Berz and Dragt,
[High-Order Computation and Normal Form Analysis of Repetitive Systems](https://www.bmtdynamics.org/pub/papers/monthnf/monthnf.pdf).
