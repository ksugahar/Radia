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

through homogeneous degree four.  The covariant longitudinal vector
potential `a_s` contains source-free normal and skew multipoles through
octupole.  `H_ref` supplies the finite-reference-speed slip term.  This gives
symmetric `H2/H3/H4` tensors and the Hamiltonian vector-field jet
`A=J H2`, `F2=J H3`, `F3=J H4`.

The runtime value kernel is C++ and is exposed by
`radia.beam.canonical_body_hamiltonian_jet` and the standalone-MEX MATLAB
entry `radia.beam.canonicalBodyHamiltonianJet`.  The topology module carries
analytic parameter tangents for all seven multipole components and checks its
values against that native kernel.  The tracked Wolfram Language script
`lie_map_symbolic_oracle.wls` independently expands the exact parent
Hamiltonian and writes the committed symbolic coefficients in
`lie_map_symbolic_reference.json`; Wolfram Language is not needed at runtime.

## RK reference orbit and Frenet--Serret frame

The reference path is acquired before the Lie expansion.  For a realized
three-dimensional field,
`radia.ffag_topopt.recover_periodic_planar_closed_orbit` integrates position
and unit tangent with DOP853 and solves the one-cell periodic closure.  Its
station values form a `PlanarDesignOrbit`.  Since the present FFAG contract is
planar, torsion is zero and the unit tangent together with the fixed bend axis
defines the straightened Frenet--Serret transverse frame.

`fourth_order_lie_map_from_tracked_orbit` is the explicit bridge from that RK
result to optics.  At every segment midpoint it constructs the moving frame,
samples the live `GridFunction` on a transverse ring, fits source-free
normal/skew multipoles through octupole, and passes those coefficients to the
canonical Lie engine.  Thus the RK pass obtains the reference trajectory and
frame; the later RK variational pass propagates deviations in that frame.  A
caller-supplied fixed `PlanarDesignOrbit` remains supported, but is not
mislabelled as an automatically recovered closed orbit.

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
HDiv-MMM magnetic field
  -> normal/skew multipoles through octupole
  -> native H2/H3/H4 and A/F2/F3
  -> forward-AD R/T/U/V and f3/f4/f5
  -> target-minus-realized map residual
  -> TSVD reachable/unreachable split
  -> ACA + thin QR + TSVD material inverse
  -> whole-HEX binary proposal
  -> complete active-system solve and Lie-map acceptance gate.
```

Skew quadrupoles control linear `x-y` coupling, normal/skew sextupoles control
second-order cross terms, and normal/skew octupoles plus lower-order cascades
control third-order cross terms.  The local TSVD residual is a certificate for
the current field basis and linearization; it is not a claim that an arbitrary
numerical transfer map is globally realizable by a binary Maxwell design.

## Deliberate boundary

This implementation completes the cascade-generated fourth-order map for
source-free, piecewise-constant body multipoles through octupole (`f3/f4/f5`).
The body Hamiltonian is still truncated at `H4`: the exact square-root
Hamiltonian's direct `H5` kinematic terms and a direct decapole term are not
silently folded into `V`.  Consequently this path has fifth-order trajectory
error against the declared `H2/H3/H4` Hamiltonian, while comparison with the
untruncated square root exposes the omitted `H5` term at fourth coordinate
degree.  Longitudinal fringe/edge vector potentials and arbitrary-order normal
forms remain separate work.  Planar closed-orbit recovery itself is available
through the RK bridge above; nonplanar torsion needs a richer frame contract.

The Hamiltonian and factorization conventions follow the published
Lie-algebraic accelerator-map formulations in Iselin,
[CERN SL-Note-2000-001](https://cds.cern.ch/record/702566/files/sl-note-2000-001.pdf),
and Berz and Dragt,
[High-Order Computation and Normal Form Analysis of Repetitive Systems](https://www.bmtdynamics.org/pub/papers/monthnf/monthnf.pdf).
