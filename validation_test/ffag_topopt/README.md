# FFAG HDiv-MMM topology-optimization validation

`validation_ffag_cell_targets.py` builds the Bell--Abell non-scaling FFAG
soft-edge one-cell target family at seven proton energies from 31 to 250 MeV.
It checks the three prerequisites for each reduced closed orbit (cell bend,
periodic position, and periodic tangent), then checks the analytic map
transform and symplectic residual.

This lane establishes the momentum-indexed target contract.  The first magnet
PoC treats each supplied entrance-to-exit design orbit as fixed and reproduces
its one-pass transfer matrix; it does not solve a ring closure problem.
Periodic closed-orbit recovery remains a separate, harder FFAG-ring validation
lane.  Vacuum observation points are not an air volume mesh.

## Ferromagnetic DUCAS lineage

The Lego proposal is a direct extension of the ferromagnetic-shape method of
Murata, Abe, Ando, and Nakayama, *IEEE Transactions on Applied
Superconductivity* 19 (2009), DOI 10.1109/TASC.2009.2018104.  Their explicit
loop solves the error-field inverse problem with DUCAS/SVD, interprets the
recovered magnetizing current as an equivalent saturated-ferromagnet volume,
updates the material mesh, and repeats the field solution.

Radia keeps that continuous step as an auditable predictor.  Candidate columns
are complete positive-material HDiv-MMM responses, the band-normalized
ACA--QR--TSVD pseudo-inverse returns signed material fractions, and multiplying
each fraction by its candidate volume gives the equivalent material-volume
change.  Positive feasible fractions propose insertion; negative feasible
fractions propose removal.  A 0--1 LP projects the continuous predictor onto
whole HEX states, and an exact Schur/full active-system solve decides whether
the binary state is accepted.  Thus the DUCAS predictor introduces neither
gray material nor a design finite difference.  Each fixed-orbit validation
iteration records singular values, normalized modal field strengths, modal
material amplitudes, and the dominant equivalent-volume changes.

Transfer-matrix error is not sent directly to those material candidates.  On
each fixed design orbit, a separate small dense TSVD first identifies the
reachable field-to-map subspace.  A Chebyshev LP inside that retained subspace
minimizes the largest engineering-band error and produces a target correction
of the sampled binormal field and normal gradient.  The default relative TSVD
cutoff is 1e-3 so near-null optics modes do not create unbounded field steps.
Only that field target enters the HDiv-MMM material inverse.  Candidate
screening therefore performs no particle tracking or closed-orbit root solve.
After each accepted binary batch, the one-pass transfer matrices are
recomputed on the same frozen paths and the optics inverse is relinearized.
The validation history stores both TSVD layers independently.

The public diagnostic entry point
`run_transfer_matrix_material_inverse_pipeline` exposes the same ordering as
five inspectable stages: magnetic-field distribution; forward-AD transfer
matrix; target-minus-realized matrix difference; band-normalized
TSVD/Chebyshev field correction; and ACA--thin-QR--TSVD binary material
screening.  The AD engine seeds the sampled `B` and `dB/dx` coordinates,
propagates scalar generator expressions automatically, differentiates each
matrix exponential with its exact Frechet derivative, and applies product
rules across longitudinal segments.  Finite differences appear only in tests
as an independent regression oracle.  The material response itself is already
an explicit linear candidate-column map, so its low-rank ACA--QR--TSVD solve
does not require a second differentiation layer.

## Second- and third-order R/T/U with x-y coupling

`radia.accelerator_taylor_topopt` extends this contract to the factorial
Taylor map `u_out = R*u + T[u,u]/2 + U[u,u,u]/6 + O(u^4)`.  Its second-order
raw vector has five component blocks per fixed-orbit segment; the third-order
contract adds normal/skew octupoles for seven blocks.  A nine-point
source-free harmonic projection turns live HDiv field evaluations into those
linear multipole rows.  Skew quadrupoles control first-order x-y blocks,
normal/skew sextupoles expose quadratic terms, and normal/skew octupoles expose
direct cubic x-y terms.

The `R/T/U` value is always rebuilt by the native C++ variational-map kernel.
Forward-mode AD differentiates the identical RK4 and Taylor-composition
algebra, including chromatic and lower-order cascade contributions, with
respect to all multipole rows.  The public
`run_second_order_taylor_material_inverse_pipeline` performs the explicit
`R/T` correction and ACA--thin-QR--TSVD material screen.
`run_third_order_taylor_material_inverse_pipeline` adds a component-wise TSVD
reachability certificate before the matching `R/T/U` material screen.
`optimize_hdiv_mmm_magnet_from_second_order_taylor_map` and
`optimize_hdiv_mmm_magnet_from_third_order_taylor_map` are the whole-element
paths.  They contract the AD Jacobian before candidate adjoints, then
completely re-solve and natively re-evaluate every committed binary state in
the original selected engineering bands.

This does not make every numerical `R/T/U` tensor realizable.  Maxwell,
fixed-orbit, multipole, and binary-material constraints still define the
reachable subspace.  `certify_taylor_map_reachability` projects the
band-normalized target error onto the local AD image and reports the residual
by dipole, `R`, `T`, and `U`; it is a local certificate, not a proof of global
nonlinear reachability.  Committed fast tests prove x-y `R` and
skew-sextupole `T_xy` control on one two-HEX problem and skew-octupole `U_xxy`
control on another.  The 42k FFAG model has not yet been rerun with either
high-order objective.

## Canonical Lie completion through f5

`radia.accelerator_lie_topopt` supplies the formal-symplectic path missing from
the raw Taylor objective.  It expands the exact curvilinear body Hamiltonian
through `H4`, builds `A/F2/F3` in the native C++ kernel, propagates analytic
rank-five forward-AD tangents, and performs the Dragt--Finn factorization
`M = R o exp(:f3:) o exp(:f4:) o exp(:f5:) + O(z^5)`.  The extraction removes
the `f3` self-cascades and ordered `f3`--`f4` cross before recovering `f5`;
the reconstructed `R/T/U/V` map is checked coefficient by coefficient for
formal symplecticity through fourth map order.
Finite-amplitude application uses implicit-midpoint Hamiltonian flows and is
symplectic as an evaluated map.

The same HDiv-MMM optimization chain remains in force: multipole response,
AD Lie map, target difference, local TSVD reachability/correction,
ACA--thin-QR--TSVD material inverse, whole-HEX proposal, full active-system
solve, and exact Lie-map acceptance.  The target itself must satisfy the
formal symplectic tolerance.  Skew quadrupole, skew sextupole, and skew
octupole rows expose first-, second-, third-, and cascade-generated
fourth-order x-y/chromatic control respectively.
The tracked Wolfram Language derivation emits an independent symbolic golden
for `H2/H3/H4`, the `f3` self-cascade, `f4`, their fourth-order cross, and an
independent `f5` kick; it is not a runtime dependency.
The complete convention, equations, and deliberate physical boundary are in
[`LIE_MAP.md`](LIE_MAP.md).

The full-field reference path is first recovered by the existing DOP853
periodic-orbit solver.  `fourth_order_lie_map_from_tracked_orbit` now carries
that `PlanarDesignOrbit` into the planar Frenet--Serret moving frame, samples
normal/skew multipoles through octupole, and builds the fourth-order Lie map.
The later variational RK therefore tracks deviations, not the reference orbit
again.

This Lie completion is deliberately bounded to source-free,
piecewise-constant body multipoles through octupole and an `H2/H3/H4`
Hamiltonian.  Direct `H5` kinematic/decapole terms, longitudinal fringe/edge
vector potentials, nonplanar torsion, and arbitrary-order normal forms remain
separate work; they are not approximated silently.

Focused regression tests cover the symbolic coefficients, native and
forward-AD map agreement, formal-symplectic gates, direct GridFunction paths,
MATLAB MEX/fallback boundaries, reachability certificates, and complete
two-HEX material solves.  The release gate owns aggregate pass counts so this
method note does not freeze a transient repository-wide test result.

The field-correction inverse is only a proposal model.  Every proposed binary
batch must also reduce the original bend-field/transfer-matrix maximum band
ratio after a complete active-system solve.  `--map-trust-region-trials`
shrinks the whole-element material-volume radius when this exact one-pass map
gate rejects a proposal; rejected states never become the incumbent.
When `--direct-map-oracle-fallback` is enabled, exhaustion of those field-target
trials triggers one bounded all-candidate contraction in the original map-error
space.  It composes the forward-mode AD field-to-map Jacobian with the HDiv-MMM
candidate response; it introduces neither a design finite difference nor an
unbounded solve loop.  After that local proposal stalls,
`--direct-map-exact-beam-width` and `--direct-map-exact-beam-depth` explicitly
enable shallow nonmonotone binary look-ahead for the direct oracle only.  The
optional `--direct-map-graph-front-proposal-limit` adds a bounded set of
connected alternatives seeded by the same global ACA--QR--TSVD modes, so beam
width corresponds to distinct binary states rather than duplicate trust-radius
proposals.  The exact fixed one-pass map gate still owns acceptance.

The original method required a suitable stochastic initial shape because a
large continuous volume update could diverge.  This lane instead uses a fixed
return-yoke seed, a physical changed-volume trust region, connectivity gates,
and exact re-solves.  Nonlocal graph/beam exploration is a separate discrete
search layer; it does not alter the DUCAS material predictor.

The all-candidate ACA--QR--TSVD screen remains global.  The
`--exact-candidate-limit` option bounds only the costly conditional
block-Schur representative front after that screen.  A positive
`--exact-beam-width` and `--exact-beam-depth` retain fully resolved,
topology-valid states across a shallow objective barrier.  The returned design
is always the best fully solved incumbent, never an intermediate worsening
state.

`validation_ffag_hdiv_mmm_poc.py` is the next, study-scale lane.  It places a
structured BDM1 HEX candidate pole slab beside the unmeshed orbit aperture,
uses native Laplace field rows at several orbit energies, and requires at least
one exact whole-element move to improve the fused bend/map objective.  The
source and one-sided slab are intentionally simple so this remains a method
PoC; a final design needs the two poles, return yoke, and CoilBuilder geometry.

The saved study-scale result uses 192 HEX / 6,912 BDM1 DoFs and three momentum
points.  Two accepted exact re-solves remove 13 whole cells and reduce the
maximum normalized bend/map residual by 5.11 %.  The target bands are not yet
met; the result is evidence for the fused optimization path and its descent,
not a completed FFAG magnet.
