# FFAG HDiv-MMM topology-optimization validation

## EarlyTimes C-type A/B route convergence

`validation_earlytimes_ctype_ab.py` isolates the field-route side of the
EarlyTimes error triangle before any Lie map is built.  It solves a compact
C-type magnet and first evaluates three representations at identical global
aperture points: the `HDiv(order=4)` B GridFunction, NGSolve-native
`curl(HCurl(order=5) A)`, and the independent axial-reflection-symmetrised
HDiv-MMM plus `rad.Fld` B source.  `curl(A)` is a boundary-consistency
diagnostic only; it is not substituted into either the canonical A Hamiltonian
or the independent B-RK route.  Only after this three-way field comparison is
acceptable may the script proceed to A-RK/B-RK or optional A-RK/Lie studies.

Both HDiv-MMM source fields are analytic.  `--a-construction exact` (the
default) integrates the equivalent-current identity
`A = mu0/(4 pi) [INT (curl M)/R dV + INT (M x n)/R dS]` in closed form, matching
the analytic charge kernels the B route already used; `--a-construction
quadrature` restores the point-dipole cloud that the recorded q14/q20/q24
studies used.  The distinction is decisive rather than cosmetic here: the beam
tube stands 1.5 mm from 30 mm iron elements, so a point-dipole cloud is deep in
their near field.  On the same 45-point set the quadrature source left
`curl(A) - source B` at 1.1e-2 T even at order 24, while the analytic
construction reaches 3.6e-7 T, which is the probe's own central-difference
truncation.  With the source exact, the residual gate error is the loft-chain
FE representation of each field, which is what the mesh/order settings control.

That representation error now converges as a finite-element error should.  On
the same 45-point set, with `HCurl(order=5)` A and `HDiv(order=4)` B, the
measured gate is

| beam loft chain | HEX | HDiv B - source B | curl(HCurl A) - source B |
|---|---|---|---|
| 1 layer, 5.0 mm, 1 subdivision | 128 | 5.26e-5 T | 1.18e-5 T |
| 3 layers, 5.0 mm, 1 subdivision | 384 | 1.40e-6 T | 2.31e-6 T |
| 5 layers, 2.5 mm, 2 subdivisions | 2560 | 1.53e-7 T | 5.06e-8 T |

against a peak source field of 0.390 T, i.e. 3.9e-7 and 1.3e-7 relative on the
refined chain.  Raising `--curve-order` from 2 to 4 changes neither column at
all: this orbit bends by roughly 1e-7 m of sagitta per longitudinal cell, so
geometric order is not the limit here.  The finite-difference curl of the raw
exact A source stays at 1.43e-6 T at a 4e-5 m step and 3.58e-7 T at 2e-5 m --
exactly second order, and independent of the beam mesh -- confirming that the
probe, not the source, sets that floor.  On the refined chain the projected
fields already agree more closely than that probe can resolve.

The three-way record retains every sampled point, all three Cartesian B
arrays, all pairwise differences, componentwise maxima, and maximum vector
norms.  The default 5-by-3-by-3 set spans longitudinal position, both x sides,
and the lower/median/upper y layers.  The parity-conditioned B
`CoefficientFunction` is also used to recover the design orbit and report the
source-to-HDiv projection error; the raw unsymmetrised source remains a
separate symmetry diagnostic.

The durable result JSON separates the zero-orbit A/B bias from the centered
first-order transfer-matrix difference for `x,px,y,py,delta`.  It also records
transverse derivative-step convergence, aperture-amplitude convergence, source
quadrature statistics, source-to-HDiv projection discrepancy, gauge residuals,
exact settings, hostname, and runtime.
Heavy mesh and aperture sweeps run on an idle compute host; the script refuses
to overwrite a result unless `--overwrite` is explicit.

The production aperture target is conservatively interpreted as
`x,y in [-20,20] mm`.  The loft-chain builder supports this without asking one
high-order element to represent the full 40 mm span: the four x macro-strips
remain fixed topology boundaries and may each be uniformly subdivided, while
vertical refinement uses an odd symmetric layer count so `y=0` stays inside
the central layer.  The initial order-6 hp smoke configuration uses two
subdivisions per x macro-strip and nine y layers, giving 72 HEX elements per
longitudinal cell and local transverse widths no larger than 5.0 mm by
4.44 mm.  This proves mesh/space support, not field accuracy; the usable
20 mm aperture is certified only after three-way field and Lie/RK convergence.
The convergence script exposes the matching physical C-yoke as
`--wide-20mm`: it sets a 50 mm full pole gap, 60 mm pole width, the
`+/-20 mm` beam mesh, order 6, 8 x strips, 9 y layers, and an 18 mm sampled
fit aperture.  The legacy small-gap defaults remain available solely so the
recorded q14/q20/q24 studies remain reproducible.

After the field gate, the script's `gauge_invariance_triangle` compares three
tracks in mechanical exit variables, converting each route's canonical momenta
with its OWN vector potential: the gauged-A fourth-order Lie map, the exact
canonical A-RK on the UNGAUGED (axial-only) field, and the HDiv B-map
Cartesian RK.  The ungauged-versus-gauged A-RK difference is a pure gauge
invariance check of the whole canonical machinery and reached 4.8e-10 on the
C-type fixture.  The two independent field routes differ through the
piecewise-constant `h(s)` orbit discretization, first order in the station
spacing (1.15e-6 at 33 stations, 5.4e-7 at 65).  The Lie-versus-A-RK leg is
dominated by the vertical plane and is the measured aperture truncation of the
declared degree-five transverse jet contract: it is insensitive to segment
count and to amplitude, halves with the fit aperture, and matches the jet fit
residual.  The maps and tracks consume a `trim_orbit` interior orbit while the
loft tube is built from an extended track (`--orbit-margin-stations`), because
a fit section on a tube boundary face samples the least-controlled one-sided
derivative content of the discrete fields; the interior fitted-Hamiltonian
floor (about 8e-6, mesh-converged) is genuine jet-contract content, which the
`--lie-reference-orbit-tolerance` gate prices at an exit-coordinate effect of
order `H1*L^2/2 ~ 3e-8 m`.

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

### Full C-yoke manufactured inverse gate

The 2026-08-22 100-machine recovery run used the 1,180-HEX, 42,480-DoF BDM1
C-yoke and an exact manufactured target obtained by adding inactive element
798.  Starting from 720 active cells, the analytic transfer-map-to-field
inverse and ACA--QR--TSVD material screen selected that element without using
the direct-map fallback.  The accepted 721-cell state was then solved in full
and reduced the exact fixed-orbit transfer-matrix maximum band ratio from 5.0
to 0.5 at all three momentum points.  The known active set was recovered, the
binary topology and coil-clearance gates passed, and the largest realized-map
symplectic residual was 1.06e-15.  Performance recording was disabled; these
are correctness data, not benchmark timings.  The compact retained evidence is
`manufactured_inverse_recovery_100_20260822.json`.

### Eight-candidate direct-map inverse gate

The 2026-08-23 follow-up used the same 1,180-HEX, 42,480-DoF model but
deliberately shuffled eight legal growth candidates, with the known target
element 798 fourth in the input list.  The original transfer-matrix error was
contracted directly with the analytic HDiv-MMM candidate response.  The global
ACA--QR--TSVD predictor retained rank 3 (ACA rank 4); a beam search was not
used.  All eight candidates were then retained in one bounded conditional
block-Schur front.  Its 37 evaluated bundles selected element 798, and the
mandatory full active-system solve recovered the known 721-cell state and
reduced the exact map maximum band ratio from 5.0 to 0.5.

This result also records an important algorithmic boundary: continuous signed
DUCAS fractions are proposal diagnostics, not the final add/remove decision.
For this case the continuous coefficient of the correct inactive element was
negative, whereas the whole-cell exact Schur response correctly selected its
insertion.  An empty/zero-rank low-rank proposal must therefore pass its
caller-bounded candidate front to the exact selector rather than silently
ending the iteration.  The initial state is reused only after the configured
H-matrix true-residual and inactive-DOF gates pass.  No finite-difference
design sensitivity, gray material, or nonmonotone beam was used.  Performance
recording was disabled; the compact correctness evidence is
`manufactured_inverse_8candidate_direct_map_100_20260823.json`.

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

### Material-aware section specification and smooth contour loop

`radia.accelerator_abe_topopt` is the continuous, smooth-interface companion
to the binary Lego path.  It preserves the four physical variables explicitly:
section transfer specification, orbit-field/multipole response, whole-cell
magnetization/fill, and iron-interface height.  An HDiv state is never bounded
DOF by DOF.  `measured_element_fill_patterns` obtains one complete solved
magnetization pattern per cell and `contract_hdiv_element_fill_response`
contracts the observation rows to one scalar fill column per cell.  Copying a
local HDiv block into an inactive cell requires either an NGSolve-owned pattern
transfer callback or an explicit assertion that a structured discontinuous
mesh has compatible local ordering.

For a section specification with AD Jacobian `J` and element field response
`E`, `compose_specification_fill_response` forms `J E`.  This does not remove
the field layer: it lets the bounded material inverse select a physically
cheap member of the many field changes that produce the requested optics, and
the implied field change is reported afterward.  Air cells have fill capacity
`[0,1]`; reference iron cells have `[-1,0]`.  The weighted Abe solver retains
its ACA--QR--TSVD mode evidence and reuses the factor during bound correction.

`bin_element_fill_to_interface_height` converts signed cell volume to a
conservative height field, and `blended_interface_displacement` keeps the
aperture and pole root fixed.  `optimize_abe_section_contour` then realizes
each absolute accumulated fill from the reference geometry, performs a
caller-owned complete HDiv solve and exact map evaluation, backtracks rejected
steps, and optionally rebuilds the analytic/AD response after every accepted
shape.  Only exact states that improve the original engineering-band objective
and pass the supplied bend/orbit/quality guard become incumbents.

The promoted 170-by-2560 saved-data regression retains two specification
modes with condition number 123.1.  Relative to the earlier research driver,
the whole-element contraction differs by `1.63e-19`, the AD specification
composition by `1.82e-17`, the bounded fills by `3.66e-16`, and the implied
field by `9.15e-20`.  This is an API-parity regression of the material stage,
not a claim that the outer smooth-contour problem is solved: the recorded
single-pass exact shape delivered 55.3 percent, so relinearized outer passes
remain the full-scale validation target.

The signed element-fill solve also has a Python-free standalone MATLAB path:
`radia.topopt.solveAbeElementFillPlan` calls
`radia_mex('topopt.abe_element_fill_plan', ...)`.  The MEX command shares the
HACApK ACA+ and QR--TSVD C++ implementation used by pybind, retains one factor
through all capacity-clipping corrections, and returns residual and volume
histories.  NGSolve DOF orientation, pattern transfer, GetTrafo deformation,
and the complete-solve outer loop deliberately remain caller-owned.  The
focused MATLAB regression checks both dense-reference and ACA paths, signed
air/iron capacities, the EarlyTimes Jacobian composition, whole-element HDiv
contraction, conservative interface height, and fixed-surface blending.

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

For a decoupled first-order target, `Symplectic2x2KAN` gives global Iwasawa
coordinates for each transverse `Sp(2,R)` block, while
`DecoupledFirstOrderTarget` restricts the static-magnetic longitudinal block
to `[[1,R56],[0,1]]`.  `certify_decoupled_first_order_reachability` measures
the six independent transverse KAN directions rather than treating eight
ABCD entries as independent.  A drift-only seed has rank three; a generic
segmented normal-quadrupole seed reaches rank six.
`solve_decoupled_first_order_continuation` follows a staged KAN homotopy and
rejects any step that loses rank six, violates the nonlinear target band, or
introduces x-y-longitudinal coupling.  It is an auditable continuation attempt,
not a global theorem: distant targets, restricted material-response bases, and
an `R56` change outside the available controls can still return unreachable.
Committed random-target tests exercise nearby targets from a rank-six seed.

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
octupole and decapole rows expose first-, second-, third-, and direct plus
cascade-generated fourth-order x-y/chromatic control respectively.
The tracked Wolfram Language derivation emits an independent symbolic golden
for `H2/H3/H4/H5`, the `f3` self-cascade, `f4`, their fourth-order cross, and an
independent `f5` kick; it is not a runtime dependency.
The complete convention, equations, and deliberate physical boundary are in
[`LIE_MAP.md`](LIE_MAP.md).

The full-field reference path is first recovered by the existing DOP853
periodic-orbit solver.  `fourth_order_lie_map_from_tracked_orbit` now carries
that `PlanarDesignOrbit` into the Bishop/RMF double-reflection moving frame,
accepts a real `HCurl(order=p)` vector-potential GridFunction, samples NGSolve's
native `curl(A)` without an HDiv projection, fits normal/skew multipoles through
decapole, and builds the fourth-order Lie map. The retained FESpace supplies
`p`; it is not duplicated as a beam-input setting.
The later variational RK therefore tracks deviations, not the reference orbit
again.

This Lie completion is deliberately bounded to source-free,
piecewise-constant body multipoles through decapole and an `H2/H3/H4/H5`
Hamiltonian.  Longitudinal fringe/edge vector potentials, nonplanar torsion,
and arbitrary-order normal forms remain separate work; they are not
approximated silently.

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

## Exact active-system performance evidence

`solver_active_hmatrix_lab_20260826.json` records a same-process paired LAB
measurement of exact inactive-leaf pruning.  The 42,480-DoF BDM1 HEX problem
keeps 25,920 active DoFs.  Pruning preserves the active principal system
exactly while reducing active charge leaves, accelerating both direct operator
application and the scalar-finishing path of mass-Riesz PCG.  The JSON keeps
kernel and complete-solve timings separate and records the remaining
mass-Riesz and block-PCG work separately.  A later entry in the same JSON
records the exact local Cholesky mass-Riesz path for broken HDiv: it retains
PARDISO for a connected conforming mass, and omits the unused Jacobi-diagonal
setup when mass-Riesz is active.  Repeat the measurement on mdx or hibino after
their cooling-system recovery before using the LAB wall times as publication
claims.
