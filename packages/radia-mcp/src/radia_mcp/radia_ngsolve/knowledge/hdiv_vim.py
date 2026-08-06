"""HDiv-type VIM knowledge for radia_mcp.radia_ngsolve.

This module describes the current Radia soft-iron direction: mesh-backed
H(div) volume-integral demagnetization using NGSolve meshes and Radia's C++
charge-Gram H-matrix.  It intentionally avoids retired solver history so MCP
answers steer agents toward the live implementation.
"""

_FAMILY = r"""
# BDM versus Raviart--Thomas: mandatory terminology

NGSolve's `HDiv` class supports two distinct H(div)-conforming families:

```python
fes_bdm = HDiv(mesh, order=p)           # BDM, the NGSolve default
fes_rt  = HDiv(mesh, order=p, RT=True)  # Raviart--Thomas
```

Radia's established `vim.Solve`, `PlanarDemagBody`, `MagnetizationSource`,
nonlinear material, IMA, charge-Gram, and persistent-field production paths
construct the first form and therefore use **BDM1/BDM2**, not RT1/RT2.
Historical Radia names and filenames containing `rt1` or `rt2` are legacy
labels; they do not change the actual NGSolve family.  Never infer
Raviart--Thomas from the class name `HDiv` alone.  Actual RT is selected only
when the code contains `RT=True`; keep it labelled as an explicit comparison
or research path until it has its own production promotion.

On simplices, the useful polynomial distinction is
`BDM_p = [P_p]^d`, with divergence in `P_(p-1)`, versus
`RT_p = [P_p]^d + x P_p`, with divergence in `P_p`.  Both have continuous
normal traces, but they are not interchangeable and have different DoF counts.
The Mathematica study implementation is executable rather than descriptive:
`simplex_ho.wls` exports `HDivTetBDM`, `HDivTetRT`, `HDivTrigBDM`,
`HDivTrigRT`, and family ledgers that prove equal `ker(div)` dimensions.
`hdiv.wls` exports the Hex/Quad RT/BDM names as aliases of NGSolve's shared
tensor-product space, matching the fact that the `RT` selector does not change
the local Hex/Quad space.
"""

_OVERVIEW = r"""
# HDiv-type VIM demag operator

Radia's soft-iron demag path is the HDiv-type volume integral method (VIM).
Its production NGSolve space is BDM: bare `HDiv(mesh, order=p)` selects BDM,
whereas Raviart--Thomas requires `RT=True`.  Never call the production
BDM1/BDM2 path RT1/RT2 merely because both families conform to H(div); request
the `family` topic for the polynomial and API distinction.
The unknown magnetization is represented in an H(div) space on an NGSolve mesh;
the demag operator is

    N = B^T G B

where B maps magnetization to volume/surface magnetic charge, and G is the
Laplace single-layer/Coulomb Gram.  The important engineering properties are:

- loop modes are field-null by construction because loops live in ker(B);
- the operator is symmetric and compatible with MINRES/CG-style solvers;
- material state, fields, and source terms live in NGSolve Mesh,
  GridFunction, CoefficientFunction, and BilinearForm vocabulary;
- this makes reduced FEM coupling cleaner than sampling a separate object
  field back into FEM.

Use this route for mesh-backed TET, HEX, and WEDGE soft iron, nonlinear BH
curves, and planar motor cross sections.  Mesh-less soft-iron solves are not a
supported Radia production path; create an NGSolve mesh and call the HDiv API.
"""

_EDDY_BUBBLE = r"""
# BDM-MMM plus eddy-bubble HCurl-VIM

The production 3-D coupled entry point is
`radia.vim.NgsolveBDMEddyBubbleVIM`.  One shared NGSolve mesh and one
`SharedMeshMaterialModel` feed two different conforming unknown spaces:

- bare `HDiv(mesh, order=1|2)` supplies the BDM-MMM magnetization branch;
- a high-order `HCurl` parent supplies `T`, with `J = curl(T)`, and is reduced
  by external-response Krylov/EVRS plus topology classes.

The production HCurl reduction is two-stage.  After the parent-space Krylov
basis is formed, `CompressHCurlResponseInCurrentGram` eigendecomposes the
sampled `J=curl(T)` Gram matrix.  Relative-null directions are independent in
the parent-T metric but carry no independent current; they are removed, and the
same whitening transform is applied to parent T and sampled J.  Only then are
bulk, conductor-cycle bridge, and air-facing SIBC blocks assembled.

`NgsolveBDMHDivMMMResponseReduction` creates the BDM parent internally and
records `parent_family="BDM"` and `parent_order` in diagnostics.  Physical
stator/rotor applied-H responses are protected before energy-POD compression;
regular-solid-harmonic or rotor-angle fields enrich the training complement.
Use the lower-level `NgsolveHDivMMMResponseReduction` only when an explicit
supplied BDM or `RT=True` comparison space is required.

The eddy-bubble reduction is face- and topology-aware.  Bulk response modes,
conductive-graph cycle/bridge modes, and surface-Omega modes are retained as
separate blocks.  Only conductor faces touching air are SIBC faces; internal
conductor-conductor faces are not SIBC.  `solve_frequency` reconstructs parent
HCurl `T`, parent BDM magnetization, sampled currents in every retained block,
average Joule loss, and residuals.  `solution.eddy_flux_density(points)` sums
the quasi-static Biot--Savart field of bulk, bridge, and SIBC currents; pass a
block name to inspect one contribution.

For the fully integrated reduction use
`mixed.solve_frequency_eddy_bubbled(frequency)`.  The production builder stores
`volume=bulk`, `volume1=bridge`, and `surface=sibc` roles from neighboring
material labels.  The convenience solve protects bridge/cycle and air-facing
SIBC blocks and eliminates only the ordinary bulk eddy-bubble block.

The transformation is built from the complete coupled HDiv/HCurl matrix.  It
therefore Schur-updates the HDiv operator and both coupling blocks, projects the
RHS, includes the affine lift for a nonzero eliminated RHS, and reconstructs
the full HCurl coefficients.  The HDiv parent demag operator continues to use
the C++ `_ChargeGramHMatrix` HACApK backend; the final mixed Schur system is a
small dense retained-mode matrix.

## Response-basis completeness gate

Do not interpret spatial h-refinement as proof that the HCurl response basis is
complete.  A single uniform-field vector-potential port can converge to the
wrong reduced subspace: increasing its Krylov depth only explores moments of
that one input.  The conducting-disk validation freezes this as a negative
control.  The single-port model remains more than 4% away from the independent
axisymmetric diffusion pole at both 8 and 16 Krylov steps.

For an axisymmetric disk driven through `A=(-y,x,0)/2`, enrich the training set
with independent spatial moments before changing the mesh:

    A, r^2 A, r^4 A, z^2 A

The four-port, three-step basis passes the public h-refinement ladder and the
HCurl p=1/2/3 ladder at less than 2% error, with p=2 reaching the recorded
sub-percent gate.  General geometries need geometry- and excitation-aware
moments; these four functions are a disk recipe, not a universal basis.

## HDiv local-response completeness gate

When mean magnetic flux converges but Joule loss is much too large, do not
immediately diagnose an unphysical HCurl star current.  First freeze the HCurl
space and replace only the reduced HDiv response basis by its full parent.  If
the Joule mismatch collapses, the missing directions are local magnetic
charge/star responses in the HDiv reduction; they are not extra current modes.
This distinction is especially important when the current branch already has
`J=curl(T)`, `nograds=True`, tangential-zero boundary conditions, and disabled
bridge cycles.

Use `NgsolveHDivLocalPolynomialTrainingPorts(mesh, max_degree=2)` as the
conservative enrichment.  It supplies 30 normalized vector-polynomial probes
through degree two to `training_fields` of
`NgsolveBDMHDivMMMResponseReduction`.  The probes deliberately include
non-source-free residual directions, so they are **training-only** and the API
rejects them as physical `external_fields`.  Keep regular-solid-harmonic ports
for physical source-free applied-H content.

The frozen nonlinear magnetic-conductor gate holds the HCurl branch fixed,
checks the reduced N2 model against the full HDiv parent, compares degree two
against degree three, then performs N4 h-refinement against an independent
A-phi solve.  Degree two to three changes Joule loss by about 0.015%; the N4
degree-two result is within 1.02% in Joule loss and 4.48% in mean B of the
independent formulation, with mixed energy-balance norm below `2e-4`.  Treat
these as one case-specific validation gate, not a universal accuracy claim.

When extracting a free-decay time constant, select the mode that is dominant
for the physical excitation by its port residue.  The largest generalized
eigenvalue can be a nearly port-invisible mode and is not automatically the
observable pole.  Always record the residue fraction, R/L passivity,
current-Gram rank, projection residual, and both h- and p-refinement.  See
`validation_test/vim_coupled/validate_hcurl_eddy_bubble_disk.py` and its
result JSON.

## Frequency-route validity gate

`EddySIBCApplicability` is part of the assembled model identity, not a label
that can be reused across an arbitrary frequency sweep.  The gate is carried
from `TopologyAwareHybridVIM` into `CoupledHDivHybridVIMSystem`.  A frequency
solve may move within the same volumetric or SIBC regime, but
`solve_frequency` rejects a point whose thickness/skin-depth and curvature
test changes the selected route.  Rebuild the topology-aware system at that
frequency so the retained volume and surface bases match the physics model.

This preflight must run before blaming a high-frequency discrepancy on BDM,
HCurl polynomial order, or response-basis rank.  For a thin magnetic conductor,
first compare an axisymmetric volume solve with an independently refined full
3-D HCurl volume solve.  A forced volumetric reduced solve above the SIBC
transition is an invalid extrapolation, even when its linear residual is at
machine precision.  Residual convergence proves the assembled algebra was
solved; it does not prove that the frequency-appropriate model was assembled.

## Nonlinear magnetic-conductor SIBC coupling

An air-facing SIBC basis does not by itself imply ESIM.  For a linear magnetic
conductor, use the ordinary frequency-dependent `SkinImpedance`.  When the
skin permeability follows a nonlinear B-H curve, use
`LocalESIMSurfaceModel`: its one-dimensional normal cell returns the local
field-dependent `Z_s(f, |H_t|)`, which is projected into a
`SurfaceImpedanceGram` rather than averaged into one scalar.

`CoupledHDivHybridVIMSystem.solve_frequency_local_esim(model, frequency)` now
runs that Karl iteration around the complete HDiv-MMM/HCurl-VIM mixed solve.
Every update reconstructs solution-shaped tangential field samples, updates
the local ESIM impedance, and resolves the magnetic, bulk-current,
conductor-cycle, and SIBC blocks together.  Nonlinear superposition is invalid,
so the method requires exactly one physical excitation.  The returned
`CoupledHDivHCurlLocalESIMSolution` retains the converged mixed solution,
surface Gram, local fields, update history, residual, and Joule loss.

The executable same-region TET gate assigns both magnetic and conductive roles
to one `iron` label.  Its 50/1000/5000 A/m ladder converges in 2--5 updates,
keeps every local Gram passive, gives positive Joule loss, and matches a
fixed-Gram replay to machine precision.  Increasing field lowers the mean
surface resistance for the supplied monotone saturating B-H curve and increases
the departure from linear SIBC.

Keep the claim boundary exact.  This API updates nonlinear **skin** impedance
while the ordinary bulk HDiv-MMM operator remains at a fixed low-field scalar
permeability.  A simultaneous constitutive iteration that updates both the
ordinary bulk nonlinear B-H operator and the local ESIM Gram is not yet
implemented.  Hysteretic or rotational skin state and multidimensional
edge/corner cells also require additional state/cell models.

## Thin magnetic-conductor adjudication

For a 20 mm diameter, 0.5 mm thick linear magnetic conductor, use an
axisymmetric Q2 volume solve as one reference and a separately refined full
3-D HCurl A-form as the second.  At the recorded 100 Hz and 10 kHz points, the
fine full 3-D results agree with the axisymmetric reference within 2%.  This
cross-formulation agreement is the reference gate before judging a reduced
coupled model.

On the static magnetic branch, the mapped-HEX BDM1 HDiv-MMM ladder decreases
from about 3.1% to 1.2% reference error as the nominal size moves from 2 mm to
0.5 mm.  The tested BDM2 HEX ladder is non-monotone and is explicitly not
promoted.  Do not infer BDM2 accuracy from polynomial order alone; require its
own h/p and observable convergence evidence.

The cause is now bounded more tightly.  A two-cell affine HEX BDM2 diagnostic
has generalized demag spectrum inside the physical interval, while the same
two-cell warped map reaches about 1.16.  Four warped cells reach about 1.57,
and an annular-sector map develops both a negative mode and an upper mode near
1.80.  Dense storage reproduces the same result, so ACA compression is not the
cause.  Volume-only and boundary-only blocks are positive but individually
large; their cross term must preserve a strong cancellation.  Raising
quadrature and widening the near region improves the spectrum but is too
expensive and does not robustly restore the annular case.  Therefore
`vim.Solve(order=2)` rejects mapped/non-affine pure HEX material models.  Use
mapped HEX BDM1, affine HEX BDM2, or pure-TET BDM2.  `ChargeGram` is still an
explicit diagnostic surface; a finite residual is not permission to promote
the rejected material lane.  Removing this gate requires a composite
mapped-HEX BDM2 charge operator that preserves the volume/surface cancellation;
raising quadrature alone is not a production remedy.  The magnetic-conductor
validation runner exposes this as ``--include-bdm2-gate`` and records the
expected rejection alongside the independent axisymmetric Q2, full 3-D HCurl,
mapped-HEX BDM1, and coupled eddy-bubble evidence.

Do not diagnose this as an intrinsic BDM2-space failure.  A separate
54-HEX air/body H1 Omega experiment assembles the discrete Hodge operator
``C.T K^-1 C`` on the same 207 active mapped-body BDM2 DoFs.  H1 orders 2 and
3 have generalized maxima about 0.992 and 0.988 with no modes outside
``[0,1]``.  The same two-cell mapped charge diagnostic reaches about 1.2175
with three out-of-range modes.  This localizes the defect to the mapped
volume/surface charge operator.  The finite H1 box is a contraction and
formulation-separation gate, not an open-boundary accuracy oracle.
The `[0,1]` contraction is specific to the standard unit H1 stiffness metric;
a weighted or Kelvin-transformed metric needs its own independent bound.

The executable API is
`H1HodgeDemagOperator(hdiv, h1, definedon=...)`.  Pass it as the
`demag_operator` of `NgsolveHDivMMMResponseReduction`; the generic NGSolve CG
path evaluates solver expressions into concrete vectors and uses
`hdiv.FreeDofs()` when the magnetic space is restricted to one material.  On
the same 54-HEX body, the resulting BDM2 response reduction feeds a shared-mesh
HCurl eddy-bubble solve with snapshot residual below `3e-12`, mixed residual
below `3e-16`, and positive Joule loss.  This verifies mixed mechanics and the
BDM2 approximation space.  It does not repair the open-boundary charge kernel,
and its sampled HCurl interaction is not an accuracy oracle.

For the one-call mixed API, pass both
`hdiv_definedon=mesh.Materials("body")` and
`demag_operator_factory=lambda hdiv: H1HodgeDemagOperator(hdiv, h1,
definedon="body")` to `NgsolveBDMEddyBubbleVIM`.  The factory receives the
exact restricted BDM space constructed by the builder.  A restricted space
without a factory, or simultaneous `demag_operator` and factory arguments, is
rejected rather than relying on coincident DOF numbering.

Accuracy is a separate gate.  For the static circular disk, the independently
refined axisymmetric Q2 reference is `1.0916977441`.  The 3-D finite-domain
H1-Hodge sequence reduces relative error from 4.65% (coarse BDM1/H1-p2) to
2.24% (coarse BDM2/H1-p3), 1.44% (fine BDM2/H1-p3), and 0.98% (fine
BDM2/H1-p4).  This strict h/p ladder supports the BDM2/H1 formulation on that
disk; it is not a universal open-boundary or solver-superiority claim.

The coarse coupled HDiv-MMM plus HCurl eddy-bubble run is a mechanics gate: it
checks the mixed solve, non-negative Joule loss, small algebraic residual, and
rejection of a stale volumetric-to-SIBC frequency route.  Its explicitly
sampled interaction is not an accuracy oracle by itself.  The subsequent
sampled h ladder is stronger: 32, 96, and 384 mapped HEX reduce the independent
reference error from 5.38% to 3.34% to 1.90%, with non-negative Joule loss and
small residual at every level.  This promotes the explicit sampled-interaction
route as a bounded h-convergent validation lane.

For a warped order-1 pure HEX mesh, the production HCurl diagonal now uses the
native Q2 isoparametric HEX Gram directly.  It projects the contravariant-Piola
reference density `K=curl(T)*abs(det(dX/dxi))`, which is tensor-Q2 to machine
precision on this map.  Do not try to meet a geometry tolerance by repeatedly
splitting a trilinear HEX into affine sub-tetrahedra: that residual decreases
too slowly and the charge count explodes.  `hex_geometry_backend="auto"`
selects the direct path for a warped order-1 pure HEX,
`hex_geometry_backend="direct"` forces it, and `"affine-subtet"` retains the
diagnostic legacy route.  Affine HEX keeps the established analytic subtet
path by default.

The production direct-Q2 lane now has its own strict 32/96/384-HEX h ladder.
Its independent-reference error decreases from 5.38% to 3.34% to 1.90%.  At
the fine level it agrees with the sampled observable to 8.03e-7 and Joule loss
to 2.93e-6; the largest direct algebraic residual is 3.90e-16.  This establishes
h convergence for this thin magnetic-conductor disk lane, not for every
geometry, material, or observable, and it is not a universal
solver-superiority claim.

Do not form the sampled diagonal as an `Nsample x Nsample` dense pair matrix
when replacing it by the direct-Q2 diagonal.  Build the stable global sampled
H-matrix once, subtract its diagonal through a matrix-free restricted operator,
and add the direct-Q2 operator on the same reduced interval.  Projection
coefficients below the established roundoff threshold may be pruned before the
native HEX Gram is built; retain unpruned/pruned charge counts in diagnostics.
On the recorded same-machine ladder these two structural changes reduced total
runtime from 1819.7 s to 339.8 s, but that observation is not a portable speed
guarantee.

Reference fidelity is part of that gate.  The 384-HEX result is 1.90% from the
18432-element fine axisymmetric Q2 reference, but 2.19% from the 6656-element
quick reference.  The latter is a frozen reject.  A coupled h-ladder accuracy
claim must use the fine reference; a coarse reference is suitable only for a
mechanics smoke.
"""

_IMPLEMENTATION = r"""
# Implementation

Primary Python entry points:

- `radia.vim.MeshSoftIron(mesh, mu_r=... | bh_table=...)`
- `radia.Solve(model, prec, maxiter, method, demag_backend="hdiv")`
- `radia.vim.Solve(mesh, mu_r=... | bh_table=..., H_ext=..., image=...)`
- `radia.vim.HDivSolver(mesh, order=1|2, ...)` for repeated 3D loads,
  history, and coupled bodies.  It owns the HDiv space and geometry-only
  ChargeGram; result mappings do not contain a reusable-operator token.
- Permanent magnets use one four-level model ladder, all on the HDiv charge and
  field machinery:
  1. fixed/given distribution:
     `MagnetizationSource(mesh, M_given)`;
  2. linear recoil/demagnetization:
     `Solve(mesh, mu_r=mu_rec, B_r=B_r, H_ext=...)` with
     `B = mu0*mu_rec*H + B_r`;
  3. simplified Play:
     `PlayHysteresisMaterial(...)` with `SolveHysteresis(...)`;
  4. full B-input EnergyStop:
     `EnergyStopMaterial(...)` with `SolveHysteresis(...)`.
  `B_r` is in tesla and may be a constant vector or a spatial NGSolve
  CoefficientFunction.  The linear-recoil path is the exact right-hand-side
  shift of the symmetric C++ HDiv solve, requires scalar `mu_rec > 1`, and
  defaults to zero applied field.  Use level 1 for the rigid `mu_rec=1` limit.
  A spatial `B_r` belongs to one conforming body.  A jump in normal
  magnetization requires separate body spaces to retain interface charge;
  fixed segments already use separate MagnetizationSource objects, while
  mutually coupled recoil segments require the multi-body block formulation.
- `radia.vim.MagnetizationSource(pm_mesh, M_given, order=1|2)` followed by
  `radia.vim.Solve(iron_mesh, ..., magnetization_sources=[source])` for fixed,
  spatially distributed permanent/given magnetization.  The source owns a
  separate HDiv space, is L2-projected once, and exposes `source.field_cf` as
  a native NGSolve CoefficientFunction plus `source.Field(points)` for batch H.
  Its charge geometry is materialized in C++ without building a Gram H-matrix.
  Separate PM and iron spaces preserve their normal jump even at a touching
  interface.  Multiple sources superpose; their coefficients are not solve
  unknowns.  This is a fixed-M source, not a recoil/nonlinear PM law.  Planar
  2D continues to use `magnets=[(mesh, M), ...]`.  Split physically distinct
  magnet segments into separate sources when `M_given` has an internal normal
  jump; one conforming source space enforces normal continuity within itself.
- `radia.vim.EnergyStopMaterial(eta, g_tables, alpha=..., gamma=...,
  b_max=...)` followed by
  `solver = radia.vim.HDivSolver(pm_mesh, order=1|2)` and
  `solver.SolveHysteresis(h_steps, material=...)` for an evolving
  permanent-magnet state.  Reuse the same solver with the returned `state`
  for continuation.  This is the C++
  isotropic vector B-input Stop law: branch states live in fixed balls,
  monotone radial `g_k` tables define convex branch energies, trial evaluation
  is pure, and commit occurs only after the HDiv step converges.  `gamma=0`
  gives the hard Stop projection; positive gamma uses the radial variational
  proximal update.  `h_steps` accepts uniform vectors or 3D NGSolve
  CoefficientFunctions.  Use `initial_b_path` for an explicit manufacturing
  history and pass the returned `state` back as `initial_state` to continue.
  Public Radia provides the generic kernel, not fitted proprietary magnet-grade
  tables.  The returned final state owns a persistent C++ field evaluator, so
  `FieldFromSolution(result, points)` provides its demagnetizing field.  This
  also feeds `SolveCoupledHysteresis([history_pm, ...], [iron, ...], h_steps)`:
  each PM trial restarts from its own committed state and all states commit
  together only after the coupled fixed point converges.  Do not confuse this
  with fixed-M `MagnetizationSource`.
- `radia.vim.PlayHysteresisMaterial(K, eta, f_k_tables)` implements the
  simplified engineering Play level through the same `SolveHysteresis`
  stepping protocol.  It carries branch history, but it is not the full
  EnergyStop claim.  Level 4 adds fixed-domain vector Stop states, convex
  branch energy/proximal updates, irreversible-demagnetization gates, and
  explicit manufacturing/restart state.
- `radia.vim.DemagOperator(HDiv(mesh, order=1|2), ...)` for the NGSolve-style
  diagnostic operator, or `radia.vim.ChargeGram(...)` for its charge map and
  C++ H-matrix components
- `radia.vim.H1HodgeDemagOperator(hdiv, h1, definedon=...)` for an independent
  finite-domain `C.T K^-1 C` contraction/reference operator.  It can be passed
  to `NgsolveHDivMMMResponseReduction(..., demag_operator=...)`; it is not an
  open-boundary substitute for `DemagOperator`.
- `radia.vim.FieldFromSolution(res, points)` -- batch demagnetizing H (A/m) at
  points from the full BDM1/BDM2 TET/HEX/WEDGE solution or the planar BDM1/BDM2
  solution.  `rad.Fld` on a solved mesh-backed
  object dispatches to this same evaluator; per-element constant-M write-back
  is metadata/visualization only, not the field oracle.  Solve materializes an
  immutable C++ source evaluator once.  TET keeps analytic volume/triangle
  near kernels; HEX/WEDGE/curved sources keep an NGSolve quadrature cloud.  Calls
  pass contiguous NumPy target arrays without rebuilding source lists, and all
  IMA terms are accumulated in one TaskManager region.  Ordinary work uses the
  exact direct sum.  Very large non-IMA maps may use a quadrupole source tree
  only after representative direct probes satisfy the configured error bound
  and show a measured speed benefit.  IMA auto remains direct to preserve the
  reduced/full roundoff contract.  `algorithm="direct"` is available on
  `FieldFromSolution` for strict validation.

Core pieces:

- `src/radia/vim/` declares NGSolve spaces/forms, handles dispatch and material
  laws, and prepares the one-time sparse charge topology.
- `src/core/rad_hdiv_vim.*` contains structured and unstructured HDiv assembly
  helpers.
- `src/core/rad_hdiv_hysteresis.*` contains the TaskManager-parallel energy
  Stop trial/commit kernel and its explicit restart state.
- `src/core/rad_hdiv_field_evaluator.*` contains the persistent direct/tree BDM1
  field source and NumPy-facing batch evaluator.
- `src/core/rad_hacapk_hdiv.*` contains `_ChargeGramHMatrix`, the C++ H-matrix
  backend for the Coulomb Gram.
- `src/radia/planar_geometry.py`, `planar_materials.py`, `planar_charges.py`,
  `planar_hysteresis.py`, and `planar_aniso.py` provide 2D planar shared
  geometry/material helpers.

TaskManager is assumed.  NGSolve assembly should run under
`with ngsolve.TaskManager():`, and the C++ kernels use parallel loops for
charge gather, dot products, preconditioner/vector updates, and sparse scatters.
The assembled NGSolve mass matrices are extracted directly in pybind; the
persistent C++ operator owns B/BT, mass CSR, Krylov iterations, and the immutable
field source.  NumPy appears only at vector/target API boundaries, not inside the
iteration loop.  This is the same Python-declaration/C++-execution split used by
NGSolve rather than an attempt to reimplement FESpace construction in Radia.
"""

_SCALING = r"""
# Scaling

The costly object is the charge Gram G.  Radia builds it through HACApK as a
charge H-matrix, then applies the material operator as B^T G B without
materializing a dense N for production runs.

Record these quantities in validation and benchmark artifacts:

- number of magnetic elements;
- H(div) unknown count and charge count;
- H-matrix build time and compression;
- solve iterations and residual;
- peak memory when available;
- field source count, evaluator build time, observation count, selected
  direct/tree route, and direct-reference field error;
- machine label (`LAB` smoke vs `mdx` validation).

Small problems are allowed to be simply "interactive".  The scaling question
matters at engineering size, where charge count and matrix build dominate.
Timing claims should be taken on mdx when it is idle.
"""

_VERIFICATION = r"""
# Verification

Fast tests should cover API contracts and small deterministic checks:

- backend selection rejects unsupported mesh-less soft iron;
- pure TET/HEX/WEDGE mesh-backed soft iron dispatches to HDiv;
- `rad.Fld` after `rad.Solve(..., image=...)` matches an explicitly mirrored
  full model for truly symmetric meshes to near roundoff;
- solve results already own `_field_evaluator`, repeated field calls reuse it,
  and large non-IMA auto-tree output stays within its direct-probe contract;
- 2D planar helpers preserve material labels and PM source regions;
- public solver names and config keys match the current API.
- BDM2 TET/WEDGE and affine HEX linear/nonlinear solves, IMA, and persistent
  field reconstruction remain consistent with their analytic and image gates;
  mapped/non-affine HEX BDM2 material solves must fail loudly.  Planar BDM2 is
  supported through Q3 geometry.
- generic H1-Hodge response reduction must solve both whole-domain and
  material-restricted HDiv spaces, preserve the `[0,1]` contraction, and feed a
  finite-residual HCurl eddy-bubble mixed solve without promoting finite-domain
  H1 or sampled interaction as an open-boundary accuracy oracle.
- local-polynomial HDiv residual probes must remain training-only, reject use as
  physical excitation, retain the 30-probe degree-two contract, and preserve
  the frozen reduced/full-parent, p-saturation, A-phi, and energy-balance gates.

Validation-class tests live under `validation_test/feec/` and should cover:

- sphere/cube demag factors and convergence trends;
- nonlinear BH curves with convergence metadata;
- IMA/image symmetry for TET/HEX/WEDGE;
- curved and high-order geometry where analytic demag truth exists;
- reduced-FEM handoff through NGSolve fields/CoefficientFunctions.
- prescribed-M source L2 projection, immutable coefficients, direct/native-CF
  field equality, multiple-source superposition, and iron response against the
  same field passed explicitly as `H_ext`.
- energy-Stop table convexity guards, hard projection, positive-gamma proximal
  stationarity, non-negative vector-loop dissipation, reverse-field remanence
  loss, arbitrary CoefficientFunction loading, and split-run restart parity.
"""

_NONLINEAR = r"""
# Nonlinear Material Solve

For BH curves, the HDiv route updates a material state on the mesh and solves
the demag equation with a robust nonlinear iteration.  Engineering defaults:

- use tolerances that match observable mesh/discretization error;
- record convergence status, iteration count, and max update;
- fail loudly on non-convergence;
- keep the solve under TaskManager.

Deep saturation can require safeguarded nonlinear steps.  Do not judge the
method from a single scalar residual; inspect the material-state update and the
field observable used by the application.
"""

_CURVED = r"""
# Curved And High-Order Geometry

Curved geometry is one of the main reasons Radia keeps soft iron in the HDiv /
NGSolve lane.  `mesh.Curve(p)`, Piola mappings, curved boundary normals, and
high-order integration are shared with the reduced FEM side.

Good validation cases:

- sphere and spheroid demag factors against analytic values;
- curved sphere external field against the exact dipole;
- curved high-order element convergence compared with flat low-order faceting;
- GMSH/Netgen export checks when the mesh originates from Cubit.
"""

_SYMMETRY = r"""
# Image Symmetry

Image symmetry is part of the HDiv field contract.  A reduced model and an
explicit full model should agree below `10 eps` relative error when the mesh is
geometrically and topologically symmetric.  Percent-level agreement is a warning sign for
asymmetric mesh cuts, incorrect image signs, wrong materialization of images, or
quadrature/charge-basis mismatch.

For reduced models, record:

- image string and reflected axes;
- real charge count and image count;
- whether `rad.Fld` was evaluated through the reduced model or a materialized
  full model;
- max/mean field difference at probes.
"""

_CROSS_METHOD = r"""
# Cross-Method Checks

Prefer analytic truth first: ellipsoid demag factors, cuboid permanent-magnet
fields, dipole limits, and closed-form thin/axisymmetric cases.  When analytic
truth is unavailable, use independent formulations:

For the linear-recoil permanent-magnet level, a sphere has `N=I/3` and the
exact vector load line

    M = (B_r/mu0 + (mu_rec-1)*H_ext) / (1 + (mu_rec-1)/3).

The curved-sphere validation must remain green before documenting a new recoil
material or spatial-remanence workflow.

- HDiv VIM on the NGSolve mesh;
- volume FEM A/phi or reduced-potential FEM where appropriate;
- boundary-element single-layer checks for surface-charge problems;
- direct full-model image materialization for symmetry tests.

Do not put local third-party provenance into public artifacts.  Public docs
should state the analytic convention and the reproduced number, not internal
comparison file names.
"""

_REFERENCE_AUDIT = r"""
# Reference Audit Ladder

When a disagreement appears:

1. Inspect mesh materials, boundaries, and finite-element spaces before solving.
2. Verify image signs and whether the mesh cut is exactly symmetric.
3. Compare charge maps and field evaluation before nonlinear iteration.
4. Check the same observable through two evaluators (`M_avg`, `rad.Fld`, probe
   grid, energy) before changing solver tolerances.
5. Move heavy sweeps to mdx and label the result as validation, not LAB smoke.
"""

_STATUS = r"""
# Status

Current direction:

- Radia soft iron: HDiv-VIM.
- BDM1/BDM2: pure TET/HEX/WEDGE operator, IMA, and persistent-field paths.
  Material solves gate mapped/non-affine HEX BDM2 while retaining mapped HEX
  BDM1, affine HEX BDM2, and pure-TET BDM2.  Planar BDM1 is supported through Q2 and planar BDM2
  through Q3.  `radia.vim.hdiv_capabilities()` is the sole order-pair table;
  geometry order is not inferred from one global p+1 rule.
- Fixed/given 3D magnetization: source-owned HDiv projection and native C++
  field coupling for BDM1/BDM2 TET/HEX/WEDGE, including Curve(2)
  and IMA; planar 2D uses its existing `magnets=` source path.
- Linear-recoil permanent magnet: scalar recoil permeability plus constant or
  spatial `B_r`, solved by the symmetric C++ HDiv path in 3D and planar 2D.
- Simplified history level: PlayHysteresisMaterial plus SolveHysteresis.
- Evolving 3D permanent magnet: C++ vector B-input EnergyStopMaterial plus
  persistent HDivSolver.SolveHysteresis, explicit manufacturing/restart state,
  and arbitrary applied CoefficientFunction steps on BDM1/BDM2 TET/HEX/WEDGE.
- Coupled evolving magnets: one or more CoupledHistoryBody objects plus
  independent recoil/linear/nonlinear bodies, each with one persistent Gram.
- Planar 2D support: HDiv/planar shared geometry and material helpers.
- Public docs: result-bearing HDiv notebooks plus synchronized JSON.
- MCP: teach the live HDiv API and reduced-FEM coupling path.

Open work:

- extend 2D and 3D validation coverage around `rad.Fld`;
- keep the image-symmetry roundoff contract green as hex/wedge coverage grows;
- implement a composite mapped-HEX BDM2 charge operator that preserves
  volume/surface cancellation before removing the material-solve gate;
- continue BDM1/BDM2 TET/HEX/WEDGE accuracy, memory, and timing measurements on mdx;
- continue charge-Gram H-matrix performance checks on mdx;
- run `validation_test/feec/bench_hdiv_field_evaluator_scaling.py` after a
  normal release to measure public `rad.Fld` on mdx/hibino;
- keep Cubit/GMSH mesh-export artifacts aligned with the HDiv API.
- extend application-level force validation while keeping the existing
  curved motor Maxwell/volume/coenergy torque agreement green.
"""

_SECTIONS = {
    "family": _FAMILY,
    "overview": _OVERVIEW,
    "eddy_bubble": _EDDY_BUBBLE,
    "implementation": _IMPLEMENTATION,
    "scaling": _SCALING,
    "verification": _VERIFICATION,
    "nonlinear": _NONLINEAR,
    "curved": _CURVED,
    "symmetry": _SYMMETRY,
    "cross_method": _CROSS_METHOD,
    "reference_audit": _REFERENCE_AUDIT,
    "status": _STATUS,
}


def get_hdiv_vim_documentation(topic: str = "overview") -> str:
    """Return HDiv-VIM knowledge for a topic; use 'all' for every section."""
    t = (topic or "overview").strip().lower()
    if t == "all":
        return "\n\n".join(_SECTIONS[key] for key in _SECTIONS)
    if t in _SECTIONS:
        return _SECTIONS[t]
    return (
        f"Unknown topic '{topic}'. Options: "
        + ", ".join(_SECTIONS.keys())
        + ", all.\n\n"
        + _OVERVIEW
    )
