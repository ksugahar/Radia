# HDiv-VIM Production Note

Radia soft iron is standardized on HDiv-VIM.  This note is the current
production checklist, not a migration archive.

NGSolve's bare `HDiv(mesh, order=p)` is BDM.  Accordingly, every established
Radia production path described below uses BDM1/BDM2.  Actual Raviart--Thomas
requires the explicit `RT=True` flag and is not implied by the `HDiv` name.

## Done

- Mesh-backed `radia.vim.MeshSoftIron` integrates with `rad.Solve`.
- The charge-Gram backend is C++/HACApK based and TaskManager-aware.
- `vim.HDivSolver` is the persistent 3D API for load sweeps, history, and
  coupled bodies.  It owns the HDiv space and ChargeGram; result mappings do
  not expose or transport geometry-cache dictionaries.
- TET/HEX/WEDGE and 2D planar paths have validation coverage under
  `validation_test/feec/`.
- `rad.Fld` is part of the public contract after HDiv write-back.
- BDM1/BDM2 solve results own a persistent C++ field evaluator: NumPy target buffers,
  one-pass IMA, TaskManager observation parallelism, analytic tet near kernels,
  and a guarded large-map quadrupole tree.  IMA auto evaluation remains direct
  to preserve the reduced/full roundoff contract.
- On geometrically and topologically symmetric full/reduced hex meshes,
  `rad.Fld` image parity is locked to the explicit full solve at the
  roundoff-level contract (`< 10 eps` relative error).
- BDM1 and BDM2 are public for pure TET/HEX/WEDGE operators, IMA, persistent
  field evaluation, and the
  NGSolve `ChargeGram`/`DemagOperator` surface.  Planar BDM1 is public through
  Q2 geometry and planar BDM2 through Q3, including IMA and the persistent
  planar field evaluator.
- Material solves reject mapped/non-affine HEX BDM2 until one composite
  volume/surface operator preserves the discrete cancellation.  The production
  alternatives are mapped HEX BDM1, affine HEX BDM2, or pure-TET BDM2.
- An independent finite-domain H1 Omega gate uses `N = C.T K^-1 C` on the
  same 207 active mapped-body BDM2 DoFs.  Its order-2/order-3 Hodge spectra are
  contractions, while the current charge diagnostic has modes above one.
  Therefore the gate is localized to the mapped charge operator; it is not a
  rejection of the BDM2 approximation space.
- `H1HodgeDemagOperator` exposes that finite-domain operator as an NGSolve
  `BaseMatrix` path.  `NgsolveHDivMMMResponseReduction` accepts it directly,
  evaluates generic CG results into concrete vectors, and restricts the mass
  preconditioner with `fes.FreeDofs()` for body-only spaces.  On the 54-HEX
  air/body gate, the resulting BDM2 reduction feeds HCurl eddy-bubble with
  snapshot residual below `3e-12`, mixed residual below `3e-16`, and positive
  Joule loss.  This is a mixed-mechanics/contraction gate, not an open-boundary
  replacement for the guarded mapped-HEX BDM2 charge solve.  The one-call path
  uses `NgsolveBDMEddyBubbleVIM(..., hdiv_definedon=...,
  demag_operator_factory=...)`, ensuring the operator is built on the exact
  restricted BDM space rather than a separately reconstructed space.
- A separate heavy static-disk accuracy artifact uses an 18432-element
  axisymmetric Q2 reference (`1.0916977441`).  The 3-D finite-domain H1-Hodge
  error decreases through 4.65% (coarse BDM1/H1-p2), 2.24% (coarse
  BDM2/H1-p3), 1.44% (fine BDM2/H1-p3), and 0.98% (fine BDM2/H1-p4).  This
  promotes a strict h/p accuracy lane for this disk without claiming exact
  open-boundary H1 behavior or universal superiority.
- `CoupledHDivHybridVIMSystem.solve_frequency_local_esim` wraps the local
  nonlinear ESIM Karl iteration around the complete HDiv-MMM/HCurl-VIM mixed
  solve.  On a regenerated TET cube whose single `iron` label is both magnetic
  and conductive, the 50/1000/5000 A/m ladder converges in 2--5 updates,
  preserves a passive surface Gram and positive Joule loss, and reproduces the
  converged mixed solution under fixed-Gram replay to machine precision.
  This promotes nonlinear skin coupling, not simultaneous bulk nonlinear B-H;
  the latter remains an explicit open item.
- `radia.vim.hdiv_capabilities()` is the sole field/geometry-order table.  Do
  not derive geometry order from HDiv order with one cross-dimensional p+1
  rule; Policy 8 rejects duplicated arithmetic guards.
- NGSolve spaces, coefficient functions, and bilinear forms are declared in
  Python, as in NGSolve itself.  Their assembled sparse matrices pass directly
  to C++ without SciPy/list materialization.  The persistent C++ operator owns
  B/BT, the geometric and material mass matrices, Krylov iterations, and field
  source evaluation; Python is not in the per-iteration solve/field path.
- NGSolve owns element orientation, Piola maps, curved geometry, local/global
  DOF transforms, quadrature, and weak assembly.  BDM2 hysteresis therefore uses
  `IntegrationRuleSpace` interpolation and mixed weak forms instead of a
  Python reconstruction of the physical high-order basis.
- Production FAR Gram blocks average both directed quadrature evaluations.
  Upper-triangle mirroring alone makes the stored matrix symmetric, but a
  one-sided finite quadrature rule is not invariant under an explicit reflected
  mesh.  One-sided environment switches are diagnostic-only.
- MCP `hdiv_vim` documents the live API and reduced-FEM coupling policy.
- `vim.MagnetizationSource(mesh, M_given, order=1|2)` L2-projects fixed/given
  3D magnetization into an independent HDiv space, builds only its C++ charge
  geometry/evaluator (no Gram H-matrix), and exposes a native NGSolve H-field
  CoefficientFunction.  `vim.Solve(..., magnetization_sources=[...])` assembles
  that field into the iron weak form; source coefficients stay immutable and
  PM/iron normal jumps remain representable because their spaces are separate.
- Prescribed sources cover BDM1/BDM2 TET/HEX/WEDGE, Curve(2), and IMA.
  Planar 2D continues to use `magnets=[(mesh, M), ...]`.
- The permanent-magnet material ladder is one canonical four-level contract:
  (1) fixed/given `MagnetizationSource`; (2) linear recoil
  `vim.Solve(mu_r=mu_rec, B_r=B_r, ...)`; (3) simplified
  `PlayHysteresisMaterial`; and (4) full B-input `EnergyStopMaterial`.  Do not
  add parallel backend names for these levels.
- The linear-recoil path accepts a constant or spatial NGSolve `B_r` in tesla
  and applies `M = (mu_rec-1) H + B_r/mu0` through the existing symmetric C++
  HDiv solve.  It supports production 3D element/curve/IMA routes and planar
  2D at their normal BDM order.  A scalar `mu_rec > 1` is required; fixed-M is
  the explicit `mu_rec=1` limit.
- A spatial `B_r` is valid within one conforming magnet body.  Do not merge
  normal-discontinuous recoil segments into one HDiv space: that removes their
  internal surface charge.  Fixed segments use separate
  `MagnetizationSource` spaces; mutually coupled recoil segments belong to the
  multi-body block-coupling increment.
- `vim.PlayHysteresisMaterial` remains the simplified engineering history
  level.  It must not be documented as the full energy-based irreversible-PM
  model.
- `vim.EnergyStopMaterial` provides the C++/TaskManager vector B-input Stop law
  for an evolving permanent-magnet state.  Shape tables are checked for the
  monotonicity needed by the convex branch energy; trial evaluation is pure,
  state advances only after convergence, and `vim.SolveHysteresis` returns an
  explicit restart state.  Applied steps may be uniform vectors or arbitrary
  three-dimensional NGSolve CoefficientFunctions.  The final state also owns
  the persistent C++ field evaluator used by `vim.FieldFromSolution`.
- The irreversible-PM validation locks hard projection, positive-gamma
  proximal stationarity, non-negative vector-loop dissipation, reverse-field
  remanence loss, and split-run restart agreement.
- `validation_test/hysteresis/test_linear_recoil_permanent_magnet.py` locks the
  level-2 curved-sphere load line against its analytic `N=I/3` solution.
- `vim.SolveHysteresis(..., order=2)` keeps constitutive state at NGSolve
  integration-rule points and records that layout/order in restart state.  BDM1
  retains the established element-average state contract.
- `vim.CoupledBody` and `vim.SolveCoupled` couple independent linear-recoil PM,
  linear iron, and nonlinear iron spaces.  Persistent C++ field
  CoefficientFunctions exchange body fields, each ChargeGram is built once,
  and nonconvergence is fail-loud.  This is also the production segmented-PM
  path for normal-discontinuous `B_r`.
- A `rad.Solve` container with multiple eligible `vim.MeshSoftIron` bodies
  dispatches to the same coupled solve.  Per-body magnetization is written back
  independently, while top-level `rad.Fld` sums every persistent BDM1/BDM2 field
  and the ordinary Radia source objects.
- `vim.CoupledHistoryBody` and `vim.SolveCoupledHysteresis` couple one or more
  stateful EnergyStop/Play PMs to independent linear/nonlinear HDiv bodies.
  Every outer trial restarts each PM from its own committed history state; all
  states commit together only after global convergence.  Every body's
  `HDivSolver` builds its ChargeGram once across outer iterations and physical
  history steps.
- The planar motor application locks curved Q2+BDM1 and Q3+BDM2 against the
  analytic ellipse torque.  Maxwell stress, magnetization-volume torque, and
  fixed-current coenergy agree, and Q3+BDM2 reduces the same-mesh-family error.

## Release Gate

Before release or `mdx`/`hibino` deployment:

- run focused HDiv smoke tests on LAB/100号機;
- run heavy validation/benchmark sweeps on an idle `mdx` or `hibino` host
  (`mdx` by default, `hibino` for MATLAB, large-memory, long-running, or
  mdx-occupied jobs);
- record the actual validation host in the result JSON/log;
- record charge count, HDiv DoF, H-matrix compression, build time, solve time,
  iteration count, and machine label;
- record field source count, evaluator build time, selected direct/tree route,
  observation count, direct-reference error, and public `rad.Fld` wall time;
- for prescribed magnetization, record its L2 projection residual, confirm no
  Gram H-matrix was built, compare native-CF and direct fields, and verify the
  source coefficients are unchanged after the iron solve;
- for linear-recoil PM runs, record the recoil `mu_r`, the `B_r` distribution,
  and the applied-field path; check the analytic ellipsoid/sphere load line
  before using a design geometry;
- for history-dependent PM runs, record the constitutive table/version,
  initialization path, reverse-field path, final explicit state, and unloaded
  remanence; never label a synthetic table as a calibrated magnet grade;
- verify image symmetry with an explicit full model when the mesh is truly
  symmetric and enforce the roundoff contract;
- keep public docs free of obsolete backend names and local validation
  provenance.

## Open Work

- add the outer constitutive update that changes the ordinary nonlinear bulk
  HDiv B-H operator and local ESIM Gram together, with energy and line-search
  safeguards;
- extend force tests around motor workflows; torque/coenergy are locked;
- keep Cubit/GMSH mesh export aligned with the HDiv API;
- continue mdx scaling measurements for charge-Gram build and solve time.
