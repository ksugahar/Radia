# HDiv-VIM Production Note

Radia soft iron is standardized on HDiv-VIM.  This note is the current
production checklist, not a migration archive.

## Done

- Mesh-backed `radia.vim.MeshSoftIron` integrates with `rad.Solve`.
- The charge-Gram backend is C++/HACApK based and TaskManager-aware.
- TET/HEX/WEDGE and 2D planar paths have validation coverage under
  `validation_test/feec/`.
- `rad.Fld` is part of the public contract after HDiv write-back.
- RT1/RT2 solve results own a persistent C++ field evaluator: NumPy target buffers,
  one-pass IMA, TaskManager observation parallelism, analytic tet near kernels,
  and a guarded large-map quadrupole tree.  IMA auto evaluation remains direct
  to preserve the reduced/full roundoff contract.
- On geometrically and topologically symmetric full/reduced hex meshes,
  `rad.Fld` image parity is locked to the explicit full solve at the
  roundoff-level contract (`< 10 eps` relative error).
- RT1 and RT2 are public for flat and isoparametric-P2 pure TET/HEX/WEDGE
  linear/nonlinear material solves, IMA, persistent field evaluation, and the
  NGSolve `ChargeGram`/`DemagOperator` surface.  Planar RT1 is public through
  Q2 geometry and planar RT2 through Q3, including IMA and the persistent
  planar field evaluator.
- `radia.vim.hdiv_capabilities()` is the sole field/geometry-order table.  Do
  not derive geometry order from HDiv order with one cross-dimensional p+1
  rule; Policy 8 rejects duplicated arithmetic guards.
- NGSolve spaces, coefficient functions, and bilinear forms are declared in
  Python, as in NGSolve itself.  Their assembled sparse matrices pass directly
  to C++ without SciPy/list materialization.  The persistent C++ operator owns
  B/BT, the geometric and material mass matrices, Krylov iterations, and field
  source evaluation; Python is not in the per-iteration solve/field path.
- NGSolve owns element orientation, Piola maps, curved geometry, local/global
  DOF transforms, quadrature, and weak assembly.  RT2 hysteresis therefore uses
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
- Prescribed sources cover RT1/RT2 TET/HEX/WEDGE, Curve(2), and IMA.
  Planar 2D continues to use `magnets=[(mesh, M), ...]`.
- The permanent-magnet material ladder is one canonical four-level contract:
  (1) fixed/given `MagnetizationSource`; (2) linear recoil
  `vim.Solve(mu_r=mu_rec, B_r=B_r, ...)`; (3) simplified
  `PlayHysteresisMaterial`; and (4) full B-input `EnergyStopMaterial`.  Do not
  add parallel backend names for these levels.
- The linear-recoil path accepts a constant or spatial NGSolve `B_r` in tesla
  and applies `M = (mu_rec-1) H + B_r/mu0` through the existing symmetric C++
  HDiv solve.  It supports production 3D element/curve/IMA routes and planar
  2D at their normal RT order.  A scalar `mu_rec > 1` is required; fixed-M is
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
  integration-rule points and records that layout/order in restart state.  RT1
  retains the established element-average state contract.
- `vim.CoupledBody` and `vim.SolveCoupled` couple independent linear-recoil PM,
  linear iron, and nonlinear iron spaces.  Persistent C++ field
  CoefficientFunctions exchange body fields, each ChargeGram is built once,
  and nonconvergence is fail-loud.  This is also the production segmented-PM
  path for normal-discontinuous `B_r`.
- A `rad.Solve` container with multiple eligible `vim.MeshSoftIron` bodies
  dispatches to the same coupled solve.  Per-body magnetization is written back
  independently, while top-level `rad.Fld` sums every persistent RT1/RT2 field
  and the ordinary Radia source objects.
- `vim.CoupledHistoryBody` and `vim.SolveCoupledHysteresis` couple one stateful
  EnergyStop/Play PM to independent linear/nonlinear HDiv bodies.  Every outer
  trial starts from the same committed history state; only the converged
  all-body trial commits.  The history and ordinary-body ChargeGrams are cached
  across outer iterations and physical history steps.

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

- extend force/energy tests around motor workflows;
- keep Cubit/GMSH mesh export aligned with the HDiv API;
- continue mdx scaling measurements for charge-Gram build and solve time.
