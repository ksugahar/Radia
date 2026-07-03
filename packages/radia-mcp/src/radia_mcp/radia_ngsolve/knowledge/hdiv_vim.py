"""HDiv-type VIM (Volume Integral Method) demag operator -- knowledge module.

SHOWCASE NOTEBOOK: docs/hdiv_vim/hdiv_curved_showcase.ipynb -- curved-geometry accuracy win: stray-field ~30-40x, nonlinear ~23x, beats shipped solver.
SHOWCASE NOTEBOOK: docs/hdiv_vim/polynomial_charge_field.ipynb -- the order>=2 charge-field kernel (H = -grad phi_M from rho=-div M + sigma=M.n), executed: degree-0/1/2 closed forms + arbitrary-degree assembler + affine-hex + curved tri/tet, each validated live vs Gauss/Duffy to machine precision (drop-rho = 90-230% error; uniform sphere -> -M/3).
VALIDATION CORPUS: examples/vim plus validation_test/feec.  Do not mirror the
full source corpus into docs; public docs carry the maintained references and
result-bearing showcase notebooks, while source history stays in git.

The HDiv-type VIM is the lab's FEEC (H(div) RT) alternative/complement to the canonical multipole-moment MMM MSC
kernel: a SYMMETRIC demag operator N = B^T G B whose loop modes are FIELD-NULL BY CONSTRUCTION, giving
mu_r-INDEPENDENT convergence with no hand-crafted loop-star.  Validated on: linear demag
(sphere/spheroid/triaxial exact vs analytic), NONLINEAR (damped Newton; cube & C-yoke <1-3% vs shipped
Radia), distorted-mesh mu_r-independence, CURVED + high-order (accuracy-per-DOF ~10-30x vs flat Radia),
and SYMMETRY models 1/2,1/4,1/8 (loops automatic + image-method demag).  Canonical reference:
docs/hdiv_vim/README.md.

Research positioning (2026-06-24): HDiv-VIM is the clean symmetric/high-order Galerkin answer, but its
charge-Coulomb Gram construction is expensive.  The production surface-charge path therefore uses the
Mathematica-derived multipole-moment MMM rows: cheaper local moment functionals for 3-DOF MMM and 5/6-DOF
MSC, with HDiv retained as the higher-order and de-Rham-exact complement.

PRIMARY (decision 2026-06-30, Sugahara -- SUPERSEDES the 2026-06-24 positioning above): HDiv-VIM is now the PRIMARY (本命) accurate soft-iron demag method, and collocation MMMM is DEMOTED to the COARSE / fast tier (optimization inner loops, mesh-less quick passes). The 2026-06-24 'production uses the multipole-moment MMM rows' framing is REVERSED: collocation MMMM gave up loop-free (its loop-free implementation was removed 2026-06-30) -- field-correct (loops field-null) but loop-polluted internal M, acceptable for coarse/optimization but NOT accurate/hysteresis; HDiv-VIM is loop-free BY CONSTRUCTION (loops = ker(B)), so it is the primary accurate route. Use HDiv-VIM (tet mesh) for production/accurate + hysteresis; collocation MMMM (hex/mesh-less) for fast coarse passes. Memory: collocation_loopfree_abandoned.

CURRENT API (2026-06-23 -- the dense Python Gram path was REMOVED): the C++ `_ChargeGramHMatrix` kernel
is the SOLE demag operator (N v = B^T (H.matvec(B v)); EXACT analytic near AND far; tet via
cell_verts/face_verts, hex/wedge via the polytope triangle soup), and folds IMA (image charges /
symmetry models) IN via image_masks / image_signs.  Production entry:
`radia.vim.hdiv_demag_solve(mesh, mu_r=/bh_table=, H_ext=, image=)` (linear or nonlinear, no Gram
switch); `build_demag(mesh, nsub)` returns only SPARSE pieces (charge map B + HDiv mass + geometry),
NO dense N/G.  Several sections below describe the OLD dense-Python call options (`analytic_gram=`,
`wilton_surface=`, `skip_dense_gram=`, `build_near_correction`, the monopole+near-correction split) --
those NO LONGER EXIST; read them as the research history behind the always-exact C++ charge Gram.

TaskManager is assumed: shared HACApK build paths and long C++ solve loops stand up or reuse NGSolve
RegionTaskManager; direct diagnostic `.matvec()` calls plus Python/NGSolve assembly follow the
caller-wraps `with ng.TaskManager():` convention.  The C++ HDiv CG/MINRES/Picard kernels use
ParallelFor/ParallelForRange for charge gather, dot products, preconditioner/vector updates, and
AtomicAdd for sparse face-vector scatters.

Exposed as the radia-ngsolve MCP tool `hdiv_vim(topic=...)`: overview, implementation, scaling,
verification, nonlinear, curved, symmetry, status, all.
"""

_OVERVIEW = r"""
# HDiv-type VIM demag operator (the FEEC alternative to MMM/MSC)

## What it is
Magnetization M lives in lowest-order H(div) (RT0): ONE signed normal-flux DOF per face.  The demag
operator is the SYMMETRIC Galerkin form

    N = B^T G B          (symmetric, [n_face x n_face])

  B = charge map   M |-> (rho = -div M per cell [P0],  sigma = M.n per boundary face [P0])
  G = Coulomb Gram (charge-charge interaction; symmetric since 1/r is symmetric)

## Why it complements multipole-moment MMM MSC
- **Loops are FIELD-NULL BY CONSTRUCTION**: loops = ker(B) (charge-free fields).  B.loop = 0 =>
  N.loop = B^T G (B loop) = 0 for ANY G.  So the loop (circulating-magnetization) null space sits at
  EXACTLY the material eigenvalue 1/chi and never pollutes the spectrum.
- => **mu_r-INDEPENDENT convergence**: plain-Jacobi MINRES iteration count is FLAT across mu_r
  (measured 80/82/78 across mu_r 10 -> 1e3 -> 1e5 on the 4^3 system).  The high-mu_r conditioning
  wall that affects surface-charge MSC iterations is ABSENT.
- **DISTORTED-ELEMENT ROBUSTNESS (the headline PRACTICAL value, 2026-06-07)**: the field-null-by-
  construction property holds on ANY mesh -- B.loop = 0 is exact on AFFINE *and* NON-AFFINE (distorted)
  hexes alike, so mu_r-independence SURVIVES distortion (measured: distort=0.18 grid, MINRES iters
  98/56/52/52 across mu_r 10 -> 1e4 -- bounded, even DECREASING).  The multipole-moment MMM constant-face
  loop patterns are field-null ONLY on affine hexes (the de-Rham defect: on distorted hexes they carry field);
  the hand-crafted surface-charge elements / the shipped MSC's installCycle retrofit (~6e-9 local-null-vector)
  exist precisely to patch this.  HDiv-VIM needs NO hand-crafted elements and NO retrofit -- it is
  EXACT (machine, ~4e-16) on distorted meshes by construction.  This is a robustness / generality /
  maintainability win on its OWN (any mesher, any distortion, provably correct; see the README
  productionization roadmap), independent of raw speed -- at mu_r<=1e4 it is at performance PARITY with the
  shipped multipole-moment MMM MSC (same demag spectrum), the win being correctness + no hand-crafting.  Golden-locked:
  tests/feec/test_hdiv_vim_symmetry_golden.py (loops field-null + PSD on distorted) and
  test_hdiv_vim_solve.py::test_minres_iters_bounded_vs_mu_r_distorted (mu_r-independence on distorted).
- **SYMMETRIC** => MINRES (symmetric indefinite Krylov); mu_r-independent, no direct factorization.

## The material system
    A = (1/chi) M_mass - N       (symmetric INDEFINITE)
  M_mass = RT0 H(div) mass.  The generalized eigenvalues of (N, M_mass) are the demag factors,
  bounded in [0,1] (basis-invariant physics).  Solved by MINRES / Jacobi-PCG (mu_r-independent).

## Origin
NGSolve FEEC prototype (examples/vim/hdiv_demag_quad_self.py) -> productionised into the Radia
C++ core.  This is the "de-Rham fix" research line: loops field-null by construction on distorted
hexes AND tets.
"""

_IMPLEMENTATION = r"""
# Current implementation (files + APIs)

[API UPDATE 2026-06-23] The dense Python Gram path -- the `analytic_gram=` / `wilton_surface=` /
`skip_dense_gram=` kwargs, `build_near_correction`, the dense `analytic_charge_gram` /
`wilton_surface_block` / `phi_tet` builders, and the monopole+near-correction split named in this section
-- was REMOVED.  The C++ `_ChargeGramHMatrix` (exact analytic near AND far, IMA via image_masks/image_signs)
is now the SOLE demag operator and `radia.vim.hdiv_demag_solve(mesh, mu_r=/bh_table=, H_ext=, image=)` the
production entry; read the removed names below as the research history behind the always-exact C++ Gram.

## C++ core
- **src/core/rad_hdiv_vim.{h,cpp}** -- structured-hex RT0 topology + charge map B + Coulomb Gram G +
  N = B^T G B (DENSE, hand-enumerated, no NGSolve).  Geometry-exact on distorted hexes via trilinear
  (hex) / bilinear (quad) sub-points.  Reusable helpers (shared by the dense path AND the on-demand
  H-matrix entry, so they agree entry-by-entry):
    - `BuildStructuredRT0(nx,ny,nz,h,distort)` -> Mesh (faces, cells, geometry)
    - `AssembleChargeMap` / `BuildChargeMapCSC` -- B (per-face <=2 charges: lo/hi cell + boundary sigma)
    - `BuildChargeQuad(m, nsub, q)` + `CoulombGramEntry(q, a, b)` -- per-pair Gram (nsub=0 centroid-
      monopole; nsub>=1 accurate sub-point); `AssembleCoulombGram` / `AssembleN` build the dense forms
    - `AssembleMass` (dense) + `BuildMassCOO` (sparse RT0 mass + per-face diagonal)
- **src/core/rad_hacapk_hdiv.{h,cpp}** -- two HACApK H-matrix managers (subclass RadHACApKBase):
    - `RadHACApKHDivManager` -- N as a FACE-based H-matrix (entry = the charge-cluster B^T G B sum);
      + `ApplySystem(x, inv_chi, y)` = inv_chi*(M_mass x) - N x; `DiagSystem(inv_chi)` (Jacobi diag).
    - `RadHACApKChargeGram` -- the CHARGE-based Coulomb Gram G as an H-matrix (the unstructured path):
      entry G[a!=b] = meas_a meas_b/(4pi r) (centroid monopole), G[a][a] = self_energy[a] (caller-
      computed, shape-aware).  N = B^T G B applied as B^T (G-Hmatvec (B m)).
- **src/core/rad_hacapk_hdiv.cpp** -- RadHACApKChargeGram::SolveLinearMaterial (Jacobi-PCG, SPD
  material system) + SolveNonlinearPicard (scalar-chi nonlinear demag), both in C++ -- the
  mu_r-independent iterative solve (the symmetric H-LDL^T was removed 2026-06-08).

## pybind (radia._radia_pybind)
- `_hdiv_vim_assemble(nx,ny,nz,nsub=0,distort=0)` -> dict{nf,n_cell,n_charge,n_bnd, N, B, M_mass}
- `_hdiv_vim_hmatrix_probe(nx,ny,nz,nsub,distort,eps,leaf,eta)` -> stats + matvec_relerr vs dense N
- class `_HDivVimHMatrix(nx,ny,nz,nsub,distort,eps,leaf,eta)` -- build-once: .ndof() / .matvec(x) (N x)
  / .apply_system(x,inv_chi) / .diag_system(inv_chi) / .stats()
- class `_ChargeGramHMatrix(...)` -- monopole ctor (centroids,measures,self_energy) OR analytic ctor
  (cell_verts,face_verts,n_el, M2); .ndof() / .matvec(q) / .stats() / .solve_linear_material(...) /
  .solve_nonlinear_picard(...) (M3, iterative solves in C++)
- `_hdiv_tri_potential(V,r)` / `_hdiv_phi_tet(V,P)` -- analytic Wilton/phi_tet charge-Gram potentials (M2)

## Python (unstructured tet ingest -- the real-geometry path)
- **examples/vim/hdiv_demag_tet.py** -- NGSolve HDiv(order=0) tet ingest.  Element-AGNOSTIC
  extraction: B = -div(u) [vol charge] + u.Trace().n [surface charge]; M_mass = HDiv mass.  Only the
  Gram self-energy geometry is element-specific -> tet/tri barycentric sub-points (c_tet~1.776 /
  c_tri~2.888).  `build_demag(mesh,nsub)` -> N (dense, monopole off-diag + sub-point self), M_mass, B,
  charge geometry.  `build_near_correction(mesh, d, nsub, near_factor)` -> sparse near-field Gram
  correction (exact sub-point MINUS monopole) -- the scalable Gram = monopole H-matrix + this.
"""

_SCALING = r"""
# Scaling: charge-Gram H-matrix

[API UPDATE 2026-06-23] The dense Python Gram path -- the `analytic_gram=` / `wilton_surface=` /
`skip_dense_gram=` kwargs, `build_near_correction`, the dense `analytic_charge_gram` /
`wilton_surface_block` / `phi_tet` builders, and the monopole+near-correction split named in this section
-- was REMOVED.  The C++ `_ChargeGramHMatrix` (exact analytic near AND far, IMA via image_masks/image_signs)
is now the SOLE demag operator and `radia.vim.hdiv_demag_solve(mesh, mu_r=/bh_table=, H_ext=, image=)` the
production entry; read the removed names below as the research history behind the always-exact C++ Gram.

The scalable Gram is the standard H-matrix split:
  - FAR field: compressed CENTROID-MONOPOLE charge H-matrix (RadHACApKChargeGram / _ChargeGramHMatrix),
    O(N log N) matvec.  ACA compresses far blocks (n_lowrank grows with N).
  - NEAR field: a SPARSE correction (exact sub-point MINUS monopole for near charge pairs, O(N) nnz),
    computed in Python (build_near_correction) and applied alongside the H-matvec:
        N m = B^T ( G_mono_Hmatvec(B m)  +  corr @ (B m) )

This lifts the sphere demag from the monopole UNDER-estimate (~0.31) to the Gram-EXACT value (~0.33),
converging to the analytic 1/3 (0.3256 / 0.3279 / 0.3291 at h = 0.5 / 0.35 / 0.25), and the scalable
(H-matvec + sparse corr) result matches the dense reference to 4 decimals.

OPEN (production refinement, not yet done): the near correction is currently computed in PYTHON
(build_near_correction); moving it INTO the C++ RadHACApKChargeGram entry needs CELL/FACE GEOMETRY
(vertices) passed to the manager (it currently holds only centroids/measures/self_energy).  The
exact-near could also use the analytic Wilton single-layer (rad_poly_analytical
RadScalarPotentialFromTriangleFaceGlobal) + the volume-potential reduction INT_V 1/r =
(1/2) SUM_faces d_face(x) * Wilton_face(x) (since lap(r)=2/r) -- exact, no nsub.
"""

_VERIFICATION = r"""
# Verification (golden tests, tests/feec/)

[API UPDATE 2026-06-23] The dense Python Gram path -- the `analytic_gram=` / `wilton_surface=` /
`skip_dense_gram=` kwargs, `build_near_correction`, the dense `analytic_charge_gram` /
`wilton_surface_block` / `phi_tet` builders, and the monopole+near-correction split named in this section
-- was REMOVED.  The C++ `_ChargeGramHMatrix` (exact analytic near AND far, IMA via image_masks/image_signs)
is now the SOLE demag operator and `radia.vim.hdiv_demag_solve(mesh, mu_r=/bh_table=, H_ext=, image=)` the
production entry; read the removed names below as the research history behind the always-exact C++ Gram.
The dense-only tests named below (test_hdiv_vim_tet_hmatrix / tet_nearcorr) were deleted; the structural
goldens (symmetry_golden, solve) now build N from the C++ kernel (conftest.hdiv_vim_dense_N_and_loops).

- test_hdiv_vim_symmetry_golden.py -- N symmetric + loops field-null (1^3..5^3); Gauss rank(B);
  PSD demag factors (accurate Gram); distorted-hex trilinear PSD.
- test_hdiv_vim_hmatrix.py -- N as HACApK H-matrix: EXACT when all-dense (6e-16), matches dense within
  ACA tol on regular+distorted, symmetric, compression grows with N.
- test_hdiv_vim_solve.py -- scalable MINRES of A=(1/chi)M_mass-N: exact-when-dense, matches dense,
  mu_r-INDEPENDENT iters (flat across mu_r 10->1e5).
- test_hdiv_vim_tet.py -- unstructured tet ingest: symmetry+loop-null exact on tets; sphere demag
  converging to 1/3.
- test_hdiv_vim_tet_hmatrix.py -- charge-Gram H-matrix on tets matches dense + demag survives compression.
- test_hdiv_vim_tet_nearcorr.py -- near-field-corrected scalable demag == dense, lifts toward 1/3.

Bug caught by verify-first (2026-06-07, commit 5d7a9823): the Python tet/tri barycentric sub-point
lattice summed to 1+c/nsub (sub-points OUTSIDE the simplex).  Fixed to proper lattices (sum to 1,
nsub=1 -> centroid).  The C++ hex/quad path was unaffected.  The self constants c_tet/c_tri are
robust to the scheme (cross-sum extrapolation -> same INT INT 1/r).
"""

_NONLINEAR = r"""
# NONLINEAR HDiv-type VIM -- ROBUST SOLVER = DAMPED NEWTON-RAPHSON

[API UPDATE 2026-06-23] The dense Python Gram path -- the `analytic_gram=` / `wilton_surface=` /
`skip_dense_gram=` kwargs, `build_near_correction`, the dense `analytic_charge_gram` /
`wilton_surface_block` / `phi_tet` builders, and the monopole+near-correction split named in this section
-- was REMOVED.  The nonlinear solve now ALWAYS uses the exact C++ `_ChargeGramHMatrix` (analytic near AND
far -- so the "needs analytic_gram for div M != 0" trap is GONE, the volume Gram is always exact); the
production entry is `radia.vim.hdiv_demag_solve(mesh, bh_table=, H_ext=, image=)` (and
`solve_nonlinear_newton(mesh, chi0, Msat, H0, ...)`), with NO Gram kwarg.  Read the removed names below as
the research history behind the always-exact C++ charge Gram.

Goal: make the HDiv-type VIM work for NONLINEAR soft-magnetic materials (BH curve / saturation).

## HEADLINE (2026-06-07): the robust strong-saturation solver is DAMPED NEWTON-RAPHSON.
The earlier sessions' simple per-element Picard / Hantila FAILED at deep saturation (NaN, wrong root);
damped Newton-Raphson with the CONSISTENT TENSOR tangent SOLVES it -- fast (4-11 iters) and accurate
(rel 1e-4 to 1e-6 vs the analytic uniform sphere) at every saturation level.  This is the Newton
counterpart of Radia's existing MMM/MSC newton_damping=True path: the OUTER Newton machinery (tangent
susceptibility + line-search damping) is SHARED; only the demag operator N = B^T G B differs.
Implemented as solve_nonlinear_newton (examples/vim/hdiv_demag_tet_nonlinear.py), golden-locked in
tests/feec/test_hdiv_vim_tet_newton.py (3 tests, feec suite 50/50).

## DONE + golden-tested (examples/vim/hdiv_demag_tet_nonlinear.py, tests/feec/test_hdiv_vim_tet_nonlinear.py)
- **Applied-field formulation RESOLVED (verify-first)**: the eigenvalue framing A_eig = (1/chi)M_mass
  - N is NOT the applied-field system.  The physical applied-field weak form is
        A+ m = M_mass h_ext ,   A+ = (1/chi) M_mass + N        (PLUS N)
  (since M = chi(H_ext + H_demag), H_demag,weak = -N m).  For a sphere this reproduces the analytic
  M/H = chi/(1 + chi D) -- VERIFIED <=2.5% for mu_r<=100; the MINUS system gives nonsense (negative/
  divergent).  This is the key gotcha: do NOT reuse the demag-factor A_eig sign for an applied-field solve.
- **BH-curve Picard works**: secant susceptibility chi^{k+1} = M(H_int)/H_int with H_int = H0 - D M_avg;
  CONVERGES (15-25 iters) and SATURATES (M -> M_sat) on the tet sphere.
- **Near-correction closes the accuracy gap**: the centroid-monopole UNDER-estimates D (0.315 vs 1/3),
  so the high-chi nonlinear M is ~12% off analytic; folding in the #3 near-field correction (N_eff =
  B^T (G + corr) B) raises D to 0.328 -> nonlinear M within ~2% of the analytic uniform-sphere answer
  (H0=0.1: monopole 0.334, near-corrected 0.305, analytic(1/3) 0.299).  Residual ~2% is RT0/mesh
  discretization (finer mesh / proper distorted M_mass -> closer).

## PER-ELEMENT machinery DONE + validated at MODERATE drive (2026-06-07)
- Per-element chi_e(H_e): weighted RT0 mass M_chi = BilinearForm(HDiv) with a piecewise-constant
  (1/chi_e) CoefficientFunction, A+ = M_chi + N, RHS = M_mass h_ext.  Per-cell |M_e| via an
  L2(order=0) projection of |gfM|; with the BH curve M(H)=chi0 H/(1+chi0|H|/Msat) the inverse gives
  the clean update chi_e = chi0(1 - |M_e|/Msat).
- VALIDATED: run on the SPHERE (uniform M) the per-element solve reproduces the scalar answer to
  ~0.3% at moderate drive (H0 <= 1e-2).  So the per-element machinery (weighted mass + per-cell field
  + per-element chi) is CORRECT.

## ROBUST SOLVER -- DAMPED NEWTON-RAPHSON (solve_nonlinear_newton, golden-locked)
Newton on the constitutive residual in RT0 coefficient form:
    F(m) = M_mass m - b_M(H) ,   H = h_ext - D_op m ,  D_op = M_mass^{-1} N ,
    b_M(H) = INT M(H).v dx   (L2 projection of the constitutive M onto RT0) ,
with the CONSISTENT TENSOR Jacobian  J = M_mass + T D_op ,  T = INT (dM/dH) u.v dx .

THREE ingredients are all required (each isolated verify-first):
  (a) CONSISTENT TENSOR tangent  dM/dH = chi_diff Hhat(x)Hhat + chi_sec (I - Hhat(x)Hhat)
      (slope along H, secant perpendicular).  The naive SCALAR tangent chi_diff*I STALLS at moderate
      drive: e.g. at |H|~4e-4 the saturating curve has chi_sec=725 but chi_diff=510, so omitting the
      perpendicular secant term gives a badly wrong Jacobian (Newton crawls / converges to M=0.08).
  (b) the #3 NEAR-FIELD CORRECTION on N.  Without it the PER-ELEMENT Newton converges to a WRONG
      root (M~0.09 vs 0.30 on the sphere at H0=0.1): the centroid-monopole N is accurate in AVERAGE
      (far field) but POOR per-element (local field), and Newton consumes the local field per cell.
      Folding in N_eff = B^T (G + corr) B restores the local fields -> correct root (M=0.29901 vs
      analytic 0.29872).  [Picard's SCALAR field hid this -- it only ever used the AVERAGE demag.]
  (c) line-search DAMPING (Armijo) + a scalar-chi PICARD WARMSTART.  Damping gives global robustness;
      the warmstart lands inside the Newton basin at the stiff BH-knee.

RESULTS (sphere, chi0=1000, Msat=1, near-corrected D=0.328):
    H0     Newton M    iters    analytic(1/3)   rel
    0.01   0.03052      3        0.02991        2.0e-2   (rel = operator demag-D systematic, not Newton)
    0.10   0.30413    ~slow      0.29872        1.8e-2
    0.30   0.88669    ~slow      0.87834        9.5e-3   (BH knee, stiff)
    1.00   0.99840     11        0.99850        1.1e-4   <- per-element Picard/Hantila FAILED here
    5.00   0.99978      4        0.99979        4.2e-6   <- deep saturation, Newton trivial
The ~2% at LOW drive is the operator's demag-factor accuracy (monopole+near-corr D=0.328 vs 1/3),
SHARED with Picard -- NOT a Newton failure (finer/analytic Gram -> closer; the #3 / Wilton path).

## Honest limitation: the BH knee is stiff (slow, not wrong)
At the knee (H0~0.1-0.3, chi_sec/chi_diff ~ 8x) Newton converges to the CORRECT answer but slowly
(line-search takes tiny steps; per-element field noise x high chi = genuine stiffness).  The scalar
Picard is well-conditioned there, so the practical recipe is: scalar Picard handles the knee, Newton
is the tool for the SATURATION regime where per-element accuracy matters and Picard fails.  The golden
test therefore locks the saturation WIN (fast + accurate at H0>=1) and the near-correction ingredient.

## CROSS-CHECK vs Radia MMM/MSC -- DONE + golden-locked
The HDiv-VIM Newton is cross-validated against a COMPLETELY INDEPENDENT codebase: Radia's trusted
MMM/MSC C++ tetrahedral solver (rad.Solve + rad.MatSatIsoTab) with the SAME saturating BH curve,
chi0/Msat, sphere, and applied field B0 = mu0 H0.  At deep saturation (chi0=1000, Msat=1e6 A/m):
    H0(A/m)   HDiv-VIM M    Radia M     analytic   HDiv rel   Radia rel
    5.0e5     993620        994080      994105     4.9e-4     2.5e-5
    1.0e6     998384        998472      998503     1.2e-4     3.1e-5
    3.0e6     999612        999621      999625     1.3e-5     4.3e-6
Two independent solvers agree with each other AND the analytic to <0.05%.  (At saturation M->Msat
regardless of the demag factor, so the operator's ~2%-at-low-drive demag-D systematic is irrelevant
-- this isolates the nonlinear SATURATION behaviour, exactly what Newton had to get right.)
Golden test: tests/feec/test_hdiv_vim_newton_vs_radia.py.

## NON-UNIFORM body (cube) cross-check -- characterized open boundary (honest, NOT a clean win)
Cube in a uniform field, HDiv-VIM Newton vs Radia MMM/MSC (same chi0=1000, Msat=1e6, same BH curve):
    H0(A/m)   HDiv Mavg   Radia Mavg   Radia Mstd   rel(H-R)
    2.0e5     642207      734557       162824       1.3e-1     <- moderate drive: 13% gap
    5.0e5     935369      985670         9573       5.1e-2
    1.0e6     986986      996728         1539       9.8e-3     <- saturation: agree to ~1%
At SATURATION the cube agrees with Radia to ~1% (M->Msat, demag-D-independent -- consistent with the
sphere result).  At MODERATE drive there is a ~13% gap, and Radia shows the cube M is GENUINELY very
non-uniform there (std = 162824 = 22% of mean: corners saturate, centre does not).
CAUSE (diagnosed, NOT convergence / NOT knee):
  - converged: maxit 120 vs 400 give the SAME M_avg (642207) -- the Newton is at its fixed point;
  - near-corr-insensitive: cube monopole+nearcorr D = 0.3125 for near_factor 2/3/4 (analytic cube
    D_uniform = 1/3 = 0.3333); the near correction helps local fields but does NOT lift the cube's
    uniform-mode D;
  - the HDiv-VIM cube is DEMAG-LIMITED at M_avg ~ H0/D_monopole = 2e5/0.3125 = 640000 (= the 642207
    observed), while Radia's non-uniform MMM/MSC solution redistributes M to a LOWER effective demag
    (2e5/734557 = 0.272) -> higher M_avg.
=> the gap is the centroid-MONOPOLE operator's PER-ELEMENT field accuracy on a SHARP-CORNERED body at
   moderate drive (where M is highly non-uniform).  The same operator-accuracy limit as the linear
   demag, exposed by non-smooth geometry; the FIX is higher-order Gram accuracy (analytic Wilton
   face-integral) -- NOT a solver change.

## WILTON analytic SURFACE Gram -- DONE (production #1, 2026-06-07, golden-locked)
build_demag(..., wilton_surface=True) replaces the surface-surface (boundary-triangle) OFF-diagonal
Gram block AND its self with the exact Wilton/Graglia analytic potential of a uniformly-charged flat
triangle (tet.tri_potential + tet.wilton_surface_block; Dunavant-5 outer x analytic inner).  Because a
UNIFORM M has div M = 0 (zero volume charge), the demag factor is PURELY surface charge, so the Wilton
surface Gram makes it EXACT:
    body     monopole+nearsub   Wilton-surface   analytic
    sphere   0.3135             0.3329 (0.12%)   0.3333
    cube     0.3108             0.3334 (0.02%)   0.3333    <- the cube the near-correction could NOT fix
The Wilton SELF is shape-exact (any triangle); the old tri_self_energy assumes equilateral (fixed
C_TRI), so on real meshes the Wilton self is the better diagonal -- DO NOT overwrite it with
tri_self_energy (that regressed the demag to 0.3256).  Golden: tests/feec/test_hdiv_vim_wilton_gram.py
(tri_potential vs numerical + cube/sphere demag <0.5%).  Ref: Wilton IEEE TAP 32(3):276 (1984);
Graglia IEEE TAP 41(10):1448 (1993).

For the NONLINEAR solve, wilton_surface=True is combined with a VOLUME-only near-correction
(build_near_correction(skip_surface_surface=True)): the surface block is exact (Wilton) and the
volume-involving (cell-cell, cell-face) near pairs keep the sub-point correction the per-element Newton
needs (skipping it -> wrong root, M=0.83 vs 0.99 at saturation -- the same wrong-root failure as the
sphere without near-corr).  On the nonlinear cube this IMPROVES the moderate-drive gap (13% -> 8.7%)
and saturation (0.98% -> 0.75%), but does NOT fully close the moderate-drive gap: the residual is the
VOLUME Gram (cell-cell, cell-face still monopole+sub-point).  HONEST: the Wilton SURFACE Gram is a
clean win for the demag factor (surface-charge-dominated, exact); the full sharp-body NON-UNIFORM
accuracy additionally needs the analytic VOLUME (tet) Gram -- the next refinement.

tri_potential is VECTORIZED over observation points (r as (M,3) -> (M,)), and wilton_surface_block
loops O(n_bf) (one batched potential eval per SOURCE triangle) not O(n_bf^2) -> the wilton_gram golden
test went 73s -> 2s.

## NON-SPHERE smooth body (production #3) -- ellipsoid, DONE (golden-locked)
A prolate spheroid has an ANALYTIC demag factor != 1/3, so it checks the Wilton Gram on a non-1/3
smooth body.  2:1 prolate (long axis z), maxh=0.6:
    monopole N_z = 0.1693 (2.5%)   Wilton N_z = 0.1742 (0.3%)   analytic N_z = 0.1736
-> the Wilton surface integral is correct for a general curved-surface triangulation, not just the
isotropic 1/3.  Golden: tests/feec/test_hdiv_vim_ellipsoid.py.  (Remaining #3 options: a real
rad.MatSatIsoTab steel TABLE -- needs table-based M(H)/M'(H) in the Newton -- and a C-yoke.)

## SCALABLE nonlinear (production #2) -- DONE (2026-06-07, golden-locked)
solve_nonlinear_newton_scalable (examples/vim/hdiv_demag_tet_nonlinear.py) replaces the dense
O(N^3)/O(N^2) demag of the dense Newton with the C++ HACApK charge-Gram H-matrix
(radia._radia_pybind._ChargeGramHMatrix, O(N log N) apply) + a sparse near-correction, and solves each
Newton step ITERATIVELY (GMRES, M_mass-preconditioned) -- NO dense factorization anywhere:
    N v = B^T ( H.matvec(B v) + corr (B v) )        # textbook H-matrix far + sparse near split
    J v = M_mass v + T M_mass^{-1} N v               # matrix-free Jacobian (M_mass factored once)
It reproduces the dense solve_nonlinear_newton (same monopole+near-corr operator) to ~MACHINE PRECISION
at saturation (sphere, rel D-vs-S: 1.1e-12 @ H0=5e5, 5.2e-14 @ 1e6, 0 @ 5e6; ~1.8e-3 at moderate drive
= the ACA/GMRES tolerance) and matches the analytic to <0.06% at saturation.  Golden:
tests/feec/test_hdiv_vim_newton_scalable.py.  The heavy cost (the Gram apply) is the C++ H-matrix
(O(N log N), established for the linear case in production #1/#2); the Newton outer loop is Python (an
O(N)-per-iteration overhead dominated by the C++ H-matvec).

## VOLUME (tet) Gram -- DONE (#B, closes the nonlinear sharp-body residual, 2026-06-07, golden-locked)
The surface Wilton Gram makes the demag FACTOR exact (uniform M = pure surface charge), but the
NONLINEAR sharp body (cube/C-yoke) left a residual from the VOLUME-involving (cell-cell, cell-face)
Gram blocks (div M != 0 -> volume charge).  The analytic tet Newtonian potential closes it WITHOUT new
singular quadrature, by reusing the Wilton triangle potential via the divergence theorem:

    phi_tet(V, P) = INT_tet 1/|P-r'| dV' = (1/2) sum_{4 faces} d_face * tri_potential(face, P)
                    (nabla'^2 R = 2/R ; d_face = (r'-P).n_hat const on a flat face)

build_demag(analytic_gram=True) (and solve_nonlinear_newton(analytic_gram=True)) builds the FULL
analytic charge Gram (tet.analytic_charge_gram: cell sources via phi_tet, face sources via
tri_potential, outer tet/Dunavant quadrature) -- every block exact, no monopole/near-correction.
RESULT: cube nonlinear vs Radia 13% -> 6.2% (H0=2e5), 5.1% -> 1.5% (H0=5e5); linear demag stays exact
(sphere 0.3329, cube 0.3335).  Substantially closes the gap (the residual ~6% at the stiff moderate
KNEE is the coarse mesh + outer-quadrature order + knee stiffness, not the Gram).  phi_tet verified vs
fine volume quadrature; golden tests/feec/test_hdiv_vim_volume_gram.py.

## HIGH-ORDER VIM -- STRUCTURE CONFIRMED at order 1; accuracy win is the continuation (2026-06-07)
The FEEC structure is ORDER-AGNOSTIC: HDiv(order>=1) assembles, the charge map B = (rho=-div M tested
on L2(p), sigma=M.n tested on SurfaceL2(p)) builds, and loops = ker(B) are FIELD-NULL BY CONSTRUCTION
at EVERY order (B.loop=0 => N.loop=B^T G(B loop)=0, any G) -- no element engineering, automatic.
DOF/charge/loop growth (tet 2x2x2): order 0 -> 120/96/25, order 1 -> 360/336/169, order 2 ->
1008/768/529.  CONFIRMED end-to-end: the order=1 demag operator BUILDS (charge field of each HDiv basis
-- -div in cells, gfu.n at boundary points -- sub-point Coulomb Gram) and its loops are field-null to
2.1e-15 (sphere coarse: ndof 288->864 order 0->1, loops 55->391).  So the high-order CAPABILITY is real
(multipole-moment MMM MSC is lowest-order only).

TWO things remain for the p-convergence SPEED WIN:
  (1) the accurate higher-order CHARGE GRAM: at order>=1 rho/sigma are POLYNOMIALS per cell/face, so
      G[i][j] = INT INT charge_i charge_j/|x-x'| needs singular quadrature for POLYNOMIAL densities --
      the higher-order analog of the Wilton/Graglia surface Gram (tri_potential / phi_tet built for the
      CONSTANT case).  The crude sub-point Gram limits accuracy (caps the structural probe).
  (2) a NON-uniform benchmark: the demag FACTOR (uniform M) is order-INSENSITIVE -- M.n is constant per
      face regardless of order, so order 0 and 1 give the SAME demag (0.3044/0.3044 on the probe).  The
      high-order benefit shows in a FIELD / non-uniform-M error, not the uniform-M demag factor.
Both are the major continuation (the p-convergence accuracy-per-DOF beats lowest-order multipole-moment MMM MSC).

## CURVED MESHES (mesh.Curve) -- a big accuracy-per-DOF win HDiv has and flat multipole-moment MMM cannot (2026-06-07)
NGSolve HDiv supports mesh.Curve(p) (isoparametric/curved elements via the Piola transform).  Multipole-moment MMM
CANNOT: Radia's ObjHexahedron/ObjTetrahedron are FLAT-faced.  CONFIRMED dramatic: a COARSE sphere
(maxh=0.8, ndof=258) flat has surface area 11.85 (-5.7%) / volume 3.76 (-10.3%) vs the true 4pi/4pi/3;
mesh.Curve(3) gives area 12.5683 (-0.015%) / volume 4.1897 (+0.02%) at the SAME ndof.  So the ~6-10%
faceting error VANISHES to ~0.02% at no DOF cost -- exact geometry.  For curved bodies (magnets/coils/
poles = most of them) this is a large accuracy-per-DOF advantage, ORTHOGONAL to the polynomial
high-order win, and unavailable to flat multipole-moment MMM MSC (which needs many more elements to resolve a curve).
FOUNDATION confirmed (HDiv + mesh.Curve assembles; B and M_mass are NGSolve curved-aware).

CURVED SAMPLING BUILT + the win MEASURED vs TRUTH (2026-06-08, examples/vim/hdiv_demag_curved.py,
golden tests/feec/test_hdiv_vim_curved.py):
  - Reusable primitive _trafo_sample(mesh, i_bnd, xi, eta, center) -> curved position / surface |J| /
    OUTWARD normal-z, via mesh.GetTrafo (same pattern as bem/sibc_hacapk.py::_trafo_eval; works on a
    linear mesh = flat AND a Curve(p) mesh = curved -- the SAME code, only mesh.Curve toggled, so a
    flat-vs-curved A/B at FIXED ndof isolates the geometry).  (Beware: it is mip.JACOBI / mip.MEASURE,
    NOT .jacobian/.weight -- the GetTrafo API gotcha that cost two probe iterations.)
  - ELEMENTARY DISCRIMINATOR (vs ANALYTIC truth): the EXTERNAL field of a uniform-M sphere = the EXACT
    point dipole phi=(1/4pi)V cos(theta)/r^2.  A surface-charge integral at an EXTERNAL point -> NO
    singular quadrature -> the ONLY error is the geometry.  MEASURED at (0,0,2) [exact 1/12]: h=0.8 FLAT
    -10.00%, Curve(3) -0.26% at the SAME ndof (~38x); h=0.5 FLAT -8.28%, curved -0.25%.  The flat error
    IS the volume error (V -10.25% flat -> +0.02% curved) inherited by m = M V.  THIS is the curved win.
  - CAVEAT (self-corrected, do NOT repeat the mistake): with the CRUDE sub-point Gram the demag FACTOR
    does NOT cleanly discriminate -- but the reason is the crude Gram's ~2% quadrature BIAS masking the
    ~0.25% geometry signal, NOT "near-isotropy" (an earlier write-up wrongly said near-isotropy).  With
    the PROPER Gram the demag factor DOES discriminate + p-CONVERGES (next bullet).  So for the crude
    elementary method use the external field; the demag factor needs the proper single-layer Gram.
  - PROPER curved + HIGH-ORDER demag Gram = SOLVED by ngsolve.bem (2026-06-08, the architectural unlock,
    examples/vim/hdiv_demag_bem_singlelayer.py + golden test_hdiv_vim_bem_demag.py): the surface
    demag Gram (uniform M -> pure surface charge sigma=M.n) IS the Laplace SINGLE-LAYER V; NGSolve
    6.2.2604 ngsolve.bem ships it high-order + CURVED-aware + FMM-accelerated.  `V =
    SingleLayerPotentialOperator(SurfaceL2(mesh,order=p), intorder=...)`; demag D_z = <sigma,V sigma>/V_vol
    with sigma=GridFunction.Set(specialcf.normal(3)[2]).  Kernel is 1/(4pi r) (no extra factor -- gives
    1/3).  MEASURED sphere demag [exact 1/3], h=0.6 coarse: FLAT +0.247% at order 0/1/2 IDENTICAL
    (order-insensitive: sigma=M.n is constant per FLAT face) -> faceting floor, only mesh-refinement helps;
    CURVE(3) order 0 -1.89% (piecewise-const sigma under-resolves the n_z that VARIES on a curved face) ->
    order 1 -0.06% -> order 2 **-0.0002% EXACT**; mesh+intorder converged.  NON-ISOTROPIC shape check (2:1 prolate spheroid,
    analytic N_z=0.17356 != 1/3): curved+order2 0.17356 (-0.001% EXACT) vs flat 0.17415 (+0.34%, order-
    insensitive) -> the single-layer gets the anisotropic SHAPE right, not just isotropy.  FULL demag
    TENSOR (curved+o2, axis=0/2 -> sigma=n_x/n_z): prolate(2:1) N_par 0.17356 + N_perp 0.41322; oblate(1:2)
    N_par 0.52720 + N_perp 0.23640 -- all EXACT vs analytic (Osborn 1945), and the formula-independent SUM
    RULE N_x+N_y+N_z=1 holds to ~1e-6.  GENERAL TRIAXIAL ELLIPSOID (a!=b!=c, the canonical fully-
    anisotropic benchmark): (a,b,c)=(1,1.5,2) -> three DISTINCT factors N_a=0.48373/N_b=0.30501/
    N_c=0.21127 all EXACT vs the analytic Osborn integral (numerically integrated, no elliptic-integral
    formula), sum=1.000000 -> the single-layer handles FULL anisotropy, not just spheroid symmetry
    (golden test_hdiv_vim_bem_demag.py::{test_prolate_nonisotropic_shape_exact,
    test_spheroid_full_tensor_and_sum_rule,test_triaxial_ellipsoid_full_anisotropy}).  => curved+high-order converges
    the demag to exactness AT COARSE MESH, fixed small ndof = the accuracy-per-DOF win over flat
    lowest-order multipole-moment MMM MSC, ON THE DEMAG FACTOR (corrects the crude-method "doesn't discriminate").  Reuses
    NGSolve, NO hand-rolled singular quadrature -- supersedes the Wilton/phi_tet SURFACE block for the
    curved/high-order/scalable path.
  - STILL TO BUILD: (1) the full operator N=B^T V B + a self-consistent linear solve on curved/high-order
    geometry (demag-factor proof done; operator+solve is next); (2) the VOLUME charge (div M != 0,
    nonlinear) on curved/high-order -- ngsolve.bem is boundary-only so this still needs the Newtonian
    volume potential phi_tet (built) on curved cells; (3) C++ maturity (single-layer surface + phi_tet
    volume Gram + Newton loop behind a Radia API) = the productionization lift.

## CURVED x NONLINEAR -- the honest magnitude (verify-first, 2026-06-08; golden test_hdiv_vim_curved_nonlinear.py)
The curved win on the nonlinear MAGNETIZATION is MODEST (~0.3%), NOT a dramatic differentiator -- an
important honest finding that redirects where to expect the curved payoff.
- WHY: a spheroid in a uniform field magnetises UNIFORMLY for ANY M-H law (linear or nonlinear), so the
  nonlinear M is the scalar fixed point M = M(H_ext - N_par M) with N_par the demag factor.  The demag
  factor is a VOLUME-NORMALISED RATIO, so it CANCELS the ~10% volume faceting error -> flat nonlinear M
  is only ~0.3% off.  curved+o2 N_par is exact -> M matches the analytic-demag fixed point <0.05%.
- It does NOT grow with severe curvature: thin oblate disks facet BETTER on the demag (flat M err 1:2
  0.20% -> 1:5 0.024% -> 1:10 0.004%; their surface is dominated by near-flat caps).  Probed prolate
  2:1/5:1 + oblate 1:2/1:5/1:10, chi0=5000/Msat=1.6e6, H_ext 1e3-1e5.
- WHERE the big curved win (~10%) IS: the FIELD output (the volume error is NOT cancelled there --
  external field flat -10% vs curved -0.26%, hdiv_demag_curved.py); the nonlinearity merely SCALES that
  field by M (no amplification).
- VALIDATION CONSTRAINT (important): RADIA IS NOT A CLEAN REFERENCE FOR CURVED geometry -- Radia's
  ObjHexahedron/ObjTetrahedron are FLAT, so Radia facets the body; comparing curved-HDiv-VIM vs flat-
  Radia confounds the geometry win with any bug.  So curved nonlinear is validatable ONLY vs ANALYTIC =
  spheroids (uniform M).  The genuinely NON-UNIFORM (volume-charge, div M != 0) curved nonlinear case
  has NO clean reference -- a mesh-convergence / self-consistency study at best, NOT a verified win.
- CONCLUSION: curved is a real flat-multipole-moment MMM-impossible CAPABILITY, but its accuracy benefit on the
  MAGNETIZATION (the primary solve output) is modest; the benefit is large on FIELD outputs.  Do NOT
  claim a dramatic curved-nonlinear magnetization win.
- (A) THE FIELD WIN -- BUILT + validated (examples/vim/hdiv_curved_nonlinear_field.py, golden
  test_hdiv_vim_curved_nonlinear_field.py): the external H field of a nonlinear soft-iron sphere
  (M_s=29982 A/m from the fixed point at H_ext=1e4) reconstructed from the curved-aware surface charge
  H(r)=(1/4pi) INT sigma (r-r')/|r-r'|^3, sigma=M.n, vs the ANALYTIC dipole (m=M V).  At 5 external
  points: FLAT ~+8.8% at EVERY point (the dipole moment inherits the ~9% volume faceting error) vs
  Curve(3) <0.4% -- a ~23x field win at the SAME ndof, vs analytic truth.  The nonlinearity SETS the
  field magnitude (physical M) but does NOT amplify the ~9% geometry error (scales it).  So the
  engineering deliverable (stray field around a nonlinear part) is ~9% wrong with flat multipole-moment MMM MSC and
  exact with curved -- the genuine curved x nonlinear payoff is HERE, on the field.

## (B) HEAD-TO-HEAD vs the SHIPPED Radia solver -- accuracy-per-resolution (2026-06-08)
The FIRST quantitative head-to-head of HDiv-VIM (curved) vs the PRODUCTION code (examples/vim/
compare_curved_vs_radia_field.py, golden test_curved_vs_radia_field.py).  Shipped Radia (rad.Fld on a
flat-tet uniform-M sphere built via netgen_mesh_to_radia) vs HDiv-VIM curved surface-charge field, both
vs the ANALYTIC dipole.  Max external-field relative error:
    h=0.6:  Radia FLAT 114 tets  8.92%  |  HDiv CURVED 120 tris  0.386%
    h=0.4:  Radia FLAT 260 tets  5.87%  |  HDiv CURVED 192 tris  0.229%
    h=0.3:  Radia FLAT 477 tets  3.48%  |  HDiv CURVED 318 tris  0.147%
    h=0.2:  Radia FLAT 2042 tets 1.71%  |  HDiv CURVED 658 tris  0.072%
=> HDiv CURVED at the COARSEST mesh (0.386%) BEATS shipped-Radia-FLAT at the FINEST (1.71%, 2042 tets);
~10-30x accuracy-per-resolution at every h.  HONEST SCOPE: ACCURACY-PER-DOF (geometry-driven, fair),
NOT wall-clock -- the HDiv-VIM here is a Python dense surface-charge prototype, not time-optimized; a
fair SPEED comparison needs the C++ productionization (NOT done).  Radia's ObjTet are FLAT = the
accessible stand-in for the also-flat multipole-moment MMM MSC; reference = the analytic dipole (Radia cannot referee
curved -- it facets).  This is the quantitative basis for the curved accuracy-per-DOF advantage over the
flat production solver; the remaining lift to a TOTAL win (incl. speed) is the C++ productionization.

## SYMMETRY MODELS (1/2, 1/4, 1/8) -- the loop handling is AUTOMATIC (2026-06-08)
The painful part of symmetry models in multipole-moment MMM MSC is the LOOP handling: the loop-star basis must be
built with a cohomology-aware installCycle on the CUT domain (the symmetry planes introduce new cycles).
In HDiv-VIM the loops are simply ker(B) (B = charge map), so on ANY cut/reduced mesh they are field-null
BY CONSTRUCTION (N = B^T G B => N.loop = 0 for loop in ker B) -- NO hand-crafted loop-star basis, NO
cohomology bookkeeping.  VERIFIED (golden test_hdiv_vim_symmetry_loops.py): a sphere as full / 1/2 / 1/4
/ 1/8 -> the loop count adapts AUTOMATICALLY (58 / 54 / 18 / 6) and the loop field-null residual
||N.loop||/||N|| ~ 4e-16 (machine zero) on EVERY reduced mesh.  => the "loop-jokyo (loop removal) is
metsuky-doi / mendokusai" problem is ELIMINATED -- ker(B) handles the cut topology for free.
SCOPE (honest 2-part split): (1) the LOOP machinery on cut meshes = automatic (verified here).  (2) the
demag VALUE of a symmetry MODEL = the IMAGE method -- NOW BUILT + verified (examples/vim/
hdiv_demag_symmetry_image.py, golden test_hdiv_vim_symmetry_image.py).  Only the REAL surface (spherical
cap) carries sigma = M.n = n_z; the flat cut faces are symmetry planes (no real charge).  Reflecting the
cap charge over the reduction planes -- sign = (-1)^(#z-reflections), since sigma = n_z flips under a
z-mirror (the IMA sign rule: field-PARALLEL mirror x=0/y=0 keeps sign, field-PERPENDICULAR z=0 flips it)
-- reconstructs the full sphere's sigma = cos(theta).  RESULT (M=z_hat, vs the full-sphere demag from the
same crude Gram): 1/2 +0.08%, 1/4 +0.11%, 1/8 -0.32% -- i.e. the reduced models reproduce the FULL demag
from ~1/2, 1/4, 1/8 the surface DOF (108 / 52 / 20 cap-tris vs the full 192).  So 1/4 and 1/8 models are
SUPPORTED: loops automatic (ker B) + demag via the image method, no hand-crafted loop-star.  (Production
note: this uses the elementary sub-point Gram; the ngsolve.bem single-layer with image kernels is the
high-order/curved/scalable production version of the same image method.)

## NON-UNIFORM NONLINEAR needs analytic_gram; C-YOKE VERIFIED vs Radia (2026-06-08, the 1/8-gate audit)
The gate before symmetry models: confirm nonlinear + C-yoke + distorted-mesh are solid.  Outcome:
- GRAM REQUIREMENT (a trap, now fail-loud): a NON-uniform-M nonlinear solve (div M != 0: cube, C-yoke,
  any non-ellipsoid) REQUIRES analytic_gram=True (the full volume Gram via phi_tet).  The surface-only
  wilton_surface Gram leaves the volume (cell) blocks crude -> wrong per-element fields -> Newton does
  NOT converge (stalls at maxit, M_avg drifts).  Only UNIFORM-M nonlinear (sphere/spheroid) converges
  with wilton_surface.  solve_nonlinear_newton now RAISES on non-convergence (require_convergence=True
  default) with a message pointing at analytic_gram -- NO silent wrong result (this trap previously made
  the C-yoke look "mesh-dependent / non-converging" when it was just the wrong Gram).
- C-YOKE VERIFIED (examples/vim/hdiv_cyoke_nonlinear.py, golden test_hdiv_vim_cyoke_nonlinear.py):
  non-convex reentrant-corner C-yoke, nonlinear (chi0=1000/Msat=1e6/H0=2e5), analytic_gram -> converges
  in 6 Newton iters, volume-avg M_z MESH-STABLE (572062/576970/580981 over maxh 0.020/0.016/0.013) and
  matches shipped Radia MMM (same flat mesh / M-H law / applied field) to <1% at EVERY mesh
  (-0.25%/+0.71%/-0.37%).  Cube likewise -0.08%/-0.15% (H0=2e5/5e5).  C-yoke is FLAT so Radia is a valid
  cross-reference.  CLOSES the C-yoke accuracy gate.
- ITERATION COUNT also dropped: analytic_gram converges in 5-6 iters where wilton_surface stalled at
  maxit (120-400) for non-uniform M -- the right Gram is both correct AND fast.
- DISTORTED-MESH: mu_r-independence (bounded iters vs mu_r 10->1e4) golden-locked
  (test_hdiv_vim_solve.py::test_minres_iters_bounded_vs_mu_r_distorted); the nonlinear ACCURACY on the
  non-convex C-yoke is now also verified (above).  => all three 1/8-gate items pass.

## REFERENCE HONESTY for the accuracy numbers (verify-first, 2026-06-07)
What the quoted accuracies are measured against, precisely:
- SPHERE / ELLIPSOID demag + nonlinear: vs ANALYTIC truth (D=1/3 or the prolate N_z; the scalar fixed
  point M = M(H0 - D M)).  These are REAL verified errors (sphere nonlinear <0.2%, demag <0.15%).
- CUBE / C-YOKE NONLINEAR: vs RADIA MMM/MSC on the SAME flat mesh (no analytic truth; Radia is the
  mature TRUSTED solver, valid here because BOTH are flat).  **CORRECTED 2026-06-08 (mesh-convergence
  study DONE)**: with analytic_gram + VOLUME-AVG M_z the agreement is <1% at every mesh -- cube -0.08%
  (H0=2e5) / -0.15% (H0=5e5) at maxh 0.28; C-yoke -0.25%/+0.71%/-0.37% over maxh 0.020/0.016/0.013, both
  converging in 5-6 Newton iters and mesh-stable.  The earlier "13% / 6.2% / ~4%" were STALE -- either
  the wrong (surface-only wilton_surface) Gram that does NOT converge for non-uniform M, or a stricter
  per-element metric (M std), NOT the volume-avg with the volume Gram.  So the honest current statement
  is "volume-avg M_z agrees with shipped Radia to <1% (cross-method, both flat)".  (A stricter
  per-element/pointwise field match would be larger and is a separate, harder claim -- not made.)

## HDiv on PYRAMIDS -- a genuine mathematical difficulty (verify-first, 2026-06-07)
NGSolve 6.2.2604 supports pyramids in HCurl + H1 (confirmed: ndof 8/20/57 and 5/5/14 across orders 0-2)
but NOT HDiv ("HDivHighOrderFESpace: Pyramid elements not implemented yet!").  This is NOT an oversight:
a pyramid is a DEGENERATE element (image of a collapsed hex, top quad -> apex), so the Piola map
(normal-flux-preserving, for H(div)) has a RATIONAL/singular Jacobian at the apex, and the
H(div)-conforming (Raviart-Thomas) pyramid space CANNOT be polynomial -- it needs RATIONAL basis
functions (Nigam-Phillips 2012; Bergot-Cohen-Durufle 2010).  H1/HCurl pyramids need rational functions
too but were implemented first; HDiv is the HARDEST of the three because the normal-flux continuity must
match adjacent tet(RT)/hex(RT) traces AND keep the de Rham sequence exact.  => HDiv-VIM pyramids are an
NGSolve-version-upgrade away. Radia now has `ObjPyramid` on the multipole-moment MMM MSC path, so pyramid meshes
should use multipole-moment MMM today; this remains an NGSolve HDiv-space limitation, not a VIM design flaw.

## NEXT (open): higher-order charge Gram (high-order build); Wilton/phi_tet-in-C++
- HIGH-ORDER charge Gram (the speed-win build above): polynomial-density singular quadrature.
- analytic VOLUME Gram: DONE (#B above) -- analytic phi_tet, cube residual 13%->6.2%; the residual ~6%
  at the moderate knee is coarse mesh + outer-quadrature order, refine if a tighter cube is needed.
- the C++ Gram H-matrix is MONOPOLE far field, so the scalable path matches the dense MONOPOLE+near-corr
  Newton, NOT the dense Wilton path; putting the Wilton surface integral into the C++ near-field
  correction gives scalable + Wilton-accurate together.
- a full C++ Newton OUTER loop (the tangent + line search in C++ too) removes the Python-loop overhead
  (minor: each iteration is dominated by the C++ H-matvec); the tangent solve is the iterative GMRES
  on the analytic H-matvec (no direct factorization).

## Why the operator is reusable (the structural win)
The demag operator N = B^T G B is GEOMETRY-ONLY (constant, mu_r-independent) -- it is assembled ONCE.
The ONLY nonlinearity is in the material term (1/chi) M_mass of A = (1/chi) M_mass - N, where chi
becomes chi(H) (or M(H) from a BH curve), evaluated PER ELEMENT PER ITERATION.  So the nonlinear
solve is an OUTER material iteration around the SAME (reusable) linear HDiv solve:

  repeat:
    1. given M (current), compute H_demag = N-related field per element (or H = H_ext + demag)
    2. update the per-element constitutive: chi_e = chi(H_e) from the BH curve (Picard), or the
       Jacobian dM/dH (Newton)
    3. re-solve the (updated-diagonal) symmetric system for M
  until ||Delta M|| small

## Solvers evaluated (verify-first outcomes)
- **Newton-Raphson (CHOSEN, SOLVED)**: per-element CONSISTENT TENSOR tangent dM/dH, line-searched +
  Picard-warmstarted.  Robust + fast at saturation; the winner.  See the HEADLINE / ingredients above.
- **Scalar Picard**: cheapest, well-conditioned, but assumes uniform M (uses the AVERAGE demag only).
  Valid + useful at MODERATE drive / the BH knee; the committed scalar foundation
  (solve_nonlinear, test_hdiv_vim_tet_nonlinear.py).  Hid the per-element near-field issue.
- **Per-element Picard**: diverges at saturation (transient |M_e|>Msat -> chi_e->0 -> ill-conditioned).
- **Hantila polarization** (constant (M_mass + alpha N) factored once -- structurally attractive,
  would need a factored solve): DIVERGED in the RT0 face metric (alpha=chi0/2 -> NaN; under-relaxed -> ~1e88).
  Two causes: chi_min->0 at saturation drives the contraction ->1, and the naive RT0 polarization
  projection r = Set((chi_e-alpha) H) amplifies.  Could be revisited with a consistent RT0 projection,
  but Newton already solves the problem -- not pursued.

## What exists to build on
- The per-element field/charge machinery (B, N, M_mass) is all in place + golden-tested (linear).
- DiagSystem(inv_chi) already produces the Jacobi diagonal for a given chi -- the per-element-chi
  update hook.  ApplySystem(x, inv_chi, y) is the per-iteration operator apply.
- Radia already has BH-curve materials (rad.MatSatIsoTab) + the Hantila solver (MMM 3DOF only today;
  the HDiv RT0 face DOF is the new target) + Play/Energy hysteresis -- reuse the constitutive layer.
- Validate against: a known nonlinear demag (e.g. a saturating sphere/cube in a strong applied field),
  and/or the existing MMM/MSC nonlinear solve on the same geometry.

## Next increment (the Newton path is done; what is left)
The Python operator Newton is solved + golden-locked.  Remaining: (1) a NON-UNIFORM saturating body
(cube / C-yoke) cross-checked against rad.Solve + rad.MatSatIsoTab; (2) the scalable C++ path -- the
per-element tensor tangent + the iterative GMRES tangent solve; (3) optionally a consistent-RT0
Hantila for the factor-once speed-up at the knee (Newton already gives correctness).
"""

_STATUS = r"""
# Status summary

[UPDATE 2026-07-03 -- hex RT1 + 2D planar SHIPPED; HDiv-VIM is the PRIMARY accurate route]
  * ROLE (Sugahara 2026-06-30): HDiv-VIM = the primary ACCURATE soft-iron demag route (loop-free by
    construction, ~6 mesh/mu_r-independent Newton iters); collocation MMMM = the COARSE/fast tier.
  * PURE-HEX RT1 charge Gram SHIPPED: Q1 volume charges (div(HDiv order-1) on a hex is Q1, NOT P0)
    + quad-face surface charges, extracted PIOLA-EXACTLY on the 27-node Q2 geometry; flat + curved.
    H-matrix build got symmetric leaf fill + far-inner Keast-15 + STATIC-SITE radial inner quadrature
    (~10x day total: cylinder benches 166/164 s -> 18.9/16.3 s flat/curved; goldens 223 s -> 10.5 s).
    The ~20k-charge use-after-free crash (thread_local cache capacity-clear while references were
    held) is FIXED (commit 20e6e9e2).  Curved-hex open item: cylinder max eig 1.0078 (> the [0,1]
    bound; halved from 1.0166) -- self/touching curved quadrature refinement pending.  `rad.Solve`
    hex dispatch STILL routes to collocation MMMM (dispatch-flip = open policy decision); reach the
    hex Gram via `radia.vim._vim.build_charge_gram(HDiv(hexmesh, order=1))`.
  * 2D PLANAR tri/quad Gram SHIPPED (commit a9999dd7, motor cross-sections): log kernel -ln(r)/2pi,
    charges = -div M on cells (P0 tri / Q1 quad -- the 2D twin of the hex gotcha) + M.n on boundary
    edges.  Closed-form gated: disk demag 1/2 EXACT (0.50000), ellipse 2:1 -> 0.33438/0.66562,
    2D Clausius-Mossotti M/H0 = chi/(1+chi/2) to 2-3e-4.  Same build_charge_gram auto-routing.
    Quadrature lessons locked in the golden: the outer MUST be product-Gauss (Dunavant-7 leaked the
    quad spectrum to 1.072 while entries agreed to 3e-5); the edge inner split-grades at the
    kernel-peak PARAMETER.  Nonlinear is tet-only for now (2D/hex nonlinear = open; wedge = open).
  * EXECUTED SHOWCASE: docs/hdiv_vim/hex_rt1_and_2d_showcase.ipynb (+ _result.json sidecar) -- hex
    spectrum/cube-1/3 gates incl. a genuine-warp real-mesh hex, fresh H-matrix build timings
    (160 -> 5632 charges: 1.1 -> 57 s on LAB), all 2D closed-form gates, and the production
    hdiv_demag_solve one-call on a tet sphere.  Goldens: validation_test/feec/
    test_hdiv_vim_hex_rt1_wiring.py + test_hdiv_vim_2d_wiring.py.

[API UPDATE 2026-06-23] The dense Python Gram path -- the `analytic_gram=` / `wilton_surface=` /
`skip_dense_gram=` kwargs, `build_near_correction`, the dense `analytic_charge_gram` /
`wilton_surface_block` / `phi_tet` builders, and the monopole+near-correction split named in this section
-- was REMOVED.  The C++ `_ChargeGramHMatrix` (exact analytic near AND far, IMA via image_masks/image_signs)
is now the SOLE demag operator and `radia.vim.hdiv_demag_solve(mesh, mu_r=/bh_table=, H_ext=, image=)` the
production entry; read the removed names below as the research history behind the always-exact C++ Gram.

## Snapshot (2026-06-08, historical -- see the API UPDATE above for current call signatures)

DONE + golden-locked (feec 85/85):
  #1  scalable mu_r-independent HDiv-VIM demag solver on REAL tet meshes (Layer A/A.5 + tet ingest)
  #2  analytic Wilton/phi_tet charge Gram in the C++ scalable path (_ChargeGramHMatrix analytic mode, M2)
  #3  bug-fixed exact Gram via near-field correction -> demag -> analytic 1/3
  NONLINEAR  damped Newton-Raphson (consistent tensor tangent + near-corr + line search + Picard
             warmstart) -- robust + fast at deep saturation where per-element Picard/Hantila failed
             (examples/vim/hdiv_demag_tet_nonlinear.py::solve_nonlinear_newton,
             tests/feec/test_hdiv_vim_tet_newton.py, +the scalar-Picard moderate-drive foundation).
             CROSS-VALIDATED vs Radia MMM/MSC (rad.Solve+MatSatIsoTab) to <0.05% at saturation
             (tests/feec/test_hdiv_vim_newton_vs_radia.py).
  WILTON     analytic SURFACE Gram (build_demag(wilton_surface=True)): exact triangle-triangle Coulomb
  GRAM (#1   integral -> demag factor 1/3 to <0.15% on cube AND sphere (the cube the near-correction
  of new     could NOT fix); the new sequence's #1.  tet.tri_potential / tet.wilton_surface_block;
  sequence)  golden tests/feec/test_hdiv_vim_wilton_gram.py.  Combined w/ a VOLUME-only near-correction
             for the nonlinear path (cube moderate gap 13% -> 8.7%).
  SCALABLE   scalable nonlinear Newton via the C++ Gram H-matrix (O(N log N) apply) + GMRES, NO dense
  (#2 of     factorization; reproduces the dense Newton to ~machine precision at saturation
  new seq)   (examples/...::solve_nonlinear_newton_scalable, tests/feec/test_hdiv_vim_newton_scalable.py).
  #3 of new  ellipsoid (non-1/3 smooth demag): 2:1 prolate N_z Wilton 0.1742 (0.3%) vs analytic 0.1736
  sequence   (monopole 2.5%) -- tests/feec/test_hdiv_vim_ellipsoid.py.  + Wilton vectorized (golden 73s->2s).
  REAL BH    solve_nonlinear_newton(..., bh_table=(Harr,Barr)) drives the operator Newton from a REAL
  TABLE      [[H,B]] table (PCHIP M(H)=B/mu0-H, saturation-clamped beyond Hmax -- the physical mu_r->1
  (real      extension Radia's MatSatIsoTab also uses; raw PCHIP extrapolation blows up).  Matches the
  material)  analytic uniform-sphere fixed point (same table, Wilton D->1/3) to <0.2% across
             H0 10^4->3e6, saturating to the table Msat.  tests/feec/test_hdiv_vim_newton_table.py.
             (Also fixed _scalar_fixed_point's hi<1e6 bisection cap -> 1e12 so a real Msat~1.5e6 is
             reachable.)  NOTE: a Radia cross-check on a COARSE sphere disagreed at LOW field -- traced
             to Radia not converging at low drive on the coarse mesh (M/H0=11.8 is impossible for a
             sphere D=1/3); HDiv matched the analytic, so HDiv was correct.  Use the analytic fixed
             point (not a coarse-mesh Radia run) as the low-field reference.
  C-YOKE     real non-convex engineering geometry (OCC box-box-box C-shape) -- VERIFIED + golden-locked
  (real      (2026-06-08, the 1/8-gate audit): with analytic_gram (the REQUIRED volume Gram for
  geometry)  div M != 0) the nonlinear Newton converges in 6 iters, the volume-averaged Mz is mesh-stable
             (572062/576970/580981 over maxh 0.020/0.016/0.013) and matches shipped Radia MMM to <1% at
             EVERY mesh (-0.25%/+0.71%/-0.37%); cube likewise -0.08%/-0.15% (H0=2e5/5e5).  The OLD
             "~4%@0.02 / 13%@0.03 / not-golden-locked" was the WRONG (surface-only wilton_surface) Gram
             that does NOT converge for non-uniform M -- now FAIL-LOUD (solve_nonlinear_newton raises,
             pointing at analytic_gram).  GOTCHA kept: compare the volume-averaged Mz (not point-sample).
             examples/vim/hdiv_cyoke_nonlinear.py, golden test_hdiv_vim_cyoke_nonlinear.py.
  CURVED+HO  curved + high-order demag via the ngsolve.bem Laplace single-layer: sphere/spheroid/triaxial
             EXACT vs analytic; field accuracy-per-DOF ~10-30x vs the shipped flat Radia solver
             (compare_curved_vs_radia_field.py).  See topic "curved".
  SYMMETRY   1/2, 1/4, 1/8 models: loops AUTOMATIC (ker B, no loop-star) + image-method demag value
             (reproduces the full demag from ~1/N DOF).  See topic "symmetry".

CLOSED since the 2026-06-07 status (do NOT re-open):
  - analytic VOLUME (tet) Gram: DONE (analytic_gram / phi_tet) -- the cube/C-yoke non-uniform nonlinear
    gap closed (cube -0.08%, C-yoke <1% vs Radia, volume-avg); the old "~8.7%/13%" were the wrong Gram.
  - real BH table + C-yoke validation: DONE (see DONE list).

OPEN (honest boundaries / next increments) -- the lift to productionize HDiv-VIM alongside multipole-moment MMM:
  - C++ PRODUCTIONIZATION (the big one): the charge Gram (Wilton surface / phi_tet volume / ngsolve.bem
    single-layer) + the Newton loop in C++ behind a Radia API.  This also enables a fair WALL-CLOCK
    comparison -- all present wins are accuracy-per-DOF (geometry-driven); the Python prototype is not
    time-optimized.  Until this lands, HDiv-VIM is a validated method but not production-sealed.
  - CURVED nonlinear VOLUME charge: ngsolve.bem is boundary-only, so non-uniform nonlinear on CURVED
    cells still needs the Newtonian volume potential (phi_tet) on curved geometry.
  - symmetry-model demag at HIGH ORDER / CURVED via the ngsolve.bem single-layer + image kernels (the
    elementary sub-point image method is the prototype).
  - BH-knee stiffness; near-field/Wilton Gram in the C++ ChargeGram entry; proper distorted RT0 M_mass
    for exact distorted demag VALUES (non-negativity already holds).
DELETED (2026-06-08): the direct symmetric H-LDL^T factorization was REMOVED from the codebase
  -- HDiv-VIM is mu_r-INDEPENDENT by construction, so the ITERATIVE solve (MINRES linear / GMRES/Picard
  nonlinear) is cheap + well-conditioned; a direct factorization bought little.  The C++ symmetric
  factorization + its self-tests are gone; H-LU/H-ILU kept (MSC solver A_SS preconditioner) (see
  docs/hdiv_vim/PRODUCTIONIZATION.md M3).  PRODUCTIONIZATION = M1 (DONE: core promoted to the
  radia.vim package + public API, feec 85/85) -> M2 (Wilton/phi_tet into the C++ ChargeGram) ->
  M3 (nonlinear Newton in C++, iterative tangent solve) -> M4 (curved/high-order + symmetry + curved
  nonlinear volume) -> M5 (the seal, AFTER the summer 静止器・回転機 meeting).
"""

_CURVED = r"""
# Curved + high-order geometry -- a win HDiv has and flat multipole-moment MMM cannot (2026-06-08)
Canonical reference: docs/hdiv_vim/README.md sec.2-3.  HDiv lives natively on curved (isoparametric)
meshes (mesh.Curve(p), Piola map); flat ObjHexahedron/ObjTetrahedron cannot represent a curved boundary.

## The production Gram = the ngsolve.bem Laplace single-layer
The uniform-M surface demag Gram IS the Laplace single-layer of sigma = M.n; NGSolve 6.2.2604
`ngsolve.bem` supplies it HIGH-ORDER + CURVED + FMM, with NO hand-rolled singular quadrature.
D_axis = <sigma, V.mat sigma>/V_vol, sigma = GridFunction.Set(specialcf.normal(3)[axis],
definedon=mesh.Boundaries(".*")); kernel 1/(4pi r) (no extra factor -> sphere gives 1/3).
Example hdiv_demag_bem_singlelayer.py; golden test_hdiv_vim_bem_demag.py.

## Verified vs ANALYTIC truth: curved + order-2 EXACT, flat floored
- SPHERE demag: FLAT +0.25% at order 0/1/2 (order-insensitive faceting floor) -> CURVE(3) o2 ~1e-4%.
- SPHEROID full tensor: prolate(2:1) N_par 0.17356 + N_perp 0.41322; oblate(1:2) 0.52720 + 0.23640 --
  all EXACT vs Osborn 1945; sum rule N_x+N_y+N_z=1 to ~1e-6.
- GENERAL TRIAXIAL ellipsoid (a!=b!=c): three DISTINCT factors all exact vs the Osborn integral.

## The FIELD win + head-to-head vs the SHIPPED Radia solver
- External field of a uniform / nonlinear soft-iron sphere = the exact dipole; FLAT ~-10% (the volume
  faceting error enters m=M V) -> CURVE(3) <0.4% (~23-38x).  hdiv_demag_curved.py,
  hdiv_curved_nonlinear_field.py.
- HEAD-TO-HEAD (compare_curved_vs_radia_field.py): HDiv CURVED at the COARSEST mesh beats shipped
  Radia-FLAT at the FINEST -> ~10-30x ACCURACY-PER-RESOLUTION.  HONEST SCOPE: accuracy-per-DOF
  (geometry-driven, fair), NOT wall-clock -- the Python prototype is not time-optimized; a fair speed
  comparison needs the C++ productionization.

## CURVED x NONLINEAR (honest magnitude)
The curved win on the MAGNETIZATION is MODEST (~0.3% -- the demag factor is a volume-normalized ratio
that cancels the ~10% volume faceting error).  The win is LARGE on the FIELD (~23x).  Radia cannot
referee curved geometry (it facets), so curved nonlinear is validatable ONLY vs analytic = spheroids.
"""

_SYMMETRY = r"""
# Symmetry models 1/2, 1/4, 1/8 (2026-06-08) -- canonical: docs/hdiv_vim/README.md sec.5
Two pieces make a symmetry model; the HDiv-VIM gets BOTH cheaply (so 1/4 and 1/8 are SUPPORTED):
- LOOPS: AUTOMATIC.  loops = ker(B) on the cut mesh, field-null ~4e-16 BY CONSTRUCTION, count adapts to
  the cut topology (sphere full / 1/2 / 1/4 / 1/8 -> 58 / 54 / 18 / 6).  NO cohomology-aware loop-star
  installCycle -- the "loop removal is painful" problem of multipole-moment MMM MSC is ELIMINATED.
  Golden: test_hdiv_vim_symmetry_loops.py.
- DEMAG VALUE: the IMAGE method.  Only the real surface (spherical cap) carries sigma = M.n = n_z; the
  flat cut faces are symmetry planes (no real charge).  Reflect the cap charge over the reduction planes
  with sign = (-1)^(#z-reflections) (sigma = n_z flips under a z-mirror = the IMA sign rule:
  field-PARALLEL mirror keeps sign, field-PERPENDICULAR flips) -> reconstructs the full sphere.  The
  reduced models reproduce the FULL demag from ~1/N the surface DOF: 1/2 +0.08%, 1/4 +0.11%, 1/8 -0.32%
  (108 / 52 / 20 cap-tris vs the full 192).  Example hdiv_demag_symmetry_image.py; golden
  test_hdiv_vim_symmetry_image.py.  Production version = the ngsolve.bem single-layer + image kernels
  (high-order / curved / scalable).
"""

_CROSS_METHOD = r"""
# Cross-method validation of the demag factor -- three INDEPENDENT discretizations agree (2026-06-10)

The demag TENSOR is the HDiv-type VIM's headline deliverable: the generalized eigenvalues of
(N, M_mass) with N = B^T G B.  The strongest evidence the operator is RIGHT is that THREE
structurally-independent ways of discretizing the SAME magnetostatic demag physics all land on the
SAME analytic (Osborn 1945) numbers -- a cross-METHOD agreement is more decisive than matching the
analytic alone, because three different discretizations cannot share the same bug.

## The three method families (each computes the ellipsoid demag factor N)
1. **FEEC SURFACE-CHARGE (this VIM)** -- M in H(div) RT0, charge map B (sigma=M.n), Coulomb/single-layer
   Gram.  EXACT to ~1e-3 via the ngsolve.bem single-layer (curved + order 2): sphere 1/3, prolate(2:1)
   N_z 0.17356, oblate(1:2) N_par 0.52720, full triaxial with sum rule sum N_i = 1 to ~1e-6 (topic
   "curved").  Loops field-null by construction => mu_r-independent.
2. **VOLUME finite-element (vector-potential A formulation)** -- the body is a RIGID permanent magnet
   (recoil mu_r = 1, so M = Hc is fixed, no induced part); inside an ellipsoid the demag field is
   uniform H = -N M, hence the interior B gives N = 1 - <B_axial>/Br (Br = mu0 Hc).  This is the
   INDEPENDENT volume leg: a totally different formulation (curl-curl A on a volume mesh, not surface
   charge) computing the same N.
3. **BEM surface-charge (boundary integral)** -- the same surface-charge sigma = M.n as (1) but solved
   by a boundary-element method; this is RADIA'S OWN formulation family (surface charge / single layer),
   so it is the most direct cousin.  (The 3D cube / C-yoke head-to-head is the next cross-method block.)

## VERIFIED (volume-FE leg, axisymmetric -- bodies of revolution)
An independent A-formulation volume-FE solve (axisymmetric: sphere / prolate / oblate spheroid)
reproduces the demag factors and AGREES with this VIM and with Osborn:
    shape      c/a   N_volumeFE   N_HDiv-VIM   N_Osborn    err
    sphere     1.00  0.333963     0.333333     0.333333    0.19 %
    prolate    2.00  0.174170     0.173560     0.173564    0.35 %  (axial N_z)
    oblate     0.50  0.528145     0.527200     0.527200    0.18 %
    isotropy sum-rule 3*N_sphere = 1.0019 (exact 1).
=> volume-FE == HDiv-VIM == analytic on the SAME three bodies, all < 0.35 %.

## The ACCURACY RECIPE that transfers (a reusable FE insight)
The demag factor is a B.n SURFACE-CHARGE jump, so it must be resolved on BOTH sides of the body
boundary.  Refining only the magnet mesh leaves the volume-FE error STUCK at ~1.2 % (flat across 16x
body-mesh refinement -- it is NOT body-side discretization); adding a FINE exterior shell (a scaled
co-ellipsoid of the surrounding medium at the same fine mesh) resolves the air side too and drops the
error to < 0.35 %.  This mirrors why the VIM Gram needs an accurate single-layer (the surface-charge
self/near term): the demag factor lives in the charge layer, not the bulk.

## DIMENSIONAL NOTE
The volume-FE leg is 2D/axisymmetric, so its native bodies are surfaces of revolution
(sphere/spheroid).  The 2D-PLANAR analog is the infinite ELLIPTIC CYLINDER magnetized transverse:
N_x = b/(a+b), N_y = a/(a+b), sum N_x + N_y = 1 (the 2D sum rule), circle -> 1/2.  The FEEC VIM and
the BEM leg are native 3D; the analytic Osborn family ties all dimensions together.

## NON-ELLIPSOIDAL bodies -- where the FULL operator beats a single factor
The ellipsoid is special: its interior demag field is UNIFORM, so one number N per axis describes it.
The practical magnet shape is the FINITE CYLINDER, which is NOT an ellipsoid -- its interior field is
NON-uniform, so the demag factor VARIES with position.  The CENTRAL (on-axis-centre) factor still has a
rigid-magnet closed form, N_c(gamma) = 1 - gamma/sqrt(gamma^2+1) with gamma = L/D (the equivalent-solenoid
centre field), verified by an independent volume-FE solve to <0.5% in the typical-aspect regime.  But the
single central (or volume-averaged) number is only a summary: the cylinder has SHARP EDGES where the
surface charge sigma = M.n is singular, and the demag field varies from centre to edge.  This is exactly
where the H(div)-type VIM earns its keep over a tabulated factor -- it returns the full position-dependent
demag OPERATOR (and, with curved/high-order elements, resolves the edge charge that a flat low-order
method smears).  So the demag-factor cross-method check (ellipsoid, exact) validates the operator's
average; the non-ellipsoidal cylinder is where the operator's spatial detail matters.

## WHY THIS IS THE RIGHT VALIDATION
analytic-only agreement can hide a shared analytic-mapping error; agreement of a SURFACE-charge FEEC
method, a VOLUME A-formulation method, and (next) a BEM surface-charge method -- three different
discretizations -- isolates the demag PHYSICS.  The demag factor (Osborn 1945) is textbook and
mu0-independent (a pure geometric ratio), so the public statement is analytic-gated; the specific
cross-method provenance is recorded internally.
"""

_SECTIONS = {
    "overview": _OVERVIEW,
    "implementation": _IMPLEMENTATION,
    "scaling": _SCALING,
    "verification": _VERIFICATION,
    "nonlinear": _NONLINEAR,
    "curved": _CURVED,
    "symmetry": _SYMMETRY,
    "cross_method": _CROSS_METHOD,
    "status": _STATUS,
}


def get_hdiv_vim_documentation(topic: str = "overview") -> str:
    """Return HDiv-type VIM knowledge for the requested topic (see _SECTIONS keys; 'all' = everything)."""
    t = (topic or "overview").strip().lower()
    if t == "all":
        return "\n\n".join(_SECTIONS[k] for k in
                           ("overview", "implementation", "scaling",
                            "verification", "nonlinear", "curved", "symmetry",
                            "cross_method", "status"))
    if t in _SECTIONS:
        return _SECTIONS[t]
    return (f"Unknown topic '{topic}'. Options: " + ", ".join(_SECTIONS.keys()) + ", all.\n\n"
            + _OVERVIEW)
