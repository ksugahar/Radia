"""HDiv-type VIM (Volume Integral Method) demag operator -- knowledge module.

The HDiv-type VIM is the lab's FEEC (H(div) RT0) alternative to the collocation MMM/MSC kernel: a
SYMMETRIC demag operator N = B^T G B whose loop modes are FIELD-NULL BY CONSTRUCTION, giving
mu_r-INDEPENDENT iterative convergence (no near-null blow-up).  This module records the CURRENT
implementation (as of 2026-06-07, on main @ feaade25) so a later session can extend it -- in
particular to NONLINEAR materials -- without re-discovering the architecture.

Exposed as the radia-ngsolve MCP tool `hdiv_vim(topic=...)`.
"""

_OVERVIEW = r"""
# HDiv-type VIM demag operator (the FEEC alternative to MMM/MSC)

## What it is
Magnetization M lives in lowest-order H(div) (RT0): ONE signed normal-flux DOF per face.  The demag
operator is the SYMMETRIC Galerkin form

    N = B^T G B          (symmetric, [n_face x n_face])

  B = charge map   M |-> (rho = -div M per cell [P0],  sigma = M.n per boundary face [P0])
  G = Coulomb Gram (charge-charge interaction; symmetric since 1/r is symmetric)

## Why it is better than the collocation MSC kernel (the whole point)
- **Loops are FIELD-NULL BY CONSTRUCTION**: loops = ker(B) (charge-free fields).  B.loop = 0 =>
  N.loop = B^T G (B loop) = 0 for ANY G.  So the loop (circulating-magnetization) null space sits at
  EXACTLY the material eigenvalue 1/chi and never pollutes the spectrum.
- => **mu_r-INDEPENDENT convergence**: plain-Jacobi MINRES iteration count is FLAT across mu_r
  (measured 80/82/78 across mu_r 10 -> 1e3 -> 1e5 on the 4^3 system).  The high-mu_r conditioning
  wall that caps the collocation MSC loop-star solver (BiCGSTAB iters ~mu_r^1.5) is ABSENT.
- **DISTORTED-ELEMENT ROBUSTNESS (the headline PRACTICAL value, 2026-06-07)**: the field-null-by-
  construction property holds on ANY mesh -- B.loop = 0 is exact on AFFINE *and* NON-AFFINE (distorted)
  hexes alike, so mu_r-independence SURVIVES distortion (measured: distort=0.18 grid, MINRES iters
  98/56/52/52 across mu_r 10 -> 1e4 -- bounded, even DECREASING).  The yano-type combinatorial +/-1
  loops are field-null ONLY on affine hexes (the de-Rham defect: on distorted hexes they carry field);
  the hand-crafted Yano elements / the shipped MSC's installCycle retrofit (~6e-9 local-null-vector)
  exist precisely to patch this.  HDiv-VIM needs NO hand-crafted elements and NO retrofit -- it is
  EXACT (machine, ~4e-16) on distorted meshes by construction.  This is a robustness / generality /
  maintainability win on its OWN (any mesher, any distortion, provably correct; the README "retire the
  yano-type" goal), independent of raw speed -- at mu_r<=1e4 it is at performance PARITY with the
  shipped MSC (same demag spectrum), the win being correctness + no hand-crafting.  Golden-locked:
  tests/feec/test_hdiv_vim_symmetry_golden.py (loops field-null + PSD on distorted) and
  test_hdiv_vim_solve.py::test_minres_iters_bounded_vs_mu_r_distorted (mu_r-independence on distorted).
- **SYMMETRIC** => MINRES (symmetric indefinite Krylov) + symmetric H-LDL^T factorization.

## The material system
    A = (1/chi) M_mass - N       (symmetric INDEFINITE)
  M_mass = RT0 H(div) mass.  The generalized eigenvalues of (N, M_mass) are the demag factors,
  bounded in [0,1] (basis-invariant physics).  Solved by MINRES (Jacobi precond) or H-LDL^T.

## Origin
NGSolve FEEC prototype (examples/feec_vim/hdiv_demag_quad_self.py) -> productionised into the Radia
C++ core.  This is the "de-Rham fix" research line: loops field-null by construction on distorted
hexes AND tets.
"""

_IMPLEMENTATION = r"""
# Current implementation (files + APIs), main @ feaade25

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
- **src/ext/HACApK/cHACApK_harith.{c,h}** -- rk-aware symmetric **H-LDL^T** (factor compressed
  H-matrices): cHACApK_hldlt_decomp / _solve_vec, the cHACApK_hldlt_factor_leafmtxp / _apply / _free
  driver (symmetric mirror of the H-LU driver).

## pybind (radia._radia_pybind)
- `_hdiv_vim_assemble(nx,ny,nz,nsub=0,distort=0)` -> dict{nf,n_cell,n_charge,n_bnd, N, B, M_mass}
- `_hdiv_vim_hmatrix_probe(nx,ny,nz,nsub,distort,eps,leaf,eta)` -> stats + matvec_relerr vs dense N
- class `_HDivVimHMatrix(nx,ny,nz,nsub,distort,eps,leaf,eta)` -- build-once: .ndof() / .matvec(x) (N x)
  / .apply_system(x,inv_chi) / .diag_system(inv_chi) / .stats()
- class `_ChargeGramHMatrix(centroids, measures, self_energy, eps,leaf,eta)` -- .ndof() / .matvec(q)
  (G q) / .stats() / .factor_solve_hldlt(b) (factor G with H-LDL^T + solve G x=b)
- `_hldlt_self_test(depth,nb)` / `_hldlt_self_test_rk(nb,rk,depth)` / `_hldlt_self_test_rk_mixed(...)`

## Python (unstructured tet ingest -- the real-geometry path)
- **examples/feec_vim/hdiv_demag_tet.py** -- NGSolve HDiv(order=0) tet ingest.  Element-AGNOSTIC
  extraction: B = -div(u) [vol charge] + u.Trace().n [surface charge]; M_mass = HDiv mass.  Only the
  Gram self-energy geometry is element-specific -> tet/tri barycentric sub-points (c_tet~1.776 /
  c_tri~2.888).  `build_demag(mesh,nsub)` -> N (dense, monopole off-diag + sub-point self), M_mass, B,
  charge geometry.  `build_near_correction(mesh, d, nsub, near_factor)` -> sparse near-field Gram
  correction (exact sub-point MINUS monopole) -- the scalable Gram = monopole H-matrix + this.
"""

_SCALING = r"""
# Scaling: charge-Gram H-matrix + near-field correction

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

_HLDLT = r"""
# rk-aware symmetric H-LDL^T (factor the compressed H-matrix)

cHACApK_hldlt_* factors a SYMMETRIC (indefinite) H-matrix as A = (I+W) D (I+W)^T (Bunch-Kaufman 2x2
pivots, lower triangle only -> ~half the off-diagonal storage + FLOPs of the non-symmetric H-LU).

rk-aware off-diagonal LEAVES (the ACA-compressed far blocks):
  (A) rk offdiag solve W_ji = U (A_ii^{-1} V)^T  (apply A_ii^{-1} to V in place; U unchanged)
  (B) snapshot V before the solve overwrites it (U read live)
  (C) trailing update -W_ji A_ki^T as a rank-kinc increment via add_lowrank_to_node, 4 operand combos
      (dense*dense rank ndiag, rk*dense AVw=A_ki Vw, dense*rk WVa=W_ji Va, rk*rk M=Vw^T Va then Uw M)
Self-tests: MACHINE precision (rk depth-1 6.4e-16, depth-2 7.0e-16, mixed-flat all-4-combos 2.4e-15).

Driver: cHACApK_hldlt_factor_leafmtxp(leafmtxp, control, nffc, *out_rc) -> opaque root;
cHACApK_hldlt_apply(root, control, r, z, nd); cHACApK_hldlt_free_factors.  Exposed as
_ChargeGramHMatrix.factor_solve_hldlt(b).

VERIFIED: factors a real HACApK Gram H-matrix end-to-end (single leaf + 7-leaf diagonal-refined tree),
solve G x=b rel err ~1.6e-15.

HONEST BOUNDARY (the remaining increment): a DEEP tree (small leaf -> the natural HACApK build creates
INTERNAL off-diagonal blocks) returns NEED_RECURSIVE (-5), fail-loud no-fallback.  The lower-only
H-LDL^T does NOT yet recurse internal off-diagonal blocks -- the symmetric analog of the H-LU's
recursive htrsm/h_addmul is future work for fully-compressed deep trees.  (Not blocking: MINRES
already solves the system mu_r-independently; the direct factor is an optimization.)
"""

_VERIFICATION = r"""
# Verification (golden tests, tests/feec/, 45/45)

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
- test_hldlt_factorization.py -- dense + rk-aware H-LDL^T self-tests (machine precision, all 4 combos).
- test_hldlt_real_gram.py -- H-LDL^T factors a real Gram H-matrix; deep tree -> NEED_RECURSIVE (boundary).

Bug caught by verify-first (2026-06-07, commit 5d7a9823): the Python tet/tri barycentric sub-point
lattice summed to 1+c/nsub (sub-points OUTSIDE the simplex).  Fixed to proper lattices (sum to 1,
nsub=1 -> centroid).  The C++ hex/quad path was unaffected.  The self constants c_tet/c_tri are
robust to the scheme (cross-sum extrapolation -> same INT INT 1/r).
"""

_NONLINEAR = r"""
# NONLINEAR HDiv-type VIM -- ROBUST SOLVER = DAMPED NEWTON-RAPHSON (SOLVED 2026-06-07)

Goal: make the HDiv-type VIM work for NONLINEAR soft-magnetic materials (BH curve / saturation).

## HEADLINE (2026-06-07): the robust strong-saturation solver is DAMPED NEWTON-RAPHSON.
The earlier sessions' simple per-element Picard / Hantila FAILED at deep saturation (NaN, wrong root);
damped Newton-Raphson with the CONSISTENT TENSOR tangent SOLVES it -- fast (4-11 iters) and accurate
(rel 1e-4 to 1e-6 vs the analytic uniform sphere) at every saturation level.  This is the Newton
counterpart of Radia's existing MMM/MSC newton_damping=True path: the OUTER Newton machinery (tangent
susceptibility + line-search damping) is SHARED; only the demag operator N = B^T G B differs.
Implemented as solve_nonlinear_newton (examples/feec_vim/hdiv_demag_tet_nonlinear.py), golden-locked in
tests/feec/test_hdiv_vim_tet_newton.py (3 tests, feec suite 50/50).

## DONE + golden-tested (examples/feec_vim/hdiv_demag_tet_nonlinear.py, tests/feec/test_hdiv_vim_tet_nonlinear.py)
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
solve_nonlinear_newton_scalable (examples/feec_vim/hdiv_demag_tet_nonlinear.py) replaces the dense
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

## NEXT (open): analytic VOLUME (tet) Gram; full-C++ Newton loop; Wilton-in-C++
- analytic VOLUME Gram (uniform-tet potential: cell-cell + cell-face) -> closes the residual nonlinear
  sharp-body non-uniform gap (the surface Wilton is done).
- the C++ Gram H-matrix is MONOPOLE far field, so the scalable path matches the dense MONOPOLE+near-corr
  Newton, NOT the dense Wilton path; putting the Wilton surface integral into the C++ near-field
  correction gives scalable + Wilton-accurate together.
- a full C++ Newton OUTER loop (the tangent + line search in C++ too) removes the Python-loop overhead
  (minor: each iteration is dominated by the C++ H-matvec); reuse the #2 H-LDL^T for a factor-once
  variant of the per-iteration solve.

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
  reuses the #2 H-LDL^T): DIVERGED in the RT0 face metric (alpha=chi0/2 -> NaN; under-relaxed -> ~1e88).
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
per-element tensor tangent + factor-once M_mass reuses the #2 H-LDL^T; (3) optionally a consistent-RT0
Hantila for the factor-once speed-up at the knee (Newton already gives correctness).
"""

_STATUS = r"""
# Status summary (2026-06-07)

DONE + golden-locked (feec 54/54):
  #1  scalable mu_r-independent HDiv-VIM demag solver on REAL tet meshes (Layer A/A.5 + tet ingest)
  #2  rk-aware symmetric H-LDL^T factoring real compressed H-matrices (+ driver)
  #3  bug-fixed exact Gram via near-field correction -> demag -> analytic 1/3
  NONLINEAR  damped Newton-Raphson (consistent tensor tangent + near-corr + line search + Picard
             warmstart) -- robust + fast at deep saturation where per-element Picard/Hantila failed
             (examples/feec_vim/hdiv_demag_tet_nonlinear.py::solve_nonlinear_newton,
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

OPEN (honest boundaries / next increments):
  - analytic VOLUME (tet) Gram: the surface Wilton is done; the residual nonlinear SHARP-body
    non-uniform gap (cube moderate ~8.7% vs Radia) is the cell-cell / cell-face Gram (still
    monopole+sub-point).  The next refinement.
  - further #3 validation options: a real rad.MatSatIsoTab steel TABLE (needs table-based M(H)/M'(H)
    in the Newton) and a C-yoke.
  - BH-knee stiffness: Newton converges to the correct answer but slowly there (scalar Picard is the
    practical tool at the knee); operator-accuracy-limited (finer/analytic Gram reduces it).
  - the scalable C++ nonlinear path (reuse #2 H-LDL^T for the factor-once tangent solve) -- new #2.
  - H-LDL^T on DEEP trees (internal off-diagonal recursion) -- currently NEED_RECURSIVE.
  - near-field / Wilton Gram in the C++ ChargeGram entry -- currently the Python overlay.
  - proper distorted RT0 M_mass for exact distorted demag VALUES (non-negativity already holds).
"""

_SECTIONS = {
    "overview": _OVERVIEW,
    "implementation": _IMPLEMENTATION,
    "scaling": _SCALING,
    "hldlt": _HLDLT,
    "verification": _VERIFICATION,
    "nonlinear": _NONLINEAR,
    "status": _STATUS,
}


def get_hdiv_vim_documentation(topic: str = "overview") -> str:
    """Return HDiv-type VIM knowledge for the requested topic (see _SECTIONS keys; 'all' = everything)."""
    t = (topic or "overview").strip().lower()
    if t == "all":
        return "\n\n".join(_SECTIONS[k] for k in
                           ("overview", "implementation", "scaling", "hldlt",
                            "verification", "nonlinear", "status"))
    if t in _SECTIONS:
        return _SECTIONS[t]
    return (f"Unknown topic '{topic}'. Options: " + ", ".join(_SECTIONS.keys()) + ", all.\n\n"
            + _OVERVIEW)
