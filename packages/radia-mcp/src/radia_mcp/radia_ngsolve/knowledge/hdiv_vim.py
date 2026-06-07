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
# NONLINEAR HDiv-type VIM -- the NEXT step (plan, NOT yet implemented)

Goal: make the HDiv-type VIM work for NONLINEAR soft-magnetic materials (BH curve / saturation),
not just linear demag.

## Why it should be clean
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

## Concrete options to evaluate (verify-first before committing)
- **Picard / fixed-point**: cheapest.  Per iteration, update (1/chi_e) on the M_mass diagonal and
  re-solve A m = b.  N (and its H-matrix) is UNCHANGED -> only the diagonal term changes; MINRES with
  the updated Jacobi diag (DiagSystem(inv_chi_e)) re-solves.  The H-matrix is built ONCE.
- **Hantila polarization** (already in Radia for MMM, src/radia/hantila_solver.py): split B =
  mu_0(1+alpha)H + mu_0 R; the LHS (I - alpha N) is CONSTANT -> factor ONCE (LU/H-LDL^T), back-
  substitute per iteration.  The HDiv symmetric N + H-LDL^T factor is a natural fit (factor once,
  apply many) -- this is likely the strongest path and reuses the #2 H-LDL^T payload.
- **Newton**: per-element tangent dM/dH; faster convergence, needs the Jacobian assembly.

## What exists to build on
- The per-element field/charge machinery (B, N, M_mass) is all in place + golden-tested (linear).
- DiagSystem(inv_chi) already produces the Jacobi diagonal for a given chi -- the per-element-chi
  update hook.  ApplySystem(x, inv_chi, y) is the per-iteration operator apply.
- Radia already has BH-curve materials (rad.MatSatIsoTab) + the Hantila solver (MMM 3DOF only today;
  the HDiv RT0 face DOF is the new target) + Play/Energy hysteresis -- reuse the constitutive layer.
- Validate against: a known nonlinear demag (e.g. a saturating sphere/cube in a strong applied field),
  and/or the existing MMM/MSC nonlinear solve on the same geometry.

## First increment suggestion
Pick ONE material model (Picard + rad.MatSatIsoTab BH curve), wire the per-element chi(H) update around
the existing _HDivVimHMatrix.apply_system / diag_system, on the STRUCTURED hex path first (golden-
testable), verify the saturated demag against MMM/MSC; THEN the tet path; THEN Hantila/H-LDL^T factor-
once for speed.
"""

_STATUS = r"""
# Status summary (main @ feaade25, 2026-06-07)

DONE + golden-locked (feec 45/45), the production sequence #1 -> #2 -> #3:
  #1  scalable mu_r-independent HDiv-VIM demag solver on REAL tet meshes (Layer A/A.5 + tet ingest)
  #2  rk-aware symmetric H-LDL^T factoring real compressed H-matrices (+ driver)
  #3  bug-fixed exact Gram via near-field correction -> demag -> analytic 1/3

OPEN (honest boundaries / next increments):
  - NONLINEAR materials (see topic "nonlinear") -- the user's next direction.
  - H-LDL^T on DEEP trees (internal off-diagonal recursion) -- currently NEED_RECURSIVE.
  - near-field correction in the C++ ChargeGram entry (needs cell/face geometry in the manager;
    or analytic Wilton) -- currently the Python build_near_correction overlay.
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
