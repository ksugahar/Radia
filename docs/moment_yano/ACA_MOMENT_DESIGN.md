# Scalable moment-yano: ACA + H-LU design (Phase 2 of the EIEM2 full-deletion track)

**Status (2026-06-23):** Phase 2 (scalable method 2) is COMPLETE and Phase 3 has landed. moment-yano is the
canonical surface-charge soft-iron demag formula for hex, wedge, and pyramid polyhedra; the EIEM2 collocation
kernels and their `SolverConfig` research opt-outs are removed. `rad.Solve(..., method=2)` solves the moment
system via the HACApK moment H-matrix + block-Jacobi BiCGSTAB for linear and nonlinear (per-element chi)
materials, with **O(N log N) storage** (the three dense O(N^2) buffers -- interaction N, BaseMatrix, dgesv
SystemMatrix -- are all removed from the method-2 path). Measured: at dof=12288 method 2 uses 673 MB vs
method 0's 5461 MB (ratio 0.12, sub-quadratic; `examples/vim/bench_moment_storage_scaling.py` + `.json`),
external B matches method 0 to ~1e-10, and the nonlinear C-yoke saturates identically (same Picard iters).

## What was built

A scalable solver for the moment-yano linear system so that `rad.Solve(..., method=2)` (HACApK) on pure
6-DOF hex soft iron solves the **moment** system (not EIEM2) with `O(N log N)` storage + matvec. The shipped
preconditioner is block-Jacobi, so iterations still grow with the high-mu demag conditioning wall; a pivoted
H-factor or symmetrized moment formulation remains future work. The nonlinear outer loop (Picard / Anderson)
wraps this linear solve unchanged.

TaskManager parallelism is part of the shipped method-2 contract. `RadHACApKBase::BuildHMatrix` and
`SolveMomentHACApK` stand up or reuse an NGSolve `RegionTaskManager`, and the HACApK C callbacks enter
NGSolve `ParallelFor` through `rad_hacapk_parallel.cpp`. Direct diagnostic `MatVec` calls are
TaskManager-preconditioned and should be made under a caller `TaskManager` scope.

## The validated picture (from the lab prototypes, all on `main`)

| Prototype (`examples/vim/`) | Result |
|---|---|
| `yano_moment_hmatrix_compressibility.py` | Gate 1 PASS: the nonlocal moment kernel `D` (centroid field+grad coupling) has **bounded ACA rank** (field ~13, grad ~16-24), `N`-independent -> H-compressible. |
| `yano_moment_matfree_solve.py` | Gate 2 PASS: matrix-free moment matvec reproduces the dense direct solve to `<1e-6` -> swapping the dense matvec for the HACApK matvec keeps the same answer. Block-Jacobi does NOT bound iters (grow `~dof^1.06`, mu_r-contrast driven). |
| `yano_moment_scalable_path.py` | The A-build kernel is **cheap** (`~0.9 us`/(elem,face), single centroid->face integral) -- lighter than HDiv's face-face charge-Gram. Cheap preconditioners (block-Jacobi, sparse-LU top-k, deflation) all FAIL in the prototype, which motivated trying H-LU; the later no-pivot HACApK diagnostic below showed that H-LU is not shippable for this non-symmetric moment matrix. |

Net: **SCALABLE moment-yano = HACApK A-build (cheap entries) + H-matvec BiCGSTAB + block-Jacobi today.**
Bounded iterations remain future work (pivoted H-factor or a symmetrized moment formulation).

## The system

Per moment element, the rows are (`BuildMomentSystemCore`, `rad_interaction.cpp`):
3 dipole, 1 monopole, and the residual quadrupole rows needed to make the block
square.  Hex has 6 face-charge DOF (2 quadrupole rows); wedge/pyramid have 5 face-charge
DOF (1 residual quadrupole row). The dense system is

    A_raw = L  -  chi * R * C

- `L`  : the per-element LOCAL geometric-moment block (block-diagonal; 6x6 for hex,
         5x5 for wedge/pyramid; cheap).
- `C`  : the centroid field+grad coupling `C[e,k,g]` = field/grad component `k`
         (`k<3` = H, `3..8` = gradH) at element `e`'s centroid from unit charge on face
         DOF `g` (`BuildCentroidFieldGrad`; the cheap kernel, now IMA-aware).
- `R`  : the per-row linear combination (dipole = combo of `C[e,0:3]`, quad = combo of
         `C[e,3:9]`); sparse + local to element `e`.
- `chi`: susceptibility (per-element under Picard).

`A_raw`'s off-diagonal block for well-separated element clusters is `-chi*R*C`, a smooth
field/grad kernel folded by local moment functionals -> **low-rank** (Gate 1). So `A_raw`
is an H-matrix with the cluster tree over element centroids.

### KEY INSIGHT -- no row normalization needed for the H-LU path

`BuildMomentSystemCore` 2-norm-normalizes each row (for the dense GMRES/block-Jacobi
conditioning). **Row normalization is a diagonal row-scaling `D` and does NOT change the
exact solution of a direct solve:** `A_norm x = b_norm` with `A_norm=D A_raw`,
`b_norm=D b_raw`, is `A_raw x = b_raw` (multiply both sides by `D^{-1}`). So the H-LU
path builds the **un-normalized `A_raw`** and gets the SAME `x` as the normalized dense LU
(method 0) -- to solver tolerance. This removes the `O(N^2)` exact-row-norm precompute
(the row norm needs the dense field/grad part of every row); the H-matrix entry
`A_raw[i][j]` is then computable **on demand** from element geometry + the on-demand
`C[e_i, k, j]` (a single centroid->face evaluation).

## The C++ template to follow

`RadHACApKHDivSystemTet` (`rad_hacapk_hdiv.{h,cpp}`) already does exactly this shape for
the HDiv-VIM: it builds the system `A = M_mass + chi*N` as a HACApK H-matrix
(`RadHACApKBase` subclass with an on-demand `ComputeSystemEntry(i,j)`) and applies the
HACApK **H-LU** (`cHACApK_hlu_*`) as a scalable direct solve / strong preconditioner.
The moment manager mirrors it for the current scalable method-2 pure-hex path:

```
class RadHACApKMomentSystem : public RadHACApKBase {
  // ExtractCoordinates(): cluster-tree points = element centroids, expanded to the
  //   6 row-DOF (and 6 col-DOF) per hex co-located at the centroid (dof = 6*nHex).
  // ComputeSystemEntry(i,j) = A_raw[i][j]:
  //   decode i -> (element e_i, row-type in {dip_x,dip_y,dip_z,mono,quad0,quad1});
  //   decode j -> face DOF g (element e_j, local face);
  //   local part: if g in e_i's faces, add the L geometric-moment coefficient;
  //   nonlocal part: subtract chi * (row-type functional applied to C[e_i, :, g]),
  //     C[e_i,k,g] evaluated ON DEMAND (single centroid->face 8x8 Gauss, IMA mirrors).
  //   NO row normalization (see KEY INSIGHT).
  // SetSystemMode(chi) + cHACApK_hlu_* -> H-LU factor of A_raw -> scalable solve.
};
```

## Increments + verification gates

1. **On-demand entry.** A C++ `MomentSystemEntry(e_i, row_type, g)` (extract the
   single-(target,source) field/grad from `BuildCentroidFieldGrad`'s inner loop + the
   `BuildMomentSystemCore` row construction, un-normalized). **Gate:** entry-by-entry ==
   the dense `A_raw` from `BuildMomentSystemCore` (drop its row-norm) to machine precision.
2. **H-matrix build + matvec.** `RadHACApKMomentSystem` over element-centroid clusters.
   **Gate:** H-matvec `A_raw @ x` == dense `A_raw @ x` to ACA tolerance; ACA rank bounded
   in `N` (re-confirms Gate 1 in C++); build sub-cubic (Benchmark Policy JSON).
3. **~~H-LU solve~~ -> BLOCKED; use a Krylov solve on the H-matvec (REVISED 2026-06-22).**
   **FINDING (Increment 2.5 de-risk + the no-pivot diagnostic):** the HACApK H-LU
   (`cHACApK_hlu_*`) is **NO-PIVOT** (the HDiv template's `A = M_mass + chi*N` is SPD, so
   no-pivot is stable there).  The moment `A_raw` is **NON-symmetric** and, in the natural
   (element) ordering, the no-pivot factorization hits a **near-zero pivot**
   (`min|U_ii|/max|U_ii| ~ 2e-16`) even though the full matrix is well-conditioned
   (`cond(A_norm) ~ 1.6e3`) -- i.e. it is a PIVOTING (ordering) problem, NOT a scaling one.
   Measured: the no-pivot H-LU round trip is 8e-7 at dof=336 (all-dense, no truncation) but
   degrades to 6.7e-2 (1080) and DIVERGES to 1.4e+4 (2760) once ACA low-rank truncation
   compounds the near-zero pivot.  **Row/col equilibration does NOT fix it** (the near-zero
   pivot persists; a Python no-pivot LU on A_raw / A_norm / Ruiz-equilibrated all keep
   `min|U|/max|U| ~ 2e-16`).  dense LU (method 0) is fine because `dgesv` PIVOTS.  The
   prototype's "moment scales via H-LU" was an inference (cheap iterative preconditioners
   fail) that was never checked against the actual no-pivot H-LU -- it does not hold.
   **REVISED Increment 3:** solve method 2 with **GMRES/BiCGSTAB on the EXACT moment
   H-matvec (Increment 2, scalable storage) + a block-Jacobi preconditioner** (invert each
   element's local 6x6 `A_raw` block).  This avoids the factorization entirely; the
   matfree prototype (`yano_moment_matfree_solve.py`) validated it converges (== dense to
   <1e-6) with iters that GROW ~dof^1.06 (the high-mu_r demag conditioning wall, shared with
   yano-MSC/HDiv-VIM -- a documented caveat, not a moment defect, NOT bounded like a true
   H-LU would give).  Route `radTRelaxationMethNo_2` (moment-eligible) to it; drop `Error204`.
   **Gate:** method-2 moment `x` == method-0 moment `x` to solver tol; storage scales (the
   dense matrix is gone); iter-growth documented.  Bounded-iter (a PIVOTED H-matrix factor,
   or a symmetrized moment formulation) stays FUTURE work.
4. **Nonlinear + storage decoupling. -> DONE (2026-06-22).** Two parts:
   - *Nonlinear (per-element chi):* `RadHACApKMomentSystem` gained a per-element-chi ctor; `SolveMomentHACApK`
     takes the `chiPerHex` vector (RHS `b[6h+t]=chi_h*Hext_h[t]`, per-element block-Jacobi), and the moment
     branch dropped its uniform-chi guard.  The Picard outer loop re-solves the H-system each iteration with
     the current chi -- the entry `MomentSystemEntry` already folds the row element's chi.  **Gate MET:**
     nonlinear C-yoke (MatSatIsoTab, driven to ~94-95% of Msat) -> method 2 == method 0 (external B ~1e-10,
     same Picard iteration count).  `tests/feec/test_moment_yano.py::test_method2_nonlinear_matches_method0`,
     `examples/vim/verify_moment_nonlinear.py`.
   - *Storage decoupling (closes the Increment-3 storage gate that was a caveat):* the method-2 path no longer
     builds ANY dense O(N^2) buffer.  `SolveGen` sets `skipDenseMatrix=1` for all method 2 except B-input
     Newton/Hantila (no dense interaction N); `radTRelaxationMethNo_0::NeedsDenseMatrix()` returns false when
     `g_yano_moment_hacapk` (no BaseMatrix); `SolveLinearStep` lazy-allocates the dgesv SystemMatrix (never
     reached on the H-BiCGSTAB happy path).  `Setup(skipDenseMatrix)` calls `PrecomputeHexaGeometry()` so the
     moment solve still sees the hexes (the index map normally built inside the skipped dense assembly).
     **Gate MET:** `bench_moment_storage_scaling.py` -- method2/method0 peak memory 0.69 -> 0.32 -> 0.18 ->
     0.12 across dof 1536..12288 (method 0 grows ~N^2, method 2 sub-quadratic ~N log N).
   - *Residual caveat (unchanged from Increment 3):* iters still GROW with N (block-Jacobi only; the high-mu_r
     demag conditioning wall shared with yano-MSC / HDiv-VIM).  Bounded-iter (a pivoted H-factor or a
     symmetrized moment formulation) stays FUTURE work; it does not block Phase 3.

## Phase 3 (DONE, 2026-06-23): EIEM2 deleted -- moment-yano is canonical

With method 0/1/2 + IMA + nonlinear all on moment, the EIEM2 surface-charge collocation
kernel was REMOVED (live/dead refactor, commits bf4424d9/99556872/15f17022/d8d4ef99):
`radTInteraction::Compute6x6/5x5/MixedBlockFast`, the former `RadHACApKMSCManager` MSC machinery
(renamed `RadHACApKMMMManager` and now MMM-3x3-only), the dead IMA-mirror block, the per-face eval-point
caches, and the
`g_yano_eval_alpha`/`g_yano_no_center_charge`/`g_yano_pyramid_cloud` research flags (+ their
`SolverConfig` kwargs).  `MscEvalPoint` (alpha=0.5 midpoint) was ALSO deleted (final polish
`6b617a91`): the per-face MSC external-field branches (`dof>=5`/`dof==6`) in `SetupExternFieldArray`
/ `AddExternFieldFromMoreExtSource` are gone (only the tet `dof==3` fill remains), so `MscEvalPoint`
appears in NO source file.  The moment formulation samples the applied field at the element CENTROID
(the moment RHS via `BuildCentroidFieldGrad`), not per-face -- the per-face MSC fill was immaterial
to the converged moment result (it fed only the initial-H guess).  moment-yano (`BuildMomentSystemCore`
dense / `RadHACApKMomentSystem` H-matrix) is the SOLE surface-charge demag = the canonical
radia MMM for hex/wedge/pyramid soft iron (tet stays 3-DOF MMM).  The quadrupole rows are
the per-element residual eigenmodes (see `examples/vim/eigenmode_quadrupole_derivation.wls`).

## References

- `examples/vim/yano_moment_{hmatrix_compressibility,matfree_solve,scalable_path}.py` (+ `.json`)
- `src/core/rad_hacapk_hdiv.{h,cpp}` (`RadHACApKHDivSystemTet`, the H-LU template)
- `src/core/rad_interaction.cpp` (`BuildCentroidFieldGrad`, `BuildMomentSystemCore`)
- `src/ext/HACApK_LH-Cimplm/` (`cHACApK_hlu_*` H-LU machinery)
