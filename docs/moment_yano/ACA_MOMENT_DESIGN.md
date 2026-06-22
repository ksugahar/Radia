# Scalable moment-yano: ACA + H-LU design (Phase 2 of the EIEM2 full-deletion track)

**Status (2026-06-22):** design locked, implementation pending. moment-yano is the
DEFAULT 6-DOF hex soft-iron demag formula (Steps 3-4) and now covers IMA (Phase 1).
The last parity item before EIEM2 can be deleted is **method 2 (HACApK) at large N** --
this document is the validated plan for it. Decision (user, 2026-06-22): *defer the
EIEM2 deletion until this lands* (keep EIEM2 as the `yano_moment=False` opt-out
meanwhile).

## What we are building

A scalable solver for the moment-yano linear system so that `rad.Solve(..., method=2)`
(HACApK) on pure 6-DOF hex soft iron solves the **moment** system (not EIEM2) with
`O(N log N)` storage + matvec and **bounded iterations**. The nonlinear outer loop
(Picard / Anderson) wraps this linear solve unchanged.

## The validated picture (from the lab prototypes, all on `main`)

| Prototype (`examples/vim/`) | Result |
|---|---|
| `yano_moment_hmatrix_compressibility.py` | Gate 1 PASS: the nonlocal moment kernel `D` (centroid field+grad coupling) has **bounded ACA rank** (field ~13, grad ~16-24), `N`-independent -> H-compressible. |
| `yano_moment_matfree_solve.py` | Gate 2 PASS: matrix-free moment matvec reproduces the dense direct solve to `<1e-6` -> swapping the dense matvec for the HACApK matvec keeps the same answer. Block-Jacobi does NOT bound iters (grow `~dof^1.06`, mu_r-contrast driven). |
| `yano_moment_scalable_path.py` | The A-build kernel is **cheap** (`~0.9 us`/(elem,face), single centroid->face integral) -- lighter than HDiv's face-face charge-Gram. Cheap preconditioners (block-Jacobi, sparse-LU top-k, deflation) all FAIL; the long-range demag coupling MUST stay in the preconditioner. **H-LU is necessary, sufficient, and affordable** (cheap factor + cheap build). |

Net: **SCALABLE moment-yano = HACApK A-build (cheap entries) + H-LU factor -> bounded iters.**

## The system

Per hex (6 face-charge DOF sigma), the moment rows are (`BuildMomentSystemCore`,
`rad_interaction.cpp`): 3 dipole, 1 monopole, 2 diagonal-quadrupole. The dense system is

    A_raw = L  -  chi * R * C

- `L`  : the per-element LOCAL geometric-moment block (block-diagonal, 6x6 per hex; cheap).
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
The moment manager mirrors it:

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
4. **Nonlinear.** Picard/Anderson outer loop reuses the linear Krylov solve per step
   (chi per-element). **Gate:** C-yoke saturation == method-0 moment.

## After Phase 2 (-> Phase 3, EIEM2 deletion)

With method 0/1/2 + IMA + nonlinear all on moment, EIEM2 (the `Use6DOF_MSC` eval-point
collocation kernel + `g_yano_eval_alpha`/`g_yano_no_center_charge`/`g_yano_pyramid_cloud`
research flags + the `yano_moment=False` opt-out) can be removed and the yano goldens
re-locked to moment-only. See the full-delete task list.

## References

- `examples/vim/yano_moment_{hmatrix_compressibility,matfree_solve,scalable_path}.py` (+ `.json`)
- `src/core/rad_hacapk_hdiv.{h,cpp}` (`RadHACApKHDivSystemTet`, the H-LU template)
- `src/core/rad_interaction.cpp` (`BuildCentroidFieldGrad`, `BuildMomentSystemCore`)
- `src/ext/HACApK_LH-Cimplm/` (`cHACApK_hlu_*` H-LU machinery)
