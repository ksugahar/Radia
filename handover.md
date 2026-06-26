# Handover: Multipole-Moment MMM Iterative Solver Design

Date policy: 2026_06_26

This handover is for Codex / Claude Code continuation.  Claude Code may be
touching source files, so treat this file as the coordination note.  Do not
edit source just to "clean up" while another agent is active.

## Current Solver Map

- `method 0`: dense direct LU.
- `method 1`: matrix-free multipole-moment BiCGSTAB / Krylov with element
  block-Jacobi.
- `method 2`: HACApK H-matrix matvec with Krylov and element block-Jacobi.

The current convergence problem is not mainly a moment-kernel speed problem.
The raw multipole-moment operator is harder for Krylov than the old yano-type
face-collocation operator.  The local element block-Jacobi preconditioner is
the best local block-diagonal choice, but it does not bound the iteration count.

## Key Finding: Loop Modes

Loop modes are part of the convergence problem:

- The centroid field/gradient demag operator is essentially blind to loop modes:
  `D @ Lb ~= 0`.
- With block-Jacobi, the preconditioned operator collapses loop modes toward
  small eigenvalues of size roughly `1 / chi`.
- Therefore the block-Jacobi preconditioned spectrum develops a near-zero loop
  cluster at high `mu_r`.

However, loop modes are not the whole story.  A previous full loop deflation
test improved GMRES only modestly (`174 -> 160` on the larger C-yoke case), so
the remaining bottleneck is the non-normal demag-complement coupling.

## ILU Experiment: What Was Proven

Prototype:

- `C:\temp\mmmm_iterative_no_hacapk.py`
- Uses dense/stored matrix matvec, GMRES, and SciPy `spilu` as the dense-side
  analog of an approximate factorization / H-LU preconditioner.
- C-yoke case: `nxy=12`, `nz=2`, `dof=648`, `nLoop=111`.

Observed GMRES iterations:

| `mu_r` | Preconditioner | Iterations | Residual |
|---:|---|---:|---:|
| 100 | block-Jacobi | 82 | `6.5e-11` |
| 100 | ILU, threshold `1e-4` | 5 | `6.2e-12` |
| 1000 | block-Jacobi | 117 | `7.8e-11` |
| 1000 | ILU, threshold `1e-4` | 8 | `4.1e-11` |

Follow-up eigen / loop-overlap check:

| `mu_r` | Preconditioner | `#(|lambda| < 0.05)` | `min |lambda|` | Loop overlap |
|---:|---|---:|---:|---:|
| 100 | block-Jacobi | 87 | `1.778e-2` | `1.000` |
| 100 | ILU, threshold `1e-4` | 0 | `9.878e-1` | n/a |
| 1000 | block-Jacobi | 111 | `1.792e-3` | `1.000` |
| 1000 | ILU, threshold `1e-4` | 0 | `8.883e-1` | n/a |

Interpretation:

- ILU does not physically remove loop modes.
- ILU makes loop modes non-bad for Krylov by moving the near-zero loop cluster
  back near eigenvalue 1 in `M^{-1} A`.
- The good iteration count comes from approximate factorization treating both
  loop modes and the non-normal demag-complement, not from loop deflation alone.

Caveat:

- This ILU was strong.  The input sparsified matrix was about `56%` nonzero, and
  `L+U` was about `88%` dense-equivalent.  This proves the factorization idea,
  not yet that the preconditioner is cheap or scalable.

## Design Direction

The convergence fix is not to keep strengthening matrix-free `method 1` with
local tricks.  The design direction is:

```text
linear solver:   GMRES
matvec:          stored matrix product
                 dense for prototype, HACApK H-matrix for scalable path
preconditioner:  approximate factorization
                 ILU for dense prototype, H-LU / robust H-factor for HACApK
```

Important separation:

- Main operator `A_H`: the operator whose action defines the linear system.
- Preconditioner `M_H`: an approximation to `A_H` used only through
  `M_H^{-1}`.

They must agree on:

- DOF ordering.
- row/column cluster ordering or an explicit permutation map.
- geometry, chi values, and moment formulation.
- matrix dimensions and row scaling convention.

They may differ in:

- ACA+ tolerance.
- max rank.
- leaf size.
- admissibility parameter `eta`.
- far-block approximation accuracy.

Use separate H-matrix instances for `A_H` and `M_H`.  H-LU factorization may
modify the tree/storage in place, so it must not destroy the matvec operator.

If `M_H` is fixed during one linear solve, ordinary GMRES is valid.  If the
preconditioner changes inside a linear solve, use flexible GMRES instead.

## Method-1 Boundary

Pure matrix-free `method 1` cannot use an ILU/H-LU preconditioner without
storing some approximation of the matrix.  A hybrid method could keep the main
matvec matrix-free while storing a sparse/H-matrix preconditioner, but that is
a different memory contract from the current method-1 design.

Therefore, do not present ILU/H-LU as a small tweak to matrix-free method 1.
It is a new stored-matrix-plus-factorization-preconditioner path.

## What Is Not Yet Proven

Do not claim that HACApK H-LU is fully validated yet.  What is proven is:

- approximate factorization fixes the loop near-zero cluster;
- approximate factorization can bound GMRES iterations on the dense prototype;
- block-Jacobi and loop-only deflation are insufficient.

Still to prove:

- the same works on `A_raw`, not only the row-normalized prototype matrix;
- the required fill/rank stays acceptable as `N` grows;
- HACApK H-LU is robust for the non-symmetric moment matrix;
- factorization time and memory beat dense LU at the target sizes.

There is an older warning in `docs/multipole_moment_mmm/ACA_MOMENT_DESIGN.ipynb`:
no-pivot H-LU had trouble on non-symmetric `A_raw`.  Treat that as a real risk.
The H-LU preconditioner implementation may need pivoting, row scaling, or a
different robust H-factorization strategy.

## Next Experiments

1. Dense ILU scaling study before HACApK implementation:
   - vary problem size;
   - test both row-normalized `A_norm` and raw `A_raw`;
   - record GMRES iterations, residual, input nnz, `L+U` nnz, factor time,
     solve time, and memory;
   - check whether the loop near-zero cluster reappears as fill is reduced.

2. If dense ILU remains convincing:
   - implement a separate preconditioner H-matrix instance;
   - factor it with H-LU / robust H-factorization;
   - use GMRES with `A_H` matvec and `M_H^{-1}` apply;
   - sweep `A_H` ACA tolerance and `M_H` ACA tolerance independently.

3. Keep policy:
   - no silent dense-LU fallback;
   - no scalar/identity substitute preconditioner;
   - fail loud if factorization/preconditioner construction fails;
   - no hard-coded thread count.

4. For MDX:
   - develop and smoke on LAB first;
   - release through the normal package path before MDX benchmarking;
   - run compute on MDX temp storage;
   - recover JSON/JSONL results into repository docs/ipynb, not only `C:\temp`.

