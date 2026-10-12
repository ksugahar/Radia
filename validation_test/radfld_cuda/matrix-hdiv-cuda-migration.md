# Matrix Build and HDiv-MMM CUDA migration boundary

This is an implementation inventory, not GPU acceptance or a performance result.

## Existing numerical owners

- `src/core/rad_hacapk_hdiv_entry.cpp`: `GetEntryStrategy` selects sampled Laplace,
  planar, hex, wedge, high-order tet, analytic and monopole entry evaluation.
  `EntryStrategy::EvaluateHostBlock` retains symmetric integrated blocks and image
  contributions. `HighOrderTetEntryStrategy::Evaluate` preserves near/far selection,
  curved-touch handling and polynomial combinations. A sampled Laplace CUDA
  evaluator cannot replace these physical element routes.
- `src/core/rad_hacapk.h`: `RadHACApKBase::BuildHMatrix` owns preparation, clustering,
  leaf creation, ACA+ fill and diagonal cache. Existing `RadHACApKStats` separates
  `t_prep`, `t_cluster`, `t_leafgen`, `t_fill` and `t_diag`; use those before choosing
  matrix-build acceleration. `GetLeafmtxp` is an internal opaque native boundary,
  not a MATLAB raw-pointer API.
- `src/ext/HACApK/cHACApK_cpp_impl.c`: `matvec_sym_thread_func` defines the symmetric
  application. Strictly lower leaves are skipped. An upper leaf is applied forward
  and transposed; a diagonal leaf is applied exactly once. Low-rank storage uses
  `a1[ndt,rank]` and `a2[ndl,rank]` in column-major external HACApK layout, so the
  forward block is `a2*a1^T`. Dense storage is row-major `[ndl][ndt]`.
- `src/core/rad_hacapk_hdiv.cpp`: `RadHACApKChargeGram::MatVecSym` applies physical
  scaling `G*x = S*(Ghat*(S*x))`. Exact dense normalized Gram has its own branch.
  `MatVecSymManyConfigured` also owns component/mask selection. These contracts
  must remain attached to a resident GPU operator.

## Narrow next implementation

1. Pack accepted HACApK upper/diagonal leaves and permutation into a checked
   mexcuda-owned context, preserving external storage layout explicitly. Construct
   it from native assembly inside MATLAB; Python snapshots are not a user dependency.
2. Apply low-rank factors and dense leaves on CUDA, including transpose mirrors,
   original-order permutation and both physical normalization factors. Keep vector
   state resident between repeated applications. Never launch NGSolve TaskManager.
3. Compare against the existing public CPU ChargeGram application for independent
   dense, low-rank, diagonal, permutation, normalization, image and element cases.
   A synthetic packed-block check is a lower-level prerequisite only.
4. Connect `B^T*G*B`, material terms, preconditioner and nonlinear solver without
   changing NGSolve FE semantics. Keep `sparsecholesky` FE subsolves and true residual
   gates. Test actual HDiv-MMM physics before release acceptance.
5. Profile entry fill separately. Batch the existing admissible physical entry
   evaluation for CUDA while keeping exact singular/near integrated blocks and ACA+
   compression. Do not replace an entry strategy with centroids or point clouds.

Release status: no Matrix Build or HDiv-MMM GPU physical acceptance is established
by this inventory. RadFld CUDA projection has its own bounded numerical evidence.
A100 is excluded from all current build, connection and execution operations.
