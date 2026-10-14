# SparseCholesky storage limits

Radia explicitly selects NGSolve SparseCholesky for supported direct FE systems.
It does not fall back to PARDISO. A successful true-residual check establishes
the accuracy of a returned solution; it cannot protect the process from a
native allocation or indexing failure during factorization.

## Windows counter risk in 6.2.2607

The [tagged upstream implementation](https://github.com/NGSolve/ngsolve/blob/v6.2.2607/linalg/sparsecholesky.cpp#L336-L390)
accumulates total factor entries and master-column entries in `long int` before
assigning them to wider storage indices. On Windows x64, `long` is 32-bit while
`size_t` is 64-bit. Therefore counts above `LONG_MAX` (2,147,483,647) cannot be
represented by those intermediate counters. This is a source-level risk,
distinct from establishing the cause of any particular process crash.

## Measured evidence and its limits

The [2606 mesh-family record](../../validation_test/esrf_three_engine/results/case3_sparsecholesky_size_limit_20260929.json)
contains successful factorizations at 324,992 and 596,748 free DOFs, with measured
factor entry counts of 596,702,979 and 1,659,756,518. Larger cases failed, but their
factor entry counts were not measured. Extrapolated counts are not observations.

A 2607 rerun did not establish a clean size threshold: available Windows commit
headroom was much smaller than free physical RAM, and the smallest previously
successful case also aborted. This memory-confounded trial is not proof of the
counter hypothesis, successful migration, or a universal DOF cutoff.

Fill depends on the graph and ordering, not just the number of unknowns. Free
physical RAM alone is also insufficient: the system must have enough commit
headroom for the factor and temporary allocations.

## Operational status

- Select a validated iterative route explicitly when direct factorization is
  unsupported at the required scale. For supported HCurl spaces this may be
  AMS or BDDC with an AMS coarse solver; that is not an H1 solver prescription.
- Do not automatically change the operator, solver or acceptance tolerance
  after an allocation failure.
- A production pre-factorization storage guard is **not yet implemented**.
  Adequate-headroom 2607 reproduction and checked symbolic-storage assessment
  remain release work. This document does not close those requirements.
- An upstream report is prepared privately and has not been submitted.
