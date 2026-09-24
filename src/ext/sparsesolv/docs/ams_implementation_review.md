# AMS implementation review

This review covers the real lowest-order HCurl path. It is not a proof of
optimal complexity or of equivalence for every material and mesh.

## Changes verified

- Galerkin products and AMG hierarchy construction own bounded TaskManager
  regions. Caller-owned setup regions remain rejected. Coarse direct solves
  are not included in the internally owned construction region.
- `setup_workers` measures participation in strength-row construction, not
  the configured thread count or utilization throughout setup.
- A fused `x=b/D` first-smoother candidate was tested but not adopted: no
  reproducible performance gain was established. The existing path is retained.

## Remaining review targets

1. PMIS coarsening is a sequential greedy update despite its name: it writes
   `cf[i]` immediately and later vertices read the modified markers. Merely
   replacing this loop with ParallelFor introduces races and changes the
   hierarchy. A parallel replacement needs snapshot/proposal/commit phases,
   deterministic tie-breaking and convergence/complexity validation.
2. SetupGeometry and some BC/norm construction loops still run outside an
   active TaskManager. A ParallelFor spelling alone does not prove parallel
   execution. Profile these phases before expanding their owned regions.
3. Coarse factorization and the historical Windows fastfail require more
   isolation. Successful bounded regions do not justify deleting the caller
   guard or running several factorizations concurrently.
4. Mult reuses mutable work vectors and profiling accumulators. A single
   preconditioner instance is not reentrant; parallel row operations within
   one application do not authorize concurrent applications to that instance.
5. Several full-vector passes and clock reads remain per application. Measure
   memory traffic and coarse correction costs before fusing more operations.
6. Large pure curl-curl cases need an explicit repeatability test on a fixed
   matrix and fixed RHS, independently of final field convergence. Compare
   repeated Mult calls as well as hierarchy rebuilds, and distinguish numerical
   nullspace sensitivity from data races before changing parallelism.

Retain separate timings for hierarchy construction, application, and the whole
solve. A faster preconditioner that increases Krylov iterations is not
automatically a faster solver. Verify true residuals and field observables,
including updates, multiple smoother sweeps, and serial/parallel execution.
