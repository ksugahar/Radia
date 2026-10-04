# ICCG Hermitian and throughput study

The original study is an isolated, serial C++ experiment. It
builds four temporary pybind modules from the headers pinned at `a0700dac3`:
baseline, natural-row triangular solves, Hermitian IC, and both changes.
It does not install modules. Production now supports explicit Hermitian ICCG;
the default remains complex symmetric. COCR's contract is unchanged.
`run_production.py` exercises the production Python binary, including
unstructured complex-symmetric systems and the original 48 AMS conditions.

## Reproduce

Use an idle compute host with Python, NumPy, SciPy, pybind11, CMake and MSVC.
The pinned baseline commit must be available locally; shallow clones need its history.
The AMS comparison additionally uses the installed Radia/NGSolve native pair.
Run in a Developer PowerShell with CMake on PATH, from the repository root:

```powershell
python validation_test/sparsesolv/iccg_study/build_variants.py --output C:/temp/iccg-variants
python validation_test/sparsesolv/iccg_study/run_study.py --build C:/temp/iccg-variants --output C:/temp/iccg-core.json --source-sha (git rev-parse HEAD) --part core
python validation_test/sparsesolv/iccg_study/run_study.py --build C:/temp/iccg-variants --output C:/temp/iccg-ams.json --source-sha (git rev-parse HEAD) --part ams
```

Use a fresh output directory. Source transformations assert their expected
matches and save exact generated patches and header hashes. Recover results,
logs and source/binary hashes before removing job-owned compute scratch.
Prototype builds use serial SparseSolv primitives, MSVC Release flags, and no
NGSolve TaskManager or OpenMP. The AMS comparison separately selects 1 or 4
NGSolve threads. These timings cannot establish parallel prototype performance.

## Correctness and measurement

- Hermitian input is **positive definite**, not merely Hermitian. The dense
  oracle is `B B^H + 2 I` with genuinely non-real off-diagonal entries; NumPy's
  direct solution is independent of the native iterative implementation.
- The IC trial conjugates the elimination products and the stored transpose
  used by backward substitution. Changing the CG dot product alone is a
  negative control. A full-pattern unshifted IC must solve the dense case in
  one iteration, including ABMC reordering and scaling on/off.
- Complex-symmetric solves with `conjugate=False` must retain bit-identical
  solutions and residual histories. Sparse Hermitian cases are unitary
  diagonal transformations of Dirichlet Laplacians, so positive definiteness
  is known by construction. RHS magnitudes include `1e-100` and `1e100`.
- Natural traversal is restricted to the serial path at 49,152 or more rows.
  It preserves each row's arithmetic and dependency ordering. This threshold
  is an experimental candidate, not a universal tuning constant.
- Core timings include scaling, setup, iteration and final true-residual
  computation. Each pair is warmed, then alternated for three measured runs.
  An additional profiled solve measures preconditioner-plus-dot time; isolated
  SpMV timing is not claimed to be the fused Krylov SpMV cost. RSS is sampled
  for the whole process, not attributed to one factorization.
- The AMS comparison uses the **same NativePCG and stopping tolerance** for
  both preconditioners. IC here is the existing standalone fixed-shift 1.05
  preconditioner, not the factory ICCG's default auto-shift/scaling contract.
  Setup is outside the caller-owned solve TaskManager as required by AMS.
- Every curl-curl RHS is manufactured as `A*x` on free DOFs. For zero mass
  this guarantees compatibility; solutions are compared by curl, not by
  coefficient vectors modulo the gradient nullspace. With small positive
  mass this RHS weakly excites gradient modes, so the results do not cover
  arbitrary loads, incompatible singular RHSs or indefinite Maxwell systems.
- All warm and measured AMS runs retain their residual, curl error and repeated
  application checks. A repeatability failure exits nonzero even when the
  physical solve passes. There is no tolerance relaxation or automatic retry.

## Results, 2026-10-04

Header baseline: `6cca23ef604f9c488e150bedcf131157f37bbfae`.
Measured on mdx2, MSVC 19.44, Python 3.12.10, NGSolve 6.2.2607, Radia 5.2.2.
The 16 core checks pass. All 24 benchmark rows pass the original-system
relative residual limit `1e-8`; paired solutions, histories and shifts are
bit-identical. An independent second timing pass confirms the large-case
improvement. The table uses its median total times:

| Matrix | Rows | Real baseline → natural (s) | Speed ratio | Hermitian trial → combined (s) | Speed ratio |
|---|---:|---:|---:|---:|---:|
| 2D Laplacian | 65,536 | 0.820 → 0.530 | 1.55 | 2.396 → 0.987 | 2.43 |
| 2D Laplacian | 147,456 | 2.857 → 1.577 | 1.81 | 8.010 → 3.271 | 2.45 |
| 3D Laplacian | 64,000 | 0.200 → 0.131 | 1.53 | 0.551 → 0.293 | 1.88 |
| 3D Laplacian | 125,000 | 0.543 → 0.332 | 1.64 | 1.362 → 0.677 | 2.01 |

The first profiled run attributes about 70–81% of iteration time to IC
application plus its dot product on these large baselines. This supports
prioritizing triangular access locality. It does not establish a speedup on
unstructured application matrices, parallel runs, Python/MEX wrappers, or a
released Radia build. Small-case timing changes are not attributed to the
natural-order branch, which is inactive below the threshold.

The AMS comparison covers 48 combinations (three meshes up to 12,217 DOFs,
four coefficient/mass settings, two thread counts and two preconditioners),
with 192 solves including warmups. All solve residuals are at most
`9.63e-11` and relative curl errors at most `1.47e-8`. At 12,217 DOFs and
coefficient contrast 1,000 with unit mass, IC takes 2,134 iterations versus
AMS's 28; serial total times are approximately 2.08 s and 0.16 s. On the
compatible pure curl-curl load, IC takes about 0.06 s versus AMS's 0.11 s.
This supports retaining both methods, not an unconditional AMS default.

**The AMS study's strict repeatability gate fails: 42/48 combinations pass.**
The six failures are four-thread AMS with coefficient contrast 1,000 and
positive mass (1 or `1e-6`), across all three meshes. Across all repetitions,
the largest relative difference between identical-input preconditioner
applications is `2.03e-4`. The difference's relative curl and matrix action
are at most `2.76e-14` and `2.12e-16`. This is consistent with amplification
in gradient directions, but does not prove the absence of a concurrency
defect. Serial applications and the tested zero-mass route pass the same
gate. Retain this failed diagnostic rather than claiming AMS is uniformly
repeatable or weakening the check.

Machine-readable evidence is in [results/core.json](results/core.json) and
[results/ams.json](results/ams.json); the latter intentionally has `pass:false`.
Operational paths/logs and prototype binaries are retained only privately.

## Original study decision and promotion requirements

These were the original experiment's open requirements. See
[production validation](PRODUCTION.md) for the subsequent implementation and
measured results; the original failed AMS evidence above remains historical.

1. Hermitian ICCG is viable for known HPD matrices. Keep this trial separate
   until structural Hermitian validation, positive-curvature/pivot failure
   behavior, reused setup and parallel ABMC are tested. The prototype assumes
   HPD inputs and is not safe as an arbitrary-matrix public entry point.
2. Natural serial row traversal is a promising measured optimization. Before
   production adoption, test unstructured application matrices around the
   threshold, native TaskManager dispatch and the Python/MEX paths. Do not
   promote the structured serial speed ratio to an application claim.
3. Retain AMS for appropriate H(curl) problems; investigate the measured
   repeated-application sensitivity before considering a broader default.
4. No public API, numerical gate, release tag or installed runtime changed.

Method background: [hypre AMS](https://hypre.readthedocs.io/en/latest/solvers-ams.html)
and [Eigen incomplete Cholesky](https://eigen.tuxfamily.org/dox/IncompleteCholesky_8h_source.html).
