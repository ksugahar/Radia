# Production ICCG and AMS validation, 2026-10-04

Native functional source: `1497e31b215fd3a68ebd8381c04cf0ea820d03d7`.
Subsequent documentation and evidence changes do not change solver arithmetic.
The default remains **complex symmetric**, `A^T = A`, `conjugate=False`.
Hermitian positive-definite input, `A^H = A`, requires explicit
`conjugate=True` (MATLAB `Conjugate=true`). COCR rejects that option.
Standalone `ICPreconditioner` / MATLAB `IC` still uses transpose IC;
use the ICCG solver option for Hermitian systems.

## Implementation

- Serial IC triangular application traverses natural rows for at least 49,152
  rows. Parallel level scheduling remains in place. Conjugation dispatch occurs
  outside factorization loops. This threshold is workload-specific.
- Hermitian IC conjugates elimination products and backward factors. Structural
  checks reject non-Hermitian matrices, including before zero-RHS returns.
  Positive pivots and encountered Krylov curvature are checked. These checks
  are not a global positive-definiteness certificate: callers must supply HPD
  matrices. Complex symmetry alone does not guarantee convergence either.
- Reused IC setup checks CSR structure and ordering configuration. The composed
  original-to-RCM-to-ABMC permutation fixes reordered SpMV, and failed setup
  invalidates prior factors.
- AMS coarse corrections reuse the existing NGSolve Cholesky factors with a
  fixed sequential block application order, including direct gradient/Pi
  subspace solves. This removes variation from concurrent accumulation without
  changing the coarse factorization or relaxing acceptance thresholds.
  See the block solve implementation in
  [NGSolve 6.2.2607](https://github.com/NGSolve/ngsolve/blob/v6.2.2607/linalg/sparsecholesky.cpp).

## Throughput

Windows x64, Python 3.12.10, NGSolve 6.2.2607; unstructured first-order HCurl,
58,298 DOFs, `nograds=True`, all-boundary Dirichlet, `K + (1+100i) M`.
The installed Radia 5.2.2 baseline and candidate used identical mesh, matrix
values and manufactured RHS hashes. The baseline JSON's checkout commit records
the harness checkout, **not** the installed baseline binary's source revision.
Each process ran one warmup and three measured solves; the table gives medians.
Timing includes solver setup/factorization and iteration inside `Solve`;
matrix assembly and the independent residual check are outside it.
NGSolve TaskManager used the listed threads; BLAS thread counts were one.
Tolerance was `1e-9`, maximum iterations 2,000, with divergence checking disabled
and default scaling/automatic shift retained (selected shift 1.0).

| Threads | Baseline (s) | Candidate (s) | Baseline / candidate |
|---:|---:|---:|---:|
| 1 | 1.90843 | 1.12748 | 1.693 |
| 4 | 0.89451 | 0.93247 | 0.959 |

The serial case improves; **the four-thread case is about 4.2% slower**.
No parallel speedup or MEX speedup is established. Both binaries take 106
iterations and produce original-system relative residuals near `9.42e-10`.
Candidate manufactured-solution relative error is `2.93e-8`. Process RSS is
approximately 172 MiB for both, not an isolated factorization-memory measure.
Residual histories and solutions are not bit-identical between binaries;
iteration count, original-system residual and the independent known solution
are the acceptance evidence. A single mesh is not a universal performance claim.

## Correctness

- Relevant Python solver/caller suite: **313 passed**, 0 failed (68.39 s).
- Full native build and MEX build with required provenance succeeded at the
  source revision above. MATLAB parity: **63 passed**, 0 failed/incomplete.
- The original AMS study's **48 configurations / 192 solves** pass. Maximum
  original-system relative residual is `9.62297e-11`; repeated application
  differences are zero. Separate regression tests require exact repeated
  application equality on the same factor at one and four threads. This does
  not promise bitwise identity across compilers or separate factorizations.
- Hermitian tests include genuinely non-real off-diagonal HPD matrices with a
  NumPy direct oracle, scaling on/off, natural/ABMC/RCM ordering, extreme
  magnitudes and invalid inputs. The MEX Hermitian test uses real coefficients
  in a complex space and a complex RHS, compared with MATLAB sparse backslash;
  it does not establish arbitrary complex tensor assembly support.
- The AMS manufactured loads and curl error checks retain the limitations of
  the original study: compatible pure-curl loads are covered, arbitrary
  incompatible singular loads and indefinite Maxwell systems are not.

Machine-readable results are in `results/production-*.json`. Operational paths
are replaced with portable descriptions; raw artifact hashes are retained.
Historical `results/ams.json` records the earlier failure and is not overwritten.
This is local validation, not release-quad or deployment certification.

## Reproduction

Use an idle compute host and the same saved mesh for both processes. Run from
the repository root. First select the installed baseline by clearing PYTHONPATH;
then select the built candidate. Neither command installs a package.

```powershell
$env:PYTHONPATH = ''
python validation_test/sparsesolv/iccg_study/run_production.py --part throughput --kind symmetric --space hcurl --mesh C:/temp/iccg-production.vol --output C:/temp/iccg-baseline.json
$env:PYTHONPATH = (Resolve-Path src).Path
python validation_test/sparsesolv/iccg_study/run_production.py --part throughput --kind symmetric --space hcurl --mesh C:/temp/iccg-production.vol --output C:/temp/iccg-candidate.json
python validation_test/sparsesolv/iccg_study/run_production.py --part ams --output C:/temp/ams-production.json
```

Check the recorded module identity and binary hash before comparing results.
The retained mesh hash identifies this measurement; regenerating an unstructured
mesh may change its connectivity and timings.
