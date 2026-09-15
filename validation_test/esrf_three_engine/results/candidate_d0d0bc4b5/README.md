# WEDGE nonlinear reproducibility investigation

Status: **HOLD**, not production acceptance. No tolerances or iteration budgets
were relaxed. Data were measured on 2026-09-15 in disposable installed-wheel
environments, not in editable solver installations.

Candidate source: `d0d0bc4b5b65b5c7e62943f1d603f7c814a8bd49`.
CI: 34938079994, successful. Radia 5.0.0 wheel SHA-256:
`6f7561812574acd550b831a91582b738451e00c48808fc9876529f51e52ca9ca`.
Old comparison wheel SHA-256:
`e9e6b9f09d105552330ac9b2b6d870187d825d11bf7a78f407553c8c7a866e25`.
The old wheel is the original locally built 4d85e72cc candidate, not the later
official 59b094d8 CI wheel. These identities must not be interchanged.

## Results

All single-case runs use two structured straight WEDGE elements, field order 2,
141 unknowns, 300000 A/m excitation, four NGSolve threads, gram_eps=1e-9,
linear tol=1e-9, nonlinear residual tol=2e-5 and nonlinear maxit=80. The generated
mesh recipe and BH table are recorded. Vertex/connectivity hashes alone are
not full .vol identities. There is no random mesh generation or explicit seed.

| Wheel / host / route | Result |
| --- | --- |
| Old / LAB / default auto-Jacobi, full six-case sequence | FAIL after five completed cases |
| Old / hibino / same six-case sequence | PASS, six cases |
| Old / hibino / WEDGE2 alone | PASS, same WEDGE2 statistics as sequence |
| Old / LAB / explicit mass-Riesz | PASS, residual 1.02027e-8 |
| Old / LAB / auto-Jacobi with MKL_CBWR=COMPATIBLE | PASS, residual 1.69512e-5 |
| v5 / LAB / default auto-Jacobi | FAIL, Newton 20, inner true residual 0.00626886 versus 0.0002745 |
| v5 / hibino / default auto-Jacobi | FAIL, Newton 17, inner true residual 0.0288164 versus 0.001 |

The original tracked `candidate_4d85e72cc/native_smoke_strong.json` is a separate
completed six-case experiment (SHA-256 a65df8e1d4b5eb7d32d90ffc574cc94ca2dfdf676f694d512a1e678ec654e2ff).
The new LAB old-wheel trial (776fe73b3ce5058cd0e43f79bb68bc27cc92b7a40773c574a206a9d3033802bb)
is incomplete. The prior PASS was not inferred from that incomplete file.
Repeating the original sequence on hibino reproduced WEDGE2's 19 Newton steps,
14318 inner iterations, 12 backtracks and residual 1.5806713897439883e-5 exactly.
Sequence: TET1, TET2, HEX1, HEX2, WEDGE1, WEDGE2. Each case runs a linear solve,
linear-BH parity solve, then the strong nonlinear solve; no explicit preconditioner
override. The single-case probe omits the two preceding parity solves.

## What this establishes

There is host/build/arithmetic sensitivity in a difficult Jacobi path, not merely
an old checkpoint mix-up. The MKL_CBWR child-process experiment is a diagnostic,
not a default configuration change or a substitute acceptance pass. Actual
dispatch instructions are not established solely by that environment variable.

The failing old LAB frozen Newton matrix was symmetric to 5.96e-16 relative;
its symmetric part had positive computed eigenvalues. Jacobi scaling left a
computed condition number about 5.46e7. These are observations at that one frozen
stage, not a certificate for every tangent or proof of root cause.

Both v5 hosts loaded MKL 2026.1-Product from their isolated venv, with DLL SHA-256
0ade2b44b4786fa3e4ed19c1e81fa2e2ac975093e373674db52f5db57d6e54b2.
Runtime JSON records CPU, DLL paths/hashes and threadpool information. Threadpool
counts are post-solve observations, not proof of the threads used inside a kernel.
LAB is i7-9700K; hibino is dual Xeon Platinum 8368. NGSolve thread count is four;
the BLAS inventories differ and must not be represented as globally four-threaded.

## Solver review and next experiment

The scalar C++ PCG verifies true residual before accepting convergence and
refreshes/restarts every 1000 iterations. Independent W+N application confirms
the failed returned iterate's residual, so this is not just a recursive-residual
reporting error. Python rejects `iters >= maxit`; that can reject an iterate
converged exactly at the cap, but these measured failures are genuinely above
their requested tolerance and are not explained by that boundary issue.

Next, compare more accurate accumulation and reliable residual replacement on
the *same frozen system*, while retaining the true-residual contract and budget.
Separately evaluate an explicit small-system mass-Riesz auto-selection policy;
it must not silently fall back after Jacobi fails, and needs both-host regression
and current-wheel tests before adoption. Large HEX/WEDGE performance evidence
must remain separate. No production solver/default was changed in this report.

Recovery owner: HDiv_MMM. The remote reproduction root is
`C:/temp/hdiv-wedge-repro-20260915`. Recover scripts, both wheels, install report,
results and logs with hash verification before deleting the disposable venv/root.
Other historical hibino directories are outside this investigation's cleanup.

Recovered archive:
`S:/Radia/validation_artifacts/hdiv_cleanup_20260915/hdiv-wedge-repro-20260915-recovery.tar.gz`.
SHA-256: `35c8dbb671895361d95613c6127054249c65995aba099bfbfd63c24d5cd29e81`.
All 15 archived files match the remote per-file hashes and sizes in
`wedge-recovery-manifest.json`. The disposable venv is reproducible from the
included exact requirements, install report and wheels; it is not archived.
