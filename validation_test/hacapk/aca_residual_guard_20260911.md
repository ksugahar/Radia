# ACA Residual Guard Candidate

## Status

An uncommitted candidate changes `cHACApK_aca` in the isolated worktree. It is
NOT approved for production or release: residual verification increases kernel
entry evaluations substantially. No installed Radia binary, Python binding,
MATLAB MEX, MCP process, or shared checkout was replaced.

## Implemented Candidate

- Reinspect unused residual rows before accepting a zero/small selected row.
- Treat a small rank-one update as a request for residual verification, not
  proof of convergence.
- Remove the dimensioned `1e-20` pivot cutoff. Normalize rows by division,
  without forming a potentially overflowing reciprocal first.
- Compute the Frobenius norm of the current low-rank approximation including
  cross terms. The residual-row threshold is relative to that norm.

If every unused residual row has norm at most
`eps * norm(approximation) / sqrt(number_of_rows)`, the unmasked residual
Frobenius norm is bounded accordingly. Pivoted rows/columns are interpolated
in exact arithmetic; their floating-point residuals remain a verification
concern. Rank-cap termination is not a convergence certificate.

## Native Verification

Two standalone diagnostic DLLs were compiled with MSVC 18.5.1 from the actual
repository C sources, before and after the edit. `aca_probe_bridge.c` supplies
a serial matrix-entry callback. `aca_probe_unreachable.c` contains abort-only
link guards for unused BLAS/tree entry points; it is never production code.
`validate_aca_residual_guard.py` calls the DLLs through ctypes and verifies
the reconstructed matrices against their original inputs.

The adjacent JSON records 29 cases, source/DLL hashes, ranks, errors, callback
counts and single-call diagnostic timings. All candidate cases meet the
explicit `1e-9` relative-error budget with requested ACA tolerance `1e-10`.

| Case | Before error | Candidate error |
| --- | ---: | ---: |
| Hidden residual, unit scale | 0.288675 | 0 |
| Hidden residual, scale 1e-80 | 1 | 0 |
| Initial zero row, nonzero matrix | 1 | 0 |
| Smooth Laplace, epsilon 0 | 1.305e-10 | 1.664e-11 |
| Smooth Laplace, epsilon 0.03 | 1.890e-10 | 1.598e-11 |
| Smooth Laplace, epsilon 0.3 | 1.546e-10 | 1.655e-11 |

Coverage includes all six row permutations of the hidden-residual matrix at
scales 1e-80, 1 and 1e80, a zero matrix, an initially zero pivot row, signed
rank-two rectangular matrices, full-rank rectangular matrices, and smooth
216-by-216 Laplace leaves. This is a leaf-compression comparison, not the
original spherical fixture's complete H-matrix tree. The latter's expected
failure remains in place because its installed native binary is unchanged.

## Performance Blocker

Using the current approximation norm improved the first candidate, but did not
remove the cost of deterministic residual scanning:

| Leaf | Before entry calls | Candidate entry calls | Ratio |
| --- | ---: | ---: | ---: |
| Signed rank-two, 64 by 48 | 330 | 3120 | 9.45 |
| Smooth Laplace, epsilon 0 | 24222 | 113067 | 4.67 |
| Smooth Laplace, epsilon 0.03 | 24820 | 191058 | 7.70 |
| Smooth Laplace, epsilon 0.3 | 28782 | 153832 | 5.34 |

These are deterministic callback counts, not end-to-end solver timing claims.
Worst-case stop-time scans are quadratic in leaf dimensions and can repeat
after finding a remaining residual row. Thus this candidate is an accuracy
reference, not an acceptable unqualified replacement for scalable ACA.

Next, reduce repeated scans and compare robust pivot reselection/ACA+ against
this reference. Any probabilistic or bounded residual check must declare its
weaker guarantee. Before promotion, rebuild the intended native target and
rerun the spherical H-matrix fixture, HDiv/HCurl/BEM/PEEC regressions, and
build-cost measurements. Do not remove the original strict expected failure
based only on these leaf-level results.

## Reproduction

Compile `aca_probe_bridge.c`, `aca_probe_unreachable.c`,
`src/ext/HACApK/cHACApK_base.c`, and `cHACApK_lib.c` into a DLL with
`cl /O2 /Gy /LD`, include paths `src/ext/HACApK` and `src/core`, and linker
option `/OPT:REF`. Build the baseline from source commit `b7e89ae53` and the
candidate from the edited source; use distinct output DLL names under
`C:/temp`. Existing shared-source printf/prototype warnings are not fixed here.

Then run:

```text
python validation_test/hacapk/validate_aca_residual_guard.py --before C:/temp/aca_probe_before.dll --after C:/temp/aca_probe_after.dll --output validation_test/hacapk/aca_residual_guard_20260911.json
```

The driver exits nonzero if any candidate case misses the accuracy budget.
Ruff and `git diff --check` passed for this change.
