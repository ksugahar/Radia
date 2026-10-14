# Mixed Omega trace-factor reuse

The solve-local cache retains the matching trace topology and symbolic SparseCholesky factorization. Numeric factors are updated when coefficients change. Original constraints and original-system residual checks remain active. No solver fallback or tolerance change is used.

On NGSolve 6.2.2607, eight threads, three P2 systems (9,261 / 42,594 / 69,243 DOFs), subsequent factorizations took 0.075–0.080 / 0.712–0.753 / 1.308–1.393 seconds, versus 1.303–1.309 / 5.814–6.876 / 9.036–9.276 seconds without reuse. First calls were slightly slower. These are single observations per changing coefficient, not statistical speedup estimates. All original-system residuals were below 1e-10.

This does not establish the ESRF case 6 end-to-end target. That validation remains required on its specified host and inputs.

To reproduce, install a native-compatible current source build, export `src/radia/kelvin_solver.py` from baseline commit `0e9f26bd6` to a temporary file, then run `python validation_test/mixed_omega_trace_reuse/run.py --baseline-source <file> --output <result.json>`. The supplied baseline file is executed and must be trusted. The runner reuses the existing two-region test fixture, records both source hashes, varies the primal coefficient on every call, and rejects residual failures.

The associated mixed Omega/Newton test run passed 64 tests with one thread. Separate P1/P2/P3 regression cases change both primal and interface coefficients to detect stale factors. Raw logs and interrupted diagnostic trials are retained privately; they are not successful benchmark samples.
