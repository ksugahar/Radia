# Radia Lightweight Test Suite

`tests/` is the fast developer and CI gate for the Radia Python 3.12
module.  It should stay small enough to run during ordinary debugging and
pre-push checks.

Heavy solver studies, Cubit export checks, GUI/panel goldens, benchmarks,
and cross-validation cases live under `validation_test/`.

## Two-Stage Test Layout

| Directory | Purpose | Typical command |
| --- | --- | --- |
| `tests/` | Lightweight debug and CI tests | `python -m pytest tests/` |
| `validation_test/` | Heavy validation, benchmarks, GUI, Cubit, golden checks | `python -m pytest validation_test/` |

The default `pytest` configuration discovers only `tests/`.  Run
`validation_test/` explicitly when doing release validation or research-grade
checks.

## Running Lightweight Tests

```powershell
python -m pytest tests/
python -m pytest tests/ -m "not slow and not golden"
python tools/ci_preflight.py --only toplevel-collect
```

Useful single-file checks:

```powershell
python tests/test_simple.py
python tests/test_radia.py
python tests/test_advanced.py
python tests/test_radia_ngsolve.py
```

`test_radia_ngsolve.py` requires NGSolve.  If NGSolve is unavailable, the
test is skipped.

## Benchmarks

Benchmark scripts are validation artifacts and live in
`validation_test/benchmarks/`:

```powershell
python validation_test/benchmarks/benchmark_parallel.py
python validation_test/benchmarks/benchmark_field_parallel.py
python validation_test/benchmarks/benchmark_correct.py
python validation_test/benchmarks/benchmark_heavy.py
python validation_test/benchmarks/benchmark_threads.py
```

## Test Categories

- Basic smoke tests: import, version, geometry creation, field evaluation.
- Functional tests: materials, transformations, relaxation, memory handling.
- Lightweight integration tests: optional-package tests that skip cleanly when
  the dependency is missing.
- Performance smoke tests: TaskManager sanity checks that are still short
  enough for local debugging.

Long p-convergence, Cubit 2025.12 export, panel GUI, benchmark, and solver
cross-validation cases belong in `validation_test/`.

## Continuous Integration

CI runs the lightweight gate by default:

```yaml
- name: Run lightweight tests
  run: python -m pytest tests/ -m "not slow and not golden"
```

Manual release validation can add:

```yaml
- name: Run validation tests
  run: python -m pytest validation_test/ -m "not slow"
```

Use `full_slow=true` in the GitHub workflow dispatch when the slow validation
set is intentionally required.

## Test Data

Small fixtures for lightweight tests belong in `tests/fixtures/`.  Validation
fixtures belong next to their validation tests under `validation_test/`.

## Writing New Tests

- `test_*.py`: pytest-discovered tests.
- `benchmark_*.py`: manual benchmark scripts under `validation_test/benchmarks/`.
- Mark long-running tests with `@pytest.mark.slow`.
- Mark reference/golden checks with `@pytest.mark.golden`.

Keep the first version of a new regression test in `tests/` only if it is fast,
deterministic, and useful during ordinary debugging.  Promote heavier checks to
`validation_test/`.

## Troubleshooting

If `import radia` fails, rebuild the module first:

```powershell
powershell.exe -ExecutionPolicy Bypass -File Build.ps1
python -m pytest tests/test_simple.py
```

If a validation test fails after this split, check whether it still references
an old `tests/...` fixture path.  Runtime fixture paths should point at the new
`validation_test/...` location.

## References

- pytest documentation: https://docs.pytest.org/
- Radia documentation: `docs/`
