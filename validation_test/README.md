# Radia Validation Test Suite

`validation_test/` contains the heavy tests that are valuable for release
confidence and research validation but too expensive or environment-specific
for the default CI/debug loop.

## What Belongs Here

- Cubit 2025.12 export and curved-mesh checks.
- Panel GUI, notebook, and golden-output checks.
- FEM/BEM/FEEC p-convergence and cross-validation studies.
- Benchmarks and long solver regressions.
- Tests requiring special optional dependencies, licenses, or long wall time.

Legacy validation scripts that are meant to be run directly, not collected by
pytest, live in `validation_test/manual/` and use `*_validation.py` names.

## Running

```powershell
python -m pytest validation_test/ --collect-only
python -m pytest validation_test/ -m "not slow"
python -m pytest validation_test/
python tools/ci_preflight.py --validation
python tools/ci_preflight.py --validation --full
```

Tests marked `compute_host` are solver-heavy and are skipped unless the actual
hostname is `mdx` or `hibino`.  Run them over SSH on an idle compute host:

```powershell
ssh mdx python -m pytest validation_test/ -m compute_host -q
```

Manual scripts:

```powershell
python validation_test/manual/batch_evaluation_validation.py
python validation_test/manual/far_field_accuracy_validation.py
python validation_test/manual/curlA_equals_B_validation.py
```

Cubit tests are opt-in:

```powershell
$env:RADIA_RUN_CUBIT_TESTS = "1"
python -m pytest validation_test/cubit -q
```

## Relationship To CI

The normal CI path runs `tests/` only.  `validation_test/` is a manual release
or operator-triggered gate, and GitHub Actions always excludes `compute_host`.
This keeps routine CI responsive while preserving deeper mdx/hibino checks
before release or publication-quality claims.
