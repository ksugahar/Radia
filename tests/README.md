# Radia Lightweight Test Suite

`tests/` is the fast developer and CI gate for the Radia Python 3.12
module.  It should stay small enough to run during ordinary debugging and
pre-push checks.

Heavy solver studies, Cubit export checks, GUI/panel goldens, benchmarks,
and cross-validation cases live under `validation_test/`.

## Use tests as executable samples

These tests also show small, concrete API calls, expected behavior and failure
conditions. Read the fixtures and assertions together, then run the relevant
case with its declared dependencies. Mocked contract tests illustrate interface
behavior; they do not establish physical accuracy.

For larger physics and application samples, use
[validation_test/](../validation_test/README.md). The
[documentation](../docs/README.md) introduces features and saved results.
Actual CAE operation uses [MCP](../packages/radia-mcp/README.md), or human–AI
collaboration through [Cubit](../docs/cubit_mesh_export/README.md) and
[Simulink](../packages/radia-simulink/matlab/README.md).

## Select a small sample

| What to learn | Executable sample | What to inspect |
| --- | --- | --- |
| Build a finite-section coil and evaluate its field | [Coil axis field](test_coil_axis_closed_form.py) | Geometry, current sign, scale and comparison with a closed-form axis field; requires Radia, NGSolve and Netgen |
| Construct and query a periodic motor ROM | [Motor ROM](test_motor_rom.py) | Synthetic periodic tables, ports, derivatives and bundle checks; these inputs are not a solved machine model |
| Transfer a heat source with checked field artifacts | [Thermal transfer](test_ih_thermal_transfer.py) | Mesh/field pairing, metadata, conserved power and rejected inputs; requires NGSolve/Netgen |
| Understand the circuit-field Simulink adapter contract | [Adapter source contract](mcp_integration/test_circuit_field_simulink_contract.py) | Static source assertions only; this does not execute Simulink or establish numerical accuracy |

Run a selected case from the repository root with the configured Radia environment,
for example:

```powershell
python -m pytest tests/test_motor_rom.py::test_periodic_fourier_value_derivative_and_continuous_skew -q
```

Read module imports and shared `conftest.py` fixtures before running. Tests may
need a built native package even when the selected calculation looks small.
For full-model comparisons, continue to the
[application validation samples](../validation_test/README.md#select-an-application-sample).

## Two-Stage Test Layout

| Directory | Purpose | Typical command |
| --- | --- | --- |
| `tests/` | Lightweight debug and CI tests | `python tools/run_test_tier.py --profile fast-contracts` |
| `validation_test/` | Heavy validation, benchmarks, GUI, Cubit, golden checks | `python -m pytest validation_test/` |

The default `pytest` configuration discovers only `tests/`.  Run
`validation_test/` explicitly when doing release validation or research-grade
checks.

## HDiv-MMM mesh families and CPU/GPU test tiers

HDiv-MMM requires three separate mesh families:

| Required family | Lightweight coverage in `tests/` | Larger coverage in `validation_test/` |
| --- | --- | --- |
| Tetrahedra only (`tet only`) | Minimal HDiv assembly and solve, orientation, true residual and a small independent reference | Physical convergence, nonlinear BH, reuse and scaling |
| Hexahedra only (`hex only`) | The same small numerical checks for hexahedra | Physical convergence, nonlinear BH, reuse and scaling |
| Tetrahedra and hexahedra joined through pyramids (`tet-hex mixed thru pyramid`) | A connected transition patch containing all three element types; shared-face orientation and normal-flux continuity, assembly, solve and true residual | A connected mixed-mesh physical model, mesh/order convergence, nonlinear BH and scaling |

A disconnected collection of elements or a single pyramid does not satisfy the
mixed transition case. NGSolve owns the HDiv basis, Piola maps, orientation and
quadrature. An unsupported upstream operation must have a reproducible explicit
failure; record that family's status as unsupported or unaccepted, never passed.
Pyramid-specific defects are reported upstream with the exact version and a
minimal reproducer. Waiting for their fix must not block the other validated mesh
families or the AMS, CompactAMG and CSR GPU development work.

Apply the same division to GPU testing. Small, deterministic GPU numerical,
backend-selection, error and lifecycle checks belong in `tests/` (or the owning
package's lightweight tests). Use an explicit GPU lane with declared prerequisites;
CPU-only CI must not silently substitute a host backend. Large GPU physics,
scaling, transfer-inclusive timing, peak-memory, repeated-run plateau and fault
sweeps belong in `validation_test/`, with machine-readable JSON evidence.
GPU HDiv-MMM coverage must track the same three mesh families, with independent
CPU agreement and original-equation true residuals. A source contract or mocked
backend test does not replace a numerical run. A missing required fixture fails;
only explicitly optional hardware checks may report a visible skip.

These are required coverage categories, not a claim that all three families or
both execution routes are currently implemented or accepted. Large Cubit-generated
mixed meshes and their export/checking workflows belong in `validation_test/`;
lightweight tests consume small checked fixtures.

## Running Lightweight Tests

```powershell
python tools/run_test_tier.py --profile fast-contracts
python tools/run_test_tier.py --profile native-smoke
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
  run: python tools/run_test_tier.py --profile fast-contracts
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
pwsh -NoProfile -ExecutionPolicy Bypass -File .\Build.ps1
python -m pytest tests/test_simple.py
```

If a validation test fails after this split, check whether it still references
an old `tests/...` fixture path.  Runtime fixture paths should point at the new
`validation_test/...` location.

## References

- pytest documentation: https://docs.pytest.org/
- Radia documentation: `docs/`
