# Radia Validation Test Suite

`validation_test/` contains the heavy tests that are valuable for release
confidence and research validation but too expensive or environment-specific
for the default CI/debug loop.

## Use validation cases as executable samples

This directory is also the sample collection for physics, numerical methods and
application workflows. Read each case's model, dependencies, inputs, assertions
and recorded results together. Select a relevant case rather than running the
whole suite as a tutorial; optional tools and compute-host requirements still
apply. A saved result supports only its recorded configuration and criteria.

[tests/](../tests/README.md) supplies smaller usage and contract examples.
[docs/](../docs/README.md) introduces features through saved demonstrations;
the [validation guide](../docs/VALIDATION_GUIDE.md) maps questions to evidence.
For operational work, use [MCP](../packages/radia-mcp/README.md) or human–AI
collaboration through [Cubit](../docs/cubit_mesh_export/README.md) and
[Simulink](../matlab/README.md).

For detailed formulation comparisons, convergence and paper-supporting work,
see the [numerical study index](STUDY_INDEX.md).

## Select an application sample

| Question | Executable case and evidence | Scope |
| --- | --- | --- |
| Does a constructed winding produce the intended field? | [Stream-function lane](stream_function/README.md) and [independent field driver](stream_function/verify_coil_field_independent.py) | Compares contour-wire field evaluations; read the driver's geometry-dependent limitations before interpreting agreement |
| Does a motor ROM reproduce response between training angles? | [Electric-machine lane](electric_machine/README.md), [angle-periodic driver](electric_machine/validation_pmsm_angle_periodic_rom.py), [recorded summary](electric_machine/validation_pmsm_angle_periodic_rom_summary.json) | Curved 8-pole/24-slot 2D model and interlaced holdout angles; solver-heavy compute-host execution |
| Does electromagnetic power reach the thermal operator consistently? | [IH operator physics test](induction_heating/test_ih_operator_physics_golden.py) and [lane instructions](induction_heating/README.md#simulink-operator-physics-golden) | Real unit-current electromagnetic solve, power closure and thermal matrix identities; does not by itself run a complete Simulink application |

These are executable sample entry points, not a report of a fresh successful run.
Use the lane instructions for dependencies and execution location. Inspect the
recorded inputs and criteria before reusing a result. The
[feature gallery](../docs/APPLICATION_GUIDE.md#browse-results) presents related
capabilities visually; each notebook and validation case retains its own geometry
and scope.

## What Belongs Here

- Cubit 2025.12 export and curved-mesh checks.
- Panel GUI, notebook, and golden-output checks.
- FEM/BEM/FEEC p-convergence and cross-validation studies.
- Solver-backed `radia-mcp` FEM/BEM and application validation under
  `validation_test/radia_mcp/`.
- Benchmarks and long solver regressions.
- Tests requiring special optional dependencies, licenses, or long wall time.

`tests/` is reserved for short checks. A test measured at 10 seconds or more in
two comparable successful CI runs moves here and is marked `slow`; a one-off
timing spike is measured again before reclassification. Real Office/Cubit GUI
startup, licensed applications, benchmarks, long golden/reference runs, solver
convergence studies, and publication validation belong here regardless of one
short measurement. The canonical directory name is `validation_test/`
(singular).

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
hostname is `mdx` or `hibino`. Run them over SSH on hibino first:

```powershell
ssh hibino python -m pytest validation_test/ -m compute_host -q
```

Use mdx only when hibino is unavailable and the mdx CI runner and job queue are
idle. A validation job must never delay CI or preflight.

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
This keeps routine CI responsive while preserving deeper hibino-first compute
checks before release or publication-quality claims.
