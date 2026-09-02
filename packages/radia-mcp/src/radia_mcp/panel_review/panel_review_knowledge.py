"""Radia Simulink application-block review and construction knowledge."""

APPLICATION_BLOCK_REVIEW = r"""
# Radia application interface review

Production human interfaces are masked blocks in
`matlab/radia_simulink_library.slx`. Notebook workbenches are removed for every
application, including IH. Desktop PySide analysis windows remain retired.

Review these files together:

```text
src/radia/<app>_design.py
src/radia/panels/calc_<app>.py
src/radia/simulink/application.py
matlab/+radia/+simulink/runApplication.m
matlab/+radia/+simulink/applicationSFunction.m
matlab/+radia/+simulink/buildLibrary.m
tests/test_simulink_application.py
tests/matlab/test_simulink_workflow.m
```

Checks:

1. `DesignSpec` and argparse describe one settings/solver contract.
2. The block mask has an application id, config, run root, timeout, and backend
   executable without callback-only hidden state.
3. Trigger, status, primary, and elapsed ports have fixed types/dimensions and
   explicit sample-time semantics.
4. Python is launched only on an explicit rising edge, never every time step.
5. The block delegates to the tested API/CLI and contains no solver duplicate.
6. Every solver-bound `.vol` passes `check-vol` after export and before solver
   initialization. Production modes use a strict versioned label contract and
   retain the `cubit-mesh-export.vol-check.v1` report. Material constants stay
   in DesignSpec/configuration; no checker guesses them from region names.
7. `command.txt`, `run.log`, `solver_result.json`, and versioned `result.json`
   survive success, timeout, invalid config, dependency failure, and solver
   failure. A spatial solve also leaves a checked GMSH `.msh v4.1` artifact in
   the run directory and indexes it from `result.json`.
8. The block sample matches the same headless golden band.
9. MEX/ROM is not promoted merely because it compiles: require parity, error
   propagation, handle lifecycle, and long-run stability.
10. No `*_notebook.py` adapter or packaged workbench is reintroduced.
11. Cubit's embedded PySide6 remains isolated and protected.

Run:

```powershell
python -m pytest tests/test_simulink_application.py -q
python -m pytest tests/test_application_interface_manifest.py -q
matlab -batch "addpath('matlab'); r=runtests('tests/matlab/test_simulink_workflow.m'); assert(all([r.Passed]))"
```
"""


BUILD_APPLICATION_BLOCK = r"""
# Build a Radia Simulink application block

The production result is a masked Simulink block in the single Radia library.

1. Enumerate every application variable and solver switch.
2. Promote reusable code to `src/`; keep heavy truth runs in
   `validation_test/`.
3. Implement `calc_<app>.py` and `<App>DesignSpec`; golden-lock every supported
   backend.
4. Add a result-bearing `docs/<topic>/*.ipynb` only for explanation and saved
   results. Published CAE examples include saved WebGUI scenes. Field scenes
   use `Draw(field, mesh, name=..., ...)` with explicit display arguments; do
   not add control widgets or `CommandWorkbench`.
5. Add the application to `radia.simulink.buildLibrary` with a masked
   `radia_application_sfun` block.
6. Use `radia.simulink.application` for the initial explicit-trigger Python
   backend and standard artifacts. Spatial modes expose `--msh-output`; the
   runner redirects it to a checked GMSH `.msh v4.1` run artifact.
7. Add Python and MATLAB tests for mask wiring, port types, model update,
   execution, timeout, missing dependency, and failure provenance.
8. Regenerate `matlab/radia_simulink_library.slx` from the builder.

The standard configuration schema is
`radia.simulink.application_config.v1`. The standard run schema is
`radia.simulink.application_run.v1`. A rising trigger executes once; outputs
are `int32 status`, selected `double primary`, and `double elapsed_s`.

MEX is optional. Keep the Python backend until independent tests establish
numerical parity, correct errors, native-handle lifecycle, and long-run
stability. The block mask and ports remain unchanged when a backend is promoted.

IH uses the same Simulink production contract as every other application.
"""


CUBIT_BOUNDARY = r"""
# Cubit and Simulink boundary

Cubit owns ACIS-backed CAD, SAT/STEP import, labeling, meshing, and checked
`.vol` export. Radia's application interface owns the subsequent Simulink
configuration and solve. Keep the process boundary explicit:

1. Export SAT when ACIS fidelity matters and STEP for neutral interchange.
2. Apply the application's versioned label contract in Cubit.
3. Export `.vol` through `cubit-mesh-export` and require `check-vol` to pass.
4. Pass the checked `.vol` and configuration to the masked Simulink block.
5. Keep solver evidence in `validation_test/` JSON and presentation notebooks
   in `docs/`; do not create another desktop panel or examples tree.

Cubit's embedded PySide6 toolbar is allowed because it runs inside Cubit's
private Python. Radia itself must not acquire a Qt/PySide dependency.
"""


TOPICS = {
    "overview": APPLICATION_BLOCK_REVIEW,
    "build_application_block": BUILD_APPLICATION_BLOCK,
    "cubit_boundary": CUBIT_BOUNDARY,
    "workflow": BUILD_APPLICATION_BLOCK,
}


def get_panel_review_documentation(topic: str = "overview") -> str:
    """Return current Simulink application-block review guidance."""
    if topic == "all":
        return "\n\n".join(f"# Topic: {key}\n{value}" for key, value in TOPICS.items())
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic: {topic!r}\n\nAvailable topics:\n"
        + "\n".join(f"  - {key}" for key in TOPICS)
    )
