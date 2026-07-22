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

The historical `build_notebook_gui` topic name is retained for MCP clients,
but the construction target is now a masked Simulink block.

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


CUBIT_PANELS_MIGRATION = r"""
# `examples/cubit_panels` migration plan

`examples/cubit_panels` is retired. Existing historical material is routed by
durable role:

- reusable geometry/solver/parser code -> `src/radia`
- numerical checks, `verify_*.py`, and comparisons -> `validation_test`
- explanatory demonstrations -> result-saved docs notebooks + JSON
- final human operation -> the Radia Simulink library
- canonical assets -> `src/radia/panels/samples/` beside their headless owner

The induction-heating legacy corpus lives under
`validation_test/induction_heating/cubit_panels_legacy`. Rescued accelerator
fixtures live under `src/radia/panels/samples/em/c_type_dipole`. Historical
files such as `coil_dipole.py`, `create_induction_model.py`, and `verify_*.py`
must not be copied into a new examples or notebook-workbench tree.

Cubit remains the mesh producer. Its C++ export plugin and embedded PySide6
toolbar are separate from MATLAB/Simulink and Python 3.12; `.vol`/`.sol` plus
artifact files are the process boundary.
"""


TOPICS = {
    "overview": APPLICATION_BLOCK_REVIEW,
    "build_notebook_gui": BUILD_APPLICATION_BLOCK,
    "presentation_template": BUILD_APPLICATION_BLOCK,
    "cubit_panels_migration": CUBIT_PANELS_MIGRATION,
    "5_skills_chain": APPLICATION_BLOCK_REVIEW,
    "13_checks": APPLICATION_BLOCK_REVIEW,
    "bug_catalogue": APPLICATION_BLOCK_REVIEW,
    "val_checkbox_trap": APPLICATION_BLOCK_REVIEW,
    "map_value_reject": APPLICATION_BLOCK_REVIEW,
    "widget_calc_gap": APPLICATION_BLOCK_REVIEW,
    "smoke_scenarios": APPLICATION_BLOCK_REVIEW,
    "red_flags": APPLICATION_BLOCK_REVIEW,
    "workflow": BUILD_APPLICATION_BLOCK,
}


def get_panel_review_documentation(topic: str = "overview") -> str:
    """Return application-block review guidance; keep old topic names stable."""
    if topic == "all":
        return "\n\n".join(f"# Topic: {key}\n{value}" for key, value in TOPICS.items())
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic: {topic!r}\n\nAvailable topics:\n"
        + "\n".join(f"  - {key}" for key in TOPICS)
    )
