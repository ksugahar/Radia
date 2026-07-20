"""Historical standalone-panel topic redirected to Simulink application blocks."""

STANDALONE_PANELS = r"""
# Radia application interfaces after the Simulink migration

Standalone PySide analysis windows are retired. EM, PCB, Motor, and Stream
Function notebook workbenches are also retired. The final human operating
surface is a masked block in the single Radia Simulink library. Python/MCP is
the AI surface; `docs/**/*.ipynb` is result-bearing documentation.

IH is the temporary comparison exception: keep both the Induction Heating
block and `radia_ih.ipynb` over the same `IHDesignSpec` and headless scripts.
If the block passes numerical, artifact/recovery, and operational gates, an
explicit policy decision may retire the IH notebook too.

Topics: quick_start, four_panels, build_notebook_gui, cubit_panels_migration,
        vol_sources, vs_cubit, ih_methods, troubleshooting

============================================================
## quick_start -- Radia Simulink library
============================================================

From MATLAB:

```matlab
addpath("matlab")
radia.setup()
radia.simulink.buildLibrary()
```

The Library Browser entry contains:

```text
Applications/Electromagnet
Applications/PCB PEEC
Applications/Motor
Applications/Stream Function
Applications/Induction Heating
```

Each standard application block consumes a versioned configuration JSON and
runs only on a boolean rising trigger. Outputs are `int32 status` (`0` idle,
`2` passed, `-1` failed), the selected primary scalar, and elapsed seconds.
Full data remains in `command.txt`, `run.log`, `solver_result.json`, and
`result.json` under the configured run root.

Fast contract checks:

```powershell
python -m pytest tests/test_simulink_application.py -q
matlab -batch "addpath('matlab'); r=runtests('tests/matlab/test_simulink_workflow.m'); assert(all([r.Passed]))"
```

============================================================
## four_panels -- compatibility topic; five active blocks
============================================================

The historical topic name is retained for MCP clients. Current surfaces:

| Application | Simulink block | Shared DesignSpec | Notebook |
|---|---|---|---|
| Electromagnet | `Applications/Electromagnet` | `radia.em_design.EMDesignSpec` | none |
| PCB/PEEC | `Applications/PCB PEEC` | `radia.pcb_design.PCBDesignSpec` | none |
| Motor | `Applications/Motor` | `radia.motor_design.MotorDesignSpec` | none |
| Stream Function | `Applications/Stream Function` | `radia.streamfunction_design.StreamFunctionDesignSpec` | none |
| Induction Heating | `Applications/Induction Heating` | `radia.ih_design.IHDesignSpec` | temporary comparison only |

Do not restore the deleted non-IH `*_notebook.py` adapters or packaged
notebooks as compatibility aliases.

============================================================
## build_notebook_gui -- retired alias; build an application block
============================================================

This historical topic now redirects to the block recipe in
`docs/panels/ADDING_NEW_PANEL.md`.

1. Enumerate all application variables and solver switches.
2. Implement and golden-lock `calc_<topic>.py` plus `<Topic>DesignSpec`.
3. Add a result-bearing `docs/<topic>/*.ipynb` only for explanation.
4. Add a masked block to `radia.simulink.buildLibrary`.
5. Delegate to `radia.simulink.application`; do not implement solver logic in
   mask callbacks or launch Python every time step.
6. Lock the mask, ports, success/failure artifacts, and headless numerical
   parity in Python and MATLAB tests.

The initial backend may be the validated Python/headless CLI. MEX/ROM is a
later optional backend and requires independent parity, error propagation,
handle lifecycle, and long-run stability tests before promotion.

============================================================
## cubit_panels_migration -- Cubit boundary
============================================================

Cubit owns geometry/mesh export through its embedded PySide6 toolbar and C++
export plugin. It exports self-contained `.vol`/`.sol` files. Simulink and
headless Python consume those files in separate processes; neither imports
Cubit's Python runtime.

Historical `examples/cubit_panels` code is not a destination. Reusable methods
go to `src/`, heavy truth runs to `validation_test/`, explanatory results to
`docs/`, samples to `src/radia/panels/samples/`, and final human operation to
the Radia Simulink library.

============================================================
## vol_sources -- durable inputs
============================================================

Application blocks consume durable `.vol`, `.sol`, `.step`, `.msh`, `.inp`,
and material/config files. Cubit is one producer, but operation must not depend
on a transient viewer state. Human documentation may use `netgen.webgui`;
automation and block runs use durable GMSH/JSON artifacts.

============================================================
## vs_cubit -- deploy boundary
============================================================

Cubit export checks remain:

```powershell
cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
cubit-smoke-test
```

These checks do not validate the Radia Simulink library. Use the Python/MATLAB
application tests for that surface. Never delete or replace Cubit's bundled
PySide6 because normal Radia does not depend on it.

============================================================
## ih_methods -- temporary dual operation
============================================================

IH uses both `Applications/Induction Heating` and
`src/radia/panels/notebooks/radia_ih.ipynb` during the comparison. Both map to
`IHDesignSpec` and the same PEEC/BEM/FEM/thermal `calc_*.py` commands. Compare
setup effort, failure recovery, result inspection, automation, and throughput
on identical inputs. The intended direction is Simulink, but the notebook is
retired only after the explicit acceptance gate.

============================================================
## troubleshooting -- common interface issues
============================================================

| Symptom | Correct response |
|---|---|
| Normal Python lacks PySide6 | Expected; use Simulink/Python/MCP, not a desktop panel. |
| Application block reports `-1` | Inspect its `UserData`, `run.log`, and `result.json`; fix config/dependency/solver failure. |
| No `result.json` | Run `tests/test_simulink_application.py`; verify the selected Python can import the same Radia checkout. |
| Block attempts Python each step | Use the rising-edge application runner or a tested prebuilt MEX/ROM state. |
| MEX compiles but behavior is uncertain | Keep the Python backend; do not promote MEX until parity/lifecycle/long-run gates pass. |
| Cubit toolbar fails | Run Cubit plugin verification; do not alter Cubit's private PySide6 runtime. |
| Old non-IH notebook path is missing | Expected after migration; use the corresponding Radia library block. |
"""


_TOPICS = (
    "quick_start", "four_panels", "build_notebook_gui",
    "cubit_panels_migration", "vol_sources", "vs_cubit", "ih_methods",
    "troubleshooting",
)


def get_standalone_panels_documentation(topic: str = "") -> str:
    """Return legacy topic knowledge redirected to Simulink blocks."""
    topic = (topic or "").strip()
    if not topic:
        return STANDALONE_PANELS
    if topic not in _TOPICS:
        return f"Unknown topic {topic!r}.  Valid topics: " + ", ".join(_TOPICS)

    headers = []
    pos = 0
    while True:
        nxt = STANDALONE_PANELS.find("\n## ", pos)
        if nxt < 0:
            break
        line_end = STANDALONE_PANELS.find("\n", nxt + 1)
        line = STANDALONE_PANELS[nxt + 1:line_end]
        for candidate in _TOPICS:
            if line.startswith(f"## {candidate} "):
                headers.append((candidate, nxt + 1))
                break
        pos = nxt + 1

    starts = [offset for name, offset in headers if name == topic]
    if not starts:
        return f"Topic {topic!r} declared but not found in document."
    section_start = starts[0]
    separator = STANDALONE_PANELS.rfind(
        "============================================================", 0, section_start
    )
    if separator >= 0 and separator > section_start - 80:
        section_start = separator

    section_end = len(STANDALONE_PANELS)
    for name, offset in headers:
        if name != topic and offset > starts[-1]:
            separator = STANDALONE_PANELS.rfind(
                "============================================================", 0, offset
            )
            section_end = separator if separator > 0 else offset
            break
    return STANDALONE_PANELS[section_start:section_end].rstrip() + "\n"
