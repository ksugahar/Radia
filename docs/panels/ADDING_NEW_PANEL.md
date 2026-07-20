# Adding a Radia Simulink Application Block

Radia's final human operating surface is a masked block in the single Radia
Simulink library. Python and MCP remain the AI surface, and result-bearing
notebooks under `docs/` explain and reproduce methods. A documentation notebook
is not an application GUI.

IH is the temporary comparison exception: its notebook workbench and Simulink
block remain supported over the same headless contract until the block passes
the operational migration gates.

## Four Stages

```text
Stage 1  application variable inventory
Stage 2  src/radia/panels/calc_<topic>.py + src/radia/<topic>_design.py
Stage 3  docs/<topic>/*.ipynb + synchronized result JSON
Stage 4  matlab/+radia/+simulink/* + matlab/radia_simulink_library.slx
```

Do not jump from a scratch script to a block. Use the
`panel-arg-selection` skill for Stage 1 even though the final surface is now
Simulink; it classifies user knobs, solver switches, derived values, geometry
inputs, and output paths.

## Stage 2: Headless Contract

The calculation remains an argparse CLI with an import-safe parser and a
`run(args) -> dict` implementation. A UI-neutral `DesignSpec` maps settings to
that exact CLI. The block must not build a second solver or configuration
language.

```python
from dataclasses import dataclass
import sys

from radia.panel_design import calc_script, json_output


@dataclass
class TopicDesignSpec:
    vol: str = ""
    frequency: str = "50"
    solver: str = "hdiv"

    def missing_required_inputs(self) -> list[str]:
        return [] if self.vol else ["Mesh .vol"]

    def build_command(self, *, python: str | None = None) -> list[str]:
        if not self.vol:
            raise ValueError("Mesh .vol is required.")
        return [
            python or sys.executable,
            calc_script("calc_topic.py"),
            "--vol", self.vol,
            "--frequency", self.frequency,
            "--solver", self.solver,
            "--output", json_output(self.vol, "_topic"),
        ]
```

The CLI result includes element count, DoF, `t_*_s` timings, headline physical
quantities, and artifact paths. Unknown solver modes and missing dependencies
fail loudly. Heavy solver imports stay inside the execution function.

## Stage 3: Documentation

Create a result-bearing notebook only when it helps readers understand the
method, equations, inputs, and representative output. Execute it before commit
and synchronize its adjacent JSON, including runtime/version metadata and the
notebook hash. Do not add widgets or a `CommandWorkbench` adapter.

## Stage 4: Simulink Block

Add the block to `radia.simulink.buildLibrary`. The standard batch-analysis
block uses:

- `radia_application_sfun` as the Level-2 MATLAB S-function entry point;
- a fixed application id and an evaluated mask for configuration JSON, run
  root, timeout, and Python executable;
- a boolean trigger input; execution occurs only on a rising edge;
- `int32 status`, `double primary`, and `double elapsed_s` outputs;
- `radia.simulink.runApplication` and `radia.simulink.application` for the
  versioned artifact contract.

Status values are `0` idle, `2` passed, and `-1` failed. Full results are in the
artifact directory rather than squeezed into a fixed-width signal.

The configuration file is versioned and maps directly to the application's
`DesignSpec`:

```json
{
  "schema": "radia.simulink.application_config.v1",
  "application": "em",
  "settings": {
    "coil_script": "C:/models/coil.py",
    "vol": "C:/models/magnet.vol",
    "method": "Omega Reduced"
  },
  "primary_key": "B_origin_mag_T",
  "working_directory": "C:/models"
}
```

Create it from MATLAB rather than hand-editing JSON:

```matlab
settings = struct("coil_script", "C:/models/coil.py", ...
    "vol", "C:/models/magnet.vol", "method", "Omega Reduced");
radia.simulink.writeApplicationConfig("em", settings, ...
    "C:/models/em_config.json", PrimaryKey="B_origin_mag_T");
```

The Simulink/MATLAB launch path writes:

```text
<run-root>/<application>/<UTC stamp>/
  launcher_command.txt
  command.txt          # after DesignSpec command construction
  run.log
  solver_result.json   # after solver output generation
  result.json
```

`result.json` is always required, including launcher failure. It uses
`radia.simulink.application_run.v1` and records the backend, status, command,
versions, platform, elapsed time, config hash, and selected primary scalar.

## Backend Policy

The validated Python/headless CLI is the initial application-block backend.
This is deliberate: a Simulink interface does not require an unproven MEX
implementation. Python is launched once on an explicit trigger, never once per
time step.

A MEX/ROM backend may later replace the internals without changing mask or port
contracts. Promotion requires independent numerical parity, checked error
propagation, handle/lifecycle tests, and long-run stability. Do not advertise a
MEX path merely because it compiles.

## Required Gates

- `python calc_<topic>.py --help` exits successfully and lists every user knob.
- Every supported solver switch passes its headless sample/golden.
- `python -m pytest tests/test_simulink_application.py -q` passes.
- `tests/matlab/test_simulink_workflow.m` finds the masked block, updates a
  model containing it, and covers success and failure artifacts.
- The block and headless paths agree within the same documented golden band.
- Timeout, missing dependency, invalid config, and solver failure remain
  inspectable in `run.log` and `result.json`.

Long solver validation runs on an idle mdx or hibino host. LAB and 100号機 are
for build, import, mask, and fast numerical checks.

## Current References

- `src/radia/simulink/application.py`: shared Python artifact runner.
- `matlab/+radia/+simulink/runApplication.m`: MATLAB process adapter.
- `matlab/+radia/+simulink/applicationSFunction.m`: rising-edge block runtime.
- `matlab/+radia/+simulink/buildLibrary.m`: single-library builder.
- `src/radia/panels/application_interface_manifest.json`: production states.
- `src/radia/ih_notebook.py`: temporary IH-only comparison adapter.

## Anti-Patterns

| Anti-pattern | Required pattern |
|---|---|
| New notebook workbench for an application | Masked block over `DesignSpec` and the headless CLI |
| Solver code inside a mask callback or S-function | Delegate to tested Radia API/CLI/MEX backend |
| Python subprocess on every simulation step | Explicit update/rising trigger or prebuilt MEX/ROM state |
| MEX promoted because it builds | Promote only after parity, lifecycle, failure, and long-run tests |
| Silent fallback to another solver | Fail with a recorded diagnostic artifact |
| Results only in block labels or the MATLAB console | Versioned JSON plus verbatim `run.log` |

Last updated: 2026-07-20.
