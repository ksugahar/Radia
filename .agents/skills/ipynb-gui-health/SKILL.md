---
name: ipynb-gui-health
description: Thoroughly verify the Radia panel GUI -- promoted (2026-06-25) from PySide6 desktop panels to Jupyter notebook workbenches (src/radia/panels/notebooks/radia_*.ipynb + radia.*_notebook CommandWorkbench backed by radia.*_design DesignSpec) -- is healthy end to end. Confirms the notebooks do NOT backslide to PySide6/PyQt, that DesignSpec dataclasses are the canonical initial-value store (NOT JSON presets), that each workbench builds the correct calc_*.py and runs headless producing a radia_result.v2 artifact with timing, and that webgui/GMSH visualization is wired. Run after any change to *_notebook.py / *_design.py / radia_*.ipynb / notebook_workbench.py / panel_notebook_manifest.json, or when a notebook panel "does not run / shows no result". This replaces the old `pyside6-health` panel gate; Cubit deploy is checked separately by `cubit-plugin-install --verify-only` and `cubit-smoke-test`.
---

# ipynb-gui-health

Since 2026-06-25 (`radia.panel_notebook_promotion.v1`) the Radia panel GUIs are
**promoted to Jupyter notebook workbenches**.  Desktop PySide6 panels
(`radia_ih.py`, `radia_em.py`, ...) are now **legacy adapters**; the canonical,
browser-native interface is:

```
src/radia/panels/notebooks/radia_<app>.ipynb     # the panel = a light notebook
        |  cell: spec = <App>DesignSpec(...)      # canonical initial values
        |  cell: <App>Workbench(spec).display()   # ipywidgets form
        v
radia.<app>_notebook.<App>Workbench  ->  CommandWorkbench  (src/radia/notebook_workbench.py)
        |  build_command()  ==  DesignSpec.build_command()  -> argv for calc_*.py
        v
src/radia/panels/calc_<app>.py    (Layer-4 headless solver, JSON on stdout)
        v
runs/radia_<app>/<UTC>/  ->  run.log + command.txt + result.json (radia_result.v2)
```

The migration state lives in
`src/radia/panels/notebooks/panel_notebook_manifest.json`:
**5 panels are `active-local-runner`** (radia-ih, radia-em, radia-pcb,
radia-motor, radia-streamfunction); **radia-export-menu is `migration-shell`**
(its in-Cubit toolbar is still desktop-only -- that part stays under
`pyside6-health`, see "What stays under pyside6-health" below).

This skill answers three questions end to end:

1. **Is the GUI paradigm-pure?** -- the notebooks are notebook-native and do NOT
   backslide to desktop Qt; DesignSpec is the single source of initial values
   (JSON is a run *artifact*, never a preset store); the manifest states match.
2. **Do the workbenches actually run headless and capture durable results?** --
   `CommandWorkbench` constructs WITHOUT Jupyter/ipywidgets, `build_command()`
   emits the right `calc_*.py`, `run_local()` writes a `radia_result.v2`
   artifact with timing, and background-run / cancel / timeout work.
3. **Is the DesignSpec <-> calc_*.py wiring consistent?** -- each
   `DesignSpec.build_command()` emits argv the target `calc_*.py` accepts, and
   the visualization channels (webgui + GMSH `.msh v4.1`) are present.

**Key advantage over the old panels (and over `pyside6-health`):** the notebook
GUI has **no PySide6 dependency**, so it is fully testable **in-process with
plain pytest** -- there is no MKL/Qt `0xc0000139` DLL clash and therefore NO
offscreen-subprocess workaround is needed (unlike `audit_pyside6_only.py`).

## When to use

- After editing any `radia.<app>_notebook` (`*_notebook.py`),
  `radia.<app>_design` (`*_design.py`), the base
  `src/radia/notebook_workbench.py`, a `radia_*.ipynb`, or
  `panel_notebook_manifest.json`.
- After promoting a new panel to a notebook workbench (a new
  `active-local-runner` manifest entry).
- When a user reports "the notebook panel doesn't run / the Run button stays
  disabled / no result appears" -- check the DesignSpec `missing_required_inputs`
  / `build_command` FIRST (the Run button is intentionally disabled until
  required inputs are present -- that is fail-loud, not a bug).

## Layer A -- one-command contract check (any dev machine)

All notebook-GUI contracts are encoded in one pytest module; run it first:

```bash
python -m pytest tests/panels/test_notebook_workbench.py -q
```

Exit 0 = every contract below holds. The module locks:

| Contract | Test |
|---|---|
| `run_local()` writes a `radia_result` artifact (schema `radia.notebook_panel_run.v2`, status / version / executed/completed UTC / wall time) | `test_run_local_writes_radia_result_artifact` |
| top-4 numeric `t_*_s` timing stages harvested from the CLI JSON output | `test_run_local_collects_top_four_cli_timing_stages` |
| timeout is recorded as `status="timeout"` | `test_run_local_timeout_is_recorded` |
| a background run is cancellable (`status="cancelled"`) | `test_background_run_can_be_cancelled` |
| each workbench has an app-specific `run_root` ending `radia_<app>` | `test_promoted_workbenches_have_app_specific_run_roots` |
| each workbench `build_command()` emits the correct `calc_*.py` | `test_promoted_workbenches_build_headless_commands` |
| `spec_cell_source()` is canonical (`<App>DesignSpec(**...)`, **no JSON**) | `test_spec_cell_source_makes_notebook_initial_values_canonical` |
| manifest states correct (5 `active-local-runner` + 1 `migration-shell`) | `test_panel_notebooks_are_marked_as_local_runner` |
| notebooks do NOT import PySide6/PyQt and carry no `active-ipywidgets` stub | `test_panel_notebooks_do_not_import_pyside` |
| active notebooks use `<App>DesignSpec()` cells + "JSON files are run artifacts, not preset storage" | `test_active_panel_notebooks_use_designspec_cells_for_initial_values` |
| active notebooks include the panel notes (Run local, `netgen.webgui`, GMSH `.msh v4.1`) | `test_active_panel_notebooks_include_panel_notes` |
| the IH notebook carries ESIM + previous-result pointers | `test_ih_notebook_carries_esim_and_previous_result_notes` |

If a contract fails, fix the **source** (the notebook / DesignSpec / workbench),
do NOT relax the test.  In particular **never** add a `try ipywidgets except
PySide6` path or store presets in JSON -- both violate the promotion policy
(`panel_notebook_manifest.json::policy`).

## Layer B -- workbench runs headless + durable artifacts (smoke)

Confirms a workbench constructs and runs WITHOUT Jupyter, and leaves a durable
`radia_result.v2` beside the run -- the property that makes a notebook panel
auditable/reproducible (parallels the panels' "Result Output Persistence").

```bash
python - <<'PY'
from pathlib import Path
import tempfile, json
from radia.notebook_workbench import CommandWorkbench

class EchoSpec:
    def build_command(self):
        import sys
        return [sys.executable, "-c", "print('ipynb-gui-health ok')"]
    def missing_required_inputs(self):
        return []

with tempfile.TemporaryDirectory() as d:
    wb = CommandWorkbench(EchoSpec(), run_root=Path(d) / "runs", timeout_s=10)
    rec = wb.run_local()
    payload = json.loads(rec.result_path.read_text(encoding="utf-8"))["radia_result"]
    assert rec.status == "passed" and payload["schema"] == "radia.notebook_panel_run.v2"
    assert payload["runtime_radia_version"] and "timing" in payload
    print("[OK] headless run ->", rec.result_path)
PY
```

Each real workbench imports cleanly and builds its command without a display:

```bash
python - <<'PY'
from radia.em_notebook import EMWorkbench
from radia.ih_notebook import IHWorkbench
from radia.pcb_notebook import PCBWorkbench
from radia.motor_notebook import MotorWorkbench
from radia.streamfunction_notebook import StreamFunctionWorkbench
for wb in (EMWorkbench(), IHWorkbench(), PCBWorkbench(),
           MotorWorkbench(), StreamFunctionWorkbench()):
    assert str(wb.run_root).endswith(wb.run_root.name)
    print(f"[OK] {type(wb).__name__:26s} run_root={wb.run_root}")
PY
```

(`display()` is intentionally NOT called here -- it requires ipywidgets in a
live kernel and raises a fail-loud `RuntimeError` otherwise, which is correct.)

## Layer C -- DesignSpec <-> calc_*.py wiring

The notebook is a thin wrapper over `calc_*.py` CLI args mapped through
`DesignSpec`.  The argument-surface health is the **same** check the desktop
panels use -- reuse it:

```bash
python -m pytest tests/panels/ -q -k "golden or notebook"
```

- `DesignSpec.build_command()` must emit only flags the target `calc_<app>.py`
  `argparse` accepts (per panel: see `panel_notebook_manifest.json` `headless`
  lists -- e.g. IH -> `calc_inductance.py` / `calc_fem_kelvin.py` /
  `calc_heat.py`).  The `panel-cli-diff` skill covers the calc-CLI side; here
  the analog is the DesignSpec emitting a runnable argv.
- Visualization stays NGSolve-through: `netgen.webgui` for browser-native 3-D
  and GMSH `.msh v4.1` for field post (the active notebooks assert both in
  their notes).  No bespoke desktop viewer.

## Definition of healthy

- Layer A: `pytest tests/panels/test_notebook_workbench.py` exits 0.
- Layer B: `CommandWorkbench.run_local()` smoke passes; all 5 workbenches import
  and report an app-specific `run_root`; no PySide6 import is triggered.
- Layer C: panel golden + notebook tests pass; each DesignSpec builds a runnable
  `calc_*.py` argv; webgui + GMSH notes present.

## What Changed From `pyside6-health`

The notebook promotion is the active panel contract.  `pyside6-health` is now
retired/guard-only: normal Radia Python environments should not require
PySide6, while Coreform Cubit's bundled PySide6 is protected and must not be
deleted.  Cubit deploy health is checked separately with
`cubit-plugin-install --verify-only` and `cubit-smoke-test`.

## Related skills / tools

- `tests/panels/test_notebook_workbench.py` -- the committed Layer-A contract.
- `pyside6-health` -- retired/guard-only notes for the old PySide6 panel era.
- `panel-cli-diff` -- calc_*.py CLI <-> caller flag matching (the calc side of Layer C).
- `panel-review` -- deeper review of remaining desktop panels / adapters.
- `verify-deploy` -- "are my src/radia edits actually loaded?".
- `src/radia/CONVENTIONS.md` -- notebook workbench conventions.
- `panel_notebook_manifest.json` -- the promotion source of truth (states, adapters, headless calc lists, policy).
