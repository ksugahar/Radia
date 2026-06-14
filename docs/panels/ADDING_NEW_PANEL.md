# Adding a New Panel — Bug-Resistant Recipe

This is the **canonical** pattern for adding a new analysis panel to
Radia (e.g. `radia_electromagnet`, `radia_thermal_axisym`, anything
else).  The pattern below makes the panel and the underlying CLI
**share a single source of truth** so the entire class of "widget
silently drops a flag" / "argparse rejects a flag" / "summary
forgets to display the value" bugs cannot occur by construction.

> **TL;DR.**  Write **one** argparse-driven Python function.
> The GUI panel, the subprocess wiring, the `.log` capture, the
> Open-GMSH button, the JSON result file, the cross-machine deploy
> are all derived from that one function automatically.

> **Before you write File 1**, decide *which* variables from the source
> example become arguments at all -- that is **Stage 1** of the Panel
> Design Workflow (CLAUDE.md) and is covered by the
> **`panel-arg-selection`** skill
> (`.claude/skills/panel-arg-selection/SKILL.md`).  This recipe assumes
> the `build_argparser()` surface has already been chosen.

## The 3-file shape

```
src/radia/panels/calc_<topic>.py          # 1. The argparse-driven CLI.  Headless.  No PySide6.
src/radia/radia_<topic>.py                # 2. The PySide6 panel.  Wraps (1) via subprocess.
tests/panels/test_<topic>_golden.py       # 3. The golden test that locks the result band.
```

That is the **entire** new-panel footprint.  Everything else (form
widgets, JSON dump, .log persistence, GMSH button enable, deploy
hooks, cross-machine consistency) reuses code that already exists in
`radia_gui_base.AnalysisWindow` / `radia_gui_base.ModePanel`.

## File 1: `calc_<topic>.py` (CLI / headless)

```python
"""calc_<topic>.py -- one-line description of what this analyses.

CLI entry point for the radia-<topic> panel.  Pure headless: imports
NGSolve / Radia inside calc() ONLY, never at module top, so the panel
can introspect build_argparser() without paying the NGSolve import
cost (and without triggering the UnboundLocalError trap if calc() ever
does its own late `from ngsolve import TaskManager`).
"""

import argparse
import sys
import os

# Allow imports from src/radia when run as subprocess.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radia.panels.calc_common import calc_main   # standard subprocess wrapper


def build_argparser() -> argparse.ArgumentParser:
    """Return the argparse.  MUST be importable WITHOUT side effects --
    the panel calls this to introspect widgets, so it cannot trigger
    NGSolve / Radia / Cubit imports.  All heavy imports go in calc()."""
    p = argparse.ArgumentParser(
        description="One-liner that becomes the panel-mode tooltip.")
    p.add_argument("--vol", required=True, help="Workpiece .vol mesh.")
    p.add_argument("--frequency", type=float, default=50.0,
                   help="Driving frequency in Hz.")
    p.add_argument("--current", type=float, default=100.0,
                   help="Coil current in A.")
    p.add_argument("--solver", choices=["mmm", "msc", "bem"], default="mmm",
                   help="Solver backend.")
    p.add_argument("--fes-order", type=int, default=1,
                   help="FE-space polynomial order.")
    p.add_argument("--n-threads", type=int, default=4,
                   help="NGSolve TaskManager thread count.")
    # standard panel outputs -- calc_main fills these if absent
    p.add_argument("--msh-output", default=None,
                   help="GMSH .msh path to write.  Panel sets this.")
    p.add_argument("--output", default=None,
                   help="Result JSON.  calc_main writes the return dict here.")
    return p


def run(args) -> dict:
    """Do the work.  Heavy imports happen HERE, not at module top.

    Returns a dict that MUST satisfy the Result Output Policy
    (2026-05-29 + 2026-05-30):
      - 'ne' or 'wp_ne'                : element count (int)
      - 'ndof' or 'wp_ndof'            : DOF count (int)
      - 't_<step>_s' for every timed step (mesh / solve / total / ...)
      - integral quantities specific to the analysis (e.g. 'B_max_T',
        'P_wp_W', 'T_mean_C / T_max_C / T_min_C')
      - 'gmsh_file' (or 'field_gmsh_file' / 'msh_file' / 'msh_output')
        pointing to a real .msh on disk if the analysis produced one
        (used to auto-enable the panel's Open-GMSH button)

    Keep this function pure data-in / data-out.  Do NOT touch the
    filesystem outside of the .msh emit + the JSON (calc_main handles
    that)."""
    import time
    from ngsolve import Mesh, TaskManager, ...    # ALL imports stay here
    from radia.panels.calc_common import ...

    t0 = time.perf_counter()
    with TaskManager():
        mesh = Mesh(args.vol)
        # ... solve ...
        t_solve_s = time.perf_counter() - t0

    if args.msh_output:
        # write .msh -> args.msh_output
        ...

    return {
        "ne":         mesh.ne,
        "ndof":       fes.ndof,
        "t_solve_s":  round(t_solve_s, 2),
        "t_total_s":  round(time.perf_counter() - t0, 2),
        "B_max_T":    B_max,
        "gmsh_file":  args.msh_output if args.msh_output else "",
    }


def main():
    parser = build_argparser()
    calc_main(run, parser)   # run(args) -> dict; calc_main writes JSON


if __name__ == "__main__":
    main()
```

### MUST DOs

- **One** `build_argparser() -> argparse.ArgumentParser` function.
  Importable WITHOUT side effects.
- **One** `run(args) -> dict` function.  All heavy imports inside.
  (The name `run` is the convention per `calc_common.calc_main`'s
  docstring; `calc` is also accepted.)
- `calc_main(run, parser)` somewhere in the module (legacy convention
  is to wrap it in `def main():`).
- Return dict has standardized keys (see "Result Output Policy").

### MUST NOT DOs

- **NO top-level `import ngsolve / import radia / import cubit`**
  (would break panel introspection).
- **NO `from ngsolve import TaskManager` placed AFTER `with
  TaskManager():` inside the same function.**  Python compiles the
  function as a whole, so the late import promotes `TaskManager` to a
  local for the whole function and the earlier `with TaskManager():`
  raises `UnboundLocalError` (keiko 100号機 2026-05-30).  Always put
  TaskManager in the **first** `from ngsolve import ...` line of calc().
- **NO silent default substitution.**  If a parameter is required
  and missing, raise -- per the "No Fallbacks" policy.
- **NO mixed-purpose names like `t_end_s` or `t_ext_C`.**  The panel's
  `_append_standard_summary` matches `t_*_s` as compute-time.  Use
  `sim_end_time_s` / `T_ext_C` (no `t_` prefix on non-time values).

## File 2: `radia_<topic>.py` (PySide6 panel)

```python
"""radia_<topic> panel -- generator-driven (POLICY 2026-05-30)."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, json_output, run_app,
)

TITLE = "<Topic display name>"
REQUIRED_LABELS = []                 # .vol material/boundary labels we need
OPTIONAL_LABELS = ["air", "kelvin"]


class _TopicPanel(ModePanel):
    def __init__(self, parent=None):
        super().__init__(parent)
        from radia.panels.calc_<topic> import build_argparser
        self.bind_argparser(
            build_argparser(),
            # `vol` is routed via AnalysisWindow's wp-vol input, not a panel widget.
            # `output` is set by build_command_from_parser; `msh_output` by extra=.
            skip=("vol", "output", "msh_output"),
            file_browse={
                # Map widget keys -> (Label, file-filter) for QFileDialog.
                # Only if the calc takes a file path.
            },
            choice_map={
                # Map widget keys -> [(display_label, cli_value), ...]
                # Use this when CLI choices need human-readable labels.
                # "solver": [("LU","lu"), ("BiCGSTAB","bicgstab"), ...]
            },
            labels={
                # widget_key -> "Custom display label:"  (else auto-derived)
            },
        )

    def is_runnable(self):
        # Optional: extra gate beyond "vol path looks valid".
        return True

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path,
            vol_flag="--vol",
            script_path=calc_script("calc_<topic>.py"),
            output_path=json_output(vol_path, "_<mode_suffix>"),
            extra=["--msh-output", msh_output(vol_path, "_<mode_suffix>")],
        )


class TopicWindow(AnalysisWindow):
    """AnalysisWindow handles Run / Stop / Output / Open GMSH / status bar.
    Subclass just to set the title + register the panel."""
    def __init__(self, vol_path="", parent=None):
        super().__init__(title=TITLE, vol_path=vol_path, parent=parent)
        self._set_panel(_TopicPanel(self))


if __name__ == "__main__":
    run_app(TopicWindow, sys.argv)
```

### MUST DOs

- One panel = one `ModePanel` calling `bind_argparser()` exactly once.
- The **same** `<mode_suffix>` string in `json_output(vol_path,
  "_<mode_suffix>")` and `msh_output(vol_path, "_<mode_suffix>")` --
  the Result Output Persistence Policy derives the .log path from the
  --output JSON path, so the suffix anchor MUST match.
- Multi-mode panel = a `QStackedWidget` holding several
  `ModePanel` subclasses, one per mode (see `radia_em.py` for the
  canonical example).

### MUST NOT DOs

- No `import ngsolve / radia / cubit` at panel module top
  (panel module loads on every Cubit toolbar launch).
- No custom widget construction outside `bind_argparser()` /
  ModePanel helpers (avoid the "orphan widget" class of bugs:
  widget defined, never read).
- No bypassing `build_command_from_parser` to hand-assemble argv
  (defeats the single-source-of-truth invariant).

## File 3: `tests/panels/test_<topic>_golden.py`

```python
"""Golden lock for radia-<topic>: assert the canonical sample yields
the expected result band.  Runs the CLI as a subprocess (no panel)
so the test is portable -- pytest runs it on any machine that has
NGSolve installed."""

import json
import subprocess
import sys
from pathlib import Path

CALC = (Path(__file__).resolve().parents[2]
        / "src" / "radia" / "panels" / "calc_<topic>.py")
SAMPLE_VOL = Path(__file__).parent / "fixtures" / "<topic>_sample.vol"
GOLDEN = Path(__file__).parent / "golden" / "<topic>_sample.json"


def test_golden_band(tmp_path):
    out_json = tmp_path / "result.json"
    out_msh  = tmp_path / "result.msh"
    cmd = [sys.executable, str(CALC),
           "--vol", str(SAMPLE_VOL),
           "--output", str(out_json),
           "--msh-output", str(out_msh)]
    subprocess.run(cmd, check=True, timeout=600)

    got = json.loads(out_json.read_text())
    golden = json.loads(GOLDEN.read_text())

    # Hard band on the headline integral.  Tighten as accuracy improves.
    assert abs(got["B_max_T"] - golden["B_max_T"]) / golden["B_max_T"] < 0.05, (
        f"B_max_T drift: got {got['B_max_T']}, golden {golden['B_max_T']}")

    # Standard housekeeping keys MUST be present.
    for k in ("ne", "ndof", "t_total_s", "B_max_T", "gmsh_file"):
        assert k in got, f"calc output missing required key '{k}': {sorted(got)}"
```

## Automated audits already in CI

After you add the 3 files above, the following gates auto-cover them
(no extra wiring):

| Audit | What it catches |
|---|---|
| `panel-cli-diff` skill | Panel widget vs CLI flag mismatch -- bug class A |
| `tests/panels/test_taskmanager_scoping.py` | Late-import TaskManager UnboundLocalError -- bug class B (keiko 100号機 2026-05-30) |
| `tests/panels/test_panel_output_health.py` | Persistence Policy + 10pt + scroll-area + result-summary keys -- bug classes C/D |
| `panel-qt-test` skill | Widget visibility / Run-button enable -- bug class E |
| `Result Output Policy` (CLAUDE.md) | Missing `t_*_s` / `ne` / `ndof` / integral keys -- bug class F |
| `Result Output Persistence Policy` (CLAUDE.md) | `.log` not written / IH-style super-override skip -- bug class G |

## Where each runtime feature comes from

Once you have the 3 files, the user gets all of these for free:

| Feature | Source |
|---|---|
| Auto-generated form widgets | `ModePanel.bind_argparser()` reads argparse + makes widgets |
| Run button / Stop button | `AnalysisWindow` base class |
| Output window (live stdout/stderr) | `AnalysisWindow._read_stdout / _read_stderr` |
| Result JSON written to disk | `calc_main()` wraps your `calc()` return dict |
| `.log` written next to JSON | `AnalysisWindow._persist_output_log()` (Result Output Persistence Policy 2026-05-30) |
| Open-GMSH button auto-enable | `AnalysisWindow._on_finished` matches `gmsh_file` / `field_gmsh_file` / `msh_output` / `msh_file` in result |
| Compute-time table in summary | `AnalysisWindow._append_standard_summary` matches every `t_*_s` key |
| Temperature mean/max/min table | Same helper matches `T_mean_C` / `T_max_C` / `T_min_C` |
| Cross-machine deploy hooks | `release_triple.py phase8` deploys `calc_<topic>.py` like any other panel file |
| Layer-3 standalone launch | `python -m radia.radia_<topic> model.vol` (handled by `run_app()`) |

## Reference implementations

Read these for working examples of the pattern:

- `src/radia/radia_em.py` -- multi-mode panel (Omega / A-Phi / MSC /
  KelvinBench), `QStackedWidget` over 4 generator-driven sub-panels.
  ~360 lines for 4 modes.
- `src/radia/panels/calc_accel_magnet.py` -- the CLI counterpart
  with `build_argparser()` + `calc()` + `calc_main()`.
- `src/radia/panels/calc_volume.py` -- minimal single-mode example
  if your panel is one-shot (no Method combo).

## Anti-patterns observed in real bug reports

These are the actual bug classes that motivated this recipe.  Do NOT
do these things:

| Anti-pattern | Real bug | Don't |
|---|---|---|
| Hand-assemble `cmd = [python, "calc.py", "--vol", vol, "--freq", str(self._freq.value())]` | Widget keys drift from CLI flags silently | Use `build_command_from_parser()` |
| Top-level `import ngsolve` in calc | Panel introspection times out at module load | Imports go inside `calc()` |
| `from ngsolve import TaskManager` inside a `try:` block after `with TaskManager():` | keiko 100号機 `UnboundLocalError` 2026-05-30 (`calc_verify_vol.py`) | Always top-of-function import |
| Return `{"t_solve": 12.3}` instead of `{"t_solve_s": 12.3}` | Panel summary never shows the solve time | Match the `t_*_s` convention exactly |
| Return `{"output_msh": path}` instead of `{"gmsh_file": path}` | Open-GMSH button stays disabled | Use one of the four recognized keys |
| Add a widget but never read it via `self.val("key")` | Orphan widget; user typed value silently ignored | `bind_argparser` keeps widget keys in lockstep with argparse |
| `if user_passed_solver == "ams" or "pardiso"` then default to BBOX without warning | HACApK `SetClusterStrategy` silent fallback (audit 2026-05-30) | Raise on unknown value -- "No Fallbacks" policy |
| IH-style override that calls `super()._on_finished()` first then appends summary | `.log` ended before the summary -- triage requested copy-paste | Override re-calls `self._persist_output_log()` at the END (see `radia_ih.py:1784`) |

## When NOT to follow this recipe

This recipe is for the **standard analysis panel** (input .vol or
.step + parameters -> JSON result + .msh + .log).  Use a different
pattern when:

- The panel is a **dialog**, not an analysis (e.g. an Export Mesh
  config sheet that just emits a Cubit APREPRO command).  These live
  in `panels/radia_export_menu.py` and use the `.ccm` directly, no
  argparse / calc_*.py.
- The panel needs **interactive 3D widgets** (vtkRenderWindow, etc.)
  -- those don't fit `bind_argparser`.  See `coil_builder.py`
  for the precedent if you genuinely need this.

---

Last updated: 2026-05-31.  Maintained alongside `CLAUDE.md` "Panel
Design Workflow Policy" + "Result Output Policy" + "Result Output
Persistence Policy" sections.
