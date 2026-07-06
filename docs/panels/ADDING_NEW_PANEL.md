# Adding a New Panel -- Bug-Resistant Recipe

This is the current Radia pattern for promoting a calculation into a user-facing
notebook workbench. The rule is simple: the computation is a headless argparse
CLI, and the notebook/workbench is only a thin operating surface around that CLI.

> TL;DR: decide the arguments, write one `calc_<topic>.py`, lock it with a
> golden/validation test, then wrap it with a DesignSpec-backed notebook
> workbench. Do not create a separate desktop-Qt implementation.

Before writing code, use the `panel-arg-selection` skill. That Stage-1 pass
decides which variables are user knobs, solver switches, derived internals,
output paths, or geometry inputs.

## The Current Shape

```text
src/radia/panels/calc_<topic>.py                 # Stage 2: argparse CLI, headless
src/radia/<topic>_design.py                      # DesignSpec defaults for the notebook
src/radia/<topic>_notebook.py                    # CommandWorkbench wrapper
src/radia/panels/notebooks/radia_<topic>.ipynb   # Stage 3: user-facing notebook
tests/panels/test_<topic>_golden.py              # fast gate when practical
validation_test/panels/test_<topic>_*.py         # heavier gate when needed
```

The repository is migrating toward root-level `panels/`, but the existing
packaged Radia notebook workbenches still live under `src/radia/panels/`.
Follow the local pattern of the panel family you are editing and keep one live
copy of each implementation.

## Stage 2: `calc_<topic>.py`

```python
"""Headless CLI for the radia-<topic> workbench."""

from __future__ import annotations

import argparse

from radia.panels.calc_common import calc_main


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="One-line description shown by the workbench.")
    parser.add_argument("--vol", required=True, help="Input .vol mesh.")
    parser.add_argument("--frequency", type=float, default=50.0,
                        help="Driving frequency in Hz.")
    parser.add_argument("--current", type=float, default=100.0,
                        help="Source current in A.")
    parser.add_argument("--solver", choices=["hdiv", "bem", "peec"],
                        default="hdiv",
                        help="Application solver backend.")
    parser.add_argument("--fes-order", type=int, default=1,
                        help="Finite-element polynomial order.")
    parser.add_argument("--n-threads", type=int, default=4,
                        help="NGSolve TaskManager thread count.")
    parser.add_argument("--msh-output", default=None,
                        help="Gmsh .msh v4.1 artifact written by the CLI.")
    parser.add_argument("--output", default=None,
                        help="Result JSON written by calc_main().")
    return parser


def run(args) -> dict:
    import time
    from ngsolve import Mesh, TaskManager

    t0 = time.perf_counter()
    with TaskManager(numthreads=args.n_threads):
        mesh = Mesh(args.vol)
        # ... solve ...

    return {
        "ne": mesh.ne,
        "ndof": 0,
        "t_total_s": round(time.perf_counter() - t0, 3),
        "gmsh_file": args.msh_output or "",
    }


def main() -> None:
    calc_main(run, build_argparser())


if __name__ == "__main__":
    main()
```

### CLI Musts

- `build_argparser()` must be importable without side effects.
- Heavy imports (`ngsolve`, `radia`, `cubit`) happen inside `run()`, not at
  module import time.
- `TaskManager` is imported before the `with TaskManager(...):` statement in
  the same function.
- The solver switch is a real CLI argument when the application exposes more
  than one supported backend.
- Unknown values raise errors. Do not silently substitute defaults.
- Result dictionaries include `ne` or `wp_ne`, `ndof` or `wp_ndof`, `t_*_s`
  timing keys, headline physical quantities, and the generated artifact path
  (`gmsh_file`, `field_gmsh_file`, `msh_file`, or `msh_output`) when present.

## Stage 3: Notebook Workbench

Notebook panels use a small `DesignSpec` dataclass plus a `CommandWorkbench`.
They persist initial values in code cells, launch the Stage-2 CLI, and collect
`run.log`, `result.json`, timing summaries, and mesh/field artifacts.

The notebook must not re-implement the computation. It maps editable settings to
CLI flags and displays ecosystem-native outputs such as `netgen.webgui` scenes
or durable Gmsh files.

### Notebook Musts

- Store persistent defaults in `DesignSpec`, not ad hoc JSON presets.
- Build the subprocess command from the same CLI flags accepted by
  `calc_<topic>.py`.
- Record Radia version, Python/platform context, total wall time, and the main
  timing stages in `result.json`.
- Keep small domain notes in the notebook when they prevent misuse.
- Use `netgen.webgui` for human-facing notebook visualization.
- Use `.msh v4.1` / `.json` artifacts for automation and LLM validation.

## Golden And Validation Tests

Fast panel checks belong in `tests/panels/` when they are deterministic and
small. Solver-heavy sweeps, convergence studies, and timing claims belong in
`validation_test/panels/` and are mdx-first when they are too large for CI.

```python
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CALC = ROOT / "src" / "radia" / "panels" / "calc_<topic>.py"


def test_golden_band(tmp_path):
    out_json = tmp_path / "result.json"
    out_msh = tmp_path / "result.msh"
    cmd = [
        sys.executable, str(CALC),
        "--vol", str(ROOT / "tests" / "panels" / "fixtures" / "<topic>.vol"),
        "--output", str(out_json),
        "--msh-output", str(out_msh),
    ]
    subprocess.run(cmd, check=True, timeout=600)
    got = json.loads(out_json.read_text(encoding="utf-8"))

    assert "t_total_s" in got
    assert "gmsh_file" in got
```

## Automated Audits

| Audit | What it catches |
|---|---|
| `panel-arg-selection` skill | Stage-1 variable classification before widgets exist |
| `panel-cli-diff` skill | CLI/workbench flag drift |
| `ipynb-gui-health` skill | DesignSpec, notebook artifact, and no desktop-Qt regression |
| `tests/panels/test_taskmanager_scoping.py` | late-import `TaskManager` errors |
| `tests/panels/test_panel_output_health.py` | output JSON/log/summary contract |
| `validation_test/panels/test_notebook_workbench.py` | notebook workbench contract |

## Reference Implementations

- `src/radia/panels/calc_accel_hdiv.py` -- HDiv-VIM electromagnet CLI.
- `src/radia/em_design.py` and `src/radia/em_notebook.py` -- notebook
  workbench pattern for electromagnet operation.
- `src/radia/panels/notebooks/radia_em.ipynb` -- user-facing notebook surface.
- `src/radia/panels/calc_volume.py` -- small one-shot CLI example.

## Anti-Patterns

| Anti-pattern | Correct pattern |
|---|---|
| Hand-assembled argv that drifts from argparse | Generate command arguments from the CLI surface |
| Heavy top-level imports in a calc module | Import heavy solver packages inside `run()` |
| A notebook cell that re-solves with local logic | Notebook launches the CLI and displays artifacts |
| A solver option that falls back silently | Raise on unknown or unsupported values |
| Timing keys such as `t_solve` | Use `t_solve_s`, `t_total_s`, etc. |
| Artifact keys such as `output_msh` | Use recognized keys such as `gmsh_file` |

---

Last updated: 2026-07-06. Maintained alongside the Panel Design Workflow Policy
in `AGENTS.md` / `CLAUDE.md`.
