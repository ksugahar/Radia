"""
Retired standalone PySide panel knowledge.

The old radia-ih / radia-em / radia-pcb / radia-streamfunction desktop
launchers were superseded by Jupyter notebook workbenches.  Keep the MCP topic
for compatibility, but make every answer point to the notebook route.
"""

STANDALONE_PANELS = """\
# Radia panel launch after notebook migration

The PySide6 standalone panel route is retired.  The canonical panel surface is
now a result-bearing Jupyter notebook workbench:

```
panels/notebooks/radia_<app>.ipynb          # target layout
src/radia/panels/notebooks/radia_<app>.ipynb # legacy during migration
radia.<app>_design.<App>DesignSpec
radia.<app>_notebook.<App>Workbench
```

Normal Radia Python on LAB / 100号機 / mdx / hibino should not install PySide6.
The only protected PySide6 runtime is Coreform Cubit's embedded Python runtime;
do not uninstall or delete that Cubit-owned copy.

Topics: quick_start, four_panels, build_notebook_gui, cubit_panels_migration,
        vol_sources, vs_cubit, ih_methods, troubleshooting

============================================================
## quick_start -- notebook panel route
============================================================

Install the production packages without the old GUI extra:

```powershell
pip install --upgrade "radia[cubit]" radia-mcp cubit-mesh-export
```

Open or execute the notebook workbench for the application:

```
src/radia/panels/notebooks/radia_ih.ipynb
src/radia/panels/notebooks/radia_em.ipynb
src/radia/panels/notebooks/radia_pcb.ipynb
src/radia/panels/notebooks/radia_motor.ipynb
src/radia/panels/notebooks/radia_streamfunction.ipynb
```

For automated checks from a source checkout:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

The notebook saves a `radia_result.v2` JSON artifact beside the run, and result
notebooks are committed with outputs plus synchronized sidecar JSON.

============================================================
## four_panels -- active notebook workbenches
============================================================

Active notebook workbenches:

| Notebook | Workbench | DesignSpec |
|----------|-----------|------------|
| `radia_ih.ipynb` | `radia.ih_notebook.IHWorkbench` | `radia.ih_design.IHDesignSpec` |
| `radia_em.ipynb` | `radia.em_notebook.EMWorkbench` | `radia.em_design.EMDesignSpec` |
| `radia_pcb.ipynb` | `radia.pcb_notebook.PCBWorkbench` | `radia.pcb_design.PCBDesignSpec` |
| `radia_motor.ipynb` | `radia.motor_notebook.MotorWorkbench` | `radia.motor_design.MotorDesignSpec` |
| `radia_streamfunction.ipynb` | `radia.streamfunction_notebook.StreamFunctionWorkbench` | `radia.streamfunction_design.StreamFunctionDesignSpec` |

The old `radia_*.py` PySide modules were removed; do not restore them as a
compatibility alias.

============================================================
## build_notebook_gui -- construction recipe
============================================================

For the full construction checklist, call the dedicated panel-review MCP topic:

```
panel_review(topic="build_notebook_gui")
```

Short version:

- Move reusable kernels to `src/`; keep heavy numerical gates in
  `validation_test/`.
- Create `<App>DesignSpec` and `<App>Workbench` around a headless
  `panels/calc_*.py` command.  Existing `src/radia/panels/calc_*.py` scripts
  are legacy-compatible during the staged migration.
- The notebook cell imports `DesignSpec` + `Workbench`, creates `spec`, and
  calls `workbench.display()`.
- `CommandWorkbench` runs locally and saves `command.txt`, `run.log`, and
  `result.json` with `radia.notebook_panel_run.v2` metadata.
- Do not store presets in JSON; persistent defaults live in the notebook
  `DesignSpec(...)` cell.
- Presentation CSS may restyle the page but must not put a Jupyter
  cell-selection layer over run buttons.

The NGSolve User Meeting `RADIA-IH.ipynb` draft is a presentation shell for
the IH workbench: title Markdown, optional dark CSS, `IHDesignSpec` +
`IHWorkbench`, and short tips.  Keep the repository notebook and validation
tests as the canonical contract.

============================================================
## cubit_panels_migration -- examples/cubit_panels route
============================================================

`examples/cubit_panels` is not a permanent destination.  The IH inductance
scripts have moved to `validation_test/induction_heating/cubit_panels_legacy`,
and the remaining accel-magnet examples were pruned after rescuing panel
fixtures to `src/radia/panels/samples/em/c_type_dipole`.  Any future legacy copy should
move into one of these lanes:

- reusable accel magnet geometry / coil builders -> `src/radia` EM APIs
- panel-only samples, notebooks, and calc wrappers -> repo-root `panels/`
- IH validation scripts (`verify_*`, `compare_*`, `test_*`) ->
  `validation_test/induction_heating` or a specific validation subtree
- IH demonstrations (`scalar_bie_sibc.py`, `bem_sibc_workpiece.py`,
  `efie_sibc.py`, `fem_esim_*.py`, `impedance_esim.py`) -> src kernels plus
  result-saved docs notebooks
- Cubit journals, `.geo`, and BH tables -> protected assets until all
  references point at the new owner

For the detailed original 35-script routing plan, call:

```
panel_review(topic="cubit_panels_migration")
```

============================================================
## vol_sources -- mesh inputs still matter
============================================================

Notebook workbenches still consume durable mesh/input files such as `.vol`,
`.sol`, `.step`, and `.msh`.  Cubit is one producer via
`cubit-plugin-install` + `export netgen`, but notebook workflows should not
depend on a transient desktop viewer state.

Human-facing visualization uses `netgen.webgui`.  Headless/LLM validation uses
durable GMSH `.msh v4.1`, JSON, and saved notebook outputs.

============================================================
## vs_cubit -- Cubit boundary
============================================================

The in-Cubit export toolbar is a Cubit plugin surface and may use Coreform
CUBIT's embedded PySide6.  That is separate from normal Radia Python.

Production deploy checks:

```powershell
cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
cubit-smoke-test
```

Notebook panel health checks:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

============================================================
## ih_methods -- IH through the notebook workbench
============================================================

Use `radia_ih.ipynb` with `IHDesignSpec` and `IHWorkbench`.  The notebook is a
thin UI over the validated `calc_*.py` scripts; it does not re-implement the
solver.  JSON files are run artifacts, not preset storage.

============================================================
## troubleshooting -- common post-migration issues
============================================================

| Symptom | Correct response |
|---------|------------------|
| `ModuleNotFoundError: No module named 'PySide6'` in normal Python | Expected after migration. Use the notebook workbench; do not install PySide6 for production. |
| Notebook produces no `result.json` | Run `ipynb-gui-health` / `pytest validation_test/panels/test_notebook_workbench.py -q` and fix DesignSpec or workbench wiring. |
| Cubit export toolbar fails | Check `cubit-plugin-install --verify-only --all-users` and `cubit-smoke-test`; do not delete Cubit's bundled PySide6. |
| Old `radia-ih` executable is missing | Use `src/radia/panels/notebooks/radia_ih.ipynb`; the executable route is not the canonical surface. |
"""


_TOPICS = (
    "quick_start", "four_panels", "build_notebook_gui",
    "cubit_panels_migration", "vol_sources", "vs_cubit", "ih_methods",
    "troubleshooting",
)


def get_standalone_panels_documentation(topic: str = "") -> str:
    """Return retired standalone-panel knowledge, redirected to notebooks."""
    topic = (topic or "").strip()
    if not topic:
        return STANDALONE_PANELS
    if topic not in _TOPICS:
        return (f"Unknown topic {topic!r}.  Valid topics: "
                + ", ".join(_TOPICS))

    headers = []
    pos = 0
    while True:
        nxt = STANDALONE_PANELS.find("\n## ", pos)
        if nxt < 0:
            break
        line_end = STANDALONE_PANELS.find("\n", nxt + 1)
        line = STANDALONE_PANELS[nxt + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, nxt + 1))
                break
        pos = nxt + 1

    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."

    section_start = req_starts[0]
    sep_above = STANDALONE_PANELS.rfind(
        "============================================================",
        0, section_start)
    if sep_above >= 0 and sep_above > section_start - 80:
        section_start = sep_above

    last_req = req_starts[-1]
    section_end = len(STANDALONE_PANELS)
    for kw, off in headers:
        if kw != topic and off > last_req:
            sep = STANDALONE_PANELS.rfind(
                "============================================================",
                0, off)
            section_end = sep if sep > 0 else off
            break

    return STANDALONE_PANELS[section_start:section_end].rstrip() + "\n"
