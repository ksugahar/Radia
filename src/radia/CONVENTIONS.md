# Radia Panel Conventions

The canonical Radia panel surface is now a lightweight Jupyter notebook
under `src/radia/panels/notebooks/`.  A panel wraps an app-specific
headless `calc_*.py` script: argparse options become `DesignSpec`
settings, the notebook workbench launches the script, and run artifacts
are saved as `run.log` plus `result.json`.

Legacy standalone PySide6 windows remain under `src/radia/radia_*.py`
while existing Cubit desktop workflows still depend on them.  The Cubit
plugin launcher (`register_toolbar.py`) discovers those files for
backward compatibility, but new user-facing panel work should prefer the
notebook route.

## Notebook Panel Policy

- Keep computation in headless CLI/function modules.  Notebook panels
  map settings to CLI arguments; they do not re-implement solvers.
- Persistent defaults live in the notebook `DesignSpec(...)` cell.
  JSON files are run artifacts, not preset storage.
- Use concise Markdown where it teaches the application.  For IH, that
  includes the linear SIBC vs nonlinear ESIM choice: ESIM solves a 1-D
  cell problem from a BH curve to update the effective surface
  impedance `Z_s`, and the previous ESIM benchmark JSONs show
  convergence and per-panel/scalar behavior.
- Use `netgen.webgui` for human-facing notebook visualization and
  durable GMSH `.msh v4.1` artifacts for LLM/headless validation.

## Legacy PySide Module-Level Metadata (required while PySide remains)

Every `radia_*.py` must define these module-level variables:

```python
TITLE = "Induction Heating"                    # Display name in launcher combo
REQUIRED_LABELS = ["source", "sink"]            # Block/sideset names that must exist
OPTIONAL_LABELS = ["workpiece", "air"]          # Block/sideset names shown but not required
OPTIONAL_FILES = {"Coil script": "Python (*.py)"}  # Optional input file browse fields
```

| Variable | Type | Description |
|----------|------|-------------|
| `TITLE` | `str` | Shown in the launcher Analysis combo box |
| `REQUIRED_LABELS` | `list[str]` | Block/sideset names required in .vol. Missing = red text + OK disabled |
| `OPTIONAL_LABELS` | `list[str]` | Block/sideset names used if present. Missing = gray, present = green |
| `OPTIONAL_FILES` | `dict[str, str]` | Key = display name, value = file dialog filter. Input files only (no output) |

## Sample Journal Files (required)

Each `radia_*.py` must have a corresponding sample `.jou` file in `src/radia/panels/samples/`.
Only golden-test-locked samples (`tests/panels/test_*_golden.py`) are shipped in the wheel.

Current layout (post-2026-04-23 demotion):

```
radia_ih.py   ->  samples/ih_bem_sample.jou            (BEM: cubit-mesh-export smoke test)
                  samples/ih_peec_bem_coarse.jou       (peec_bem mode, golden-locked)
                  samples/ih_fem_kelvin_skin_fine.jou  (fem_coilmesh mode, golden-locked)
                  samples/ih_peec_inductance.jou       (peec_inductance mode, golden-locked)
radia_em.py   ->  samples/em_sample.jou
radia_pcb.py  ->  samples/pcb_sample.jou
```

Non-canonical / research-stage samples live under `examples/induction_heating/demoted_samples/`
(see that directory's README for the demotion history — e.g., the old
`ih_fem_sample.jou` no-Kelvin baseline and the small-mesh `ih_fem_kelvin_sample.jou`
with its misleading "auto-Kelvin" comment).

Naming convention: `{stem}_sample.jou` where `{stem}` is the part after `radia_` (e.g., `radia_em.py` -> `em`). Multiple samples per stem are allowed when distinct solver methods need different mesh strategies — radia_ih ships four shipped samples (one per panel mode).

These samples are:
- Packaged with `pip install radia` (included in wheel)
- Used as the default working folder in the launcher
- Self-contained (geometry + mesh + blocks/sidesets, no external dependencies)

## .vol is Always Required

The launcher always exports `.vol` before launching the analysis window.
Every analysis window receives `.vol` as its first argument.
There is no `REQUIRES_VOL` flag -- it is always true.

## Working Folder Memory

The launcher remembers the user's last working folder per machine:

- Stored in: `~/.cubit/radia_launcher.json`
- Key: `"last_jou_dir"`
- Fallback chain: last folder -> parent folder -> package samples folder

## Entry Point

Each `radia_*.py` must define a `main()` function for console_scripts:

```python
def main():
    run_app(MyWindow)

if __name__ == "__main__":
    main()
```

## Mesh Evaluation Policy

**Format QA (our quality guarantee)**:
- .msh, .bdf, .vtk at order 1-2 verified via GMSH API `getJacobians()`
- Volume and area compared against ACIS CAD values
- Negative Jacobian determinant = inverted element (node ordering bug)
- This is the Cubit plugin's responsibility

**p-Convergence (curving accuracy)**:
- .vol at order 1-5 read by NGSolve `Integrate(CF(1), mesh)`
- Volume, area, length compared against CAD
- Monotonic convergence expected (each order gains ~2-3 digits)
- NGSolve reading .vol correctly is NGSolve's responsibility

**Responsibility boundary**:
- We guarantee: export files (.msh/.bdf/.vtk) contain correct geometry
- NGSolve guarantees: .vol is read correctly by `Mesh("model.vol")`

## Visualization Routing

Panel visualization should ride existing viewer ecosystems.

- GUI / notebook route for humans: use `netgen.webgui` (`Draw(...)`,
  browser/Jupyter scene widgets).
- LLM / headless route for automation: use `.msh v4.1` plus `gmsh`
  / `GmshPostExport` so screenshots, field views, and validation
  artifacts are durable files.
- `.vol` and `.sol` are notebook input/output files. Their
  double-click viewer should be plain Netgen (`netgen.exe "%1"`), not
  GMSH and not the Radia-specific viewer by default.

Do not make a panel depend on a GMSH desktop window for ordinary
interactive viewing, and do not make an LLM workflow depend on a
transient `netgen.webgui` browser state.

Implementation note: pip Netgen's command-line dispatch lives in
`netgen.__main__`; package startup/DLL setup lives in `netgen.__init__`.
If `.vol` double-click does not launch, add a `.vol` handler to
`netgen.__main__` using Netgen's native `Ng_LoadMesh` sequence and keep
the Windows association pointed at `netgen.exe "%1"`. The
`radia-vol-viewer` association is a legacy/helper path for custom `.sol`
companion-mesh inference, not the default notebook IO route.

## Qt Compatibility

Analysis windows use PySide6 (system Python 3.12).
They do NOT run inside Cubit's process.
The launcher (`register_toolbar.py`) runs inside Cubit and uses PySide6 (Cubit 2025.12 ships PySide6; the Qt5/PyQt5 path was removed in radia 4.80.0).
