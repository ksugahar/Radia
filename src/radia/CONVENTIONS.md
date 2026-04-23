# Radia Analysis Window Conventions

Each analysis mode is a standalone PySide6 window defined in `src/radia/radia_*.py`.
The Cubit plugin launcher (`register_toolbar.py`) discovers these files automatically.

## Module-Level Metadata (required)

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

## Qt Compatibility

Analysis windows use PySide6 (system Python 3.12).
They do NOT run inside Cubit's process.
The launcher (`register_toolbar.py`) runs inside Cubit and uses PySide6 or PyQt5.
