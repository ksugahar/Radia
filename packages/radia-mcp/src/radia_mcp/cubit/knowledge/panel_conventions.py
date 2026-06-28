"""
Radia-NGSolve analysis window conventions for the Cubit MCP server.

Covers: radia_*.py module conventions, label requirements, sample .jou,
and the Cubit -> .vol -> NGSolve pipeline as mediated by the Cubit plugin.
"""

PANEL_CONVENTIONS = """
# Radia-NGSolve Analysis Window Conventions

## Architecture

```
Cubit Solve -> Radia-NGSolve...
  +-- Launcher dialog (Cubit PySide6 / Qt6): mode + order + folder + label check
  |     |
  |     +-- export netgen "radia_model.vol" order N overwrite
  |     +-- Launch: python radia_*.py radia_model.vol [--optional-files ...]
  |
  +-- radia_ih.py    (Induction Heating: BEM / FEM)
  +-- radia_em.py    (Electromagnet: Omega / A-Phi / MSC)
  +-- radia_pcb.py   (PCB: PEEC)
  +-- (user-defined radia_*.py auto-discovered)
```

Each analysis window is a standalone PySide6 app (Python 3.12, NOT Cubit's Python 3.10).
The .vol file is the SOLE interface between Cubit and NGSolve.
Coreform Cubit 2025.12+ is the supported in-Cubit host.  The retired
Qt5 `.ccl` launcher is not part of the current install; `cubit-plugin-install`
registers the PySide6 toolbar through a generated startup shim and
`--verify-only --all-users` checks that registration.

## Module-Level Metadata (REQUIRED in every radia_*.py)

```python
TITLE = "Induction Heating"                       # Launcher combo display name
REQUIRED_LABELS = ["source", "sink"]               # Must exist in .vol
OPTIONAL_LABELS = ["workpiece", "air"]             # Used if present
OPTIONAL_FILES = {"Coil script": "Python (*.py)"}  # Input files only
```

| Variable | Type | Purpose |
|----------|------|---------|
| TITLE | str | Shown in launcher Analysis combo |
| REQUIRED_LABELS | list[str] | Missing = red + OK disabled |
| OPTIONAL_LABELS | list[str] | Missing = gray, present = green |
| OPTIONAL_FILES | dict[str,str] | Browse fields in launcher, input only |

## Label Convention for .jou Files

### Induction Heating (radia_ih.py)
Required:
- `source` (sideset): terminal face for current injection
- `sink` (sideset): terminal face for current extraction

Optional:
- `workpiece` (block): conductive workpiece for SIBC/ESIM
- `air` (block): air domain for field calculation

### Electromagnet (radia_em.py)
Required: (none)
Optional:
- `iron` (block): ferromagnetic material
- `air` (block): air domain
- `kelvin_int` (block): interior Kelvin sphere
- `kelvin_ext` (block): exterior Kelvin sphere

### PCB (radia_pcb.py)
Required: (none)
Optional: (none) -- uses FastHenry .inp, not .vol labels

## Sample .jou Files

Each radia_*.py MUST have a corresponding sample in `panels/samples/`.
Only samples that are locked by a golden test under `tests/panels/` are
shipped in the wheel. Non-canonical IH history that is still useful for
debugging lives under `validation_test/induction_heating/demoted_samples_legacy/`
(see AGENTS.md § "Sample Promotion Ladder").

Current layout (post-2026-04-23 demotion):

```
radia_ih.py  -> panels/samples/ih_bem_sample.jou            (BEM: cubit-mesh-export smoke test)
                panels/samples/ih_peec_bem_coarse.jou        (peec_bem mode, golden)
                panels/samples/ih_fem_kelvin_skin_fine.jou   (fem_coilmesh mode, golden)
                panels/samples/ih_peec_inductance.jou        (peec_inductance mode, golden)
radia_em.py  -> panels/samples/em_sample.jou
radia_pcb.py -> panels/samples/pcb_sample.jou

Legacy validation fixtures:
  ih_fem_sample.jou              (no-Kelvin FEM baseline, truncation demo)
  ih_fem_kelvin_sample.jou/.py   (small-mesh FEM+Kelvin, misleading auto-Kelvin comment)
  ih_fem_kelvin_skin.jou/.py     (coarse variant of skin_fine)
  ih_closed_torus.jou/.py        (closed-torus research variant)
```

Naming: {stem}_sample.jou where stem is the part after radia_ in the filename.
A single radia_*.py panel may ship multiple samples when distinct solver methods
need different mesh strategies — radia_ih ships four (BEM smoke, peec_bem,
fem_coilmesh, peec_inductance).

Samples are packaged with pip install radia and serve as:
- Default working folder for the launcher
- Self-contained examples (geometry + mesh + labels, no dependencies)

## Block/Sideset Registration in .jou

```python
# Material blocks (volume elements)
block 1 add volume 1
block 1 name "iron"

# Boundary sidesets (surface elements, preferred for FEM BC)
sideset 1 add surface 1
sideset 1 name "source"

# DO NOT mix volume and surface elements in the same block
```

## .vol Export (always done by launcher)

The launcher always calls:
```
export netgen "radia_model.vol" order N overwrite
```

The analysis window receives .vol as first argument. It does NOT export.
Curve order is determined at export time in the launcher dialog.
"""

LABEL_GUIDE = """
# Label Guide for Radia-NGSolve

## How Labels Flow

```
Cubit blocks/sidesets -> export netgen -> .vol SetMaterial/SetBCName
                                          -> NGSolve mesh.Materials() / mesh.Boundaries()
```

## Naming Rules

- Use lowercase, descriptive names: "iron", "source", "air"
- Block names become material labels (mesh.Materials("iron"))
- Sideset names become boundary labels (mesh.Boundaries("source"))
- Sideset takes priority over block for boundary labels

## Checking Labels

The launcher dialog checks labels BEFORE export:
- Required labels missing: red text, OK button disabled
- Optional labels missing: gray text (informational)
- All required present: green text, OK enabled

This prevents the common mistake of forgetting to name blocks/sidesets
before running a computation that requires them.
"""
