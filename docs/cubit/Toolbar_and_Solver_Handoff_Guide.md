# Cubit Toolbar and Solver Handoff Guide

Radia no longer adds solver application panels to Cubit. Cubit owns CAD,
meshing, labels, and checked Netgen `.vol` export. Human solver operation lives
in the masked blocks of the Radia Simulink library; Python/MCP provides the AI
operation surface over the same headless application contracts.

## Supported Cubit Extension

`cubit-mesh-export` installs one Cubit-owned toolbar for mesh export and
diagnostics. Its PySide6 dependency is supplied by Cubit's private Python
runtime and must not leak into the `radia` package or an external solver
process.

The toolbar may:

- import or inspect ACIS SAT and STEP geometry;
- create and inspect Cubit blocks, sidesets, and nodesets;
- export a high-order Netgen `.vol` file;
- run `check-vol` and display its diagnostic report;
- export VTK only as a Cubit mesh-export capability.

It must not become another solver GUI, launch a Radia application dialog, or
duplicate a Simulink mask.

## Handoff Contract

```text
Cubit / cubit-mesh-export
  CAD (SAT or STEP)
  -> mesh and labels
  -> Netgen .vol
  -> check-vol --strict-labels

Radia Simulink block or Python/MCP
  validated .vol + DesignSpec/configuration
  -> canonical headless solver API
  -> run.log + result.json + spatial GMSH .msh
```

The `.vol` file carries topology, element order, material-region names, and
boundary names. Physical values such as conductivity, permeability, BH data,
frequency, current, and thermal properties belong to the application
configuration, not to Cubit UI state.

Every solver-bound `.vol` must satisfy its versioned
`radia.vol-label-contract.v1` contract before initialization. Keep the
`cubit-mesh-export.vol-check.v1` report with the run artifacts and reference it
from `result.json`.

## Labeling Rules

- Use lowercase ASCII names and underscores.
- Put volumes in material blocks and faces in boundary sidesets.
- Do not infer physical constants from a label.
- Keep application-specific requirements in a checked label-contract file.
- Treat missing, duplicate, or unexpected required labels as errors.

Common IH labels include `coil`, `workpiece`, `air`, `kelvin`, `source`,
`sink`, `sibc`, `kelvin_int`, `kelvin_ext`, `outer`, and `GND`; the selected
solver mode determines which are required.

## Development Rules

- Cubit-loaded Python may import Cubit's bundled PySide6 for the export toolbar.
- `src/radia`, `radia-mcp`, solver CLIs, and Simulink code must not depend on Qt
  or PySide6.
- Reusable solver behavior belongs to importable `radia.*` APIs. A `calc_*.py`
  entry point only parses configuration, calls that API, and writes artifacts.
- A Cubit toolbar callback must not contain electromagnetic or thermal solver
  logic.
- Production spatial results are GMSH `.msh v4.1`; VTK remains available only
  inside `cubit-mesh-export` and NGSolve's upstream API.

## Verification

```powershell
check-vol model.vol --strict-labels --contract path\to\labels.json
cubit-plugin-install --verify-only --all-users
python tools/audit_pyside6_only.py
```

Application acceptance is performed through the corresponding Simulink block
and its headless golden validation, not through a Cubit solver panel.
