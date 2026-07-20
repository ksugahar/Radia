"""Cubit-to-Radia Simulink application and mesh-label conventions."""

PANEL_CONVENTIONS = r"""
# Radia Simulink application conventions for Cubit users

## Process architecture

```text
Cubit 2025.12 embedded Python/PySide6
  -> C++ export netgen
  -> self-contained .vol/.sol files
  -> Radia Simulink application block (separate MATLAB process)
  -> DesignSpec + headless calc_*.py (separate Python 3.12 process)
  -> run.log / solver_result.json / result.json
```

Cubit is a mesh producer, not the application GUI host. The single Radia
Simulink library owns EM, PCB, Motor, Stream Function, and IH human operation.
IH temporarily also keeps `radia_ih.ipynb` for comparison; all other packaged
analysis workbenches are removed.

The source of truth is
`src/radia/panels/application_interface_manifest.json`, each application's
`DesignSpec`, and its `calc_*.py` argparse contract. Simulink must not import
Cubit's Python runtime, and Cubit must not import normal `radia`/NGSolve.

## Label conventions

### Induction Heating block / `IHDesignSpec`

For BEM-A or full FEM coil paths:
- `source` sideset: current injection terminal
- `sink` sideset: current extraction terminal

Typical materials/boundaries include `workpiece`, `coil`, `air`, `sibc`,
`kelvin_int`, and `kelvin_ext` according to the selected headless method.

### Electromagnet block / `EMDesignSpec`

Typical labels are `iron`, `air`, `kelvin_int`, and `kelvin_ext`. Exact
requirements belong to the selected method and its CLI/sample, not to Cubit's
toolbar.

### PCB PEEC block / `PCBDesignSpec`

The application consumes a FastHenry `.inp`; it does not require `.vol` block
labels.

### Stream Function block / `StreamFunctionDesignSpec`

The block consumes one or more checked surface/evaluation `.vol` files.
Material-aware modes may additionally consume iron/shield meshes.

## Canonical samples

Golden-locked assets live under `src/radia/panels/samples/`:

```text
IH block / comparison notebook -> ih_bem_sample.jou,
                                  ih_peec_bem_coarse.jou,
                                  ih_fem_kelvin_skin_fine.jou,
                                  ih_peec_inductance.jou
Electromagnet block            -> em_sample.jou
PCB PEEC block                 -> pcb_sample.jou
```

Samples are owned by the headless application contract, not by a notebook.
Research history belongs under `validation_test/`.

## Export boundary

The Cubit toolbar exports; it does not launch a non-IH notebook or run the
solver. Curve order and labels are fixed at export time. Application config
then references the durable files and the Simulink block executes explicitly.
"""


LABEL_GUIDE = r"""
# Label Guide for Radia-NGSolve

## How labels flow

```text
Cubit blocks/sidesets -> export netgen -> .vol SetMaterial/SetBCName
                                          -> NGSolve Materials/Boundaries
```

## Rules

- Use lowercase descriptive names such as `iron`, `source`, `sink`, and `air`.
- Block names become material labels; sideset names become boundary labels.
- Do not mix volume and surface entities in one block.
- Validate labels against the selected application's headless sample/CLI.
- A Simulink mask must not silently rename or substitute a missing label.

## Example

```python
block 1 add volume 1
block 1 name "iron"
sideset 1 add surface 1
sideset 1 name "source"
```

Use `check-vol`/NGSolve inspection before a physics solve. Missing required
labels fail and leave an application diagnostic artifact.
"""
