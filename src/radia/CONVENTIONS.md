# Radia Application Interface Conventions

The canonical human-facing Radia applications are masked blocks in
`matlab/radia_simulink_library.slx`. EM, PCB, Motor, and Stream Function have
no notebook workbench. IH temporarily keeps
`src/radia/panels/notebooks/radia_ih.ipynb` beside its Simulink block so the two
operating styles can be compared over the same `IHDesignSpec` and headless CLI.

Python and MCP are the first-class AI interfaces. Result-bearing notebooks
under `docs/` explain and reproduce methods; they are not production GUIs.

## Application Contract

- Computation lives in tested Radia APIs and headless
  `src/radia/panels/calc_*.py` modules.
- UI-neutral `*_design.py` dataclasses map application settings to the exact
  argparse contract.
- Simulink blocks delegate to that contract and do not reimplement solvers.
- The initial application backend may launch Python once on an explicit rising
  trigger. It must not launch Python every simulation step.
- MEX/ROM is an optional later backend after numerical parity, lifecycle,
  failure propagation, and long-run tests.
- Every run leaves `launcher_command.txt`, `run.log`, and versioned
  `result.json` provenance. Solver command/output artifacts are added once the
  run reaches those stages.

## Simulink Library

The single Radia Library Browser entry contains these application blocks:

```text
Applications/Electromagnet
Applications/PCB PEEC
Applications/Motor
Applications/Stream Function
Applications/Induction Heating
```

The standard batch block has one boolean rising-trigger input and three
outputs: `int32 status`, `double primary`, and `double elapsed_s`. Status is
`0` idle, `2` passed, or `-1` failed. Full results stay in the configured run
artifact root.

## Samples

Canonical `.jou`, `.vol`, `.step`, BH, and related assets live under
`src/radia/panels/samples/` and are owned by the headless application contract,
not by a GUI implementation.

```text
Induction Heating block / IH comparison notebook
  -> ih_bem_sample.jou
  -> ih_peec_bem_coarse.jou
  -> ih_fem_kelvin_skin_fine.jou
  -> ih_peec_inductance.jou

Electromagnet block -> em_sample.jou + em_sample_coil.py + em_sample_bh.txt
PCB PEEC block      -> pcb_sample.jou
```

Only golden-locked samples ship as canonical production inputs. Research and
superseded formulations stay in `validation_test/` or local `C:\temp`, not as
alternate interface implementations.

## Cubit Boundary

The Cubit Export Mesh toolbar remains a Cubit-owned PySide6 surface. Normal
Radia Python and Simulink must not load Cubit's private Python runtime. Cubit
exports self-contained `.vol`/`.sol` assets; applications consume those files
in separate processes.

## Mesh Evaluation

Mesh p-convergence is a documentation/validation workflow, not a Cubit toolbar
action or application block. Cubit batch export produces the files and
`src/radia/panels/calc_mesh_eval.py` evaluates them.

- GMSH/BDF/VTK format QA checks Jacobians, volume, and area.
- NGSolve `.vol` p-convergence checks order 1-5 geometry integrals.
- Cubit export owns file correctness; NGSolve owns interpretation of `.vol`.

## Visualization

- Human documentation notebooks use `netgen.webgui` when an interactive scene
  teaches the method.
- Simulink/headless/AI workflows use durable GMSH `.msh v4.1` and JSON
  artifacts.
- Application operation must not depend on a transient desktop viewer state.

## Qt Boundary

Normal Radia Python does not depend on PySide6. Coreform Cubit's embedded
PySide6 is allowed only for its export toolbar. Do not add desktop PySide6
analysis windows or new notebook workbenches.
