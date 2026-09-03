# Radia Application Interface Conventions

The canonical human-facing Radia applications are masked blocks in
`matlab/radia_simulink_library.slx`. No application, including IH, has a
notebook workbench.

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
- Every solver-bound `.vol` is checked after export and before application
  initialization. Production application/mode meshes use a versioned strict
  label contract, and the `cubit-mesh-export.vol-check.v1` report stays in the
  run artifact directory. Mesh labels identify regions; material constants are
  validated by the DesignSpec/configuration and are never guessed from names.
- A run that computes a spatial field writes a runner-owned GMSH `.msh v4.1`
  post-processing artifact in its run directory and records it in
  `result.json`. A scalar/circuit-only mode records GMSH as not applicable; it
  does not create a dummy field.
- `cubit-mesh-export` is the only supported VTK producer. Radia application,
  validation, docs, and MCP field outputs use GMSH `.msh v4.1`; they do not call
  `ngsolve.VTKOutput` or write `.vtk`, `.vtu`, or `.vts` files. NGSolve's own
  `VTKOutput` remains a valid upstream feature; this is a Radia artifact-policy
  choice, not a restriction on NGSolve.
- A Simulink user may treat the checked `.vol` as the complete geometry and
  region-label input. Materials are authored as a MATLAB `dictionary`, mapped
  explicitly to those `.vol` region names, and compiled once during model
  initialization to the fixed-width numeric `RadiaMaterialBus`. MEX
  S-Functions never perform string or dictionary lookup per simulation step.
- The same typed winding-terminal contract feeds either the native reduced
  field/circuit state-space MEX path or an LTspice interval block. Switching
  circuit backends must not change the `.vol`, material dictionary, winding
  polarity, series/parallel identity, or mechanical state convention.
- FEMM replacement is a multiphysics contract, not a magnetics-only claim.
  `Coupling/Field Study Configuration` compiles electrostatic, DC/AC current-flow,
  steady-heat, and time-harmonic eddy-current settings to `RadiaStudyBus`.
  The corresponding `Applications/Field Study` worker runs once on an
  explicit trigger and writes versioned JSON, timing, Gmsh v4.1, and launch
  companions. A solver callback is never invoked once per Simulink step.
- Harmonic eddy current uses `(K+j*omega*M_sigma)a=S*i` and is accepted only
  when frequency, linear material law, winding identity, residual, and branch
  real-power/Joule-loss closure pass. Nonlinear harmonic permeability remains
  rejected until a separately validated iteration is available.

## Simulink Library

The single Radia Library Browser entry contains these application blocks:

```text
Applications/Electromagnet
Applications/PCB PEEC
Applications/Motor
Applications/Stream Function
Applications/Induction Heating
Applications/Magnetic Levitation
Applications/Field Study
Material Models/Material Dictionary
Coupling/Winding Dictionary
Coupling/Field Study Configuration
LTspice/LTspice Circuit
LTspice/Hysteretic LTspice Plant
```

The standard batch block has one boolean rising-trigger input and three
outputs: `int32 status`, `double primary`, and `double elapsed_s`. Status is
`0` idle, `2` passed, or `-1` failed. Full results stay in the configured run
artifact root.

`Coupling/Winding Dictionary` binds winding names to exact `.vol` regions,
including one signed polarity per coil-side region, turn count, parallel paths,
resistance, and positive/negative terminal names. Initialization compiles names
to integer ids in `RadiaWindingBus`.

Dynamic electromechanical coupling uses `RadiaMachineCommandBus` and
`RadiaMachineResponseBus`. Simulink or Simscape owns angle, speed, position,
load, controller, and mechanical integration; the field backend returns
currents, flux linkage, back EMF, torque, force, and losses. A field backend
must not hide its own incompatible motion convention.

## Samples

Canonical `.jou`, `.vol`, `.step`, BH, and related assets live under
`src/radia/panels/samples/` and are owned by the headless application contract,
not by a GUI implementation.

```text
Induction Heating block
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
in separate processes. `check-vol` is the boundary gate: it checks NGSolve
loading, curved-map Jacobians, labels and label relations, then adds CAD
volume/area/total-edge and mesh-metadata comparisons when `.vol.json` exists.

## Mesh Evaluation

Mesh p-convergence is a documentation/validation workflow, not a Cubit toolbar
action or application block. Cubit batch export produces the files and
`src/radia/panels/calc_mesh_eval.py` evaluates them.

- Cubit GMSH/BDF/VTK format QA checks Jacobians, volume, and area. VTK remains
  confined to the `cubit-mesh-export` component.
- NGSolve `.vol` p-convergence checks order 1-5 geometry integrals.
- Cubit export owns file correctness; NGSolve owns interpretation of `.vol`.

## Visualization

- Every repository-published CAE example is an executed `docs/**/*.ipynb` with
  saved WebGUI output. Use `ngsolve.webgui.Draw` for meshes and computed fields,
  and `netgen.webgui.Draw` for pre-mesh CAD geometry. A field view uses an
  explicit call such as `Draw(field, mesh, name="B_magnitude", ...)`, including
  the relevant display arguments. A static image or source excerpt does not
  replace the interactive scene.
- Simulink field-producing workflows use durable, checked GMSH `.msh v4.1`
  and JSON artifacts. GMSH is the post-processing target, not the mesh
  generator or the NGSolve interchange format.
- Application operation must not depend on a transient desktop viewer state.

## Qt Boundary

Normal Radia Python does not depend on PySide6. Coreform Cubit's embedded
PySide6 is allowed only for its export toolbar. Do not add desktop PySide6
analysis windows or new notebook workbenches.
