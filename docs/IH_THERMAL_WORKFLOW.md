# IH Eddy-Thermal Simulink Workflow

Radia's production induction-heating interface is the masked **Induction
Heating** block in `matlab/radia_simulink_library.slx`. The retired `radia-ih`
PySide panel and notebook workbench are not supported interfaces.

## Runtime Structure

The block exposes the coupled field state explicitly:

```text
coil_current_A ---------> Eddy ---------> heat_density_W_per_m3
workpiece_angle_rad ----> Eddy                    |
temperature_K ----------> Eddy                    v
                                               Thermal ----> temperature_K
ambient_temperature_K --------------------------> |
workpiece_angle_rad -----------------------------> |
```

- `radia.simulink.ihEddySFunction` is a readable Level-2 MATLAB S-Function.
  It owns three input ports, one distributed heat output, sample time, and the
  native Eddy handle lifecycle.
- `radia.simulink.ihThermalSFunction` is a separate Level-2 MATLAB S-Function.
  It publishes the accepted temperature in `Outputs` and advances the thermal
  state in `Update`, giving the feedback loop an explicit one-step delay.
- `radia_mex('ih.*', ...)` owns the numerical objects behind checked `uint64`
  handles. Python is not called once per simulation step.

The fixed update order is:

```text
Eddy output -> conservative transport(theta_prev, theta_now) -> Thermal update
```

Temperature is stored in the workpiece coordinate system. In
`periodic-uniform` mode, rotation transports the accepted temperature field
conservatively before the implicit thermal solve. Eddy maps the source-frame
heating into workpiece coordinates.

## Geometry And Operators

Use `radia.simulink.assembleIHOperatorsFromGeometry` or the
`radia-ih-assemble` CLI at the explicit model-update boundary. The runtime
does not mesh or assemble finite-element operators during a time step.

```powershell
python -m pip install "radia[cubit]"
radia-ih-assemble workpiece.vol coil.step
radia-ih-assemble workpiece.vol coil.vol -o native_ih.json
```

- Workpiece geometry is a checked Netgen `.vol` or `.vol.gz` mesh.
- A coil `.step` selects PEEC assembly.
- A coil `.vol` selects BEM-A assembly.
- Physical method selection may also name FEM or BIM when the supplied
  configuration contains the corresponding validated operators.
- Every solver-bound mesh must pass its versioned
  `radia.vol-label-contract.v1` contract before initialization.
- The `cubit-mesh-export.vol-check.v1` reports, `run.log`, `result.json`, and
  GMSH `.msh v4.1` post-processing fields stay in the run directory.

The shape-only assembler supports fixed linear material operators. It declares
rotation disabled when a conservative transport operator is unavailable; it
does not guess a node ordering or silently substitute a lumped/LUT model.

## Configuration

`radia.simulink.makeIHNativeConfig` validates and converts assembled operators
to the row-major numeric ABI consumed by `radia_mex`. A production
configuration includes:

- complex Eddy matrix and right-hand side;
- distributed heat projection and heat-cell weights;
- thermal mass and stiffness matrices in CSR form;
- heat-to-temperature projection and temperature-cell weights;
- initial temperature, sample time, and rotation convention;
- physical method, linear solver, thermal solver, and checked `.vol` reports.

The native preview currently accepts `bh_mode="linear"`. If validated
temperature-dependent real and imaginary Eddy-matrix slopes are supplied, a
temperature change rebuilds the Eddy operator. For a fixed linear operator,
changing current scales the unit-current solution and heat quadratically
without rebuilding the operator.

Configuration is loaded from a MAT or JSON file by
`radia.simulink.configureIHNativeModel`. The model stores the validated config
in its model workspace, not in a process-global GUI object.

## Open And Run

From MATLAB:

```matlab
addpath("matlab");
install_radia_ih();
```

The installer opens the tracked `matlab/radia_ih.slx`. Use its Geometry Update
block to replace the workpiece and coil paths and rebuild the checked operator
artifact. Current and angle remain ordinary Simulink source signals.

## Verification

The release gates cover these contracts separately:

- standalone MEX commands: ABI, numerical behavior, errors, stale handles;
- Level-2 S-Functions: initialization, Outputs/Update ordering, termination,
  repeated runs, and closed-loop wiring;
- rotation: zero, small, large, periodic, and conservative-energy cases;
- geometry update: `.vol`, `.vol.gz`, `.step`, role normalization, label
  contracts, and persistent artifacts;
- application model: clean open, compile, simulation, and full-window visual QA.

Numerical evidence and benchmarks belong under `validation_test/` with JSON
results. Public demonstrations belong in executed `docs/**/*.ipynb` notebooks
with saved, parameterized WebGUI scenes. Neither replaces the Simulink
production interface.
