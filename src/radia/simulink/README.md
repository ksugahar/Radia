# Radia Simulink Blocks

This directory contains the Simulink boundary for production Radia reduced
models.  The numerical kernels remain in the Radia C++ API; the Simulink
files only define the block contract and marshal MATLAB arrays.

## Application blocks

The final human operating surfaces for Electromagnet, PCB PEEC, Motor, Stream
Function, and Induction Heating are built by
`matlab/+radia/+simulink/buildLibrary.m`. Electromagnet, PCB, Motor, and Stream
Function retain their explicit-trigger application boundary. Induction
Heating is separate: readable Level-2 MATLAB Eddy and Thermal S-Functions own
its ports and lifecycle, while checked `radia_mex` handles own native numerical
state. Initialization locates NGSolve's runtime through the configured Python
installation, but no Python solver or per-step fallback is launched. The first
distributable IH package is a preview for preassembled operators; native `.vol`
operator assembly remains its production gate.

For a field-producing mode, the runner overrides `--msh-output` with a path in
the run directory, requires a valid GMSH `.msh v4.1` file, and lists every
generated GMSH artifact in `result.json`. Scalar/circuit-only modes explicitly
remain non-spatial. GMSH is used for post-processing only; solver interchange
continues to use Netgen `.vol`.

MEX/ROM remains an optional backend promotion for the four batch blocks and is
mandatory for IH. Compilation alone is not a production gate; failure
propagation, handle lifecycle, repeated runs, and long-run stability must pass
while the block mask and port contract stay unchanged. Notebook workbenches are
retired for every application, including IH.

## Motor ROM block

`radia_motor_rom_sfun.cpp` is a fixed-step C-MEX S-function.  Build it from a
MATLAB session with:

```matlab
addpath('src/radia/simulink');
radia_motor_rom_sfun_build('S:/Radia/01_GitHub');
```

The S-function takes five parameters:

```text
RadiaMotorROM, sample_time_s, initial_state, max_iterations, tolerance
```

`RadiaMotorROM` is the struct loaded by the generated `*_load.m` bundle
helper.  `initial_state` is a MATLAB struct with optional fields
`time_s`, `rotor_angle_rad`, `rotor_speed_rad_s`,
`temperature_K`, and `generalized_currents_A`.

The single input vector is ordered as:

```text
[phase_voltages_V; load_torque_Nm; ambient_temperature_K]
```

The single output vector is ordered as:

```text
[phase_currents_A;
 eddy_currents_A;
 phase_flux_linkage_Wb;
 rotor_angle_rad;
 rotor_speed_rad_s;
 electromagnetic_torque_Nm;
 resistive_loss_W;
 hysteresis_loss_W;
 temperature_K;
 energy_balance_residual_W;
 nonlinear_iterations]
```

The output is the accepted state after the previous block update.  The ROM
step is performed once per fixed sample time in `mdlUpdate`, so the block is
used as a discrete Co-Simulation-style plant.  The adapter rejects bundles
that require an external hysteresis callback; that callback needs a separate
stateful bridge and must not silently fall back to a linear model.

## Production block policy

The Simulink layer exposes physical model units, not individual internal
matrix kernels:

* MagLev LTI models use the existing MATLAB `ss` / Simulink State-Space block.
* Dynamic motor ROMs use `radia_motor_rom_sfun` or an FMI Co-Simulation FMU.
* Fixed reduced IH and HCurl Eddy Bubble/CLN models can use
  `matlab/+radia/+simulink/stateSpaceMexSFunction.m`, backed by the native
  `simulink.state_space.*` MEX handle commands. Their matrices are copied once
  at `Start`; no Python process or per-step state-vector transfer is used.
* Moving height-family models remain MATLAB S-functions because their reduced
  operators are interpolated as the mechanical state changes.

The C ABI remains the canonical implementation boundary so that Python,
Simulink, FMI, and other hosts cannot silently acquire different numerical
kernels.

## Induction-heating native preview

The tracked `matlab/radia_ih.slx` model separates readable Level-2 MATLAB Eddy
and Thermal S-Functions and closes the temperature feedback explicitly:

```matlab
addpath('matlab');
install_radia_ih();
```

The wrappers call independent `radia_mex('ih.*', ...)` object-handle commands
that consume checked row-major electromagnetic operators, a
heat-to-temperature projection, thermal CSR matrices, cell weights, and an
explicit rotation mode. Thermal publishes only the accepted DWork state and
advances it in `Update`, so the closed loop has a real one-step delay.
Rotation transports the accepted workpiece temperature conservatively before
the implicit thermal solve; Eddy maps the stationary-source heat distribution
back into workpiece coordinates.

The release model contains a one-DOF diagnostic configuration. Physical cases
must supply preassembled operators through a checked MAT/JSON configuration.
The preview does not yet assemble PEEC/BEM-A/BIM/FEM operators from Cubit
`.vol` files. LUT and lumped-state-space helpers are neither reachable from
`openIH` nor included in the native IH release package.

For learning and design optimization, wrap a `sim` call or the fast waveform
function in an objective and use:

```matlab
objective = @(x) ihObjectiveFromSimulation(x);  % returns one scalar
opt = radia.simulink.optimizeIH(objective, x0, lower, upper);
```

The adapter uses `fmincon` when Optimization Toolbox is available and a
bounded `fminsearch` fallback otherwise. A Simulink objective should create a
fresh `Simulink.SimulationInput` and set model variables there; this makes the
same function compatible with `parsim` for batch learning runs.
