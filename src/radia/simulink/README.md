# Radia Simulink Blocks

This directory contains the Simulink boundary for production Radia reduced
models.  The numerical kernels remain in the Radia C++ API; the Simulink
files only define the block contract and marshal MATLAB arrays.

## Application blocks

The final human operating surfaces for Electromagnet, PCB PEEC, Motor, Stream
Function, and Induction Heating are built by
`matlab/+radia/+simulink/buildLibrary.m`. Their initial backend is
`radia.simulink.application`: a boolean rising trigger launches the validated
`DesignSpec`/headless Python CLI once and records `command.txt`, `run.log`,
`solver_result.json`, and versioned `result.json`. It is not a per-time-step
Python bridge.

MEX/ROM remains an optional backend promotion. Compilation alone is not a
production gate; numerical parity, failure propagation, handle lifecycle, and
long-run stability must pass while the block mask and port contract stay
unchanged. IH temporarily keeps its notebook workbench for comparison.

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

## Induction-heating control workflow

The MATLAB package under `matlab/+radia/+simulink` is the control-side
surface for the IH solver. It deliberately separates the fast electromagnetic
calculation from the slower thermal and mechanical loop:

```matlab
addpath('S:/Radia/01_GitHub/matlab');

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=2.4e3, ...
    ThermalConductance_W_per_K=12, ...
    SampleTime_s=1e-2, InitialTemperature_K=293.15);

% The power signal can come from Radia VIM/FEM/SIBC/ESIM, a measurement,
% or a position/drive/temperature lookup table.
t = (0:plant.sample_time_s:30).';
P = 2.0e3 * ones(size(t));
Tamb = 293.15 * ones(size(t));
waveform = radia.simulink.simulateIHWaveform(plant, t, P, Tamb);
plot(waveform.time_s, waveform.temperature_K);
xlabel('time (s)'); ylabel('workpiece temperature (K)');
```

`makeIHPlant` produces the exact zero-order-hold discrete matrices for a
lumped thermal envelope. Its inputs are `power_W` and
`ambient_temperature_K`; its outputs are temperature, heat loss, accumulated
input energy, and temperature rate. The model does not pretend that a lumped
thermal equation replaces the Radia field solve. The high-precision result
must feed `power_W`.

For moving workpieces, build the electromagnetic loss table on physical
grids and use the same table in MATLAB or an n-D Lookup Table block:

```matlab
lut = radia.simulink.makeIHPowerLUT( ...
    {position_rad, drive_A, temperature_K}, P_wp_W, ...
    InputNames=["position_rad", "drive_A", "temperature_K"]);
P_now = radia.simulink.evaluateIHPowerLUT(lut, ...
    [position_now, drive_now, temperature_now]);
```

The default LUT policy clips outside the Radia training domain. This avoids
silent extrapolation of SIBC/ESIM data; use `Extrapolation="error"` while
building a new training envelope. Position and speed are kept as explicit
ports in the waveform result so a mechanical solver can drive the same plant.

When Simulink is installed, `radia.simulink.buildIHControlModel` creates a
fixed-step model with the plant and named outputs. `IncludePID=true` adds a
temperature-setpoint PID loop with power saturation. The external power
provider can then be replaced by a Radia LUT, the motor-ROM S-function, or the
native fixed reduced state-space block selected by
`PlantBlock="radia-mex"`.

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
