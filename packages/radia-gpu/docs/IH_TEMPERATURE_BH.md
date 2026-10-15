# Element-temperature B-H / PEEC-BEM(ESIM) research workflow

This feature is a research preview. Its package/release acceptance and its
experimental electromagnetic accuracy are different decisions. The original
tetrahedron uses synthetic amplitude curves, not measured steel data. An
independent full FEM comparison and experimental agreement remain unverified.

## Installation and first run

Use only the release wheel selected after the exact-artifact acceptance gates.
Check its SHA-256 against the release record. A wheel without `ih-build.json`
and `matlab/radia_cuda_ih.mexw64` does not contain this feature. The previous
0.1.3/0.1.4 releases are not temperature-dependent IH distributions.

On Windows x64, use the declared supported MATLAB release, Parallel Computing
Toolbox, Simulink and a compatible NVIDIA GPU/driver. Unzip the verified wheel
to a user-owned folder, then in MATLAB add its `radia_gpu/matlab` folder to the
session path. This route requires no Python solver or Python-generated input.
Do not replace another user's MATLAB path or shared runtime.

Run `matlab/samples/ih-temperature-bh/runTemperatureBHResearch.m`. It creates
a new model in a fresh temporary folder and opens it. Save your research copy
to durable project storage before editing. Keep real material measurements,
real-device inputs and comparison results in authorized private storage.

The model inputs are current (A, peak-phasor convention) and ambient temperature
(K). Its outputs are accepted nodal surface temperature (K), face-integrated
heat (W), and real/imaginary element Zs (ohms). The model workspace owns
`ihConfig`. Inspect and change a copy through Model Explorer or:

```matlab
workspace = get_param(bdroot,'ModelWorkspace');
cfg = evalin(workspace,'ihConfig');
% Edit and check cfg before assigning it back.
assignin(workspace,'ihConfig',cfg);
save_system(bdroot);
```

The sample's 10,000 A ideal-filament excitation is a synthetic numerical choice
to exercise the nonlinear range, not a laboratory operating instruction.
Choose experimental excitation from the declared specimen/instrument ratings
and measurement plan rather than copying this numerical demonstration.

Stop and restart simulation after changing material, geometry or thermal
operators. The parameters are nontunable; geometry changes require reassembly
with `radia.simulink.makeCompletePEECBEMConfig`, not reuse of stale operators.

## Measured material table contract

`temperature_axis_K` is positive and strictly increasing. `bh_knots` is a
common H axis (A/m), starts at zero and strictly increases. `bh_values_T` has
dimensions element x temperature x H knot and contains B in tesla. The native
frontend constructs piecewise-linear B(H) coefficients, then interpolates
them linearly in temperature. Do not pass a signed hysteresis loop as this
single-valued amplitude law.

Each curve starts at B=0, has differential permeability >= mu0 and meets the
tail B=mu0*(H+Ms) continuously. Set `initial_permeability_H_per_m` to its first
segment slope and `saturation_magnetization_A_per_m` to B(last)/mu0-H(last).
Provide positive `resistivity_ohm_m` in element x temperature form. The explicit
`temperature_evaluation` CSR map evaluates each element temperature from the
thermal coefficients; these coefficients are not automatically nodal values.
No temperature extrapolation occurs. Invalid tables and out-of-range candidate
temperatures raise before accepted state is changed.

## Update and recovery semantics

Outputs observes accepted state. Update evaluates local material at the
accepted temperature, converges the complete electromagnetic reaction block
with ESIM, computes coherent loss, and advances the thermal equation. The
accepted electromagnetic, thermal, heat and Zs states commit together after
the checks. Zs belongs to the last electromagnetic update's starting
temperature; it is zero before the first update and after reset.

The electromagnetic-to-thermal exchange stays on device. The optional Zs
observation downloads only requested diagnostics. A failed native update
retains accepted state. The Simulink wrapper terminates its own context on
simulation failure; restarting simulation starts from the configured initial
state. It does not resume an interrupted experiment automatically.

## Supported scope and limits

- Flat P1, closed, consistently oriented, connected genus-zero body surface;
  supplied coil impedance and supported finite-filament coil geometry.
- Bounded GPU block preconditioning (currently <=512 coupled unknowns), bounded
  material tables, and explicit device-memory reserve. No silent CPU fallback.
- Single-valued amplitude B(H,T). No remanence, coercivity, minor-loop history,
  or hysteresis-loop loss. Frequency-dependent data must match the declared run.
- The example uses a surface thermal-layer approximation; it does not resolve
  the body's bulk temperature or certify an actual specimen's heat transport.
- Independent CPU cell/direct-equation comparison shares native BEM geometry
  assembly. Small true residuals do not certify BEM/SIBC physical accuracy.
- No general speedup, process-cold or global peak-memory acceptance follows
  from the small demonstration. Numerical, interface and release gates do not
  replace measured B-H, power/phase and temperature comparison.

Follow `IH_EXPERIMENT.md` for staged measurements, uncertainty and held-out
conditions. Read the release's exact artifact and numerical evidence before
claiming a particular version has passed its distribution acceptance.
