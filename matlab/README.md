# Radia MATLAB MEX

`radia_mex` exposes Radia C++ kernels directly to MATLAB without routing
solver calls through Python. NGSolve still owns TaskManager, mesh loading,
finite-element orientation, curved mappings, and HCurl/HDiv space updates.

## Build

Configure the normal MSVC build with MATLAB enabled, then build the target:

```powershell
cmake -S . -B build-msvc `
  -DRADIA_BUILD_MATLAB_MEX=ON `
  -DMatlab_ROOT_DIR="C:\Program Files\MATLAB\R2026a"
cmake --build build-msvc --config Release --target radia_mex -j 8
```

The build copies `radia_mex.mexw64` into this directory. On Windows it also
places the required oneMKL dispatch DLLs beside the MEX; these local binary
artifacts are intentionally excluded from Git. MATLAB's existing OpenMP
runtime is reused, so `libiomp5md.dll` is not copied beside the MEX.

In MATLAB:

```matlab
addpath("<radia-repository>\matlab")
radia.setup()
info = radia.spaceInfo("model.vol", 6);
```

The regression suite is run with:

```powershell
matlab -batch "addpath('matlab'); results = runtests('tests/matlab'); assert(all([results.Passed]))"
```

`radia.setup` discovers the NGSolve runtime associated with the selected
Python installation and places its DLL directories before older Netgen DLLs.
This setup operation is cached; numerical calls enter C++ directly.

## Current production slice

- NGSolve TaskManager execution and HCurl/HDiv p-space construction
- B-input `radia.EnergyStopMaterial` with checked native-handle lifetime
- complex dense solve and mixed-Galerkin Schur complement
- skin impedance and SIBC tail/termination kernels
- CLN Lanczos reduction, tridiagonal construction, impedance sweeps,
  loop-star coupling transforms, and ACA/SVD star compression
- PEEC filament geometry with a stateful HACApK inductance manager
- monopole HDiv charge-Gram H-matrix with regular, transpose, and symmetric matvecs
- EVRS/T-method algebra and tetrahedral HCurl reduced Gram integration
- complex filament/surface Biot-Savart kernels and P1/P2 SLDL Galerkin assembly
- stateful HACApK scalar BEM construction, build, matvec, and statistics
- first legacy/pybind overlap slice: tetrahedron, hexahedron, wedge, pyramid,
  container, linear/nonlinear material, solve, field, and object utilities

The MEX boundary uses MATLAB column-major arrays externally and converts them
to Radia's row-major `[target][source]` convention internally.

## Simulink and induction heating

The MATLAB package also contains the Simulink-facing control layer:

```matlab
addpath("<radia-repository>\matlab")
plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=2.4e3, ...
    ThermalConductance_W_per_K=12, SampleTime_s=1e-2);
```

The plant accepts Radia-computed `power_W` and ambient temperature and
returns temperature, heat loss, input energy, and temperature rate. A
position/drive/temperature power LUT is created with
`radia.simulink.makeIHPowerLUT`, so motion coupling remains an explicit
mechanical input rather than a hidden thermal approximation. Use
`radia.simulink.buildIHControlModel` when Simulink is installed, and
`radia.simulink.optimizeIH` to wrap either `sim`/`parsim` or a fast MATLAB
waveform objective for controller and process optimization.

## Binding policy

The pybind surface contains 174 binding entries with 165 unique names. The
legacy `radentry` C ABI contains 118 functions, and 62 names directly overlap
the pybind surface. These families should be ported as follows:

| Family | MATLAB representation |
|---|---|
| Legacy object, transform, material, solve, and field calls | Thin numeric wrappers over `radentry` |
| Pure array kernels such as SIBC, Schur, CLN, EVRS/T, Biot-Savart, and Galerkin assembly | Direct MEX commands |
| Stateful EnergyStop, HACApK, and field evaluators | Typed `uint64` handles with `mexLock` lifetime |
| NGSolve mesh and finite-element work | Mesh paths and numeric contracts; do not emulate Python NGSolve objects |
| Python callbacks and `CoefficientFunction` adapters | Keep in Python; expose the underlying numeric method instead |
| Internal self-tests and underscore-only diagnostics | Test surface, not public MATLAB API |

This is deliberately not a mechanical one-MEX-file-per-`.def` translation.
One gateway avoids duplicating the Radia/NGSolve runtime and gives every
stateful class the same checked lifetime and exception boundary.
