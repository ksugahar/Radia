# MATLAB and Python API Parity

Radia uses two complementary parity checks. The low-level binding audit proves
that the declared pybind11 numerical surface has a MEX command. The module
audit proves that every production Python module has a named MATLAB route,
including explicit Python fallback where a native port is not practical.

## Current inventory

The 2026-07-21 audit classifies all 209 Python files under `src/radia`:

| MATLAB route | Modules |
|---|---:|
| Native MEX | 5 |
| Native MATLAB | 4 |
| Explicit Python fallback | 157 |
| Private or not applicable | 43 |

The MEX gateway currently exposes 324 commands and covers all 94 public
top-level pybind11 numerical names, 27 internal numerical kernels, and 119
stateful class members declared by the binding audit. This broad low-level
coverage is distinct from high-level module parity: many Python modules compose
NGSolve objects, CAD topology, application settings, and artifact workflows.

Run both checks after changing a Python, pybind11, MATLAB, or MEX API:

```powershell
python tools/check_matlab_python_parity.py
python -m pytest packages/radia-mcp/tests/test_matlab_radia_mex_contract.py -q
```

The checked source of truth for module classification is
`matlab/python_api_parity_manifest.json`.

## Remaining native gaps

The largest MATLAB-native gaps are:

| Family | Current MATLAB route | Native promotion target |
|---|---|---|
| Acoustic CQ-BEM and FSI | `radia.python.acoustic` | Keep NGSolve FE/BEM plumbing in Python; promote reusable numeric kernels only |
| Axisymmetric FEM (`axifem.pyd`) | Native Q1 element matrices plus Python/NGSolve composition | Promote Q2/heat arrays, then checked sparse assembly |
| VIM, ESIM, and IH composition | Family fallback plus focused MEX kernels | Stable numeric/handle contracts around validated headless workflows |
| CoilBuilder and CAD topology | Python fallback | File/value contracts where native MATLAB adds operational value |
| Kelvin, DtN, and open-boundary composition | Family fallback plus native primitives | Promote repeated array kernels, retain NGSolve-owned object work in Python |
| High-level BEM, PEEC, SIBC, motor, and MagLev design | Family fallback plus focused MEX kernels | Application-specific MATLAB APIs over shared artifacts |
| Pure-Python analytical formula collections | Python fallback | Native `.m` functions when they are independently testable and frequently used |

Fallback is deliberate, visible, and testable. `radia.internal.callPython`
requires Python 3.12 through MATLAB `pyenv` in InProcess mode, returns the
loaded executable and DLL in its metadata, and labels the backend
`python-fallback`. A fallback may run during initialization, explicit update,
artifact generation, or a batch solve. It may not run at each Simulink sample.

## Native promotion order

The checked backlog in `matlab/python_api_parity_manifest.json` fixes the next
promotion sequence and the owning tests:

1. axifem Q2/heat element arrays and checked numeric assembly
2. VIM/ESIM/IH artifact-to-handle and ROM initialization
3. high-level BEM/PEEC/SIBC numeric factories and artifact loaders
4. Kelvin/DtN factor arrays and boundary-operator handles
5. batched acoustic CQ transfer and FSI interface-coupling kernels
6. Motor/MagLev model and operating-point artifact initialization
7. coil centerline/cross-section/terminal value and file contracts

This order targets repeated numerical work and Simulink execution first.
NGSolve spaces/forms and build123d/OCC object identity remain with their
owning ecosystems unless a stable numeric, sparse-matrix, handle, or file
boundary provides a measurable operational benefit.

## axifem promotion

`radia.axifem.q1_magnetic_element_matrices` / `q2_magnetic_element_matrices`
and their MATLAB adapters now expose the same shared C++ Q1/Q2 Henrotte
kernels through pybind11 and MEX. The returned 4-by-4 or 9-by-9 stiffness and
sigma-mass matrices use nodal `A_phi` values. The production NGSolve Q1/Q2
bilinear-form integrators call that same source, so these are not separate
MATLAB formulas.

The focused gates compare the native array API with the independent Python
formula, the NGSolve BFI assembly, and the MATLAB adapter, including
axis-touching elements and fail-loud invalid geometry. The Q2 matrices also
feed `makeAxiEddyElementModel`, whose exact-ZOH state model runs through the
native Simulink state-space MEX with no Python-per-step path. NGSolve continues
to own the mesh, FESpace, CoefficientFunction, BilinearForm, and global
assembly.

## Acoustic promotion

The acoustic analytic sphere models and CQ primitives now use one C++ source
for Python and MATLAB. The native surface includes:

- soft, rigid, fluid, and elastic sphere scattering
- soft-sphere scattering at a complex wavenumber
- BDF1/BDF2 delta evaluation
- convolution-quadrature grid construction

Python exposes these through pybind11 and MATLAB through seven `radia_mex`
commands plus `+radia/+acoustic` wrappers. Focused tests compare both adapters
with independent formulas and with each other. The full CQ-BEM and FSI solvers
remain Python/NGSolve-owned because their value lies in NGSolve's finite-element
and boundary-element object model, not in duplicating that plumbing.

## Measured performance

The committed validation driver compares the same shared C++ kernel on an idle
`mdx` compute host. It uses 20,000 scattering points, 28 partial-wave terms,
five warmups, and 31 measured repetitions. Process startup is reported
separately from warm calls.

| Measurement | Python/pybind11 | MATLAB/MEX |
|---|---:|---:|
| First scattering call | 41.703 ms | 119.784 ms |
| Warm scattering median | 31.684 ms | 22.319 ms |
| BDF transfer median | 0.307 ms | 0.392 ms |
| Cold process end to end | 2.904 s | 10.833 s |

The scattering checksums agree to a relative error of `8.1e-16`. In this run,
MEX delivered about 1.42 times the warm scattering throughput, while pybind11
had about 1.28 times lower latency for the transfer-dominated BDF call. The
result supports measuring each workload: neither adapter is inherently faster.
Full metadata and unrounded values are stored in
`validation_test/acoustics/acoustic_mex_python_benchmark.json`.
