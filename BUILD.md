# Radia Build Guide

Radia's production native build targets Windows x64 and Python 3.12. The build
uses MSVC, the exact NGSolve/Netgen versions pinned in `pyproject.toml`, and the
Intel MKL files installed into the selected Python environment by `mkl-devel`.

## Quick Start

Create or activate a Python 3.12 virtual environment, then install the native
build inputs into that environment:

```powershell
python -m pip install --upgrade pip
python -m pip install ngsolve==6.2.2606 netgen-mesher==6.2.2606 `
  mkl-devel pybind11==3.0.2 ninja cmake build wheel
pwsh -NoProfile -ExecutionPolicy Bypass -File .\Build.ps1 -Rebuild -Test
```

`Build.ps1` uses the `python` selected on `PATH`. Do not point CMake at one
Python installation while importing the result from another.

## Build Targets

| Command | Purpose |
|---|---|
| `.\Build.ps1` | Incremental native build |
| `.\Build.ps1 -Rebuild` | Clean configure and build |
| `.\Build.ps1 -RadiaOnly` | Main `_radia_pybind` extension only |
| `.\Build.ps1 -AxiFemOnly` | Fast axifem C++ iteration |
| `.\Build.ps1 -MatlabMexOnly` | Standalone Radia MEX gateway and native S-Function targets |
| `.\Build.ps1 -OptunaMexOnly` | Lightweight standalone `optuna_mex` target |
| `.\Build.ps1 -Test` | Build, source import check, and focused pytest |

The four `*Only` switches are mutually exclusive. `-InstallToSitePackages` is
for an intentional local installation test; editable development normally
loads the extensions from this checkout.

Coreform Cubit is optional for the Radia core. If installed, `Build.ps1`
discovers it through `tools/find_cubit.ps1`. The independent
`cubit-mesh-export` package owns the tracked `.ccm` plugin and its deployment.

## Toolchain Contract

- Visual Studio Build Tools with the x64 MSVC toolchain are required.
- Python 3.12 is the supported wheel ABI.
- NGSolve and Netgen must match the exact pins in `pyproject.toml`.
- `pybind11==3.0.2` is the checked adapter toolchain version.
- `mkl-devel>=2026,<2027` provides headers, import libraries, and runtime DLLs.
- MATLAB is required only for MEX and Simulink native targets.

The selected Python environment's `Library` directory is the first MKL search
location. An explicitly set `MKLROOT` is a fallback for a controlled external
installation. There is no hard-coded oneAPI installation path.

NGSolve 6.2.2606 uses its own OpenBLAS distribution. Radia's HACApK, PARDISO,
and dense native kernels intentionally use MKL. These are separate native
dependencies; Radia does not replace or bundle the user's NumPy implementation.

## MATLAB and MEX

Reusable numerical kernels are exported as independently callable MEX
functions. Level-2 MATLAB S-Functions provide readable Simulink ports,
lifecycle, sample-time behavior, and diagnostics on top of those kernels.
A complete C/C++ MEX S-Function is reserved for a measured real-time,
zero-copy, generated-code, or native-lifecycle requirement.

Use checked `uint64` handles for long-lived native objects. The registry must
validate type, generation, ownership, and liveness; public APIs never expose a
raw pointer. Test the standalone MEX ABI independently before the Simulink
lifecycle tests.

## Wheel Build

Build a local wheel without publishing it:

```powershell
pwsh -NoProfile -ExecutionPolicy Bypass -File .\Build_Wheel.ps1 -DryRun
```

Wheel creation fails on the first packaging error and preserves the log tail.
Retries belong at the infrastructure/job level only after the failure has been
classified as transient.

Native outputs are not committed or copied directly to another machine. CI and
release workflows retain artifacts built from the exact checked commit; wheels
and versioned GitHub Release packages are the distribution boundary.

## CI Build

Routine native CI runs on the self-hosted `mdx` runner. Each job creates an
isolated virtual environment, installs the repository pins and `mkl-devel`, and
builds from the checked commit. LAB and 100号機 remain editable development
hosts and are not routine CI runners.

Normal source CI runs only compact regression tests affected by changed paths.
Solver benchmarks, machine comparisons, convergence studies, and publication
evidence belong to `validation_test/` and run explicitly on mdx or hibino.

## Troubleshooting

### Wrong Python ABI

```powershell
python --version
python -c "import sys,struct; print(sys.executable); print(struct.calcsize('P')*8)"
```

Confirm Python 3.12 x64, remove `build-msvc`, and rebuild with that interpreter
active.

### NGSolve or Netgen mismatch

```powershell
python -c "import importlib.metadata as m, ngsolve; print(ngsolve.__version__); print(m.version('netgen-mesher'))"
```

Both versions must match `pyproject.toml`. Rebuild every extension that links
their C++ ABI after either pin changes.

### MKL not found

```powershell
python -m pip install --upgrade --force-reinstall "mkl-devel>=2026,<2027"
python -c "import pathlib,sys; print(pathlib.Path(sys.prefix) / 'Library')"
```

That `Library` directory must contain `lib\mkl_rt.lib` and
`bin\mkl_rt.3.dll`. Set `MKLROOT` only when deliberately using another checked
MKL installation.

### Native import failure

Inspect the imported path and dependency versions before copying files:

```powershell
python -c "import radia; print(radia.__file__)"
python -m pip show radia ngsolve netgen-mesher mkl
```

For an editable-source mismatch, use the `verify-deploy` skill. For a release,
use `release-quad`; do not repair an ABI mismatch with a direct `.pyd` drop.
