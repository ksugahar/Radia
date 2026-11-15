# NGSolve fork lanes

User decision, 2026-11-16. The Netgen/NGSolve forks are validation lanes for
Radia's GPU and MATLAB work, not Radia source. They follow upstream NGSolve
and carry small, reviewable fixes.

## Supported lanes

| Lane | CPU | GPU | Python |
|---|---|---|---|
| pybind11 | yes | yes | one wheel per Python version |
| MEX | yes | yes | none: no Python version, DLL or import |

- **pybind11** is the Python package of the fork, built per Python ABI.
- **MEX** uses a Python-free C++ SDK: netgen and ngsolve built with
  `USE_PYTHON=OFF`, installed with relocatable CMake config files, and
  consumed through `find_package(NGSolve CONFIG)`. No SDK DLL and no MEX file
  may import a `python*.dll`.

## MEX + CUDA build

The MEX GPU lane is a **MEX + CUDA build**. MSVC builds the MEX through
`matlab_add_mex`. The MEX links the SDK's `ngscudalib`, which links the CUDA
runtime bundled with MATLAB (`cudart`, `nvrtc`, the driver). Kernels are
compiled at run time by nvrtc. MATLAB's `mexcuda` command is not used.

A MEX initialises the GPU backend itself, without Python:
`ngs_cuda::InitCUDA`, `ngs_cuda::InitCudaDevice` and `ngla::InitCuLinalg`.

## TaskManager

The CPU lanes use NGSolve's TaskManager, so the SDK links `libngsolve` with
TaskManager by design. Radia-owned MEX kernels still never use TaskManager;
this exception covers only routes that call NGSolve through the fork SDK.

## BLAS and LAPACK

The MEX SDK links MATLAB's own MKL through `libmwblas` and `libmwlapack`,
which use 64-bit integers (ILP64). The SDK is built with `LAPACK_ILP64=ON`
and explicit `LAPACK_LIBRARIES`; it bundles no separate MKL or OpenBLAS.
The libraries are exported by name, so the consumer supplies MATLAB's
`extern/lib` directory, as `mex` does. A consumer that omits the ILP64
definition fails at link time rather than at run time.

## Evidence

Build receipts, hashes and consumer checks stay in private storage. A GPU
result counts only after a run on the validation host; a build or link check
is not GPU evidence.
