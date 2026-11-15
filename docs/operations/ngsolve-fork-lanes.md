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
runtime and NVRTC bundled with MATLAB (`cudart`, `nvrtc`) and the NVIDIA
driver library (`nvcuda`), which comes with the driver, not with MATLAB.
Kernels are compiled at run time by NVRTC. MATLAB's `mexcuda` command is not
used.

A MEX initialises the GPU backend itself, without Python:
`ngs_cuda::InitCUDA`, `ngs_cuda::InitCudaDevice` and `ngla::InitCuLinalg`.

## TaskManager

Radia-owned MEX kernels never include, link or use TaskManager directly.

The four fork lanes above are the only exception. Their CPU lanes need
NGSolve's TaskManager, so `libngsolve`, `ngscudalib` and every MEX linked to
them carry it as an inherited dependency (TaskManager symbols come from
`ngcore`). This does not allow a Radia-owned kernel to call TaskManager.

## BLAS and LAPACK

The MEX SDK links MATLAB's own MKL through `libmwblas` and `libmwlapack`,
which use 64-bit integers (ILP64). The SDK is built with `LAPACK_ILP64=ON`
and explicit `LAPACK_LIBRARIES`; it bundles no separate MKL or OpenBLAS.
The libraries are exported by name, so the consumer supplies MATLAB's
`extern/lib` directory, as `mex` does. The imported `ngsolve` target carries
the `MKL_ILP64` definition, and consumers must compile with it. A consumer that
drops it fails to link the tested C++ wrappers, whose mangled names encode the
integer width; a direct call to an unmangled LAPACK routine with a 32-bit
prototype would still link and fail only at run time.

## Evidence

Build receipts, hashes and consumer checks stay in private storage. A GPU
result counts only after a run on the validation host; a build or link check
is not GPU evidence.
