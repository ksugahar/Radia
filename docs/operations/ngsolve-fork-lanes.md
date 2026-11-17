# NGSolve fork lanes

User decision, 2026-11-16. Radia depends on the Sugahara forks of Netgen and
NGSolve (`ksugahar/netgen`, `ksugahar/ngsolve`), not on the official wheels.

## Maintenance model

The forks are not a development line. Each maintenance round:

1. merges the latest upstream `master` into the fork's `master`, keeping the
   fork as close to upstream as possible;
2. adds only bug fixes that Radia needs, as small labeled commits; general
   fixes may be offered upstream only after the user approves the exact text;
3. builds the pybind11 wheels and the MEX SDK (CPU and GPU) and runs the
   fork checks;
4. tags netgen and ngsolve with the same tag on `master` and publishes the
   wheels and SDK on the forks' GitHub Releases;
5. pins Radia's metadata to the exact published wheel versions of `ngsolve`
   and `netgen-mesher` from that release (a PEP 440 version such as
   `6.2.2608.post100.dev0+radia.1`, not the tag string), then rebuilds and
   revalidates Radia against them.

The fork has one long-lived branch, `master`. Radia's own binaries link the
NGSolve C++ ABI, so every pin move rebuilds Radia; upstream changes are not
assumed to be ABI compatible. The ngsolve fork's `external_dependencies/netgen`
submodule points at the netgen fork commit with the same tag. Users install
with `pip install --find-links <wheelhouse>`, where the wheelhouse is a page
or directory that lists the actual wheel files of that release; an ordinary
release page is not assumed to work as one. PyPI carries no fork wheels.
Every Radia CI lane that installs NGSolve receives the same fork wheels and
checks the installed distribution versions and native ABI.

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

## GPU development order

Establish the CPU reference through pybind11 and MATLAB on the build host,
develop the Python GPU lane on Flow2, then develop the MATLAB GPU lane through
plain MEX on the validation host. Follow the shared
[CUDA development policy](../induction_heating/CUDA_VALIDATION.md#development-order)
for parity and host-specific evidence. Keep the four supported fork lanes and
the maintenance release checks for each applicable artifact.

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
integer width. A direct call to an unmangled LAPACK routine with a 32-bit
prototype can link despite the ABI mismatch and may produce silent wrong
results or memory corruption.

## Evidence

Build receipts, hashes and consumer checks stay in private storage. A GPU
result is scoped to the binding, source and GPU that actually ran it. Flow2
qualifies Python development cases; the validation host checks the Windows
MATLAB lane. Record build and link checks as separate evidence.
