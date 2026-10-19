# Python GPU: NGSolve precedent and Radia deviations

Status: R2 design evidence, not a built or distributable Python GPU implementation.
The reference is NGSolve commit `736417e0624fa555b5db82281382860ef4a8c9bc`;
a moving branch or a CPU wheel is not evidence of Windows GPU support.

## Verified reference architecture

NGSolve builds a shared `ngscudalib` (or standalone `ngscudalib_local`) and
Python extension `ngsolve._ngscuda`, exposed by `ngsolve.ngscuda`.
At this revision CMake uses a host compiler and `find_package(CUDAToolkit)`;
kernels are compiled at runtime with NVRTC. Direct link targets are
`CUDA::cudart`, `CUDA::nvrtc`, `CUDA::cuda_driver`, NGSolve and the Netgen Python
build interface. This is not the older nvcc-compiled implementation or evidence
that cuBLAS/cuSPARSE are direct dependencies of this exact target.
The Python module initializes CUDA during extension import.
[Build source](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/ngscuda/CMakeLists.txt),
[binding](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/ngscuda/python_ngscuda.cpp).

The Python loader explicitly preloads `libcudart.so.12` and
`libnvJitLink.so.12` from NVIDIA Python package locations. These are Linux paths,
not Windows DLL discovery. Packaging uses scikit-build-core and pins Netgen to
the build ABI. The Linux override enables CUDA; the Windows override does not
establish CUDA-wheel support. The metadata provider inspected here does not
declare NVIDIA runtime requirements, so loader code alone does not establish a
complete dependency installation contract.
[Loader](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/python/ngscuda.py),
[packaging](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/pyproject.toml),
[dependency provider](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/tests/ngsolve_metadata.py).

## Proposed Radia contract and explicit differences

| Aspect | Radia decision / evidence still needed |
|---|---|
| Native layout | Mirror thin Python extension plus domain CUDA library, NGSolve/Netgen ABI matching; no MATLAB Engine required by users. Do not copy a MEX ABI into a Python process. |
| Compilation | Retain MATLAB-bundled CUDA 12.8 and MSVC 14.44. NGSolve's generic CUDAToolkit discovery does not prove this bundle supplies a supported standalone build. A bounded compile/link/import experiment must record actual compiler, headers, import libraries and DLL identities. No external toolkit or compiler override is authorized. |
| NVRTC | NGSolve's runtime-compilation mechanism is a reference, not authorization to add NVRTC or JIT requirements. Prefer sharing Radia's already compiled kernels; if NVRTC is chosen, record that deviation from the existing Radia build and separately qualify its dependency, cache, architecture and runtime-compilation limits. |
| Windows runtime | Expected reference roles are CUDA runtime, NVRTC if used, and NVIDIA driver. Resolve actual versioned Windows DLL names and transitive dependencies from the resulting artifact; do not translate Linux filenames by assumption. Driver remains a host prerequisite, not a bundled DLL. |
| Backend / import | Require explicit GPU selection and lazy GPU loading; CPU import must not initialize CUDA. Do not inherit implicit host fallback or eager CUDA import. Failure raises. |
| Accuracy | Python CPU validation remains independent of shared CUDA kernels and is not a prerequisite for GPU users. Preserve FE ownership, caller-owned TaskManager for CPU and original-equation residual gates. |
| Ownership | Preserve checked lifetime/generation/stream/error contracts; a native module layout alone does not validate array ownership. |

## Licensing and distribution decision

NGSolve declares LGPL-2.1-only. Linking or distributing its components requires
the applicable notices, license and corresponding-source/relink obligations;
using its architecture does not transfer a CUDA redistribution license.
[NGSolve license](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/LICENSE).

NVIDIA CUDA 12.8 terms designate redistributable runtime components, including
Windows CUDA runtime and NVRTC components, subject to the license's conditions.
They do not make the driver or every toolkit file freely redistributable.
Use an explicit component/version/license/notice allowlist and preserve the
project's prohibition on casually bundling third-party DLLs.
[NVIDIA CUDA 12.8 EULA, sections 1.1 and Attachment A](https://docs.nvidia.com/cuda/archive/12.8.0/eula/index.html).

The MATLAB bundle's presence and mexcuda documentation do not demonstrate
permission or support for building and distributing a standalone Python GPU
extension with that bundle. Record the applicable MathWorks/bundled-component
terms and a supported build route before distribution; no categorical license
prohibition or permission is inferred here.
[MathWorks R2026a mexcuda](https://www.mathworks.com/help/releases/R2026a/parallel-computing/mexcuda.html).

**Disposition:** the NGSolve architectural precedent is verified. Compatibility
of the standalone Windows Python build with the mandated MATLAB toolchain, and
the exact bundled-component distribution terms, remain unverified. No external
CUDA exception is assumed or requested by this document. If the required route
cannot pass under the approved toolchain, return the precise missing capability
and a proposed exception for an explicit decision; do not silently install a
toolkit or claim a distributable implementation. F14 evidence is therefore
partially resolved, with build/redistribution closure still required.
