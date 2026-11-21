---
name: linux-gpu-build
description: Prepare and verify the Radia _hdiv_cuda extension from source on Linux using the documented aarch64 sm_100 development route.
---

# Build the Radia HDiv CUDA extension on Linux

Read the [Linux/aarch64 CUDA procedure](../../../docs/compilation/python-hdiv-gpu.md#isolated-linuxaarch64-build)
and the [GPU builder](../../../tools/build_python_gpu.py) at the selected
revision. Use an isolated clean checkout and job environment. Build the CPU
modules at the same commit using the [Linux CPU build guide](../../../docs/compilation/linux-hpc.md).

## Prepare and build

- Use nvcc from an explicit CUDA 12.8 toolkit, an explicit GCC 11.5 host
  compiler, and sm_100 via the builder's architecture argument. Use the
  procedure's CMake and Ninja prerequisites. Unset the compiler environment
  variables listed in the builder's toolchain check, including NVCC_CCBIN,
  CFLAGS and CUDAFLAGS. The builder rejects any nonempty value.
- Resolve the toolkit and compiler paths, scheduler allocation, and scratch
  directory on the target system. Keep build output outside the checkout.
  Match NGSolve/Netgen and the pip CUDA runtime to the selected checkout's
  package metadata. The toolkit provides compile tools and headers. The
  pinned pip libraries provide the runtime.
- Build only the documented component with --components _hdiv_cuda. Use an
  approved resource allocation and admission JSON with the exact source_sha,
  a nonempty owner and gpu_build_approved set to true. Supply the independently
  recorded SHA-256 of that JSON. Keep build outputs private except the module
  and its minimal public receipt. Logs, commands and staging data can contain
  site-specific paths.

## Verify and report

Record the source and binary identities, toolchain, target architecture and
checks actually run. A build only establishes artifact creation. Follow the
guide's checked-mesh requirements before attempting the GPU smoke test. Report
runtime import, kernel execution, CPU/GPU agreement and performance separately.
Recheck each target environment and source revision. The documented route is
experimental. The toolchain check does not establish host architecture support. Linux wheel
publication and qualification of other GPU architectures are outside this skill.
