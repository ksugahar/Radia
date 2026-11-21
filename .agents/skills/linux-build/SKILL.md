---
name: linux-build
description: Build and verify Radia CPU extensions from source on Linux using tools/build_linux.py. Use for x86_64 Intel LLVM with oneMKL or aarch64 GCC with NVPL.
---

# Build Radia from source on Linux

Read the [Linux HPC build guide](../../../docs/compilation/linux-hpc.md) and
check [build_linux.py](../../../tools/build_linux.py) at the selected revision.
For a user installation, start from a release tag. For development, record the
exact source revision being tested.

## Prepare the target environment

- Select the recipe by architecture. The guide specifies Intel LLVM and
  oneMKL on x86_64, and GCC 11.5 with LP64 sequential NVPL on aarch64.
- Resolve compiler, library and module locations on the target host. Site
  module names, library paths and scheduler allocations depend on the site.
  Use a fresh job environment and a separate build directory. Set the job
  count within the allocated resources.
- Match NGSolve and Netgen to the selected checkout's pyproject.toml. On x86_64
  the guide assumes an already reviewed MKL dependency environment. The
  repository has no hash-locked recipe for preparing it. On aarch64, prepare
  NGSolve/Netgen with OCC, NumPy and SciPy against the same NVPL provider
  before building Radia. Stop if the required
  provider or ABI cannot be established. Follow the guide's architecture
  commands and inspect the relevant flags with python tools/build_linux.py --help.

## Verify and report

Run the build on the target, import the rebuilt extensions in a fresh Python
process, and run the relevant numerical checks from the guide. Record source
revision, architecture, dependency versions, compiler, provider and check
results. Report skipped checks explicitly. Keep linux-build.json and any
.so.build.json sidecars private because they can include site-specific paths
(for example, in ldd output). The sidecars require --require-clean-source.

The earlier trial described in the guide does not establish support for the
versions pinned in the selected checkout's pyproject.toml or for another site stack. Verify each target environment before describing it as supported.
This skill covers CPU source builds. GPU builds and Linux wheel publication
are outside this procedure.
