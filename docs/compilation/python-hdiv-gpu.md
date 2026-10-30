# Experimental Python HDiv-MMM CUDA build preparation

`radia.hdiv_gpu.PreparedHDivNonlinear` owns the Python GPU workflow. The shared
numerical sources live once in `src/core/gpu/hdiv`; MATLAB/Simulink gateways
consume those sources through mexcuda. Python never imports `radia_simulink`.
The MATLAB host publisher only converts Radia-owned arrays to its ABI layout.

This is a development build route, not a release or numerical acceptance.
HDiv-MMM development changes require a new build of the exact source commit;
previous binaries and acceptance records do not qualify a new revision.
No performance testing is part of this preparation.

## Isolated Linux/aarch64 build

Use an isolated clean checkout and job environment with the package metadata's
exact NGSolve/Netgen pins and CUDA 12.8 pip runtime dependencies. Build the CPU
modules for that same commit first with `tools/build_linux.py`. The optional
GPU builder uses the installed SDK's CMake package and ABI for its **host-only**
material staging; nvcc compiles the device translation unit separately, without
NGSolve TaskManager or MEX headers.

Only an explicitly admitted resource frame may execute these commands. The
admission JSON names the exact `source_sha`, nonempty `owner`, and
`gpu_build_approved: true`; supply its independently recorded SHA256.

```sh
python tools/build_python_gpu.py --components _hdiv_cuda \
  --cuda-root "$CUDA_ROOT" --host-compiler "$HOST_CXX" --architecture 100 \
  --output "$JOB_ROOT/hdiv-gpu" --admission "$ADMISSION" \
  --admission-sha256 "$ADMISSION_SHA256"
```

For the currently prepared aarch64 route the reviewed host compiler is GCC
11.5, nvcc is CUDA 12.8, and the target is sm_100. Use the installed Netgen
wheel's headers rather than another pybind11 SDK. Missing or ambiguous SDK
metadata fails. CMake and Ninja are required; nothing installs or downloads
packages. Runtime libraries are exact-pin pip-owned files; the toolkit supplies
compile headers/tools, not a runtime fallback.

Outputs include the ABI-suffixed `_hdiv_cuda` module, a minimal `.build.json`,
and private command/SDK receipts. Stage the module and public receipt under the
isolated checkout's `src/radia` for an admitted smoke. CPU and GPU provenance
must identify the same clean source and hash-correct binary. No shared editable
installation is replaced. A build receipt is not evidence of GPU execution.

## Before a later basic GPU smoke

Supply a checked `.vol` and its label contract. The existing pure-Python
`cubit_mesh_export.check` checker must be available in that job for strict label
checking; it does not launch Cubit. No checker bypass is introduced here.
NGSolve owns FE preparation under the caller's `TaskManager`. Select
`backend="gpu"` explicitly; unsupported state or CUDA failures raise without a
CPU fallback. Reuse the existing small component input, independent CPU oracle
and unchanged thresholds. Kernel execution, CPU/GPU numerical parity, full GPU,
whole-package HDiv-MMM, topology optimization and performance remain separate
unverified steps for this source revision.

See the [targeted API ownership audit](../api_inventory/API_SURFACE_2026-10-29_HDIV_GPU.md).
