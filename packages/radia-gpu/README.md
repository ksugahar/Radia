# radia-gpu

Independent GPU distribution from the Radia monorepo: MATLAB CUDA MEX
kernels, readable Simulink blocks, and LTspice integration code.
CPU electromagnetic execution belongs to `radia` (Python/pybind11/NGSolve).
`radia-mcp` is the AI operation interface. Neither package is required merely
to locate and use the installed GPU MATLAB files.

This alpha contains FP64 CUDA CSR Jacobi-PCG, a Simulink block for the same
kernel, and the selected standalone LTspice workflow. It also contains experimental FP64 assembled-operator AMS and its versioned
HDF5 exchange. This does not certify production AMS, ICCG, HDiv-MMM, or
migration of the older CPU MEX application models.
There is no CPU solver fallback.

## Install

```powershell
python -m pip install --pre radia-gpu
radia-gpu-doctor --json
radia-gpu-path
```

Add the reported directory to MATLAB's path. MATLAB R2026a, Parallel Computing
Toolbox, an NVIDIA GPU and a compatible NVIDIA driver are external requirements.
Python delivers and verifies the MATLAB payload; it does not execute GPU solves.
The installed CUDA solvers do not require the Python `radia` or `ngsolve`
packages. CSR-PCG consumes caller-provided matrix arrays; AMS consumes an
already-built hierarchy. The current CPU hierarchy-generation route uses
Radia/NGSolve in a separate preparation environment. Optional
`radia-gpu[coupling]` installs Radia for workflows that explicitly need it.
The Windows x64 binary contains sm_75/sm_80 code and compute_80 PTX.
Develop and build `radia-gpu` on 100. Final numerical acceptance compares against
Python/NGSolve on SSH host `NGSolve` (A100). Formal GPU performance
evaluation uses that A100 when it is online; record source/runtime versions,
artifact hashes and machine-readable results. CPU-only mdx execution verifies
the reference calculation, not CUDA execution. RTX 5000 installation acceptance
does not certify A100 scaling.
Cubit mesh export is tested separately on 100 and 102. NGSolve validation
starts by checking and loading an exported `.vol` with its versioned label
contract; it does not require Cubit on the NGSolve host.

Deployment targets are 100 (build and RTX 5000 acceptance) and A100, excluding
102. The 102 host tests `radia` and `radia-mcp` without installing `radia-gpu`.
A100 is not always running: record its acceptance as pending while offline,
then validate the same artifact hash when available. Do not treat an offline
target as a successful test.

Call `radia.gpu.csrPCG` with zero-based int32 CSR row pointers/columns and FP64
GPU values, RHS and positive inverse diagonal for a symmetric positive-definite
system. The caller owns FE assembly
and the GPU arrays. The solver verifies the true residual and raises on invalid
inputs, unsupported state, breakdown or nonconvergence.

`radia.simulink.buildCudaPCGBlock` builds a typed Level-2 MATLAB S-Function
block. Matrix upload occurs at Start; each RHS upload and solution download
is explicit in the block. Simulink is optional for standalone MEX use.

LTspice integration includes MATLAB `radia.ltspice`, Python
`radia_gpu.ltspice`, and the selected sampled-data
Simulink workflow. LTspice executes in its own installed runtime; it is not
a GPU circuit solver. LTspice, MATLAB, Simulink, CUDA Toolkit and NVIDIA
drivers are not redistributed. Install external products separately under
their own licenses. This distribution includes Radia integration code only.

## Build and release

Build on 100 from a clean exact commit with `Build.ps1 -CudaMexOnly
-RequireNativeProvenance` and `Build.ps1 -AmsCudaOnly
-RequireNativeProvenance`, then run `python -m build --wheel` in this directory.
mexcuda compilation and GPU execution never use NGSolve `TaskManager`.
NGSolve preprocessing in this route does not use TaskManager. AMS setup may
use CUDA or OpenMP when validated for the workload; compare numerical parity,
transfer-inclusive timing and memory before selecting either backend. The
current snapshot adapter uses a serial reference. CUDA kernels and cuBLAS
provide GPU solve parallelism. Independent CPU Radia workflows retain their
own threading policy.

NGSolve's [device tutorial](https://docu.ngsolve.org/latest/i-tutorials/unit-5.5-cuda/poisson_cuda.html)
uses `TaskManager` around `ngsolve.gpu` calls, but device operators are provided
by `ngscuda`; [TaskManager itself](https://docu.ngsolve.org/latest/i-tutorials/appendix-taskmanager/taskmanager.html)
manages host threads. Reconsider it for this route only with a demonstrated
CUDA integration benefit and measured parity/performance, rather than assuming
host parallel loops execute on the GPU.

The wheel build requires the matching CUDA build receipt and binary hash.
The independent CPU Radia reference follows its own Python threading policy.
Canonical CUDA and MATLAB sources live in this package, under `native/` and
`matlab/`. Packaging stages those sources; it does not read a root MATLAB tree.
Use only the resulting verified wheel for deployment and PyPI publication.
The initial alpha requires independent review before a stable release.

`radia-gpu.yml` builds and accepts the installed wheel on the build host.
`radia-gpu-release.yml` publishes that successful main CI artifact without
rebuilding, requiring its exact SHA-256 and installed acceptance provenance.
Stable versions additionally require the exact-source independent review and
same-wheel A100/CPU parity record described in
[`validation_test/radia_gpu/releases`](../../validation_test/radia_gpu/releases/README.md).
Configure PyPI Trusted Publishing for project `radia-gpu`, GitHub owner
`ksugahar`, repository `Radia`, workflow `radia-gpu-release.yml`, environment
`pypi`. A new project needs a pending publisher; checking availability does
not reserve its name. Publish only the Windows wheel, not a source archive
that would require an end user's proprietary build toolchain.

GPU AMS is developed separately. HDiv-MMM is a subsequent evaluation target,
after solver/preconditioner parity: preserve NGSolve's FE ownership and validate
the electromagnetic formulation, nonlinear convergence and optimization outcome
against the Python CPU route. Report transfer-inclusive timing and peak memory;
neither successful PCG nor AMS alone certifies HDiv-MMM support.

## Experimental AMS and large numeric exchange

AMS CUDA and its snapshot adapter are package-owned native sources, built by
`tools/build_ams_cuda.py` through the root build dispatcher. The CPU setup
retains NGSolve and Radia ownership; the validation adapter is not installed
as a CPU solver in this wheel. `radia.gpu.amsPCG` accepts a checked numeric
hierarchy and RHS already uploaded to GPU memory, with no CPU fallback.

All AMS file exchange uses `radia.ams-hdf5.v1` plain HDF5. Install optional
`radia-gpu[exchange]` dependencies for the Python writer, and read the result
with `radia.gpu.readAMSExchange`. Values stay FP64 and CSR indices stay
zero-based int32. Shapes and array dimension order are explicit. HDF5 supports
the project's >2 GB exchange requirement; it does not bypass int32 kernel
capacity checks or the GPU working-memory budget. See
[AMS validation](../../validation_test/ams_cuda/README.md) for numerical gates,
transfer-inclusive timing and current limitations. Small-case A100 numerical
and timing evidence is recorded in that validation; it does not certify
A100 scaling or the exact artifact of a future release.

Python circuit converters use `radia-ltspice`; install `radia-gpu[ltspice]`
for the optional schemdraw route. MCP entry points, including
`mcp-server-radia-ltspice`, belong to `radia-mcp` and load the GPU distribution
only when circuit execution is requested. No circuit import alias is left in
the CPU Radia distribution.

The optional `coupling` extra requires Radia 5.3.1 or later with the CPU thermal-operator API. Existing Radia 5.3.0 wheels do not satisfy this contract; until that CPU release is available, use the standalone CUDA PCG, Simulink PCG and LTspice interfaces.
