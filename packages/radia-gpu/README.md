# radia-gpu

Independent GPU distribution from the Radia monorepo: MATLAB CUDA MEX
kernels, readable Simulink blocks, and LTspice integration code.
CPU electromagnetic execution belongs to `radia` (Python/pybind11/NGSolve).
`radia-mcp` is the AI operation interface. Neither package is required merely
to locate and use the installed GPU MATLAB files.

Version 0.1.0 contains FP64 CUDA CSR Jacobi-PCG, a Simulink block for the same
kernel, and the selected standalone LTspice workflow. It also contains experimental FP64 assembled-operator AMS and its versioned
HDF5 exchange. This does not certify production AMS, ICCG, HDiv-MMM, or
migration of the older CPU MEX application models.
There is no CPU solver fallback.

## Install

```powershell
python -m pip install radia-gpu
radia-gpu-doctor --json
radia-gpu-path
```

Add the reported directory to MATLAB's path. MATLAB R2026a, Parallel Computing
Toolbox, an NVIDIA GPU and a compatible NVIDIA driver are external requirements.
Python delivers and verifies the MATLAB payload; it does not execute GPU solves.
CSR-PCG consumes caller-provided matrix arrays and requires no Python Radia or
NGSolve runtime. The experimental native AMS candidate constructs its hierarchy
inside mexcuda and requires compatible NGSolve/Netgen native libraries, including
the permitted Python DLL dependency. GPU use requires no Python adapter,
interpreter initialization, or Python-generated solver snapshot. Native-library
deployment/packaging acceptance is separate from numerical candidate validation.
Optional
`radia-gpu[coupling]` installs Radia for workflows that explicitly need it.
The Windows x64 binary contains sm_75/sm_80 code and compute_80 PTX.
Develop and build `radia-gpu` on 100. Required release acceptance compares each
published CUDA solver against Python/NGSolve on 100 (Quadro RTX 5000), and
requires an independent review of the exact source. Record source/runtime
versions, artifact hashes and machine-readable results. A100 validation and
scaling are optional additional evidence, not publication requirements.
Formal timing remains subject to the idle-mdx host policy; 100 timings are diagnostic. CPU-only mdx execution verifies
the reference calculation, not CUDA execution. RTX 5000 installation acceptance
does not certify A100 scaling.
Cubit mesh export is tested separately on 100 and 102. NGSolve validation
starts by checking and loading an exported `.vol` with its versioned label
contract; it does not require Cubit on the NGSolve host.

Deployment targets are 100 (build and RTX 5000 acceptance) and A100, excluding
102. The 102 host tests `radia` and `radia-mcp` without installing `radia-gpu`.
A100 is optional and never blocks development or publication. It is not always running: record its acceptance as pending while offline,
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
same-wheel 100 RTX 5000/CPU parity record described in
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

`radia.gpu.prepareAMS(model)` constructs the shared native AMS directly in
mexcuda host code and uploads its numeric hierarchy by default. GPU users
supply numeric FEM operators, discrete gradient, coordinates and free-DOF mask;
no Python adapter or Python-generated solver snapshot is required. The legacy
MEX constructor is reused, with NGSolve sparsecholesky kept for coarse setup.
The default host setup uses one thread; `SetupThreads` selects OpenMP row/chunk
work. No TaskManager region is activated. Indirect NGSolve header/DLL and Python
DLL dependencies are permitted; no Python interpreter initialization is needed.

Use `validation_test/ams_cuda/native_setup_validation.py` for independent Python
CPU parity and MATLAB-owned serial/OpenMP setup comparison, including transfer,
solve and download times. The Python snapshot adapter remains a validation
oracle only. Build-host timings are diagnostics; formal performance acceptance
requires an idle GPU-equipped timing host.

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

### CUDA Graph experiment (2026-10-12)

`radia.gpu.amsPCG(hierarchy,rhs,tol,maxIterations,false,true)` enables
experimental CUDA Graph replay of the existing AMS cycle. The default remains
direct CUDA launches. Graph capture, instantiation and cleanup belong to each
solve; the graph does not outlive its arrays or workspace. No CPU fallback is
introduced. NGSolve sparsecholesky setup remains unchanged.

The diagnostic direct/Graph ABBA experiment is recorded in
`validation_test/ams_cuda/results/cuda-graph-dfeec57af.json`. At 702073 DOFs,
transfer-inclusive solve time changed from 0.40324 to 0.38472 seconds, while
setup-inclusive time changed from 4.74394 to 4.74596 seconds. This does not
establish an end-to-end improvement. Host hierarchy construction remains the
next acceleration target. CPU parity was verified independently; CPU timing
was single-shot and does not establish repeated CPU/GPU speedup. These RTX
5000 timings are diagnostic, not idle-mdx performance acceptance.

### CUDA AMS matrix generation

`radia.gpu.prepareAMS(model,SetupThreads=4)` now generates the auxiliary and
AMG coarse Galerkin matrices on CUDA. `SetupOnGPU=false` explicitly selects
the previous host matrix-generation path for validation. Both paths retain
NGSolve AMG coarsening/interpolation and serial `sparsecholesky` coarse
factorization. Python is not a prerequisite and TaskManager is not activated.

CUDA computes sorted CSR `AP` followed by `Pt(AP)`, preserving the current
NGSolve sparse-product accumulation order and separate FP64 multiply/add
rounding. Each row has one owner; floating-point atomics are not used. The
matrix cache and all scratch allocations belong to one setup command and are
freed on normal/error exit. Allocation checks reserve at least 2 GiB or 20%
of device memory. Unsupported capacity or device state raises an error.

The shared CPU code leaves its optional accelerator callback disabled.
`graph_validation.py --comparison generation --rounds 2` compares both
matrix-generation paths in ABBA order, excluding backend warmups. Full time
includes host setup, device generation transfers, hierarchy upload, solve and
solution download. Setup memory records the exact maximum concurrently owned
CUDA allocation bytes; solve memory records a live allocation delta, not a
whole-process peak. RTX 5000/100 measurements remain diagnostic.
