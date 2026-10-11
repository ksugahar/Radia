# radia-gpu

Independent GPU distribution from the Radia monorepo: MATLAB CUDA MEX
kernels, readable Simulink blocks, and LTspice integration code.
CPU electromagnetic execution belongs to `radia` (Python/pybind11/NGSolve).
`radia-mcp` is the AI operation interface. Neither package is required merely
to locate and use the installed GPU MATLAB files.

This alpha contains FP64 CUDA CSR Jacobi-PCG, a Simulink block for the same
kernel, and the selected standalone LTspice workflow. It does not claim
GPU AMS, ICCG, HDiv-MMM, or migration of the older CPU MEX application models.
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
The Windows x64 binary contains sm_75/sm_80 code and compute_80 PTX.
RTX 5000 acceptance does not certify A100 scaling.

Deployment targets are 100 (build and RTX 5000 acceptance) and A100, excluding
102. A100 is not always running: record its acceptance as pending while offline,
then validate the same artifact hash when available. Do not treat an offline
target as a successful test.

Call `radia.gpu.csrPCG` with zero-based int32 CSR row pointers/columns and FP64
GPU values, RHS and positive inverse diagonal. The caller owns FE assembly
and the GPU arrays. The solver verifies the true residual and raises on invalid
inputs, unsupported state, breakdown or nonconvergence.

`radia.simulink.buildCudaPCGBlock` builds a typed Level-2 MATLAB S-Function
block. Matrix upload occurs at Start; each RHS upload and solution download
is explicit in the block. Simulink is optional for standalone MEX use.

LTspice integration includes `radia.ltspice` and the selected sampled-data
Simulink workflow. LTspice executes in its own installed runtime; it is not
a GPU circuit solver. LTspice, MATLAB, Simulink, CUDA Toolkit and NVIDIA
drivers are not redistributed. Install external products separately under
their own licenses. This distribution includes Radia integration code only.

## Build and release

Build on 100 from a clean exact commit with `Build.ps1 -CudaMexOnly
-RequireNativeProvenance`, then run `python -m build --wheel` in this directory.
The wheel build requires the matching CUDA build receipt and binary hash.
It stages canonical MATLAB sources rather than maintaining a package-local fork.
Use only the resulting verified wheel for deployment and PyPI publication.
The initial alpha requires independent review before a stable release.

`radia-gpu.yml` builds and accepts the installed wheel on the build host.
`radia-gpu-release.yml` publishes that successful main CI artifact without
rebuilding, requiring its exact SHA-256 and installed acceptance provenance.
The initial publisher rejects stable versions pending independent review.
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
