# radia-simulink

Breaking source migration, not a published or accepted new-name wheel.
Historical 0.1.x releases and measurements below belong to the retired
distribution and do not certify this distribution. CPU mexcuda execution and direct Python GPU
bindings are separate implementation/acceptance work, not completed by renaming.
The experimental CSR adapter now has explicit CPU/GPU source routes and bounded
component tests; this does not certify a release wheel or the other solvers.
No old-name package, import or CLI compatibility is provided.

Independent GPU distribution from the Radia monorepo: MATLAB CUDA MEX
kernels, readable Simulink blocks, and LTspice integration code.
Python CPU electromagnetic execution remains in `radia` (pybind11/NGSolve).
Explicit CPU/GPU routes and optional direct Python GPU bindings now have
bounded component validation. Each distributed native payload still requires
its own qualification; renaming alone does not establish acceptance.
`radia-mcp` is the AI operation interface. Neither package is required merely
to locate and use the installed GPU MATLAB files.

The 0.1.5 release adds experimental GPU IC(0)-PCG with retained factors
and `radia.gpu.HCurlICSystem`. CUDA performs the fixed-shift numerical
factorization; cuSPARSE performs the forward and transpose triangular solves.
Host code owns symbolic analysis and iteration control. The published 0.1.4
wheel remains unchanged. A100 comparison is pending while that host is shut
down; build-host measurements are numerical diagnostics, not formal timing.
See `validation_test/ic_cuda` for accuracy, lifecycle and scaling evidence.

Version 0.1.4 contains FP64 CUDA CSR Jacobi-PCG, a Simulink block for the same
kernel, and the selected standalone LTspice workflow. It also contains experimental FP64 assembled-operator AMS and its versioned
HDF5 exchange, bounded HDiv-MMM and Laplace FMM operators. This does not certify
arbitrary production AMS, ICCG, general HDiv-MMM, or migration of the older CPU
MEX application models.
The experimental AMS also includes sparse nonlinear affine-tetrahedron assembly
and a persistent `radia.gpu.AMSContext` with reusable GPU PCG workspace. Its
manufactured monotone-material verification does not certify measured steel or
arbitrary nonlinear devices. See the sparse study in validation_test/ams_cuda.
There is no CPU solver fallback.

`radia.gpu.HCurlSystem` retains the affine linear HCurl matrix, source vector,
AMS hierarchy and GPU workspace across solves. Geometry, boundary labels,
reluctivity and gauge mass remain fixed; `CurrentScale` changes only the
source amplitude. Construct a new system when any fixed input changes.
Nonlinear materials continue to use `radia.gpu.solveHCurl` and its adaptive
hierarchy updates. Supported HCurl geometry is first-order affine tetrahedra.
The AMS iteration keeps curvature control on the GPU and caches bounded
small coarse actions. Independent CPU comparisons, true residuals and
transfer-inclusive tuning evidence are documented in
[the FEM tuning study](../../validation_test/ams_cuda/FEM_TUNING.md).

FMM accuracy depends on the geometry and selected settings. The candidate's
independent scalar/dipole check uses `Order=8, Separation=2`; its measured
relative error is 1.64e-9. The same case at separation 1.5 gave 2.42e-8 and
does not meet the stable publication gate of 1e-8. Defaults are unchanged.
These bounds apply to the declared case, not arbitrary inputs; the separation-2
setting has no formal speed acceptance yet.

## Python solver ABI environments

The CPU `radia` distribution retains its exact NGSolve and Netgen 6.2.2607
pins. The optional `radia-simulink[gpu]` Python HDiv route requires NGSolve
6.2.2608 and native payloads built for that ABI. These two dependency sets
cannot coexist in one Python environment: use separate virtual environments
and explicit artifact boundaries for independent CPU validation. Do not force
installation with `--no-deps`, replace a solver DLL, or load both solver ABIs
in one process. The GPU user route does not require a CPU `radia` solve.
The base MATLAB/MEX distribution has no Python solver dependency; selecting
`[coupling]` adds `radia` and does not make it compatible with `[gpu]`.

This HDiv candidate does not migrate the CPU `radia` release to 2608.
Existing CPU validation artifacts built with 2608 remain evidence of their
recorded ABI only; they do not certify a CPU release with 2607 metadata.
A CPU ABI migration requires its own dependency, exporter compatibility,
native rebuild and release validation review.

## Install

```powershell
python -m pip install radia-simulink
radia-simulink-doctor --json
radia-simulink-path
```

Add the reported directory to MATLAB's path. MATLAB R2026a, Parallel Computing
Toolbox, an NVIDIA GPU and a compatible NVIDIA driver are external requirements.
Those requirements apply to the MATLAB/MEX interface. Qualified optional
Python GPU bindings use separately installed NVIDIA runtime packages and do
not require MATLAB or MATLAB Engine at execution time:

```powershell
python -m pip install "radia-simulink[gpu,coupling]"
```

The wheel must include the qualified Python native binding for the requested
solver; missing bindings or runtime dependencies raise without CPU fallback.
See [Python GPU runtime dependencies](docs/python-gpu-runtime.md). HDiv workflows
also require a checked `.vol` and its label contract; `check-vol` is supplied
by the separately installed mesh-export checker and needs no CAD application.
Its independent CPU accuracy validation is not a GPU-user prerequisite.
CSR-PCG consumes caller-provided matrix arrays and requires no Python Radia or
NGSolve runtime. The experimental native AMS candidate constructs its hierarchy
inside mexcuda and requires compatible NGSolve/Netgen native libraries, including
the permitted Python DLL dependency. GPU use requires no Python adapter,
interpreter initialization, or Python-generated solver snapshot. Native-library
deployment/packaging acceptance is separate from numerical candidate validation.
Optional
`radia-simulink[coupling]` installs Radia for workflows that explicitly need it.
The Windows x64 binary contains sm_75/sm_80 code and compute_80 PTX.
Develop and build `radia-simulink` on 100. Required release acceptance compares each
published CUDA solver against Python/NGSolve on 100 (Quadro RTX 5000), and
requires an independent review of the exact source. Record source/runtime
versions, artifact hashes and machine-readable results. A100 validation and
scaling are optional additional evidence, not publication requirements.
Formal GPU timing uses the on-demand A100 host when idle; 100 timings are diagnostic.
RTX 5000 installation acceptance does not certify A100 scaling.
Cubit mesh export is tested separately on 100. NGSolve validation
starts by checking and loading an exported `.vol` with its versioned label
contract; it does not require Cubit on the NGSolve host.

Deployment targets are 100 (build and RTX 5000 acceptance) and A100, excluding
102. The 102 host tests `radia` and `radia-mcp` without installing `radia-simulink`.
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
`radia_simulink.ltspice`, and the selected sampled-data
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

`radia-simulink.yml` builds and accepts the installed wheel on the build host.
`radia-simulink-release.yml` publishes that successful main CI artifact without
rebuilding, requiring its exact SHA-256 and installed acceptance provenance.
Stable versions additionally require the exact-source independent review and
same-wheel 100 RTX 5000/CPU parity record described in
[`validation_test/radia_simulink/releases`](../../validation_test/radia_simulink/releases/README.md).
Every CUDA module in the wheel needs independent CPU parity. CSR-PCG, AMS and
the HDiv material solve also require an original-equation true residual. FMM
and projection apply operators and must explicitly record that a solve residual
is not applicable, with a null residual and the measured operator error. A zero
invented residual cannot qualify these operator checks.
Configure PyPI Trusted Publishing for project `radia-simulink`, GitHub owner
`ksugahar`, repository `Radia`, workflow `radia-simulink-release.yml`, environment
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
`radia-simulink[exchange]` dependencies for the Python writer, and read the result
with `radia.gpu.readAMSExchange`. Values stay FP64 and CSR indices stay
zero-based int32. Shapes and array dimension order are explicit. HDF5 supports
the project's >2 GB exchange requirement; it does not bypass int32 kernel
capacity checks or the GPU working-memory budget. See
[AMS validation](../../validation_test/ams_cuda/README.md) for numerical gates,
transfer-inclusive timing and current limitations. Small-case A100 numerical
and timing evidence is recorded in that validation; it does not certify
A100 scaling or the exact artifact of a future release.

Python circuit converters use `radia-ltspice`; install `radia-simulink[ltspice]`
for the optional schemdraw route. MCP entry points, including
`mcp-server-radia-ltspice`, belong to `radia-mcp` and load the GPU distribution
only when circuit execution is requested. No circuit import alias is left in
the CPU Radia distribution.

The optional `coupling` extra requires Radia 5.3.1 or later with the CPU thermal-operator API.

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
5000 timings are diagnostic, not formal A100 performance acceptance.

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

The clean candidate `73b33291d` passed bitwise comparison of every hierarchy
array and sparsecholesky factor against the independent CPU oracle. The
continuously timed public API pipeline (four samples per backend, two ABBA
rounds after warmup) produced:

| HCurl DOFs | Before, seconds | CUDA generation, seconds | Speedup |
| --- | ---: | ---: | ---: |
| 86735 | 0.61350 | 0.43760 | 1.40 |
| 702073 | 4.99129 | 2.68588 | 1.86 |

At 702073 DOFs, maximum CPU solution relative difference was 4.12e-8 and
independently checked true relative residual was 7.64e-10. CUDA setup's peak
owned allocation was 1317974332 bytes (1.23 GiB); the public pipeline's live
GPU allocation delta was approximately 464841076 bytes (0.43 GiB), unchanged
between backends. These are different stages, not a whole-process peak.
Small 5326-DOF stage-sum timing regressed from 0.10574 to 0.11374 seconds.
The explicitly selected backend does not silently switch based on size.

Full JSON, memory definitions, artifact/source hashes, numerical gates and
limitations: `validation_test/ams_cuda/results/cuda-galerkin-73b33291d.json`.
The candidate is not a deployed release; formal performance acceptance remains
pending. No A100 operations were performed.
# 0.1.3 common scope

This candidate combines bounded HDiv-MMM, native affine lowest-order HCurl
FEM/Newton/AMS, and experimental Laplace FMM in one exact-source wheel.
Set `RADIA_HDIV_HMATRIX_BUILD`, `RADIA_HDIV_PROJECTION_BUILD` and
`RADIA_FMM_BUILD` to accepted clean builds of the wheel commit. Stale artifacts
are rejected even when their earlier numerical evidence passed.

HDiv supports tetrahedra with affine or P2 geometry, orders 1/2 and linear
isotropic material. Host-cache is opt-in and consumes additional memory.
General adaptive high-order HACApK, images, internal interfaces, constraints
and nonlinear HDiv are not accepted. FMM boundary maps are flat P1; the
independently measured direct-sum speedup is not a comparison with CPU FMM.
Small workloads can be slower on GPU. IH preview modules are not distributed.
Compatible user-installed native dependencies remain required; no implicit
CPU/Python solver fallback is provided. Separate historical performance
results do not substitute for installed common-wheel acceptance.


## Exact-wheel acceptance and performance

The main CI artifact from `88049e809` has SHA-256
`76a1b64ed290194a95aba755eb83ad6781538839012948c5adee1495259554ba`.
The same wheel passed isolated installed acceptance on the build GPU and an
A100, including independent CPU comparisons, original-equation residuals,
actual loaded MEX identity, repeated Simulink runs and failure cleanup.
The build-machine check also passed LTspice and the official Toolkit model check.

Fresh idle same-machine measurements used two ABBA cycles. GPU whole time
includes construction, all transfers, solve/application, result recovery and
context destruction. HDiv includes mesh-file loading; reference-fixture reads
and output serialization are excluded. FMM input-file reads and output
serialization are excluded.

| Bounded public case | CPU median | GPU median | CPU/GPU |
| --- | ---: | ---: | ---: |
| P2 sphere, HDiv order 2, host-cache; CPU 1 thread | 5.1781 s | 0.26823 s | 19.30 |
| Same HDiv case; CPU 16 threads, BLAS 1 thread | 0.43456 s | 0.25845 s | 1.68 |
| FMM, 8192 sources, 2048 targets, 3 channels, order 8, separation 2 | 4.4744 s | 0.99416 s | 4.50 |

The FMM CPU route is a single-threaded NumPy direct sum, not CPU FMM.
HDiv direct/cache coefficients were identical; independent true relative
residual was `7.96e-12` or lower. Its exact accounted owned GPU peak was
2286980 bytes with cache, compared with 450288 bytes for direct assembly.
Cache improves this workload's time while consuming more memory.
The FMM measured relative difference from direct sum remained below `1e-8`.

These cases do not establish large-model scaling or a universal GPU speedup.
The first numerical payload and MATLAB/GPU initialization are recorded
separately; complete process-cold solver time is not certified. HDiv host and
global CUDA/runtime peaks remain unmeasured. FMM native buffer counters and
50 ms process samples are distinct measures; samples include runtime/baseline
memory and may miss short allocation spikes. The HDiv CPU native file hash is
verified, but its original build-source receipt remains unverified.

Machine-readable qualification is in
`validation_test/radia_simulink/releases/0.1.3.json`; measured samples, input/source/
binary hashes and memory limits are in
`validation_test/radia_simulink/results/0.1.3-performance.json`.


## Experimental explicit CSR execution

`radia.csrPCG(..., "cpu")` accepts full host int32 CSR indices and FP64 arrays;
`radia.csrPCG(..., "gpu")` accepts the corresponding gpuArrays. Both use the
same mexcuda-built module. The current CPU kernel is serial and reports one
host thread. No thread pool is created and no route fallback is performed.
The dedicated `radia.gpu.csrPCG` entry remains an explicit GPU request.

For an assembled NGSolve real SPD matrix with constrained degrees of freedom
already handled by the caller, the optional direct Python adapter is:

```python
from radia_simulink.csr import solve_spd
solution, info = solve_spd(matrix, rhs, backend="gpu")
```

The GPU binding shares the CSR CUDA implementation with MEX, but uses a pybind11
ABI with explicit host/device transfers. GPU calls do not import MATLAB Engine.
A configured CUDA runtime and a separately qualified native extension are
required. Missing GPU support raises. The independent Python CPU selection is
`backend="cpu"`; it uses NGSolve `sparsecholesky` inside the caller's TaskManager
region. CPU accuracy validation is a release responsibility, not a GPU-user
initialization step.

Build the candidate with `Build.ps1 -CudaMexOnly -CsrPythonBindings` under an
explicit build admission. The binding is compiled with NGSolve/Netgen's pybind11
headers and the approved MATLAB-bundled CUDA/MSVC toolchain. A wheel containing
the extension carries its CPython ABI tag. No CUDA DLL is bundled. Runtime
redistribution and complete distribution qualification remain separate gates.
These CSR results do not establish HDiv, AMS, nonlinear or large-model parity.

The Python HDiv GPU solver and pybind build belong to `radia.hdiv_gpu`; see the Radia [development build preparation](../../docs/compilation/python-hdiv-gpu.md). This distribution owns mexcuda and the thin MATLAB host publisher only. The Python workflow requires the optional pure-Python `cubit-mesh-export` checker for strict `check-vol` validation; Cubit itself is not required.

Affine linear P2 GPU results obtained with the former degree-5 quadrature do not certify the degree-10 route. Independent CPU re-certification passed for one 30-DoF affine tet P2 linear isotropic case (inverse susceptibility 1, uniform applied H): coefficient error 4.11e-13 and original-equation residual 5.53e-15. The exact wrapper and native identities are in `validation_test/hmatrix_cuda/results/hdiv_linear_p2_degree10.json`. This bounded result does not certify other meshes or final release artifacts; nonlinear P2 evidence has its separately recorded scope.

The experimental numeric HMatrix boundary accepts explicit `assembly="normalized_leaves"` with `schema_version=2`. It preserves normalized leaf factors, scale and permutation through the same shared packer as Python; schema 1 remains the legacy HACApK struct conversion route, which also preserves factor values. This distinction declares the payload contract and bounded admission; it does not fix a re-normalization in the current schema 1 adapter. Offsets and permutations are one-based in MATLAB. Only upper/diagonal leaves are supplied; the packer mirrors off-diagonal leaves. Dense coefficients are flat row-major, low-rank a1/a2 are flat column-major columns-by-rank/rows-by-rank. This is a host-prepared numeric boundary, not CUDA geometry assembly or mixed-mesh certification.

`radia_simulink.hdiv.prepare_hdiv_host_setup(..., backend="gpu", label_contract=..., order=...)` publishes the setup produced by `radia.hdiv_gpu.prepare_hdiv_host_setup` as MATLAB numeric data under the caller's NGSolve TaskManager, without a CPU solution or resident GPU operator. MATLAB can pass the in-memory numeric setup to `PreparedHDivNonlinear('', order, H, bhTable, Assembly='ngsolve-host-normalized', HostSetup=setup)`. This explicit route preserves Python's NGSolve FE and Radia charge preparation and uses the shared GPU Newton; it requires neither MATLAB Engine in Python nor a file exchange. `native-tet` remains the default geometry route and rejects HostSetup. Host setup is CPU preparation, not CUDA geometry assembly; unsupported mixed/pyramid charge construction still fails. Mesh/label hashes describe the bytes checked at preparation, not a cryptographic attestation of an arbitrary caller-supplied struct.

Schema-2 normalized leaves must partition the upper matrix exactly once. A bounded per-row bitmap rejects duplicate, overlapping or missing blocks before resident allocation. Diagonal squares cover their upper triangle; off-diagonal blocks are mirrored. Legacy schema-1 and internal additive operator packing retain their existing contracts. Prefer the checked NGSolve/Radia snapshot publisher for HDiv workflows. Matrix symmetry/physical validity also remains a caller contract. `HostSetup` FE dimensions/types are checked before GPU initialization; native leaf/key/charge-map validation runs after GPU selection, before resident allocation.

## License

Owner-authored code is licensed under Apache-2.0 from the next releases.
See LICENSE, NOTICE and THIRD_PARTY_NOTICES.md. Third-party code and assets
retain their original licenses; previously published artifacts are unchanged.
