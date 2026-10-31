# CUDA development and validation

Develop the reusable CPU and shared C++/CUDA kernels in `radia`, first exposed
through Python/pybind11. `radia-simulink` exposes them through thin MATLAB
plain-MEX gateways (mex/matlab_add_mex) and Simulink blocks. Follow the repository's current compiler,
runtime, explicit-backend, TaskManager and package-ownership policies.

## Shared numerical source

Actively share reusable C++/CUDA numerical kernels, domain data structures,
validation/error logic, lifecycle and host-control code between Python/pybind11
and MATLAB/plain MEX. Do not maintain a second implementation of the same
numerical algorithm merely for another binding. Shared electromagnetic domain source belongs
to `radia`; `radia-simulink` depends on it one-way. Independently released
packages retain ownership of their own kernels. `radia` must not require
MATLAB/MEX or `radia_simulink` imports, linking or runtime dependencies.

Bindings remain thin adapters for ABI, types, memory ownership and checked
handles. Share the reusable state/validation/lifetime semantics; keep Python
and MATLAB ownership/conversion details in their adapters. Preserve provider-
owned allocation/destruction where applicable. Common host control does not
mean common host scheduling: MATLAB scheduling/S-Function lifecycle, compiler/
ABI build settings and package integration remain package-owned. Never bring
NGSolve TaskManager into MEX paths through a shared source dependency.

Build Radia-owned `.cu` kernels and host wrappers once as a shared CMake/nvcc
CUDA library. Develop first through pybind11; plain `mex` or `matlab_add_mex`
links the same compiled objects/library. mexcuda is retired; no fallback.
Record the common library hash and both bindings' link receipts: common source
alone does not prove common compiled code. Owner-scoped implementation migration
is separate from this policy; no existing route is certified by declaring it.

On Windows use MATLAB-bundled nvcc/CUDA 12.8 and MSVC 14.44; CPU/NGSolve C++
remains MSVC. Preserve the coordinated baseline review on MATLAB upgrades,
rather than independent binding-specific compiler overrides. This policy
performs no installation, compiler upgrade or host configuration change.

Link `libmwgpu` only when an adapter uses `mxGPUArray`. Verify CUDA runtime and
context identity/ownership across both bindings; do not assume unproven static-
`cudart` sharing is safe. Keep MATLAB headers and GPU-array adapters outside the
shared domain kernel. Explicit CPU MEX selection must still create no CUDA
context. Shared host code does not transfer MATLAB scheduling ownership.

Upstream NGSolve NVRTC is a separate, versioned/qualified route, not a fallback
for Radia-owned kernels. CUDA may use cuBLAS, cuSPARSE and cuSOLVER within its
explicitly selected, validated solver scope. Radia-owned native CPU math uses
MKL on x86 (PARDISO remains off) or NVPL on aarch64, not OpenBLAS; this statement
does not claim third-party wheel dependencies were rebuilt. CuPy is permitted
for development/validation; production uses pybind11 rather than an implicit
CuPy backend or fallback.

Common source and same-kernel agreement do not provide an independent physical
oracle. Keep pybind11/MEX bitwise parity and separately user-approved CPU/GPU
numerical parity, plus analytic/independent formulation validation. A shared
bug can reproduce in both bindings. No new implementation is certified here.

## Routine execution and GPU overflow

The validation host is the primary GPU development role; the build host is the
secondary GPU development role and owns CPU development, CI and formal release
artifact production. Validation-host GPU development uses isolated source
checkouts and job environments, separate from exact-wheel QA and the deployed
runtime. It does not expand normal deployment rights or permit Cubit/cme.
This role assignment alone authorizes no install, remote operation, host
configuration change or otherwise unapproved source-build route.

CPU jobs using at most eight physical cores each, with at most six concurrent
CPU jobs across the build host, are unconditionally pre-authorized. Count
existing CI, build, validation and research jobs by job root/process tree (compiler
children within one build are not separate jobs), bound CPU/thread use, preserve deadlines and owned-process
containment, and record resource observations. Idle cores, CPU-load thresholds
and fixed owner pools are not additional admission gates. Do not stop, repin or
overwrite another owner's jobs or runtime.

If the validation-host GPU is busy, GPU execution may overflow to a free
build-host GPU slot without another per-job grant for that routing. Count all
existing GPU jobs: at most one per host, at most two across the two hosts;
placing two on one host is not authorized by this rule. The build-host GPU job
uses a separate slot and is not also counted in the six CPU-job slots. Record
its actual host-control CPU threads/resources and avoid oversubscription. Retain the declared
artifact, architecture, compiler/ABI, runtime, ownership and feature scope.
Record the actual execution role, concurrent jobs and timings. This routing
decision does not authorize installs, configuration changes or an unverified
artifact/architecture substitution. Other remote, environment, deployment and
publication actions retain their existing authorization requirements.

## Physical accuracy and backend equivalence

Physical acceptance compares declared engineering quantities with an independent
converged reference. A tolerance on terminal resistance or inductance does not
imply the same tolerance at every field/current-density point. Check field
phase, conservation, boundary conditions and mesh convergence separately:
agreement of integrated quantities alone does not certify the entire field.

Backend equivalence compares identical declared inputs against the validated
Radia Python CPU reference, using separately user-approved numerical tolerances.
A physical 1% target is not automatically a CPU/GPU parity tolerance. On a shared
kernel, pybind11/MEX parity remains bitwise. Declare Simulink-only behavior
(ports, sample time, S-Function lifecycle or LTspice) and GPU-only kernels without
a CPU counterpart in the parity manifest, with reasons and separate validation.

GPU development can proceed alongside CPU validation; Python CPU accuracy
validation is not a prerequisite imposed on GPU users. Adoption must identify
the validated reference and its qualified scope. A backend name is not speedup
evidence: measure setup, transfers, solve, total time and memory, preserving
source/artifact identities. Do not claim an unexecuted backend is qualified.

## Fin observables

Retain the user-approved 1% targets on declared R/L observables and region-
integrated root/tip Joule heat. There is no blanket all-point 1% current-density
requirement. Freeze regions, frequency, excitation, current normalization,
ports and the inductive return/partial-inductance convention before judging.

For mean regional Joule power and peak drive-current phasors,
`R_region,eq = 2 P_region / |I_drive,peak|^2`; with RMS phasors the factor is one.
This is a loss contribution under the full-model excitation, not the intrinsic
two-terminal resistance of an isolated fin. A window-only sum is not a whole-
device resistance, and the full drive current need not pass through the fin.

The selected fin-L direction is an in-situ partial inductance including coupling
to the body, not an isolated-fin inductance. Freeze the self/mutual allocation
and matching reference convention before its 1% verdict; report both parts and
verify allocated contributions sum to the declared whole-model L, including
any region outside the measurement window. Regional magnetic energy alone
does not define a self-contained fin inductance. Do not double-count internal
reactive terms when comparing volume-current and surface-impedance models.

These criteria do not change the linear true-residual gate, reference
convergence requirements, baseline-regeneration or mutation procedures.

Thin MEX gateways may depend on NGSolve-inherited DLLs. This permission does
not permit including, linking or using NGSolve TaskManager in MEX builds;
shared-library dependency and scheduling ownership are separate contracts.

TaskManager governs CPU parallel regions only, including CPU assembly or
preconditioning inside an explicitly selected GPU workflow. The caller-owned
CPU rule and the documented p1 solver-owned exceptions remain CPU rules.
GPU device execution follows upstream NGSolve device-matrix/NVRTC conventions;
Radia does not extend TaskManager into a GPU execution model or add implicit
CPU fallback. MEX TaskManager prohibitions remain unchanged.

CPU numerical execution defaults to one thread. TaskManager, OpenMP and BLAS
parallelism are explicit opt-ins with recorded limits. Prefer independent
cases/trials for throughput, retaining the six-total-CPU-job and eight-physical-
core-per-job bounds. Compiler children within one build are not independent
jobs; numerical defaults do not change build parallelism or running jobs.
A 16-thread CPU/GPU comparison is an optional declared profile requiring the
appropriate resource scope, not the numerical default. Historical validated
records, benchmark values and baselines remain unchanged.

## Portable text contract

Source, scripts, configuration, documentation and policies are UTF-8 with LF
on Windows and Linux. Writers set encoding and LF explicitly. Git attributes
set the default while preserving binary, raw data, signed artifacts and frozen
hash-bound evidence; neither historical evidence nor another owner's WIP is
normalized. Shell/batch files have no blanket CRLF exception. The tracked-index
text contract rejects new CRLF/mixed blobs; fresh checkouts follow attributes.
No global Git configuration change is required or authorized by this policy.
