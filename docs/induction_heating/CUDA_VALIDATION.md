# CUDA development and validation

Develop the reusable CPU and shared C++/CUDA kernels in `radia`, first exposed
through Python/pybind11. `radia-simulink` exposes them through thin MATLAB
mexcuda gateways and Simulink blocks. Follow the repository's current compiler,
runtime, explicit-backend, TaskManager and package-ownership policies.

## Shared numerical source

Actively share reusable C++/CUDA numerical kernels, domain data structures,
validation/error logic, lifecycle and host-control code between Python/pybind11
and MATLAB/mexcuda. Do not maintain a second implementation of the same
numerical algorithm merely for another binding. Shared domain source belongs
to `radia`; `radia-simulink` depends on it one-way. `radia` must not require
MATLAB/MEX or `radia_simulink` imports, linking or runtime dependencies.

Bindings remain thin adapters for ABI, types, memory ownership and checked
handles. Share the reusable state/validation/lifetime semantics; keep Python
and MATLAB ownership/conversion details in their adapters. Preserve provider-
owned allocation/destruction where applicable. Common host control does not
mean common host scheduling: MATLAB scheduling/S-Function lifecycle, compiler/
ABI build settings and package integration remain package-owned. Never bring
NGSolve TaskManager into mexcuda paths through a shared source dependency.

On Windows, Python/pybind11 CUDA compilation uses nvcc directly with MSVC host
C++; CPU/NGSolve C++ remains MSVC, and MATLAB uses mexcuda over shared kernels.
Use the current reviewed Windows baseline, not independent per-binding compiler
overrides. This policy authorizes no installation or host configuration change.

Common source and same-kernel agreement do not provide an independent physical
oracle. Keep pybind11/MEX bitwise parity and separately user-approved CPU/GPU
numerical parity, plus analytic/independent formulation validation. A shared
bug can reproduce in both bindings. No new implementation is certified here.

## Routine execution and GPU overflow

CPU jobs using at most eight physical cores each, with at most six concurrent
CPU jobs across the build host, are unconditionally pre-authorized. Count
existing jobs, bound CPU/thread use, preserve deadlines and owned-process
containment, and record resource observations. Idle cores, CPU-load thresholds
and fixed owner pools are not additional admission gates. Do not stop, repin or
overwrite another owner's jobs or runtime.

If the validation-host GPU is busy, GPU execution may overflow to a free
build-host GPU slot without another per-job grant for that routing. Count all
existing GPU jobs: at most one per host, at most two across the two hosts;
placing two on one host is not authorized by this rule. Retain the declared
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
