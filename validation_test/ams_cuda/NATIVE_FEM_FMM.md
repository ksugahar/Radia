# Native FEM, AMS and FMM candidate

The GPU workflow consumes checked `.vol` geometry and numerical densities in
MATLAB. Python is used only by independent validation. There is no Python solver
snapshot prerequisite or CPU solver fallback.

`radia.gpu.assembleHCurl` uses NGSolve mapped basis, curls, orientation, quadrature,
gradient and free DOFs, followed by CUDA element integration and deterministic CSR
reduction. `radia.gpu.solveHCurl` completes resident Newton, native AMS construction,
CUDA PCG and independent true-residual checks. The supported formulation is
lowest-order HCurl on affine tetrahedra, homogeneous scalar reluctivity, positive
gauge mass and constant volume current. Curved/higher-order geometry and more
general coefficient fields are rejected rather than approximated silently.

`radia.gpu.LaplaceFMM` implements numeric multilevel scalar/dipole `1/r` evaluation.
`biotSavartFMM` maps weighted current to SI vector potential and magnetic field;
`fieldCloudFMM` evaluates retained numeric field clouds. Neither adapter infers
wire singular corrections, periodic images or analytic RadFld geometry.
`prepareBoundaryFMM` prepares exterior flat P1 Galerkin SL/DL maps inside mexcuda.
NGSolve owns the trace basis, normals and quadrature; the existing shared
Sauter–Schwab kernel supplies explicit near corrections. The boundary operator
can be used in FEM–BEM coupling; this validation does not certify a complete
coupled formulation. Core compact HDiv-MMM interactions remain HACApK.

The Chebyshev interpolation method follows Fong and Darve, JCP 228 (2009),
doi:10.1016/j.jcp.2009.08.031, independently implemented. Expansion order and
separation are explicit accuracy choices, not universal guarantees. Sixth-order
physical field tests failed the original accuracy gate; eighth order with
separation two passed. No tolerance was relaxed.

Run `native_hcurl.py`, `fmm_numeric.py`, `fmm_physics.py` and `fmm_boundary.py`
with their explicit build/output arguments. Each saves success or failure JSON.
The boundary test uses independent NGSolve Galerkin integration, the physical
test uses analytic sphere layers and a circular current, and numeric tests use
NumPy direct summation. Native live-handle/repeated-call/error checks accompany
the numerical checks.

Build on 100 using `Build.ps1 -AmsCudaOnly -RequireNativeProvenance` and
`Build.ps1 -FmmCudaOnly -RequireNativeProvenance`. Use isolated binary paths;
do not replace another session's runtime. Native NGSolve DLL deployment remains
a separate packaging gate. MATLAB owns host scheduling; no active NGSolve
TaskManager region is created by mexcuda.

100 timings are diagnostic, including upload/download and warmed median.
Native owned memory and sampled device-free memory are separately identified;
neither is a whole-process peak measurement. Formal idle-mdx GPU performance
The original mdx timing gate could not run because that host has no NVIDIA GPU;
the user subsequently authorized NGSolve(A100), as recorded below.
Numerical agreement alone is not a claim of GPU acceleration or release readiness.

## Same-artifact A100 and native DLL acceptance (2026-10-13)

`accept_payload.py stage` freezes both previously tested MEX binaries, MATLAB
wrappers, public checked inputs and recursively resolved native dependencies.
`run_payload.ps1` starts an operation-owned MATLAB batch process with only the
payload runtime, MATLAB and Windows on its child PATH. No Python Engine,
NGSolve Python package path, interpreter initialization command or solver
snapshot is used by this execution route. Dependency DLLs are allowed; this
gate does not claim the absence of Python DLL imports.

The same complete 73-file payload passed on 100 and the already running
NGSolve/A100 host. Both loaded MEX hashes and all 33 native dependency DLL hashes
and paths were checked in the actual MATLAB process. The A100 independent CPU
re-evaluation of the nonlinear operator gave true residual 4.36e-12; the largest
cross-host relative difference over FE and FMM outputs was 1.45e-14. Both native
registries were empty after execution. The two MATLAB R2026a patch levels differ
and are recorded. See `results/payload-a100-20261013.json` and the matching
100 result. This accepts the declared candidate DLL payload, not a complete
installed radia-gpu wheel or a common-commit release candidate. Third-party DLLs
in the private acceptance bundle are not authorized wheel redistribution assets.

Fresh mdx inspection found only Microsoft Remote Display Adapter and VMware
SVGA 3D, with no NVIDIA device or nvidia-smi. The original pending JSON retains
that observation and is marked superseded by the user's explicit change to
NGSolve(A100). Previous diagnostic samples retain their original classification.

## Formal NGSolve(A100) evaluation (2026-10-13)

`measure_native_cpu_gpu.py` and `measure_native_sample.m` executed fresh same-host
measurements with the already accepted native binary/DLL hashes. The HDiv owner
confirmed no running/reserved job; CPU/GPU utilization and free-memory guards
passed before every sample. CPU and GPU host setup used one thread. Each case
has first CPU/GPU calls followed by two CPU/GPU/GPU/CPU cycles: 50 samples total.
GPU API times include construction, all transfers, execution, result recovery
and context release. Input fixture loading and result-file serialization are
excluded equally from numerical timing. MATLAB startup (5.738 s), device
initialization (0.415 s) and CPU script/import setup (2.989 s, excluding
interpreter launch) are separately recorded. "First" is the first call of that
case in the owned process, not a fresh operating system or device reset.

| Case | CPU warm median, s | GPU whole warm median, s | CPU/GPU |
| --- | ---: | ---: | ---: |
| Native HCurl assembly + nonlinear Newton/AMS, 5,326 DOF | 0.84460 | 0.52494 | 1.61 |
| Scalar/dipole FMM, 512 to 128 points | 0.01117 | 0.13308 | 0.084 |
| Scalar/dipole FMM, 8,192 to 2,048 points | 4.94516 | 0.23619 | 20.94 |
| Small sphere SL/DL, circular current and field-cloud pipeline | 0.07305 | 0.29879 | 0.244 |
| Flat P1 Galerkin SL/DL construction and application | 20.51374 | 3.01749 | 6.80 |

FMM and physical-pipeline CPU references are independent NumPy direct sums,
not an optimized CPU FMM. Boundary CPU assembly uses NGSolve directly. The
small FMM and physical pipeline are slower on GPU; no general speedup claim
follows. Nonlinear GPU solutions passed independent CPU operator re-evaluation
with maximum true residual 4.36e-12. All original field/operator accuracy gates,
native registry cleanup and actual loaded MEX/DLL identity checks passed.

Native AMS maximum phase-owned upper bound was 9,901,548 bytes. The large FMM
reported 3,570,276 geometry + 3,411,072 workspace + 196,608 output bytes, excluding
MATLAB-owned density arrays/runtime allocation overhead. Per-sample GPU free
memory endpoints and 50 ms CPU/MATLAB process RSS/private-memory samples are
saved. MATLAB host samples also cover output serialization outside timed API
work. These samples do not certify transient host/global-device peaks; no memory
reduction claim is made. See `results/formal-a100-20261013.json` for all samples,
hashes, limits, first calls, setup breakdowns and version identities.

Complete installed-wheel acceptance remains with the radia-gpu release owner
after common release source/artifact integration. Broader geometry, high-order
HCurl and larger-scale or different accuracy-parameter acceptance are outside
these bounded cases.
