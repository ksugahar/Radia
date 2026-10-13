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
acceptance is pending because the probed mdx host has no NVIDIA GPU.
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
SVGA 3D, with no NVIDIA device or nvidia-smi. Under this chat's idle-mdx timing
policy, formal GPU cold/warm ABBA, matched CPU/GPU throughput, scaling and peak
memory performance acceptance cannot run there. Remaining work is to designate
an authorized GPU-equipped formal timing host, verify its idle queues/processes
and memory, then run fresh same-input CPU/GPU ABBA measurements. Existing 100
and A100 diagnostic times must not be relabelled as formal measurements. The
HDiv chat's separately authorized A100 formal measurements do not certify these
AMS/FMM operators. Complete installed-wheel acceptance remains with the
radia-gpu release owner after common release source/artifact integration.
