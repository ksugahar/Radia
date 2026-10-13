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
acceptance is pending because the probed mdx host has no NVIDIA GPU. A100
same-artifact acceptance is also pending. Numerical agreement alone is not a
claim of GPU acceleration or release readiness.
