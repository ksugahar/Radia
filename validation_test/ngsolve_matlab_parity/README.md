# NGSolve Python/MATLAB Native-MEX Validation

This validation proves that Radia's MATLAB MEX boundary uses the same NGSolve
finite-element meaning as the public Python API. MATLAB owns native NGSolve
`Mesh`, `FESpace`, `BilinearForm`, `LinearForm`, `GridFunction`, `BaseMatrix`,
and vector handles for the duration of each case. No Python callback executes
inside a MATLAB case.

This is a manual, result-bearing `validation_test` lane. The committed JSON
artifacts are checked by fast tests, but the solver-heavy run itself does not
belong in routine CI.

## Tiers

| Tier | Scope | Numerical evidence |
|---|---|---|
| Core 100 | 50 two-dimensional and 50 three-dimensional linear cases | Sparse entries, free DoFs, matvec, native solve, and residual |
| Breadth 500 | TRI, QUAD, mixed TRI/QUAD, TET, HEX, and WEDGE; real, complex, spatial-coefficient, and boundary forms | Python/MATLAB sparse matrices and matvecs |
| Scale 20 | H1, HCurl, and HDiv problems from 10,000 to fewer than 1,000,000 DoFs | Matrix-free matvec norms and energies, without transferring the sparse matrix to MATLAB |
| Manufactured 15 | Five element families at three orders | Analytic field samples, free-DoF residuals, and order convergence |

Python NGSolve is the differential oracle. MATLAB reads the exact same `.vol`
files and every tier must finish with zero live native handles. The deterministic
Netgen fixtures isolate API and numerical parity. Production-mesh validation is
a separate Cubit export -> `check-vol` -> Python/MATLAB parity lane; replacing
these fixtures with CAD-dependent meshes would weaken reproducibility rather
than improve it.

## Run

From MATLAB through the official MathWorks MCP server:

```matlab
addpath("S:\Radia\01_GitHub\matlab", "-begin");
addpath("S:\Radia\01_GitHub\validation_test\ngsolve_matlab_parity");
report = validate_ngsolve_matlab_100_cases();
extended = validate_ngsolve_matlab_extended();
```

The extended run writes each tier immediately after it completes, then writes
the aggregate result. This preserves completed evidence if a later heavy tier
fails:

- `results_ngsolve_matlab_100_case_parity.json`
- `results_ngsolve_matlab_breadth_500.json`
- `results_ngsolve_matlab_scale_20.json`
- `results_ngsolve_matlab_manufactured.json`
- `results_ngsolve_matlab_extended.json`

## HIBINO

HIBINO is the preferred host for the extended lane. Stage the exact checked MEX,
its five oneMKL sequential runtime DLLs, the MATLAB wrappers, and the generated
Python oracle. The aggregate result records the host, MATLAB/NGSolve versions,
MEX path, MEX SHA-256, duration, and final native-handle count. Copy the durable
machine result back as `results_ngsolve_matlab_extended_hibino.json` (and keep
the three tier JSON files with matching `_hibino` suffixes).

Temporary meshes and MAT oracles belong under `C:\temp`; only the result JSON
belongs in the repository.

## Verified Result

HIBINO passed the complete extended lane with MATLAB R2026a Update 4, NGSolve
6.2.2606, and MEX SHA-256
`9149f254621d30128a661e1d4159e5939459eaa49f82427204fd1ebb78e07beb`:

- Breadth: 500/500; maximum matrix relative error `3.21e-15`.
- Scale: 20/20; 21,515 to 134,130 DoFs; maximum matvec relative error
  `9.59e-16`.
- Manufactured solutions: 15/15; maximum final sample error `7.01e-16`.
- Total tier time: 1,172.31 s; final native handle count: 0.
