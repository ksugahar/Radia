# sparsesolv — Compact AMS / COCR for NGSolve

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

In-tree C++17 library that adds AMS (auxiliary-space Maxwell) and
COCR (Conjugate Orthogonal Conjugate Residual) on top of official
[NGSolve](https://ngsolve.org/), shipped as part of the Radia wheel.

- **Compact AMS** — auxiliary-space preconditioning for HCurl
  (Hiptmair–Xu 2007). HYPRE-free, header-only C++.
- **COCR** — short-recurrence Krylov solver for complex-symmetric
  systems (A^T = A) (Sogabe–Zhang 2007).
- **ICCG** — auto-shift IC(0) + ABMC parallel triangular solve.

Forked from [JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv).
Header-only C++17, supports both `double` and `std::complex<double>`.

## Distribution

This library is **not** distributed as a standalone PyPI package
(the legacy `ngsolve-sparsesolv` PyPI release was retired
2026-05-08).  It is built into the radia wheel and exposed as the
Python submodule `radia.sparsesolv_ngsolve`.

```bash
pip install radia
```

```python
import radia.sparsesolv_ngsolve as ssn
from radia.sparsesolv_ngsolve import (
    SparseSolvSolver,
    ICPreconditioner,
    CompactAMSPreconditioner,
    ComplexCompactAMSPreconditioner,
    CompactAMGPreconditioner,
    COCRSolver,
    GMRESSolver,
)
```

`import radia` (which transitively imports `ngsolve` and registers
DLL search paths on Windows) must run before the first
`import radia.sparsesolv_ngsolve`; calling the submodule import on
its own without going through `import radia` will fail with
"DLL load failed".

## Build (from the Radia repo root)

```powershell
pwsh -ExecutionPolicy Bypass -File Build.ps1
```

This drives a top-level CMake target, defined in `CMakeLists.txt` of
the Radia repository:

```cmake
add_ngsolve_python_module(sparsesolv_ngsolve
    src/ext/sparsesolv/ngsolve/python_module.cpp)
```

The output `sparsesolv_ngsolve.pyd` is copied to `src/radia/` and
ships with the radia wheel.  Do **not** run `pip install
src/ext/sparsesolv` — there is no standalone `pyproject.toml`
anymore.

## Solver-selection guidance

| Problem | FE space | Recommended solver | Reason |
|---|---|---|---|
| Poisson (H1, real) | H1 | ICCG | memory-efficient, fast |
| Curl-curl (real) | HCurl `nograds=True` | Shifted ICCG | auto-shift IC handles semi-definite |
| Magnetostatic, large | HCurl real p=1 | Compact AMS + CG | mesh-independent iteration count |
| Magnetostatic, nonlinear | HCurl real | Compact AMS + CG | `Update()` works inside Newton |
| Eddy current, complex, large | HCurl complex p=1 | Compact AMS + COCR | mesh-independent iteration count |
| Eddy current, complex, small/medium | HCurl complex | ICCG (`conjugate=False`) | memory-efficient |

## Performance

Hiruma 30 kHz Cu coil + Fe core (`mu_r=1000`), tol = 1e-10.
Intel Xeon 8-core, MSVC 2022, MKL 2024.2.

| Mesh | HCurl DOFs | Iters | Time | ms/iter | Memory |
|---|---:|---:|---:|---:|---:|
| mesh1_2.5T  |   155,527 | 144 |   4.5 s |   25.8 |   368 MB |
| mesh1_3.5T  |   197,395 | 168 |   7.3 s |   37.7 |   460 MB |
| mesh1_5.5T  |   331,595 | 249 |  16.2 s |   57.3 |   725 MB |
| mesh1_20.5T | 1,441,102 | 499 | 222.6 s |  396.2 | 2,933 MB |

ABMC-ICCG (IC only) on mesh1_3.5T diverges (17,178 iters, 438 s,
residual 2.8e-10) — IC cannot handle the curl-curl null space.
AMS resolves it via discrete-gradient + Nedelec interpolation
correction.

Reproduce: `python examples/hiruma/bench_compact_ams.py --all`.

## Quick start

### ICCG on a 2-D Poisson problem

```python
import radia                                   # registers DLL paths
from radia.sparsesolv_ngsolve import SparseSolvSolver
from ngsolve import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
fes  = H1(mesh, order=2, dirichlet="left|right|top|bottom")
u, v = fes.TnT()
a = BilinearForm(fes); a += grad(u)*grad(v)*dx; a.Assemble()
f = LinearForm(fes);   f += 1*v*dx;             f.Assemble()
gfu = GridFunction(fes)

solver = SparseSolvSolver(a.mat, method="ICCG",
                          freedofs=fes.FreeDofs(), tol=1e-10)
gfu.vec.data = solver * f.vec
```

### Compact AMS + COCR on a complex eddy-current problem

```python
import radia
import radia.sparsesolv_ngsolve as ssn
from ngsolve import *

# Complex system A = K + jw*sigma*M assembled as a (omitted).
# Real SPD auxiliary matrix on the same FES topology:
fes_real = HCurl(mesh, order=1, nograds=True,
                 dirichlet="dirichlet", complex=False)
u_r, v_r = fes_real.TnT()
a_real = BilinearForm(fes_real)
a_real += nu_cf * curl(u_r) * curl(v_r) * dx
a_real += 1e-6 * nu_cf * u_r * v_r * dx
a_real += abs(omega) * sigma_cf * u_r * v_r * dx("cond")
a_real.Assemble()

G_mat, h1_fes = fes_real.CreateGradient()
coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat, grad_mat=G_mat,
    freedofs=fes_real.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
    ndof_complex=fes.ndof, cycle_type=1, print_level=0)

with TaskManager():
    inv = ssn.COCRSolver(a.mat, pre, freedofs=fes.FreeDofs(),
                         maxiter=500, tol=1e-10)
    gfu.vec.data = inv * f.vec

print(f"Converged in {inv.iterations} iterations")
```

See [docs/compact_ams_cocr.md](docs/compact_ams_cocr.md) for the
formulation, configuration parameters and convergence theory.

### Stand-alone C++ (header-only)

```cpp
#include <sparsesolv/sparsesolv.hpp>

sparsesolv::SparseMatrixView<double> A(rows, cols, row_ptr, col_idx, values);
sparsesolv::SolverConfig config;
config.tolerance      = 1e-10;
config.max_iterations = 1000;
auto result = sparsesolv::solve_iccg(A, b, x, size, config);
```

## Directory layout

```
src/ext/sparsesolv/
├── include/sparsesolv/         # header-only library
│   ├── sparsesolv.hpp          # main header
│   ├── core/                   # types, config, matrix views, ABMC
│   ├── preconditioners/        # IC, Compact AMG, Compact AMS
│   ├── solvers/                # CG, COCR
│   └── ngsolve/                # NGSolve BaseMatrix wrappers + pybind11 bindings
├── examples/hiruma/            # complex eddy-current benchmark (30 kHz)
├── tests/                      # solver / preconditioner pytest suite
├── docs/                       # algorithm + API + tutorial markdown
└── LICENSE
```

`CMakeLists.txt` and `pyproject.toml` were removed in the 2026-05-08
in-tree merge: this directory is no longer a buildable target on its
own.  Use the top-level Radia `Build.ps1` / `CMakeLists.txt`.

## License

MPL 2.0.  See [LICENSE](LICENSE).  Forked from
[JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv) (also MPL 2.0).
