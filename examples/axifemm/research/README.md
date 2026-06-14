# radia-axifemm

FEMM/Henrotte axisymmetric finite elements as an NGSolve extension. Brings the
Henrotte (1993) / FEMM (Meeker) trick of using `{1, r², z, r²z, ...}` shape
functions in physical `(r, z)` coordinates to NGSolve, so axisymmetric solvers
get FEMM-grade smooth `B_z` (piecewise constant per element) and well-behaved
`B_r ∝ 1/r` on the axis without forking NGSolve.

Part of the **Sugahara-lab Radia public-fork family** — independent NGSolve
extensions absorbing the best techniques from worldwide EM solvers (FEMM, ELF,
COMSOL, …) and exposing them as pip-installable Python packages.

## Status

✅ **Phase A2 complete.** Q1 / Q2 axis-aligned quadrilateral and P1 triangle
elements are implemented in C++ with closed-form (Mathematica-derived) element
matrices. The 3-way Cauer-ladder cross-validation against an independent
BEM-Foster reference passes 6/6 stages on a Cu disk benchmark.

| Element             | DOFs / cell | Basis on physical (r, z)                | Status |
|---------------------|-------------|-----------------------------------------|--------|
| `Q1_AxisAligned`    | 4           | `{1, r², z, r²z}`                       | ✅     |
| `Q2_AxisAligned`    | 9           | `{1, r², z, r²z, r⁴, r⁴z, r²z², ...}`   | ✅     |
| `P1_Triangle`       | 3           | `{1, r², z}`                            | ✅     |

## Cross-validation result (Cu disk, R=10 mm, t=2 mm, σ=5.8e7 S/m)

Cauer-ladder per-rung time constants `τ_pair[k] = L_{2k+1}/R_{2k}`
(Nagamine et al. 2026 convention, R_{2k} series, L_{2k+1} shunt):

| k | BEM Cauer (ref) | Q2 fine (ne=2530) | Q1 v.fine (ne=15170) | Q2 / BEM | Q1 / BEM |
|---|----------------:|------------------:|---------------------:|---------:|---------:|
| 0 | 219.32 μs       | 218.71 μs         | 218.05 μs            | -0.28 %  | -0.58 %  |
| 1 |  78.65 μs       |  78.12 μs         |  77.77 μs            | -0.68 %  | -1.12 %  |
| 2 |  40.04 μs       |  39.54 μs         |  39.37 μs            | -1.24 %  | -1.66 %  |
| 3 |  23.74 μs       |  23.16 μs         |  23.14 μs            | -2.46 %  | -2.54 %  |
| 4 |  17.07 μs       |  16.07 μs         |  16.06 μs            | -5.86 %  | -5.91 %  |
| 5 |  14.70 μs       |  13.12 μs         |  13.01 μs            | -10.77 % | -11.50 % |

Q2 beats Q1 at every stage (basis-order convergence study within the same FE
solver). The independent BEM-Foster reference (integral equation, elliptic
kernel, 1920-element ring mesh) is run via Mathematica (`bem_disk_axisym_cauer.wls`)
and a 50-digit mpmath classical Cauer extraction (`disk_bem_cauer.py`).

## Build

axifemm is built into the radia wheel (since 2026-05-10) — no separate
package install.  Just rebuild radia:

```powershell
pwsh -ExecutionPolicy Bypass -File S:/Radia/01_GitHub/Build.ps1
```

The C++ source lives in `src/ext/axifemm/` and the top-level
`CMakeLists.txt` defines an `add_ngsolve_python_module(axifem ...)`
target.  Output is copied to `src/radia/axifem.pyd` and ships in
the radia wheel.  Requires NGSolve 6.2.2603+, CMake ≥ 3.16, MSVC.

This `examples/axifemm/research/` directory is now a **research workspace**
(tests / scripts / demos) — it is not installable as a separate
package.  The `pyproject.toml`, `CMakeLists.txt`, and `src/` were
removed when axifemm was absorbed into the radia wheel (2026-05-10
cleanup).

## Quick usage

```python
from ngsolve import *
import radia                                          # registers DLL paths
from radia.axifem import (
    H1Henrotte,
    AxiHenrotteStiffnessBFI,
    AxiHenrotteSigmaMassBFI,
)

mesh = Mesh(...)                                  # axisymmetric (r, z) mesh
fes  = H1Henrotte(mesh, order=2,                  # 1 → Q1/P1, 2 → Q2
                  dirichlet="axis|outer")
u, v = fes.TnT()

mu_cf    = mesh.MaterialCF({"air": mu0, "conductor": mu0}, default=mu0)
sigma_cf = mesh.MaterialCF({"conductor": 5.8e7}, default=0.0)

a = BilinearForm(fes); a += AxiHenrotteStiffnessBFI(mu_cf);    a.Assemble()
m = BilinearForm(fes); m += AxiHenrotteSigmaMassBFI(sigma_cf); m.Assemble()
# ... feed (a.mat, m.mat) into eigsh / Hiruma 3-term / etc.
```

Standard NGSolve `BilinearForm` weak forms (`grad(u) * grad(v) * dx`) also
work via the registered DiffOps, but the closed-form custom integrators
above bypass NGSolve Gauss quadrature and reproduce the Mathematica
analytic element matrices to machine precision.

## Tests

```bash
python tests/test_closed_form_vs_python.py        # interior 4e-16 vs Python prototype
python tests/test_hiruma_disk_q1.py               # Q1, very-fine mesh, τ_1 = 223.06 μs
python tests/test_hiruma_disk_q2.py               # Q2, fine mesh
python tests/test_3way_cauer_cross_validation.py  # BEM vs Q1 vs Q2 (6/6 stages PASS)
```

## References

- F. Henrotte, K. Hameyer, et al., "A new method for axisymmetric linear and
  nonlinear problems," *IEEE Transactions on Magnetics* **9**(2):1352–1355,
  March 1993.
- D. Meeker, FEMM 4.2 axisymmetric formulation notes
  (`https://www.femm.info/wiki/AxisymmetricFormulation`).
- Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune, Matsuo, "Verified Numerical
  Computations of the Cauer Network Representation of a Square Prism
  Conductor," 2026 (Japan J. Industrial Appl. Math. submission). Cauer-ladder
  convention `R_{2k}` series, `L_{2k+1}` shunt.
- Hiruma & Igarashi, "Eddy-Current Analysis Using Cauer Ladder Network Method,"
  *IEEE Trans. Magn.* **56**(3), 2020. Three-term recurrence basis generation.

## License

MIT. Henrotte's algorithm is published; FEMM's `prob3big.cpp` is referenced
under the FEMM Aladdin Free Public License (only the algorithm structure is
inspected, no code is copied verbatim).
