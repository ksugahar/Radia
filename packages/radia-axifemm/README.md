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

🚧 **Phase 2-A** (bootstrap). Build infrastructure is in place; the C++
finite-element classes (`AxiHenrotteFE_Q1_AxisAligned`,
`AxiHenrotteFE_P1_Triangle`) are stubbed. A Python prototype that this
extension is being built to replace lives in
`W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/` (Triangle P1 0.1%
match vs FEMM NMR; Q1 0.5% match vs BEM-Foster on Cu disk axisym).

## Reference

- F. Henrotte, K. Hameyer, et al., "A new method for axisymmetric linear and
  nonlinear problems," *IEEE Transactions on Magnetics* **9**(2):1352–1355,
  March 1993.
- D. Meeker, FEMM 4.2 axisymmetric formulation notes
  (`https://www.femm.info/wiki/AxisymmetricFormulation`).

## Build (development)

```bash
pip install --no-build-isolation -e .
```

Requires NGSolve 6.2.2405+ and a C++17 compiler with CMake ≥ 3.16.

## Quick usage (Phase 2-D target API)

```python
from ngsolve import *
from radia_axifemm import H1Henrotte

mesh = Mesh(...)                          # axisymmetric (r, z) mesh
fes  = H1Henrotte(mesh, dirichlet="axis|outer")
u, v = fes.TnT()
nu   = 1.0 / mu
a    = BilinearForm(nu * grad(u) * grad(v) * x * dx).Assemble()
# ... standard NGSolve weak form usage
```

## License

MIT. Henrotte's algorithm is published; FEMM's `prob3big.cpp` is referenced
under the FEMM Aladdin Free Public License (only the algorithm structure is
inspected, no code is copied verbatim).
