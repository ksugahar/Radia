# radia-vim

**Volume Integral Method (VIM)** for cuboid eddy currents — NGSolve extension.

Implements the Newton-kernel volume Galerkin method for HDiv div-free interior
basis on an axis-aligned hex (cuboid). Used to extract Cauer ladder networks
for eddy currents in **isolated 3D conductors** under the boundary condition
`J·n = 0` (correct vacuum formulation, no air-box pathology).

Part of the **Sugahara-lab Radia public-fork family** — independent NGSolve
extensions absorbing the best techniques from worldwide EM solvers (FEMM, ELF,
COMSOL, …) and exposing them as pip-installable Python packages.

`radia-vim` complements **`ngsolve.bem` (ngbem)**, which is surface-only:
where ngbem provides surface integral operators (Laplace/Helmholtz/Maxwell
SL/DL), `radia-vim` provides the missing **volume Galerkin** route on a single
hex cell with HDiv div-free interior basis. No fork of NGSolve required.

## Why volume integral on a single hex?

Eddy-current Cauer ladder extraction in an isolated rectangular conductor
needs eigenvalues + Foster amplitudes of the operator
```
    K[u, v] = ∫_Ω ∫_Ω  G(r, r') · u(r) · v(r') dΩ dΩ
```
where `G(r, r') = 1 / |r - r'|` is the Newton kernel and `Ω` is a single
hex cell. There is no mesh — just a single 6-face hex. Standard FEM (NGSolve
`H1`/`HDiv`/`HCurl`) requires an air box that introduces an `AIR_SCALE`
pathology (τ_0 diverges with box size). Standard surface BEM (ngbem) cannot
handle volume basis functions. `radia-vim` fills this gap.

## Status

🚧 **Phase F-1** (bootstrap, 2026-05-07). Build infrastructure + precision
plumbing. Subsequent phases:

  - F-2: HDiv div-free hex basis (any order p, Schöberl-Zaglmayr style)
  - F-3: Self-coincident Sauter-Schwab quadrature for cube × cube
  - F-4: K matrix assembler (OpenMP parallel)
  - F-5: Foster solver + Cauer pipeline (Python wrapper)
  - F-6: order 6-8 production runs + paper figures
  - F-7: PyPI publication via OIDC trusted publisher

## Precision

The K matrix entries can be assembled at three precision levels, selected
at compile time via CMake options `RADIA_VIM_WITH_QUAD` and `RADIA_VIM_WITH_MPFR`:

| Type                       | Digits | Compile flag              | Use case                |
|----------------------------|--------|---------------------------|-------------------------|
| `double`                   | 16     | (always on)               | Order 3-4 (≤25 Cauer pairs) |
| `cpp_bin_float_quad`       | 32     | `RADIA_VIM_WITH_QUAD=ON`  | Order 5-6 (≤40 pairs)   |
| `mpfr_float` (60+ digits)  | 60+    | `RADIA_VIM_WITH_MPFR=ON`  | Order 7-8 (50+ pairs)   |

The empirical study driving these targets is in
[`scripts/precision_sensitivity_study.py`](../radia-axifemm/scripts/precision_sensitivity_study.py)
(in the sibling package `radia-axifemm`). Doubles fail catastrophically at
50 modes (87% mid-stage error in Cauer extraction); quad already gives 6%
mid-stage error; MPFR 60-digit recovers full reference accuracy.

## Build (development)

```bash
cd packages/radia-vim
pip install --no-build-isolation -e .              # double precision only
# or
CMAKE_ARGS="-DRADIA_VIM_WITH_QUAD=ON" \
  pip install --no-build-isolation -e .            # add quad precision
# or
CMAKE_ARGS="-DRADIA_VIM_WITH_QUAD=ON -DRADIA_VIM_WITH_MPFR=ON" \
  pip install --no-build-isolation -e .            # full precision spectrum
```

Requires NGSolve 6.2.2405+, CMake ≥ 3.16, a C++17 compiler, and (for quad/mpfr)
Boost.Multiprecision headers, plus libquadmath (GCC/Clang) or libmpfr/libgmp.

## Quick usage (Phase F-1)

```python
from radia_vim import version, precision_info, newton_kernel_double

print(version())                                  # 0.1.0
print(precision_info())                           # <PrecisionInfo double=yes ...>
print(newton_kernel_double([0,0,0], [1,0,0]))     # 1.0
```

## License

MIT.

## References

- S. Sauter and C. Schwab, *Boundary Element Methods*, Springer 2011.
- J. Schöberl and S. Zaglmayr, "High order Nédélec elements with local
  complete sequence properties," *Int. J. Comp. Math. Electr. Electron.
  Eng.* **24**(2):374–384, 2005.
- S. Hiruma and H. Igarashi, "Eddy-Current Analysis Using Cauer Ladder
  Network Method," *IEEE Trans. Magn.* **56**(3), 2020.
- H. Nagamine, T. Yamaguchi, K. Sugahara, S. Hiruma, T. Mifune, T. Matsuo,
  "Verified Numerical Computations of the Cauer Network Representation of
  a Square Prism Conductor," 2026 (preprint).
