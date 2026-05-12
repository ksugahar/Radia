# Tanimoto canonical CLN reference implementations

Reference Jupyter notebooks from M.~Tanimoto's master's thesis (修論, 2025),
illustrating the Cauer-Ladder-Network (CLN) extraction in four formulations
on a 1 cm × 1 cm Cu cylinder ($\sigma=10^{6}$ S/m). Mirrored here from
`S:/NGSolve/谷本/修論/` on 2026-05-12 as the canonical reference against
which the radia-ngsolve CLN scripts (`../ngsolve_validation/...`) cross-validate.

## Notebooks

| Notebook | Formulation | Function spaces | Source | Cauer rung formulas |
|---|---|---|---|---|
| `2次元CLN.ipynb` | 2D scalar (axisym disk) | H1 order=1 with Dirichlet on outer boundary | volumetric $J = \sigma E$, $E=\text{const}$ | $R_k = 1/\int J^2/\sigma\,dx$, $L = \int R\,J\,A_{\text{acc}}\,dx$, $J\!\leftarrow\!J - \sigma A_{\text{acc}}/L$ |
| `CLN_AT.ipynb` | 3D vector A-T | HCurl(A, nograds, dirichlet="all") × HCurl(T, nograds, free) | $f_T^{(0)} = \int_{\partial_c}\!-(E_s\times W)\cdot n\,ds$ (Tanimoto boundary source) | $R_0 = 1/\int J_0^2/\sigma\,dx$, $A$-stage RHS $= N\cdot R J$, $L_{k+1} = \int B^2/\mu\,dx$ |
| `CLN_APhi.ipynb` | 3D A-Phi | HCurl(A) × H1(φ on conductor) | Dirichlet voltage $\phi=h$ on "in", $\phi=0$ on "out" | Same A-step as A-T; J update via $\phi$-equation $\sigma\nabla\phi\cdot\nabla\psi = \sigma\nabla\psi\cdot(-A_{\text{acc}}/L)$ |
| `CLN_T-Omega.ipynb` | 3D T-Omega | HCurl(T, nograds, free) × H1(Ω on conductor) | $f_T^{(0)} = \int_{\partial_c}\!-(E_s\times W)\cdot n\,ds$ (same as A-T) | Same R, L as A-T |
| `メッシュ.ipynb` | mesh utility | — | — | — |
| `CLN_AT-Copy1.ipynb` | early variant | — | — | development copy |

## Key Tanimoto patterns (reproduce these in any new CLN implementation)

1. **Boundary source for stage 0 in vector formulations**:
   ```python
   f += -Cross(Es, W.Trace())*n * ds("conductorBND")
   ```
   Enforces $\operatorname{curl} T = J = \sigma E_s$ on the conductor boundary via Stokes' theorem (no volumetric $J$ assignment needed).
2. **`HCurl(..., nograds=True)`** removes the grad-zero null space so the curl-curl matrix is invertible without an explicit tree-cotree gauge.
3. **No Dirichlet on the T-space**, full Dirichlet on the A-space (`in|out|conductorBND`). The asymmetric BC is part of the Cauer-CLN structure.
4. **Accumulator**: in 3D A-T / A-Phi / T-Omega, `Apot = sum_k gfA_k` (no R multiplier). In 2D scalar, `Apot = sum_k R_k * gfA_k`. The two normalisations give different absolute L, R values but the same $\tau_{\text{pair}}=L/R$.
5. **ICCG solver**: Tanimoto uses `SparseSolvPy` (JP-MARs) with `BadDivCount=10`, `BadDivVal=10.0`, `tol=1e-16`, `max_iter=200`, omega `1.1`. `Type1` direvgence detection.
6. **Closed-form validation reference** for 1 cm × 1 cm Cu cylinder:
   $R_{\text{theory}}\!=\!(2k+1)/(\pi r^2 \sigma h)$, $L_{\text{theory}}\!=\!\mu/(8(k+1)\pi h)$. These are *long-cylinder asymptotic* (1D radial diffusion) and break down on finite-aspect geometry.

## Why these live in examples/CLN/scripts/

Per the 2026-05-12 policy (memory: `feedback_no_ngsolve_py_in_cln_workdir.md`)
all radia-ngsolve consumers live in `examples/CLN/scripts/`. These six
notebooks are the canonical Tanimoto reference. The validation-stage and
debug variants from `S:/NGSolve/谷本/{定式_誤差検証,20240910_静止器回転機用}/`
were intentionally NOT mirrored here -- they are development scratch
notebooks rather than the final canonical implementation.

## What *not* to read from these as "best practice"

- The notebooks pre-date the H-H source projection that's required on
  *rectangular* / sharp-corner geometries (memory: `project_tanimoto_AT_HH_projection_breakthrough.md`).
  Tanimoto's notebooks all use the smooth cylinder, where the boundary
  source converges cleanly without an explicit H-H step.
- **T-Ω requires single-connectedness.** `CLN_T-Omega.ipynb` uses a Cu
  cylinder, which is simply connected, so the T-Ω matrix is non-singular
  out of the box. For multiply-connected conductors (plate with holes,
  torus, gear teeth ring — genus g ≥ 1), the T-Ω formulation needs g extra
  loop-current DOFs coupled to the FE system via the Hiptmair-Ostrowski
  Loop Method. See `../multiconn_loop_method/` for the reference
  implementation (EMPY-derived) and memory key
  `reference_loop_method_multiconnected_TOmega.md` for the algorithm
  summary.
- The 1 cm cylinder is in the *long-aspect* regime where the analytical
  $R_{\text{theory}}, L_{\text{theory}}$ formulas above hold. For Stoll-style
  sphere / TEAM 28 / cuboid 5×2×1 you must use the proper analytical reference
  (Stoll Bessel, BEM-Foster, etc.).
- ICCG omega=1.1 + BadDiv=10 are tuned for `order=1, nograds=True` HCurl;
  larger `order` or `complex=True` may require retuning.
