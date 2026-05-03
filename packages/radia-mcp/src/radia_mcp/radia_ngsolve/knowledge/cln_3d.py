"""3D Cauer Ladder Network (CLN) knowledge base for Radia MCP server.

Captures Tanimoto-Kameari iterative CLN methods for 3D eddy current
analysis as developed in W:/00_CAE/NGSolve/谷本/ (Tanimoto's master's
thesis + production code for static/rotating machine analysis).

Covers ~25 notebooks across three thematic groups:
  - 修論 (master's thesis): A-T, T-Ω, A-Φ formulations + 2D reference
  - 定式_誤差検証 (formulation/error verification): A-Penalty, A-gauge,
    CASE C/D/E ablation
  - 20240910_静止器回転機用 (production for static/rotating machines):
    optimized A + ICCG, mixed A-Φ (Lukas), TEAM-28 scale

Each formulation produces a Cauer-II ladder {R_n, L_n} via Kameari's
iterative orthogonalization on impressed J source.
"""

CLN_3D_OVERVIEW = """
# 3D Cauer Ladder Network (CLN) — Tanimoto-Kameari Methods

## Mathematical Foundation

The Cauer Ladder Network reduces a 3D eddy current problem to a series
of (R_n, L_n) circuit elements via iterative orthogonalization in the
σ-weighted L^2 inner product on the conductor:

  R_n = 1 / ⟨J_n, J_n/σ⟩_cond
  L_n = R_n × ⟨J_n, A_pot_n⟩_cond     (accumulated Tanimoto)
  J_{n+1} = J_n - σ·A_pot_n / L_n     (Schmidt orthogonalization)

The first Cauer admittance Y_zz(s) = -μ_0·s·Σ a_n/(1+sτ_n) is
reconstructed from the ladder; Foster modes τ_n correspond to physical
eddy current decay times.

## Three Equivalent Formulations (equal in physics, differ in FE space)

### (1) A-T Formulation — Most Common
Two HCurl spaces:
  fesA = HCurl(mesh, order, nograds=True, dirichlet="<conductor BCs>")
  fesT = HCurl(mesh, order, nograds=True)         # auxiliary T field

Curl-curl in conductor with boundary E source:
  a(T,W) = ∫_cond (1/σ) curl T · curl W dx
  f(W) = -∫_∂cond (E_s × n) · W ds                # boundary E flux
  J = curl T          # current via curl of auxiliary

Then iterate Kameari on (R_n, L_n) using J and A:
  a(A,N) = ∫ (1/μ) curl A · curl N dx
  f(N) = ∫ R·J · N dx in conductor

### (2) T-Ω Formulation
Hybrid HCurl × H1(conductor-only):
  fesT = HCurl(mesh, order, nograds=True)
  fesΩ = H1(mesh, order, definedon=mesh.Materials("conductor"))

Scalar potential Ω confines current source to conductor; coupled
solve for (T, Ω) with Lagrange multiplier or coupled bilinear form.

Distinctive feature: Ω is conductor-only, reducing global DoF.

### (3) A-Φ Formulation
HCurl + H1 with Φ in conductor:
  fesA = HCurl(mesh, order, dirichlet=...)
  fesΦ = H1(mesh, order, definedon="conductor",
            dirichlet="in|out")

Body current J = σ∇Φ (implicitly through weak form). Lukas variant
explores order=3 with mixed space `V = HCurl × H1`.

## Iterative Algorithm (Common to All Formulations)

```python
# Initial impressed J_0 from boundary E or applied A_ext
J = sigma * A_ext         # or: J from boundary integral
Apot = 0
for n in range(N_STAGES):
    # Solve curl-curl A_n = J_n with chosen formulation
    solve(a, f_with_J_source, gfA_n)
    # R_n = 1 / <J, J/sigma>
    R_n = 1 / Integrate(J*J/sigma * dx, mesh)
    Apot += R_n * gfA_n
    # L_n = R_n × <J, Apot>
    L_n = R_n * Integrate(J * Apot * dx("conductor"), mesh)
    # Schmidt-orthogonalize J for next stage
    J = J - sigma * Apot / L_n
```

## Validation Cases (Cylinder, r×h)

For a cylindrical conductor with TM modes, analytical formulas exist:
  R_theory[n] = (2n+1) / (π·r²·σ·h)
  L_theory[n] = μ_0·(2n²+2n+1) / (8·(n+1)·π·h)

(Bessel-mode decomposition; reference: Tanimoto 修論.)

## Constraint/Gauge Variants (定式_誤差検証)

### A-Penalty
Add small penalty term to stabilize near-singular curl-curl:
  a_penalty(A, N) = (1/μ) curl A · curl N + ε/μ A · N    # ε ~ 1e-6
Avoids explicit gauge correction; loses ∇·A=0 exactly but in practice
gives same Cauer ladder values within tolerance.

### A-Gauge (Coulomb explicit)
Two-step:
  1. Solve A via HCurl(nograds=True): gives A with arbitrary gradient
  2. Project to Coulomb gauge: solve H1 problem ∇²φ = ∇·A, then A := A - ∇φ
Maintains ∇·A = 0 to machine precision but adds H1 solve per stage.

## Solver Variants

| Solver | Notebook | Distinctive |
|---|---|---|
| Direct sparse | A_direct, CASE_*_direct | Reference baseline |
| NGSolve CG + local pre | A_CG, CASE_*_CG | Standard Krylov |
| SparseSolvPy ICCG | CLN_AT (修論) | JP-MARs research backend |
| accICCG | CASE_*_accICCG, A_ICCG_最新版 | Acceleration param tuning |

Production recommendation: A + ICCG (20240917_A_ICCG_最新版.ipynb,
includes inline gauge correction).

## Common Boilerplate

```python
from netgen.occ import Cylinder, Z, OCCGeometry, Pnt
from ngsolve import (Mesh, HCurl, H1, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, curl, dx, ds,
                     Integrate, TaskManager, Cross, specialcf)
import scipy.sparse as sp

# Geometry
cyl = Cylinder((0,0,0), Z, r=0.01, h=0.01)
cyl.faces.Max(Z).name = "out"
cyl.faces.Min(Z).name = "in"
cyl.faces.name = "conductorBND"
cyl.mat("sig")
cyl.maxh = 0.001
mesh = Mesh(OCCGeometry(cyl).GenerateMesh(maxh=0.001)).Curve(3)

sigma = 1e6   # or 5.8e7 for Cu
mu = 4*pi*1e-7

# Per-stage solve (A-T form)
fesA = HCurl(mesh, order=1, nograds=True, dirichlet="in|out|conductorBND")
fesT = HCurl(mesh, order=1, nograds=True)
A, N = fesA.TnT()
T, W = fesT.TnT()

Es = (0, 0, 1)  # impressed E along z
n = specialcf.normal(mesh.dim)
a_T = BilinearForm(fesT)
a_T += 1/sigma * curl(T) * curl(W) * dx
f_T = LinearForm(fesT)
f_T += -Cross(Es, W.Trace()) * n * ds("conductorBND")
```
"""


CLN_3D_NOTEBOOK_INDEX = """
## Notebook Index (W:/00_CAE/NGSolve/谷本/)

### 修論/ — Master's Thesis (canonical reference)
- CLN_AT.ipynb         : A-T formulation, primary 3D Kameari (10-stage)
- CLN_T-Omega.ipynb    : T-Ω with H1 confined to conductor
- CLN_APhi.ipynb       : A-Φ with HCurl + H1(conductor)
- 2次元CLN.ipynb        : 2D scalar reference for validation
- メッシュ.ipynb         : OCC geometry + mesh utilities

### 定式_誤差検証/ — Formulation/Error Verification
- 20231211_A_(Penalty)_CLN.ipynb  : Penalty stabilization, single-stage test
- 20231221_A_gauge_CLN.ipynb       : Explicit Coulomb gauge, 4-stage
- 20240108_gauge_test.ipynb        : Gauge fix validation diagnostics
- 20240208_accICCG_practice.ipynb  : Solver tuning sandbox
- CASE_C_*.ipynb (CG, accICCG, direct)  : A nograds=True, no penalty
- CASE_D_*.ipynb                         : A with penalty 1e-6
- CASE_E.ipynb                           : Single-stage gauge variant

### 20240910_静止器回転機用/ — Production (latest)
- A-T_formulation.ipynb            : Production A-T baseline
- 20240917_A_ICCG_最新版.ipynb      : ⭐ Latest stable A + ICCG
- A_CG.ipynb                       : CG variant
- A_direct.ipynb                   : Direct sparse reference
- A_ICCG_*.ipynb (5 ablations)     : Penalty/stage/order/type1 sweeps
- A_CG_TEAM28size.ipynb            : TEAM-28 scale benchmark
- Lukas_A-Φ_test.ipynb             : Mixed HCurl×H1, order=3
- curlT = A.ipynb                  : T-from-A reconstruction validation

## Open Problem (本研究 2026-05)

**Kameari + Kelvin transformation has not been combined to date.**
The 2D practice example (W:/30_CauerLadderNetwork/2020_11_04_線形のCLNの練習/
CLN_H1_mode_Kelvin_NG.m) is COMSOL-based with auxiliary HelmholtzEquation
fields. 3D HCurl + Kelvin pullback for Kameari iteration remains a
research direction. The 2026-05-04 attempt (cuboid 5×2×1, NGSolve A
formulation) gave τ_0 = 326 μs vs target ~14 μs — discrepancy attributed
to A_ext = (B_0/2)(-y, x, 0) being unbounded at infinity, incompatible
with Kelvin pullback. Future direction: H-formulation or T-Ω with
reduced-Ω = -H_0 z + Ω_r where Ω_r decays at infinity.
"""


CLN_3D_KEY_FORMULAS = """
## Key Formulas Quick Reference

### Cauer-II Ladder Synthesis (from Kameari iteration)
After N stages, the admittance Y_Foster(s) = Σ a_n/(1+sτ_n) is
reconstructed from {R_n, L_n} via Cauer-II continued fraction at s=0:
  Y(s) = 1/(s·L_0 + 1/(R_0 + 1/(s·L_1 + 1/(R_1 + ...))))

### Foster Pole Identification
The dominant Foster pole τ_lead corresponds to MAX τ_n where τ_n=L_n/R_n,
not necessarily the first Cauer stage. For closed PEC cuboid 5×2×1:
  Analytical TE_z(1,1,0) τ = μ_0·σ·a²·b² / (π²·(a²+b²))
  Kameari max τ_n over 12 stages = 25.33 μs (vs analytical 25.46 μs,
                                              0.5% match)

### Schmidt Orthogonality Drift Diagnostic
  drift_n = max_{m<n} |⟨J_n, J_m/σ⟩_cond| / ⟨J_n, J_n/σ⟩_cond

For 3D HCurl order=3 with bonus_intorder=8 (closed PEC):
  N ≤ 11: drift ≤ 1e-12 (machine precision)
  N = 12: drift = 3e-12 (exponential growth onset, ~×6/stage)
  N = 25: drift = 5.5% (1% breakdown threshold)
  N ≥ 26: corrupted regime

### bonus_intorder Critical Setting
For order ≥ 3 HCurl with quadratic A_ext source, default integration
order is too low. Use `bonus_intorder=8` in all `dx()` calls to keep
Schmidt drift at machine precision through 11+ stages.
"""


def get_cln_3d_notebook(name: str = "list") -> str:
    """Retrieve raw Python code from a Tanimoto CLN notebook.

    Args:
        name: One of:
            "list"       - list available notebooks
            "AT"         - CLN_AT.py (A-T formulation, primary)
            "T_Omega"    - CLN_T_Omega.py (T-Ω formulation)
            "APhi"       - CLN_APhi.py (A-Φ formulation)
            "2D"         - CLN_2D.py (2D scalar reference)
            "production" - A_ICCG_production.py (latest 2024-09-17)

    Returns:
        Python script content (or list of available notebooks).
    """
    from .cln_notebooks import get_notebook, list_notebooks
    if name == "list":
        return list_notebooks()
    return get_notebook(name)


def get_cln_3d_documentation() -> str:
    """Return comprehensive 3D CLN (Cauer Ladder Network) documentation.

    Covers Tanimoto-Kameari iterative methods for eddy current analysis:
    A-T, T-Ω, A-Φ formulations with constraint variants (penalty, gauge)
    and production-grade solver implementations (A + ICCG).

    Source: W:/00_CAE/NGSolve/谷本/ (Tanimoto's master's thesis +
    production code for static/rotating machine analysis).

    Returns:
        Multi-section markdown string suitable for LLM agent consumption.
    """
    return (
        CLN_3D_OVERVIEW
        + "\n\n"
        + CLN_3D_NOTEBOOK_INDEX
        + "\n\n"
        + CLN_3D_KEY_FORMULAS
    )
