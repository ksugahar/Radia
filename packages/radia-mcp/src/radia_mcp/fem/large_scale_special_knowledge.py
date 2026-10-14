"""Large-scale FEM, error theory, multi-scale (Hollaus link), special techniques."""

LARGE_SCALE = r"""
# Large-scale FEM (millions of DOFs)

For motor / accelerator / fusion simulations with N > 1M DOFs, special
techniques are required.

## Parallelization strategies

| Approach | Where |
|----------|-------|
| **Shared-memory (TaskManager)** | NGSolve `with TaskManager():` ★ |
| **Distributed (MPI)** | NGSolve MPI build; Hypre, MUMPS |
| **GPU** | (lab not production) |
| **H-matrix acceleration** | HACApK for BEM, see `radia_mcp.bem` |

## Solver scaling

| Solver | N scaling | Notes |
|--------|-----------|-------|
| LU direct | O(N^2) memory, O(N^{1.5}) for 3D ellipt | breaks at ~500k |
| BiCGSTAB + IC | O(N) memory, O(N log N) iter | classic baseline |
| AMS-PCG (Compact AMS) | O(N) memory, O(N log N) iter | ★ lab default for HCurl |
| BoomerAMG (HYPRE) | O(N) memory, O(N) iter | distributed (MPI) |

## Mesh strategies for large-scale

1. **Adaptive refinement**: refine only where solution is rough
2. **Hierarchical decomposition**: nested meshes, restriction/prolongation
3. **Domain decomposition**: split into sub-domains with overlap
4. **Multi-scale FEM (Hollaus)**: macro + micro problems coupled

## References (this folder)

[LOCAL] `10_FEM_定式化/32_大規模問題/` (16 files, 219 MB)
- Iterative methods, parallel solvers, AMG variants
- Distributed memory implementations
- Cross-reference: `radia_mcp.matrix_solvers` for solver theory
"""


ERROR_THEORY = r"""
# Error theory & a-posteriori error estimators

For adaptive mesh refinement (hp-adaptive), need a-posteriori error
estimators to decide WHERE to refine.

## Classic estimators

| Estimator | Type | What it detects |
|-----------|------|-----------------|
| Residual-based | local, cell-wise | element residual + jump term |
| Hierarchical | difference of orders | h vs p+1 solution gap |
| Recovery (Z-Z, SPR) | post-process gradient | smoothed vs raw gradient |
| DWR (dual-weighted residual) | goal-oriented | error in functional, not norm |

## In NGSolve

```python
# Built-in residual error estimator
err = ZienkiewiczZhuRefinement(gfu, mesh)
for el, e in zip(mesh.Elements(), err):
    if e > threshold:
        el.SetRefinementLevel(...)
mesh.Refine()
```

## Lab use

- Not yet production. Lab typically uses h-uniform refinement
  (the docs mesh-evaluation demo sweeps p=1..5 on the same mesh)
- Reference layer; could add adaptive refinement to motor / accelerator
  panels later

## References (this folder)

[LOCAL] `10_FEM_定式化/31_誤差の理論/` (7 files, 96 MB)
"""


HOLLAUS_MSFEM = r"""
# Hollaus MSFEM — Multi-Scale FEM for laminated cores (★ lab application)

Reference: Karl Hollaus et al. (TU Wien, 2014+).
See `radia_mcp.motor.hollaus_eddy` for detailed code recipe.

## Why MSFEM for laminated cores

Iron cores in motors / transformers are LAMINATED (stack of thin steel
sheets ~0.35mm with insulation). To resolve eddy currents in each
lamination directly with FEM:
- Mesh size h ≪ 0.35mm → millions of elements per cubic cm
- Skin depth at 50 Hz ≈ 0.7mm; at 1 kHz ≈ 0.1mm → even finer mesh
- Impractical for full machine analysis

## Hollaus's solution: MSFEM

Decompose A into:
- A_0: macro-scale (smooth across many laminations)
- A_1: micro-scale (1D eddy current within each lamination)

Solve macro problem on coarse mesh + cell problem in 1D representative
volume. Effective material law μ_eff(ω, B²) captures lamination physics.

## Lab implementation

`radia.calc_motor_lamination` (v3 voltage drive, v4 nonlinear ν(B²)
with DEIM, v5 Frljić 2026 anisotropic μ_eff tensor)

```python
import radia
result = radia.calc_motor_lamination(
    vol_path="laminated_motor.vol",
    frequency=50,
    effective_target="reaction",   # Schöbinger 2020
    formulation="voltage_drive",   # Hanser 2025
)
```

See `radia_mcp.motor.hollaus_eddy` for the full genealogy and code.

## References (this folder)

- Hollaus papers indirectly via `18_MOR_モデル縮約/12_Karl_Hollaus_MSFEM/`
  (34 files)
- DEIM hyperreduction → `radia_mcp.mor.systematic.hyperreduction`

## Cross-references

- `radia_mcp.motor.hollaus_eddy` — full code recipe
- `radia_mcp.mor.systematic` — DEIM theory
- `radia_mcp.matrix_solvers.preconditioners_ams` — solver for the system
"""


MISC_TECHNIQUES = r"""
# Miscellaneous FEM techniques (residual + niche)

10 files in `10_FEM_定式化/99_misc/` cover:

- **Fixed-point methods** for nonlinear (Picard variants, Polarization)
  - `Fast_FixedPoint_Hysteresis.pdf`
  - `Locally_Convergent_FixedPoint_NonLinear.pdf`
  - `JA_FixedPoint_TimePeriodic_Hysteresis.pdf`
- **Polarization method** (Hantila's α-trick)
  - `Polarization_Method_Static_Fields.pdf`
- **Newton-Raphson + Krylov coupling** (relaxed convergence)
  - `Fast_nonlinear_NR_Krylov_relaxed.pdf`
- **Helmholtz-Weyl decomposition** (3D L^r vector field)
  - `Helmholtz_Weyl分解_3次元Lr.pdf`
- **Electric vector potential F** formulation
  - `電気ベクトルポテンシャルFを用いた定式化.docx`
- **2nd-order element reordering**
  - `2次要素の並び替え.pptx`

## When to consult these

- Implementing custom nonlinear solver (Picard variants, Newton-Krylov)
- Need decomposition for divergence-free / curl-free splitting
- Special element re-ordering for memory access patterns

## Cross-references

- `radia_mcp.motor.hollaus_eddy` for Picard-DEIM coupling
- `radia_mcp.matrix_solvers.krylov.bicgstab` for Krylov methods
"""


def get_large_scale_special_knowledge(topic: str = "large_scale") -> str:
    """Dispatch large-scale + special techniques topics.

    Topics:
        large_scale       - Parallel scaling, AMS/AMG, DD (DEFAULT)
        error_theory      - A-posteriori error estimators, adaptivity
        hollaus_msfem     - ★ Multi-scale FEM for laminated cores (lab app)
        misc_techniques   - Fixed-point, polarization, Helmholtz-Weyl, F potential
        all               - Everything
    """
    topic = topic.lower().strip()
    if topic in ("large_scale", "large", "scaling", "parallel"):
        return LARGE_SCALE
    if topic in ("error_theory", "error", "adaptivity", "a_posteriori"):
        return ERROR_THEORY
    if topic in ("hollaus_msfem", "hollaus", "msfem", "multi_scale"):
        return HOLLAUS_MSFEM
    if topic in ("misc_techniques", "misc", "fixed_point", "polarization"):
        return MISC_TECHNIQUES
    if topic == "all":
        return "\n\n".join([LARGE_SCALE, ERROR_THEORY, HOLLAUS_MSFEM,
                            MISC_TECHNIQUES])
    return (f"Unknown topic '{topic}'. Available: large_scale, error_theory, "
            "hollaus_msfem, misc_techniques, all.")
