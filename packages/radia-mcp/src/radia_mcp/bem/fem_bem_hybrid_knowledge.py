"""FEM-BEM hybrid methods for open-boundary EM."""

OVERVIEW = r"""
# FEM-BEM Hybrid — open-boundary magnetostatics and eddy current

When the problem has both:
- A complex INTERIOR region (e.g. saturated iron core with nonlinear BH)
- An UNBOUNDED exterior (e.g. air around a transformer)

Pure FEM needs PML or large air mesh; pure BEM struggles with nonlinear
material.  Hybrid FEM-BEM splits the problem:

- **FEM** inside V_iron (nonlinear iron, complex geometry)
- **BEM** outside (open boundary, exact radiation condition)
- Coupling via Dirichlet-Neumann data exchange on interface Γ

## When to use what

| Strategy | When |
|----------|------|
| **Radia HDiv-VIM + Kelvin transform (FEM)** | Lab production: iron+air natural |
| **Pure FEM + PML** | When iron geometry is simple |
| **Pure FEM + Kelvin** | NGSolve native, no BEM needed |
| **Pure BEM (linear only)** | Conductors only, simple shapes |
| **FEM-BEM hybrid** | Mixed: complex iron + accurate open boundary |

## Lab's preferred alternative: Kelvin transform

The lab does NOT typically use FEM-BEM hybrid.  Instead:

```python
# calc_fem_kelvin.py (radia panels)
# Solve inside V_iron + air + V_kelvin (compact ball)
# Apply Kelvin transform on interface → open BC
# All in pure FEM (NGSolve)
```

This avoids the FEM-BEM coupling complexity and gets equivalent
accuracy.  See `radia_mcp.radia_ngsolve.kelvin_transformation`.

## When FEM-BEM hybrid IS used in the lab

- Coupling HDiv-VIM / Radia field evaluation with FEM (NGSolve) for transformer-like
  problems → via `netgen_mesh_to_radia` (M from FEM → Radia field evaluation)
- Coupling PEEC (Radia) with BEM (ngsolve.bem) for inductance
  validation → see `radia_mcp.peec.peec_inductance`

## References (this folder)

- [LOCAL] 08_fem_bem_hybrid/Iterative_FEM_BEM_Transformer.pdf
- [LOCAL] 08_fem_bem_hybrid/3D_Eddy_FEM_BEM_Coupling.pdf
- [LOCAL] 08_fem_bem_hybrid/Hybrid_FEM_BEM_open_boundary.pdf
- [LOCAL] 08_fem_bem_hybrid/Applications_Hybrid_FEM_BEM.pdf
- [LOCAL] 08_fem_bem_hybrid/HFSS_Hybrid_FE_IE.pdf — HFSS commercial
"""


COUPLING = r"""
# FEM-BEM coupling formulations

## Iterative coupling (weak coupling)

Iterate:
1. Solve FEM in V_iron with current Dirichlet/Neumann from Γ
2. Compute BEM update from FEM-given fields on Γ
3. Update interface conditions on Γ
4. Repeat until convergence

| Pros | Cons |
|------|------|
| Reuse independent solvers | Slow convergence for strong coupling |
| Black-box compatible | May not converge for high-contrast μ_r |

## Monolithic coupling (strong coupling)

Build a single block matrix:

```
[ K_FEM   C^T   ] [u_FEM]   [f_FEM]
[ C       K_BEM ] [u_BEM ] = [f_BEM]
```

where C is the interface coupling block.

| Pros | Cons |
|------|------|
| Robust for any μ_r contrast | Complex implementation |
| Single solve (no iteration) | Cannot reuse existing solvers |

## In ngsolve.bem + ngsolve.fem

NGSolve supports both via `ngsolve.bem`:
```python
from ngsolve.bem import (SingleLayerPotentialOperator,
                          DoubleLayerPotentialOperator)
# Build hybrid system manually with NGSolve BlockMatrix
```

Status in lab: documented; not yet a panel.  For Radia-NGSolve
hybrid see `netgen_mesh_to_radia` workflow in CLAUDE.md.
"""


def get_fem_bem_hybrid_knowledge(topic: str = "overview") -> str:
    """Dispatch FEM-BEM hybrid topics.

    Topics:
        overview    - When to use hybrid vs alternatives (DEFAULT)
        coupling    - Iterative vs monolithic coupling formulations
        all         - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "when"):
        return OVERVIEW
    if topic in ("coupling", "iterative", "monolithic"):
        return COUPLING
    if topic == "all":
        return "\n\n".join([OVERVIEW, COUPLING])
    return f"Unknown topic '{topic}'. Available: overview, coupling, all."
