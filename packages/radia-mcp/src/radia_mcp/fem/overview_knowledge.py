"""FEM formulations overview, decision tree, lab stack."""

LAB_STACK = r"""
# Sugahara Lab FEM stack (production, 2026-05)

★ Production FEM library = **NGSolve 6.2.2603+** (PyPI)
  with lab additions in `src/ext/sparsesolv/` and Radia for open-boundary/MMM.

| Layer | Tool | Use |
|-------|------|-----|
| Mesh generation | Cubit (hex) / Netgen (tet) | STEP -> .vol |
| FE space | NGSolve `H1`, `HCurl`, `HDiv`, `HDivSurface` | standard |
| Axisym magnetic | `radia.axifem` (Henrotte basis) | `H1Henrotte` for A_phi curl-curl |
| Open boundary | Kelvin transform + radia | natural BC via Kelvin |
| Solver+precond | `radia.sparsesolv_ngsolve` | Compact AMS / COCR (HCurl) |
| Nonlinear | Picard / Newton with Hantila polarization | radia.hantila_solver |
| Hysteresis | `MatEnergyHysteresis` / `MatPlayHysteresis` | radia C++ |
| BEM coupling | `ngsolve.bem` (Weggler ST) | optional |
| Post / viz | `GmshPostExport` -> .msh v4.1 | GMSH GUI |

## Lab CORE formulation choices

| Problem class | Formulation | Why |
|---------------|-------------|-----|
| **Magnetostatic (scalar)** | **H1 + reduced Omega** | natural for nonlinear iron |
| **Magnetostatic (vector)** | **HCurl A (ungauged) + CompactAMS** | open boundary via Kelvin |
| **Eddy current MQS** | **HCurl A or A-V** with COCR + ComplexCompactAMS | complex symmetric |
| **Axisym magnetic curl-curl** | **Henrotte basis** (radia.axifem) | 1/r kernel handled exactly |
| **Axisym heat / scalar** | **standard H1 + 2*pi*r weighting** | FEMM convention |
| **Laminated stack** | **Hollaus MSFEM** (effective material) | radia.calc_motor_lamination |

## Cross-references

- `radia_mcp.matrix_solvers` — CG/COCR + preconditioners (Compact AMS theory)
- `radia_mcp.bem` — when to use BEM vs FEM (MMM/MSC) vs FEM-BEM hybrid
- `radia_mcp.motor` — Hollaus MSFEM application
- `radia_mcp.differential_forms` — exterior calculus foundation of H(curl)/H(div)
"""


DECISION_TREE = r"""
# Decision tree: which FEM formulation for which problem

```
Problem class
├── Scalar (no curl) - magnetostatic / electrostatic / heat
│   ├── Magnetostatic Phi
│   │   ├── Linear: H1 + Omega (single scalar)
│   │   └── Nonlinear iron: Reduced Omega (Omega-Omega split)
│   │       -> see "potential_formulations" topic
│   ├── Electrostatic
│   │   └── H1 + boundary Dirichlet/Neumann
│   ├── Heat / diffusion
│   │   └── H1 (standard)
│   └── Axisymmetric scalar
│       └── H1 + 2*pi*r weighting (NOT Henrotte; FEMM convention)
│         -> radia_ngsolve.axifemm_documentation
│
├── Vector curl-curl - magnetostatic vector / magnetodynamic
│   ├── Static, vector A
│   │   ├── ungauged + CompactAMS (lab default)
│   │   ├── tree-cotree gauge (small/medium)
│   │   └── Coulomb gauge penalty (legacy)
│   │      -> see "gauge_open_boundary" topic
│   ├── Eddy current MQS (frequency domain)
│   │   ├── Bulk conductor: A-V Biro-Preis
│   │   ├── Surface SIBC: HCurl A + Robin BC (lab IH workflow)
│   │   └── Thin sheet: shell formulation
│   ├── Axisymmetric magnetic A_phi
│   │   └── Henrotte basis {1, r^2, z}  ★ lab core
│   │      -> radia.axifem
│   └── Full Maxwell (high-freq, scattering)
│       └── HCurl + Helmholtz kernel (NOT lab focus per CLAUDE.md)
│
├── Mixed scalar+vector
│   ├── A-V coupled (conductor + air)
│   │   -> see "potential_formulations.av_biropreis"
│   ├── T-Omega (mixed)
│   │   -> see "potential_formulations.t_omega"
│   ├── H-formulation (HCurl with conductor body force)
│   │   -> see "potential_formulations.h_formulation"
│   └── Reduced potential (split into source + reaction)
│       -> see "potential_formulations.reduced_potential"
│
├── Time domain
│   ├── Frequency domain (single freq) → see eddy current branches
│   ├── Time-stepping → see "time_domain_axisym.time_domain_fem"
│   └── Periodic steady-state → harmonic balance
│      -> see "time_domain_axisym.harmonic_balance"
│
└── Special techniques
    ├── Open boundary
    │   ├── Kelvin transform (radia / NGSolve) ★ lab default
    │   ├── Asymptotic BC (older, simpler)
    │   └── PML (full-wave only)
    ├── XFEM (cracks, interfaces)
    ├── Isogeometric (smooth geometry)
    ├── Discontinuous Galerkin (DG)
    ├── Hierarchical / p-version
    │   -> see "elements" topic
    └── Multi-scale (laminated stacks)
        -> see "motor_hollaus_eddy" MCP (radia_mcp.motor)
```

## Common pitfalls (lab-observed)

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| HCurl A without gauge | iterative solver diverges | Use CompactAMS (handles kernel) or tree-cotree |
| Eddy current ω → 0 | system ill-conditioned | shifted preconditioner (matrix_solvers) |
| Reduced potential edge sources | wrong B at conductor surface | use total potential in conductor |
| H-formulation in air | natural BC mismatch | switch to A in air, H in conductor |
| Axisym at r=0 | 1/r singularity | Henrotte basis ★ |
| Laminated stack at h<delta | impractical mesh | Hollaus MSFEM effective material |
"""


HISTORY = r"""
# Genealogy: FEM-EM lineage 1980 -> 2025

```
1980  Bossavit edge elements (early)
1985  Nedelec edge basis (existence proof)
1989  Biro-Preis A-V gauge (multiply connected eddy current)
1990  Bossavit "Computational Electromagnetism" textbook
1992  Hiptmair tree-cotree gauge analysis
1993  Reduced scalar potential formulations
1995  Mur-style ABC, asymptotic BC
1998  Hiptmair multigrid for H(curl) (SIAM SINUM 36:204)
2000  Polarization method (alpha-trick) — Hantila revisits
2000  Biro-Preis edge element A_red (IEEE TMAG 36:3128)
2003  Hsiao-Wendland BIE textbook (BEM-FEM coupling)
2005  Isogeometric Analysis (Hughes)
2007  Hiptmair-Xu Nodal Auxiliary Space (★ basis of Compact AMS)
2009  XFEM for EM (cracks, multiscale)
2010  Sauter-Schwab BEM textbook
2013  Vico-Greengard MFIE-LF stabilization
2014  Hollaus MSFEM (multiscale FEM for laminated cores) ★
2017  Hanser Schur complement coil modeling
2020  Schöbinger effective-target homogenization
2021  Ostrowski-Hiptmair Two-Step Maxwell (frequency-stable)
2024  Hollaus-Schöbinger DEIM hyperreduction
2025  Frljic anisotropic effective μ tensor
```

## Sugahara lab specific milestones

```
2020  Compact AMS C++ implementation (HYPRE-free)
2022  Radia 4.x with NGSolve integration
2024  radia.axifem (Henrotte basis ★)
2025  calc_motor_lamination v5 (Frljic anisotropic μ_eff)
2026  radia_mcp ecosystem (this knowledge layer)
```
"""


def get_overview_knowledge(topic: str = "decision_tree") -> str:
    """Dispatch FEM overview topics.

    Topics:
        lab_stack       - Production FEM components + which to pick
        decision_tree   - Solver+formulation decision (DEFAULT)
        history         - Genealogy 1980 Nedelec -> 2025
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("lab_stack", "stack", "production"):
        return LAB_STACK
    if topic in ("decision_tree", "decision", "choose", "tree"):
        return DECISION_TREE
    if topic in ("history", "genealogy", "lineage"):
        return HISTORY
    if topic == "all":
        return "\n\n".join([LAB_STACK, DECISION_TREE, HISTORY])
    return (f"Unknown topic '{topic}'. Available: lab_stack, decision_tree, "
            "history, all.")
