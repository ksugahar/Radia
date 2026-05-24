"""Potential formulations: A-Omega, T-Omega, H-formulation, Reduced, Darwin."""

CATALOG = r"""
# Potential formulations catalog

For magnetostatic / magnetodynamic FEM, multiple potential choices exist.
Pick based on: source type, conductor regions, multiply-connectedness,
nonlinearity, and frequency range.

| Formulation | Unknowns | Where applicable | Lab use |
|-------------|----------|------------------|---------|
| **A** (vector) | A (HCurl) | Whole domain | ★ lab default for static |
| **A-V (Biro-Preis)** | A (HCurl) + V (H1, in conductors) | Bulk eddy current | for multiply-conn conductors |
| **A_red (reduced)** | A_red = A - A_source | Source-and-iron problems | reduces matrix size |
| **T-Omega** | T (HCurl, in conductors) + Omega (H1, outside) | Eddy current in conductors | ★ classical alternative |
| **H** (HCurl direct) | H (HCurl) everywhere | Superconductor / nonlinear | when sigma->inf |
| **Omega** (scalar) | Omega (H1) | Magnetostatic only, no eddy | nonlinear iron friendly |
| **Reduced Omega** | Omega_reduced (in iron) + Omega_total (in air) | Magnetostatic with current source | ★ accelerator magnets |
| **Darwin** | A + V with displacement current correction | DC to MHz transition | radia.darwin_model |

## How to choose

```
1. Static / time-stepping?
   ├── Static
   │   ├── Pure magnetostatic, no current source → Omega (scalar)
   │   ├── With current source (coils)            → Reduced Omega ★
   │   └── Vector quantity needed                 → A
   │
   ├── Frequency-domain MQS (eddy current)
   │   ├── Conductors are simply-connected         → T-Omega
   │   ├── Multiply-connected conductors           → A-V (Biro-Preis) ★
   │   └── Surface skin effect only                → A + SIBC (lab IH)
   │
   ├── Transient (TD-FEM)
   │   ├── Linear quasi-static                     → A or A-V
   │   └── Nonlinear iron                          → A + Picard/Newton
   │
   └── DC to MHz transition (capacitive + inductive)
       └── Darwin model (radia.darwin_model)
```

## Reference (this folder)

- [LOCAL] `10_FEM_定式化/01_A_Omega/` (4 files)
- [LOCAL] `10_FEM_定式化/02_T_Omega/` (1 file)
- [LOCAL] `10_FEM_定式化/03_H_formulation/` (5 files)
- [LOCAL] `10_FEM_定式化/04_Reduced_Potential/` (10 files) ← lab heavy use
- [LOCAL] `10_FEM_定式化/05_Darwin_Model/` (2 files)
"""


A_OMEGA = r"""
# A-Omega formulation (mixed vector-scalar)

Use A (HCurl) in conductor regions and Omega (H1 scalar) in air. Couple
at interfaces.

## Formulation

In conductor V_c (sigma > 0):
    curl(nu * curl A) + sigma * dA/dt + sigma * grad V = J_s

In air V_a (sigma = 0):
    H = -grad Omega
    grad Omega . n  =  (1/mu_0) * curl A . t   at interface

## When to use

- Eddy current in conductors surrounded by air
- Linear or nonlinear iron in V_a (or split into V_iron with sub-formulation)
- Lab uses this less often than pure A or T-Omega; mostly for legacy code

## References

[LOCAL] `10_FEM_定式化/01_A_Omega/` (4 files)

## Code pattern

```python
from ngsolve import *
fes_A = HCurl(mesh, order=2, definedon="conductor")
fes_O = H1(mesh, order=2, definedon="air")
fes = fes_A * fes_O
(uA, uO), (vA, vO) = fes.TnT()
a = BilinearForm(fes)
a += nu * curl(uA) * curl(vA) * dx("conductor")
a += sigma * uA * vA * dx("conductor") * timestep
a += mu_0 * grad(uO) * grad(vO) * dx("air")
a += interface_coupling_term  # n x H continuity
a.Assemble()
```
"""


T_OMEGA = r"""
# T-Omega formulation

Use T (electric vector potential, HCurl) in conductors, Omega (magnetic scalar
potential, H1) in air. T satisfies J = curl T.

## Formulation

In conductor V_c:
    H = T - grad Omega
    curl(rho * curl T) + d/dt(mu*(T - grad Omega)) = 0

In air V_a:
    H = -grad Omega    (no T)

## Advantages over A-V

- T directly encodes current (J = curl T), so current sources naturally enforced
- Scalar Omega in air → fewer DOFs than vector A
- Well-suited to simply-connected conductors

## Disadvantages

- For MULTIPLY-CONNECTED conductors: need "cuts" (cohomology) to make T well-defined
- Modern fix: use A-V (Biro-Preis) instead which handles multiply-connected naturally

## References

[LOCAL] `10_FEM_定式化/02_T_Omega/` (1 file)

## When to use

- Bulk eddy current in **simply-connected** conductor (e.g. solid torus → cut needed)
- When scalar potential in air is preferred for memory
- Legacy / historical formulation; lab production uses A or A-V more often

## Cross-reference

- For multiply-connected: see `bem.surface_ie.pmchwt` and `matrix_solvers.em_specific.tree_cotree`
- Cohomology cuts: see `radia_mcp.differential_forms.cohomology`
"""


H_FORMULATION = r"""
# H-formulation (HCurl direct)

Use H (HCurl) directly as the unknown. Suitable for superconductor /
high-conductivity / nonlinear B(H).

## Formulation

In conductor:
    curl E = -dB/dt
    where B = mu(|H|) * H (nonlinear),
          E = rho(|J|, T) * J,  J = curl H

Result:
    curl(rho * curl H) + d/dt(mu(|H|) * H) = 0    on V_c

In air:
    curl(curl H) = 0  (subject to grad of scalar Omega above)

## When to use

★ Superconductor modeling (E-J power law: rho = E_c/J_c * (J/J_c)^(n-1))
- High-frequency conductors (sigma -> infinity)
- Nonlinear B(H) directly without M reconstruction

## References (this folder)

[LOCAL] `10_FEM_定式化/03_H_formulation/` (5 files)

## Lab use

- Superconducting magnet design (in `radia_mcp.accelerator`, `radia_mcp.fusion`)
- High-Tc YBCO tape modeling
- Not yet a production lab panel; reference layer only
"""


REDUCED_POTENTIAL = r"""
# Reduced potential formulations (★ accelerator magnet lab core)

Split the field/potential into a "source part" computed analytically
(Biot-Savart) and a "reaction part" solved by FEM.

## A_reduced

Define A = A_source + A_red
- A_source = analytical Biot-Savart from coils
- A_red = FEM unknown, with reduced equation:
    curl(nu * curl A_red) = -curl(nu * curl A_source) + J_iron

Saves DOFs: A_source captures the "smooth" part, A_red only the
"reaction near iron".

## Omega_reduced

For static magnetostatic with current source:
- H = H_source + H_red = -grad Omega_reduced  (in iron)
- H = H_source                                (in air)
- Continuity of B.n at interface

## When to use

★ Accelerator magnet design (Radia/NGSolve, lab production):
- Coil source = Biot-Savart from `CoilBuilder.to_radia()`
- Iron yoke = FEM with reduced Omega
- No coil mesh needed → huge mesh savings
- See `radia_mcp.electromagnet` and `radia_mcp.accelerator`

## References

[LOCAL] `10_FEM_定式化/04_Reduced_Potential/` (10 files) ← heavy lab use
[LOCAL] `10_FEM_定式化/01_A_Omega/` for A_red discussion

## Critical pitfall

In nonlinear iron, reduced potential schemes can have **cancellation errors**:
H_source - grad Omega ≈ 0 in heavily-saturated regions, losing precision.
Switch to total Omega in heavily-saturated regions, reduced in air.
"""


DARWIN = r"""
# Darwin model (DC to MHz, capacitive + inductive)

The Darwin approximation includes magnetic-induced electric fields
(inductive coupling) AND charge-induced electric fields (capacitive
coupling) but **NEGLECTS radiation** (no dB/dt → ∇×E term).

## Where it matters

| Frequency | Regime | Equations |
|-----------|--------|-----------|
| DC | static | div E = ρ/ε, curl H = J |
| < 1 kHz | MQS | curl H = J, add Faraday curl E = -dB/dt |
| 1 kHz - 1 MHz | Darwin (lab) | MQS + capacitive currents |
| > 100 MHz | Full Maxwell | + displacement current ∂D/∂t |

The Darwin regime is ideal for **PCB power electronics** and PEEC where
both inductive (eddy current) and capacitive (parasitic C) effects matter
but radiation is not yet relevant.

## Lab use

- `radia.darwin_model` (in `rad_darwin.cpp`) — Darwin gauge implementation
- Kaimori-Mifune-Kameari-Wakao TD model (2024 paper) → `radia_mcp.motor`
- PEEC unstructured formulation (Codecasa) → `radia_mcp.peec`

## References

[LOCAL] `10_FEM_定式化/05_Darwin_Model/` (2 files)
- Kaimori-Mifune Darwin TD paper (lab knowledge in motor MCP)
- Codecasa "Unstructured PEEC formulations considering R, L, C effects"

## Cross-reference

- `radia_mcp.peec.peec_knowledge` — PEEC's natural Darwin support
- `radia_mcp.motor.hollaus_eddy` — Hollaus + Darwin coupling
"""


def get_potential_formulations_knowledge(topic: str = "catalog") -> str:
    """Dispatch potential formulation topics.

    Topics:
        catalog          - All formulations + decision tree (DEFAULT)
        a_omega          - A-Omega mixed vector-scalar
        t_omega          - T-Omega (electric vector + magnetic scalar)
        h_formulation    - H-formulation (HCurl direct) for SC / nonlinear
        reduced          - ★ Reduced potential (lab accelerator magnet core)
        darwin           - Darwin model (DC to MHz transition)
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("catalog", "overview", "compare"):
        return CATALOG
    if topic in ("a_omega", "a-omega", "aomega"):
        return A_OMEGA
    if topic in ("t_omega", "t-omega", "tomega"):
        return T_OMEGA
    if topic in ("h_formulation", "h-formulation", "h_form"):
        return H_FORMULATION
    if topic in ("reduced", "reduced_potential", "reduced_omega", "omega_reduced"):
        return REDUCED_POTENTIAL
    if topic in ("darwin", "darwin_model"):
        return DARWIN
    if topic == "all":
        return "\n\n".join([CATALOG, A_OMEGA, T_OMEGA, H_FORMULATION,
                            REDUCED_POTENTIAL, DARWIN])
    return (f"Unknown topic '{topic}'. Available: catalog, a_omega, t_omega, "
            "h_formulation, reduced, darwin, all.")
