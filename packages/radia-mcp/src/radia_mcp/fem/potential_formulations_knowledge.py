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
   ├── Static, NONLINEAR (saturable) iron + coil source
   │   ├── B-input A_red  nu(|B|)  (convex)        → ★ well-conditioned at ANY mu  (topic `nonlinear`)
   │   └── NOT reduced-Omega mu(|H|): STALLS at high mu (the saturation knee regime)
   │
   ├── Transient (TD-FEM)
   │   ├── Linear quasi-static                     → A or A-V
   │   └── Nonlinear iron                          → A_red B-input + Picard/Newton (topic `nonlinear`)
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

- Superconducting magnet design (in `radia_mcp.accelerator`, `radia_mcp.fusion_reactor`)
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


NONLINEAR_AFORM = r"""
# Nonlinear (saturable) magnetostatics: B-input A-formulation vs H-input reduced-Omega

For NONLINEAR iron (saturable B-H) the CHOICE of formulation governs the
CONDITIONING of the nonlinear iteration, not just the DOF count.  Pick the
B-input (convex) form; it converges where the H-input reduced-Omega stalls.

## The conditioning split: B-input (convex) vs H-input (non-convex)

| input | constitutive | energy | conditioning |
|-------|--------------|--------|--------------|
| **B-input** nu(\|B\|) | H = nu(\|B\|) B | INT W(\|B\|), W **convex** | ★ well-conditioned at ANY mu |
| **H-input** mu(\|H\|) | B = mu(\|H\|) H | non-convex in H | ill-conditioned at high mu |

For an ordinary saturating material nu(q) is monotone, so the B-input
co-energy INT W(\|B\|) is **convex** -> the B-input fixed point / Newton is
well-posed at any permeability.  The H-input mu(\|H\|) form is not convex in H
and is ill-conditioned at high mu.

## The reduced-Omega trap (the obvious 3-D choice that fails on saturable iron)

The reduced scalar potential `H = H_s - grad Omega`, `mu(|H|)` Picard is the
natural 3-D magnetostatic choice (scalar unknown, Biot-Savart coil source, no
coil mesh -- see the `reduced` topic).  But for SATURABLE iron it is
ill-conditioned exactly where it matters: the saturation **knee** of an
iron-dominated circuit sits in the LOW-drive (unsaturated, `mu_r ~ mu_r0`)
regime, and there the reduced-Omega Picard **STALLS** -- the high-mu scalar
potential is poorly determined inside the iron, giving a spurious `|H|` ->
spurious `mu` -> oscillation.  The points that DO converge are all
post-saturation; the knee is in the non-converged regime (a catch-22).

## The cure: the reduced VECTOR potential, B-input

`B = B_s + curl A_r`:
- `B_s = mu0 H_s` = the coil's Biot-Savart field (`curl H_s = J`), from Radia /
  `CoilBuilder.to_radia()` -- **no meshed coil, no div-J** issue.
- `A_r in H(curl)` = the iron reaction.
- `nu(|B|)` saturable in iron, `nu0 = 1/mu0` in air.  Since
  `curl(nu0 B_s) = curl H_s = J`, the governing `curl(nu(|B|) B) = J` becomes
  the source-free-looking weak form

      INT nu(|B|) curl(A_r) . curl(v)  +  INT_iron (nu(|B|) - nu0) B_s . curl(v)  =  0

  whose ONLY source is the IRON-localised `(nu - nu0) B_s` term (the coil
  drives the iron through its known `B_s`).  A tiny `eps INT A_r . v`
  regularises the gradient gauge so a direct (sparsecholesky) solve works.
  Under-relaxed Picard (or Newton on the convex energy), continuation in the
  drive.

This is well-conditioned at ANY mu -- it converges in the low-drive high-mu
regime the reduced-Omega cannot reach.

## Code pattern (reduced-A, B-input, Picard)

```python
from ngsolve import *
import radia as rad
NU0 = 1.0 / mu0
def mur_B(Bm):                                  # Froehlich (or any monotone nu)
    return 1.0 + (mur0 - 1.0) / (1.0 + (Bm / Bk)**2)
Bs = GridFunction(VectorL2(mesh, order=1))
Bs.Set(mu0 * rad.RadiaField(coils, "h"))        # Biot-Savart, project ONCE
fes = HCurl(mesh, order=1, dirichlet="outer"); u, v = fes.TnT()
Ar = GridFunction(fes); eps = 1e-6 * NU0
for it in range(niter):                          # B-input Picard
    B = Bs + curl(Ar); Bm = sqrt(InnerProduct(B, B) + 1e-30)
    nu   = mesh.MaterialCF({"iron": NU0 / mur_B(Bm)},        default=NU0)
    nu_m = mesh.MaterialCF({"iron": NU0 / mur_B(Bm) - NU0},  default=0.0)
    a = BilinearForm(nu*InnerProduct(curl(u),curl(v))*dx + eps*InnerProduct(u,v)*dx)
    f = LinearForm(-nu_m*InnerProduct(Bs, curl(v))*dx)
    a.Assemble(); f.Assemble()
    An = GridFunction(fes)
    An.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    Ar.vec.data = (1.0 - relax)*Ar.vec + relax*An.vec        # under-relax
```

## Verified (lab)

`docs/clebsch_hodograph/demos/clebsch_dipole_saturation_3d.py` (golden
`test_clebsch_dipole_saturation_3d_aform`): on an H-frame dipole, at a LOW
drive (`NI=3 kA-t`, iron `<mu_r>~1200`) the B-input A-formulation converges to
`resid 8e-6` in **15 iters**, while the reduced-Omega `mu(|H|)` Picard STALLS
at `resid ~2e-2`.  Same geometry, same B-H knee: the difference is the
formulation's conditioning.

The 2-D linear-cost companion (the iron flux path as a saturable Chaplygin
guide -> a magnetic-circuit 1-shot, `B_gap(NI)` at linear cost) is
`docs/clebsch_hodograph/demos/clebsch_dipole_saturation_2d.py`.

## Cross-reference

- Same B-input convexity that fixes the **HDiv-VIM** nonlinear demag solve
  (the de Rham dual) -- see `radia_mcp.fem.equivalence_source`.
- For the reduced-potential CANCELLATION pitfall (a separate issue:
  `H_s - grad Omega ~ 0` losing precision in saturated regions), see the
  `reduced` topic above.
- Gauge regularisation of the curl-curl null space: `fem.gauge_open_boundary`.
"""


def get_potential_formulations_knowledge(topic: str = "catalog") -> str:
    """Dispatch potential formulation topics.

    Topics:
        catalog          - All formulations + decision tree (DEFAULT)
        a_omega          - A-Omega mixed vector-scalar
        t_omega          - T-Omega (electric vector + magnetic scalar)
        h_formulation    - H-formulation (HCurl direct) for SC / nonlinear
        reduced          - ★ Reduced potential (lab accelerator magnet core)
        nonlinear        - ★ B-input A-formulation cure (saturable iron conditioning)
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
    if topic in ("nonlinear", "b_input", "b-input", "conditioning", "saturable",
                 "a_reduced_nonlinear", "aform", "a_formulation"):
        return NONLINEAR_AFORM
    if topic in ("darwin", "darwin_model"):
        return DARWIN
    if topic == "all":
        return "\n\n".join([CATALOG, A_OMEGA, T_OMEGA, H_FORMULATION,
                            REDUCED_POTENTIAL, NONLINEAR_AFORM, DARWIN])
    return (f"Unknown topic '{topic}'. Available: catalog, a_omega, t_omega, "
            "h_formulation, reduced, nonlinear, darwin, all.")
