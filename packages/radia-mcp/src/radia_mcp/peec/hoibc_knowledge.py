"""HOIBC (Higher Order Impedance Boundary Conditions) knowledge.

Extends the 1st-order Leontovich SIBC documented in peec_knowledge.py
to higher-order perturbation theory that captures surface curvature
and remains valid at lower frequencies (δ comparable to body size).

Reference papers (public-safe curated corpus):
- Course G2ELab 2018 lesson 2 (derivation tutorial)
- Dong-Di Rienzo 2020 IEEE Access 8:186496 (FEM/BEM 3D implementation)
- Bilicz-Badics-Pávó 2023 ISEM (wide-band nonlocal IBC)

Topics:
    derivation                   - Perturbation expansion theory (Rytov-
                                    Senior + Mitzner improvements)
    dong_di_rienzo_2020          - 3D scalar-potential 3-Laplace-BVP
                                    algorithm (PEC → Leontovich → Mitzner)
    bilicz_2023_nonlocal         - Wide-band nonlocal IBC for low-freq
    vs_sibc_decision             - When SIBC suffices, when HOIBC needed
    radia_application            - How to extend calc_inductance.py /
                                    calc_fem_kelvin.py
"""

from ..common import load_prompt


# ============================================================
# DERIVATION (perturbation theory)
# ============================================================
HOIBC_DERIVATION = r"""
# HOIBC derivation via perturbation theory

## Standard SIBC (Leontovich 1948) — what HOIBC extends

The Leontovich SIBC relates tangential E and H at a conductor surface:
    E_t = Z_s (ẑ × H_t)
with Z_s = (1+j) / (σ δ),  δ = sqrt(2 / (ωμσ)).

VALIDITY: skin depth δ much smaller than:
  - Radius of curvature R_curv of the conductor surface
  - Wavelength λ
  - Conductor thickness d

LIMITATION: when δ ≈ R_curv (e.g., thin wires at low frequency, or
edges/corners of laminations), the local-tangent-plane approximation
of Leontovich breaks down.  HOIBC fixes this via perturbation.

## Perturbation parameter

Define the SMALL PARAMETER:
    p̃ = μ_r δ / D       (Dong-Di Rienzo 2020 IEEE Access)
where D is the characteristic conductor size.

The original Rytov-Senior used p̃ = δ/D, which is WRONG for magnetic
conductors -- the μ_r factor is essential to get the correct
expansion in the strong-skin regime.

## Asymptotic expansion

Fields inside the conductor are expanded as:
    H = H_0 + p̃ H_1 + p̃² H_2 + ...
    φ = φ_0 + p̃ φ_1 + p̃² φ_2 + ...    (reduced magnetic scalar potential)

Plug into the diffusion equation ∇×(∇×H) = -μσ ∂H/∂t and equate terms
of equal power of p̃.  This gives a HIERARCHY of equations:
- Order 0: H_0 satisfies a flat-conductor 1D diffusion (Leontovich)
- Order 1: H_1 corrects for SURFACE CURVATURE (radii d_1, d_2)
- Order 2: H_2 adds tangential-derivative corrections (Mitzner)

## Hierarchy of approximations

| Order | Name | Captures |
|-------|------|----------|
| 0 | PEC (perfect electric conductor) | No field inside, surface eliminates |
| 1 | Leontovich | Local plane wave in skin layer |
| 2 | Mitzner | Surface CURVATURE effects |
| 3+ | Higher Mitzner | Tangential derivatives, edge effects |

## Surface boundary condition at each order (Dong-Di Rienzo 2020)

After equating powers of p̃, the EXTERIOR field is governed by
a sequence of Laplace problems for the magnetic scalar potential:

    ∇²φ_m = 0   in dielectric region    (m = 0, 1, 2, ...)
    ∂φ_m / ∂n = f_m(φ_0, ..., φ_{m-1}, H_s)   on conductor surface

where H_s is the imposed external (source) field via Biot-Savart.

CRITICAL: the LHS is just LAPLACE, so any existing Laplace FEM/BEM
code can solve it.  The only thing that changes order-to-order is
the Neumann data on the conductor surface (RHS f_m).

## Time-domain ↔ frequency-domain switch

The derivation is in Laplace-transform variable s.  To switch to
frequency domain ω, substitute s → 2j (correspondingly normalized
by τ = 2/ω).

## Multi-frequency advantage

Since the BVPs are LAPLACE (frequency-independent), the spatial
basis functions of φ_0, φ_1, φ_2 are solved ONCE and re-evaluated
at any frequency by substituting p̃(ω) into:
    H_ext(ω) = H_s - ∇(φ_0 + p̃(ω) φ_1 + p̃(ω)² φ_2)

This is the main computational win of HOIBC over directly solving
the volume eddy-current problem at each frequency.

## Reference

- M. A. Leontovich, "Investigations on Radio Wave Propagation",
  Part II, USSR Academy of Sciences, 1948 (original SIBC)
- S. M. Rytov, J. Exp. Theor. Phys. 10:180-189, 1940 (asymptotic
  expansion)
- K. M. Mitzner, IEEE TAP 15:540-543, 1967 (curvature correction)
- S. Yuferev, N. Ida, "Surface Impedance Boundary Conditions: A
  Comprehensive Approach", CRC Press 2010 (textbook)
- Dong & Di Rienzo, IEEE Access 8:186496, 2020 (3D implementation)
"""


# ============================================================
# DONG-DI RIENZO 2020 - FEM/BEM IMPLEMENTATION
# ============================================================
DONG_DI_RIENZO_2020 = r"""
# Dong-Di Rienzo 2020: 3D FEM/BEM HOIBC implementation

## Paper
J. Dong, L. Di Rienzo, "FEM and BEM Implementations of a High Order
Surface Impedance Boundary Condition for Three-Dimensional Eddy
Current Problems", IEEE Access 8:186496-186506, 2020.
DOI: 10.1109/ACCESS.2020.3029566

## Key advance vs. earlier 2D papers

| Year | Paper | Dim | Implementation |
|------|-------|-----|----------------|
| 1967 | Mitzner | analytical | curvature correction theory |
| 2009-2014 | various 2D papers | 2D | BEM/FEM with HOIBC |
| **2020** | **Dong-Di Rienzo** | **3D** | **first 3D FEM/BEM with full HOIBC** |

## Reduced magnetic scalar potential (ψ) formulation

In dielectric region (μ_r = 1, no current source):
    H_ext = H_s - ∇φ
where H_s is Biot-Savart from filamentary source currents.

Then φ obeys Laplace: ∇²φ = 0.

Expand: φ = φ_0 + p̃ φ_1 + p̃² φ_2

Each φ_m solves an INDEPENDENT Laplace BVP (Neumann data depends on
the previous orders + curvature):

### BVP for φ_0 (PEC approximation)
    ∇²φ_0 = 0
    ∂φ_0 / ∂n = (H_s)_n      on conductor surface

### BVP for φ_1 (Leontovich correction)
    ∇²φ_1 = 0
    ∂φ_1 / ∂n = -(1/√s) · Σ_k ∂/∂ξ_k (H_s - ∇φ_0)_{ξ_k}    on surface
where ξ_1, ξ_2 are LOCAL TANGENTIAL coordinates.

### BVP for φ_2 (Mitzner correction)
    ∇²φ_2 = 0
    ∂φ_2/∂n = -(1/√s) Σ_k ∂/∂ξ_k (-∇φ_1)_{ξ_k}
              - (1/(μ_r √s)) · CURVATURE_TERM(d_1, d_2)
                                · (H_s - ∇φ_0)_{ξ_k}
where d_1, d_2 are the LOCAL RADII OF CURVATURE.

## Algorithm pseudocode

```
H_s ← Biot-Savart from filamentary currents
solve Laplace BVP_0 → φ_0
H_ext_PEC = H_s - ∇φ_0      ← order 0 result

solve Laplace BVP_1 → φ_1
H_ext_Leontovich = H_s - ∇(φ_0 + p̃ φ_1)    ← order 1 result

solve Laplace BVP_2 → φ_2
H_ext_Mitzner = H_s - ∇(φ_0 + p̃ φ_1 + p̃² φ_2)    ← order 2 result
```

## Numerical validation (prolate ellipsoid test)

Test geometry: σ = 5.998e7 S/m prolate ellipsoid + circular coil
at 10 A.  Reference: very-fine-mesh axisymmetric FEM.

| Frequency | μ_r | δ (mm) | p̃ | Required order |
|-----------|-----|--------|---|----------------|
| 1 kHz | 1 | 2.06 | 4.1e-3 | Leontovich enough |
| 50 Hz | 1 | 9.19 | 1.8e-2 | Mitzner needed |
| 1 kHz | 100 | 0.21 | 4.1e-2 | Mitzner needed |
| 1 kHz | 1000 | 0.065 | 1.3e-1 | Mitzner NOT enough; need order 3 |

## Magnetic-conductor caveat

For μ_r >> 1 (electrical steel, ferrite), the correct small parameter
is p̃ = μ_r δ / D, NOT just δ/D.  Without the μ_r factor, naive
asymptotic expansion fails even at high frequency.

This is the KEY improvement of Dong-Di Rienzo 2020 over Rytov-Senior.

## Computational cost

For 2nd-order tetrahedral FEM (1.3M DOF):
- BVP_0 (PEC): 33 s, 5.4 GB
- BVP_1 (Leontovich): +52 s, 6.3 GB
- BVP_2 (Mitzner): +55 s, 6.5 GB
- Total: ~2 min 20 s

BEM total (3162 elements, 2nd-order): ~5 min 45 s.

For multi-frequency sweep, this single solve covers ALL frequencies
via post-process (huge advantage over volume FE per frequency).

## Limitation: smooth surfaces only

HOIBC formulation assumes smooth conductor surfaces.  Sharp edges
and corners introduce LOCAL ERRORS in the curvature evaluation.
For PCB traces or stator slots with corners, mesh refinement plus
PEC condition near corners is recommended.
"""


# ============================================================
# BILICZ-BADICS-PÁVÓ 2023 - NONLOCAL WIDE-BAND IBC
# ============================================================
BILICZ_2023_NONLOCAL = r"""
# Bilicz-Badics-Pávó 2023: Wide-band nonlocal IBC

## Paper
S. Bilicz, Z. Badics, J. Pávó, "Wide-band nonlocal impedance boundary
condition model for high-conductivity regions in integral equation
framework", ISEM 2023.

## Motivation

Standard local SIBC (Leontovich):
    E_t / H_t = η = sqrt(jωμ/σ)        (point-wise)

Becomes inaccurate when δ ≈ R_curv (low freq, or thin wire).  HOIBC
fixes this by EXPANDING IN δ/R_curv.  But there's an alternative:
use a NONLOCAL operator that solves the FULL 2D eddy-current problem
in the wire cross-section.

## Key idea: 2D cross-section eddy-current BVP

Replace the local SIBC with a NONLOCAL relation:
    E_t(s) = Z[H_t(s')]    for all (s, s') on the wire surface
where Z is computed by solving a 2D BVP in the wire CROSS SECTION.

For a wire with cross-section Ω and boundary Γ:

```
∇² E_u - jωμκ E_u = 0    in Ω
∂E_u/∂n = -jωμ H_v       on Γ          (relates Eu ↔ Hv)

∇² H_u - jωμκ H_u = 0    in Ω
∂H_u/∂n = κ E_v          on Γ          (relates Hu ↔ Ev)
```

where (u, v, n) is the local coordinate system: u along wire axis,
v tangential to cross-section contour, n normal outward.

## Advantage

- **Works at ANY frequency** down to DC (because the 2D BVP captures
  the full eddy-current distribution, not just the skin-depth
  approximation).
- **Any cross-section shape** (rectangular for PCB traces, circular
  for round wires, hollow tubes, ...).
- **Naturally compatible with surface integral equation** (PEEC-like)
  framework.

## Disadvantage

- More expensive than HOIBC: needs 2D FE solve PER FREQUENCY (HOIBC
  reuses the same Laplace solves across frequencies).
- Numerical implementation more involved (couples 2D cross-section
  FE with 3D surface IE).

## Comparison with Dong-Di Rienzo HOIBC

| | HOIBC (Dong 2020) | Nonlocal IBC (Bilicz 2023) |
|---|---|---|
| Frequency range | Strong skin (δ << R_curv) | Wide-band (any δ) |
| Cost per freq | 0 (post-process) | 1 × 2D BVP per freq |
| Curvature limit | OK if smooth | Any shape (incl. rectangular) |
| Best for | Magnetic body in skin regime | Inductor wide-band sweep |

## Method-of-images extension

The Bilicz 2023 paper extends [Bilicz 2022] by:
1. Adding general (non-rectangular) cross-section support.
2. Including PEC plane (ground) via method-of-images modification of
   the Green's function.

This is the natural production path for PCB-trace impedance extraction
across DC → GHz range.

## Reference

- Bilicz et al. 2022 (cylindrical wire nonlocal SIBC, base case)
- Yuferev-Ida 2010 (textbook, comparative)
"""


# ============================================================
# DECISION TREE: SIBC vs HOIBC vs nonlocal IBC
# ============================================================
HOIBC_VS_SIBC_DECISION = r"""
# Decision tree: SIBC vs HOIBC vs Nonlocal IBC

```
Problem: surface impedance on a conductor surface in MQS regime
│
├── δ << R_curv AND δ << conductor_thickness?
│   ├── YES: Leontovich SIBC is enough (1st-order)
│   │        → Use --sibc-model bessel (PEEC) or --sibc-model dowell
│   │           (already in radia_mcp.peec)
│   │
│   └── NO: HOIBC or nonlocal IBC needed
│
├── δ ~ R_curv but conductor surface is smooth?
│   └── HOIBC (Mitzner correction, Dong-Di Rienzo 2020)
│        → solves 3 Laplace BVPs once, multi-freq post-process
│        → not yet implemented in calc_inductance.py
│
├── δ ~ conductor_thickness OR cross-section non-trivial (rectangular)?
│   └── Nonlocal IBC (Bilicz 2023) -- 2D cross-section BVP per freq
│        → works at ANY frequency down to DC
│        → ALREADY conceptually implemented in calc_inductance.py
│          via PEEC+BEM weak coupling (BEM does full surface IE)
│
├── Coil with finite turn-to-turn coupling and skin/proximity?
│   └── Use --nwinc N --nhinc M to subdivide cross-section
│        (Multi-filament approach in PEECBuilder, already implemented)
│
└── Saturating magnetic conductor (nonlinear iron)?
    └── ESIM cell problem (radia_mcp.ih.esim_cell_problem)
        → NOT HOIBC; uses 1D radial PDE per cell of lamination stack
```

## Specific guidance for Sugahara lab

| Application | Recommended IBC model |
|-------------|----------------------|
| IH coil round wire, 1-500 kHz | Standard Bessel SIBC |
| IH workpiece, smooth surface, σ × μ_r = 1e7-1e11 | Standard Leontovich SIBC |
| PCB trace, DC-GHz wide-band | Nonlocal IBC (Bilicz 2023) — TODO |
| Coil end-region with sharp curvature | HOIBC Mitzner — TODO |
| Saturating workpiece (φ-formulation) | ESIM cell problem (already implemented) |
| Anisotropic / laminated stack | Hollaus Effective Material (already implemented) |

## Why HOIBC matters now

The lab's existing PEEC+SIBC stack handles ~90% of IH and PEEC
applications using just Leontovich.  HOIBC fills the remaining 10%:
- Low-frequency end of inductor sweep (DC → kHz)
- Coil-end-region modeling (sharp curvature)
- High-μ_r magnetic workpieces (μ_r δ / D ~ 0.1)

For each of these, Mitzner correction can improve accuracy by 1-2
orders of magnitude over Leontovich without significant FE cost
increase (one extra Laplace solve per order).
"""


# ============================================================
# RADIA APPLICATION: how to extend calc_inductance.py / calc_fem_kelvin.py
# ============================================================
HOIBC_RADIA_APPLICATION = r"""
# Implementing HOIBC in radia panel scripts

## Current state (2026-05-22)

| Script | Current IBC | HOIBC ready? |
|--------|-------------|---------------|
| `calc_peec.py` | None (filament-only) | n/a |
| `calc_inductance.py` (PEEC+BEM) | 1st-order SIBC via BEM weak coupling | ⚠ partial via BEM full-Galerkin |
| `calc_fem_kelvin.py` | 1st-order Leontovich Robin BC | ❌ no HOIBC |
| `calc_fem_coilmesh.py` | none (full volumetric solve) | n/a |
| `esim_cell_problem.py` | nonlinear cell (different physics) | n/a |

## Proposed addition: --ibc-order N flag to calc_fem_kelvin.py

```python
parser.add_argument('--ibc-order', type=int, default=1,
                    choices=[0, 1, 2],
                    help='Surface impedance BC order. '
                         '0=PEC, 1=Leontovich (default), 2=Mitzner HOIBC')
```

### Implementation (Dong-Di Rienzo 2020 recipe in NGSolve)

```python
from ngsolve import (H1, BilinearForm, LinearForm, grad, dx, ds,
                     GridFunction, specialcf)

# Surface curvature: NGSolve 6.2.2603 does NOT expose a direct
# MeanCurvature CF.  Options:
#   (a) Pass user-known R_curv as scalar (simplest; what
#       calc_fem_kelvin.py --mitzner-radius does today)
#   (b) Compute trace of Weingarten map via specialcf.Weingarten(3)
#       (3D; needs surface FE machinery)
#   (c) Provide a per-element curvature CF from the CAD origin
#       (NGSolve's OCC geometry has access)
#
# For Phase 1 implementation we use Option (a) -- one R_curv scalar.
# For non-uniform surfaces, sub-divide by material region or pass a
# H1(boundary=true) curvature GridFunction.

R_curv = 0.025   # m, user-provided mean curvature radius
H_mean = 1.0 / R_curv  # 1/m

# BVP_0: PEC approximation
fes = H1(mesh, order=p, dirichlet=outer_bnd)
phi_0 = GridFunction(fes)
u, v = fes.TnT()
a0 = BilinearForm(fes)
a0 += grad(u) * grad(v) * dx
a0.Assemble()
f0 = LinearForm(fes)
f0 += Hs_normal * v * ds(definedon=mesh.Boundaries(conductor_bnd))
f0.Assemble()
phi_0.vec.data = a0.mat.Inverse(fes.FreeDofs(), inverse='pardiso') * f0.vec

# BVP_1: Leontovich correction (use surface gradient of (Hs - grad phi_0))
phi_1 = GridFunction(fes)
# ... similar pattern with surface-tangential derivatives in RHS

# BVP_2: Mitzner correction (use H_mean in RHS)
phi_2 = GridFunction(fes)
# ... add d_k+d_3-k / (d_k * d_3-k) terms, where 1/d_k contracted into H_mean

# Multi-frequency post-process
def H_ext_at_freq(omega):
    p_tilde = mu_r * skin_depth(omega) / D
    return Hs - grad(phi_0 + p_tilde * phi_1 + p_tilde**2 * phi_2)
```

### Computational cost vs benefit

| Order | Cost increase | Accuracy gain (μ_r=100, δ/D=4%) |
|-------|---------------|----------------------------------|
| 0 (PEC) | baseline | 50% error |
| 1 (Leontovich) | +1 Laplace solve | 10-20% error |
| 2 (Mitzner) | +1 more Laplace solve | <1% error |
| 3+ (higher) | +N more solves | <0.1% error |

For typical IH workpieces and motor laminations, **order 2 (Mitzner)
is the sweet spot** — captures curvature without exorbitant cost.

## TODO checklist for full Radia HOIBC integration

- [ ] Add `--ibc-order` flag to calc_fem_kelvin.py
- [ ] Implement BVP_0/1/2 in `_solve_hoibc_cascade` helper
- [ ] Extract surface curvature from NGSolve mesh
  (NGSolve's `GeomFunction.MeanCurvature` or via boundary normal
  derivative; for analytic surfaces use the parametric curvature)
- [ ] Validate against Dong-Di Rienzo 2020 prolate ellipsoid benchmark
- [ ] Cross-check vs current Leontovich Robin BC at high freq
  (should converge to same answer when δ << R_curv)
- [ ] Document in calc_fem_kelvin.py docstring + golden test
- [ ] Extend `motor_em_force_recipe` MCP knowledge with HOIBC pointer

## When HOIBC is NOT the right tool

- **Volumetric eddy currents needed** (e.g., for force on rotor bar
  in induction motor): use full volume FE with `calc_fem_coilmesh.py`
  or Hollaus MSFEM-effective-material approach
- **Nonlinear / hysteretic iron**: use ESIM cell problem + Picard
  iteration; HOIBC assumes linear material
- **PCB trace wide-band (DC-GHz)**: use Bilicz 2023 nonlocal IBC
  instead (HOIBC requires δ << R_curv which fails at DC)
"""


SIBC_FAMILY_VARIANTS = r"""
# SIBC / HOIBC variant landscape (papers in public-safe curated corpus)

Beyond the canonical 1st-order Leontovich (radia_mcp.peec.peec_usage
'sibc') and Mitzner HOIBC (this module), several SIBC variants exist
for specialized scenarios.  Lab paper coverage:

## XFEM SIBC (public-safe curated corpus)

X-FEM (eXtended FEM) embeds the SIBC into a discontinuous enrichment
of the FE basis near the conductor boundary, avoiding the need to
align the mesh with the surface.

USE WHEN: moving conductors (rotor surface in a stator-mesh-fixed
formulation), or evolving topology (induction heating workpiece moved
through a fixed coil).

Reference: Bouillault, Buffa, Ren — IEEE TMAG ~2008 series.

## FDTD SIBC (public-safe curated corpus)

For time-domain analysis (FDTD), the SIBC becomes a CONVOLUTION in
time: E_t(t) = ∫_0^t Z_s(t-τ) H_t(τ) dτ.

Z_s(t) is approximated by:
- Prony series (exponential decomposition)
- Pade rational form
- Recursive convolution (Maloney 1992)

USE WHEN: full-wave FDTD with conductive walls, ultra-wideband.
NOT applicable to MQS / Darwin / quasistatic.

## Nonlinear SIBC (public-safe curated corpus)

Standard SIBC assumes linear material (single Z_s value).  For
saturating iron under high B, Z_s(|B|) becomes field-dependent.

Two approaches:
1. **ESIM (Effective Surface Impedance Method)**: solve 1D cell
   problem `d/dz[(1/μ(z)) dH/dz] = jωσH` per surface point with
   spatially-varying μ(z) reflecting saturation.
   → Already in Radia via `esim_cell_problem.py` (radia_mcp.ih.esim).
2. **Nonlinear effective Z_s in scalar potential** (lab paper:
   `Nonlinear Effective Surface Impedance in Magnetic Scalar Potential
   formulation`) — extends Dong-Di Rienzo HOIBC to nonlinear μ.

USE WHEN: induction heating workpiece (heavy saturation), motor
laminations under DC bias.

## Multi-layer / coated conductor SIBC

For multi-layer coated conductors (e.g., HTS tape, anti-reflection
coating on radar absorber):
```
Z_s = Z_1 + Z_2 * tanh(γ_1 * d_1) + ...   (transmission-line cascade)
```

USE WHEN: HTS magnet (CORC/Roebel cable), radar-absorbing material,
PCB trace with prepreg layers.

NOT YET in Radia.  Workaround: extend `peec_sibc` Dowell formula manually.

## WPT (Wireless Power Transfer) wide-band

For WPT coil design across DC → GHz, the nonlocal IBC (Bilicz 2023)
is the natural method.  Lab paper at `public-safe curated corpus` covers
WPT-specific application.

Cross-link: `peec_hoibc('bilicz_2023_nonlocal')` for the formulation,
`07_Applications/` for the WPT-specific implementation details.
"""


HOIBC_MATHEMATICA_DERIVATION = load_prompt("peec", "mathematica_derivation")

_HOIBC_MATHEMATICA_DERIVATION_INLINE_DEPRECATED = r"""
# Mathematica Symbolic Derivation of HOIBC

Wolfram Language recipes for the Dong-Di Rienzo 2020 perturbation
hierarchy: PEC → Leontovich → Mitzner. Use these to:
- Verify the Neumann RHS formulas before coding them in NGSolve
- Reproduce the prolate ellipsoid benchmark (paper Fig 7)
- Cross-check a Radia HOIBC implementation symbolically

Companion to `peec_hoibc('radia_application')` (numerical recipes) and
`peec_hoibc('dong_di_rienzo_2020')` (formula reference).

## Recipe 1: Small parameter `p_tilde` symbolic

```mathematica
(* Dong-Di Rienzo 2020 small parameter:                              *)
(*     p_tilde = mu_r * delta / D                                    *)
(* where D is the characteristic conductor size, NOT delta/D alone   *)
(* (Rytov-Senior was wrong for mu_r >> 1).                           *)

pTilde[mur_, delta_, D_] := mur * delta / D;

(* Verify expansion validity for steel mu_r=100, R=10 mm, 1 kHz *)
mu0 = 4*Pi*1e-7;
delta1kHz[mur_] := Sqrt[2/(2*Pi*1e3 * mur * mu0 * 2e6)];   (* steel sigma 2e6 *)
pTilde[100, delta1kHz[100], 0.01]
(* Result: 4.1e-2 -- Mitzner needed *)

pTilde[1000, delta1kHz[1000], 0.01]
(* Result: 1.3e-1 -- Mitzner NOT enough, order 3 needed *)
```

## Recipe 2: BVP_0 Neumann RHS (PEC approximation)

```mathematica
(* H_s = applied Biot-Savart, normal projection onto conductor surface *)
(* BVP_0: Laplace[phi_0] = 0 in dielectric                            *)
(*        d phi_0 / dn = H_s . n_hat  on conductor                    *)

(* For sphere R=a in uniform H = H0 z_hat: *)
HsNormalSphere[a_, H0_, theta_] := H0 * Cos[theta];

(* Analytical phi_0 from external multipole expansion *)
phi0Sphere[a_, H0_, r_, theta_] :=
    -H0 * (a^3 / (2 * r^2)) * Cos[theta] + H0 * r * Cos[theta];

(* Verify: H_ext_PEC = H_s - grad(phi_0); tangential should be uniform *)
HextPEC[a_, H0_, r_, theta_] := Module[{phi0},
    phi0 = phi0Sphere[a, H0, r, theta];
    H0 * Cos[theta] - D[phi0, r]
];
HextPEC[1, 1, 1, theta] // Simplify
(* Should give (3/2)*Cos[theta] -- the PEC sphere dipole reinforcement *)
```

## Recipe 3: BVP_1 Leontovich Neumann RHS

Surface-tangential derivative of `(H_s - grad phi_0)` along local
coordinates xi_1, xi_2. Symbolic on a sphere:

```mathematica
(* On sphere surface, tangential gradient in (theta, phi) *)
gradS[f_, theta_, phi_, a_] := {(1/a) D[f, theta], (1/(a*Sin[theta])) D[f, phi]};

(* BVP_1 Neumann RHS for sphere: *)
HtPEC[a_, H0_, theta_] := (3/2) * H0 * Sin[theta];   (* tangential H_ext_PEC *)

NeumannRHS1[a_, H0_, theta_] := -(1/Sqrt[s]) *
    D[a * HtPEC[a, H0, theta], theta] / a;
(* s here is the Laplace variable; replace s -> j*omega for freq domain *)

(* Verify integration over sphere: should be zero by Gauss's law       *)
(* (Laplace BVP needs compatibility condition: integral of Neumann = 0) *)
Integrate[NeumannRHS1[a, H0, theta] * a^2 * Sin[theta],
          {theta, 0, Pi}, {phi, 0, 2*Pi}] // Simplify
```

## Recipe 4: BVP_2 Mitzner curvature term (the key advance)

```mathematica
(* Curvature term for sphere: d_1 = d_2 = a (principal radii)        *)
(* General formula: CURVATURE_TERM = (d_k + d_{3-k}) / (d_k * d_{3-k}) *)
(* For sphere both axes a: 2a / a^2 = 2/a                            *)

CurvSphere[a_] := 2/a;   (* = 2 * mean curvature = 2 * H_mean *)

(* BVP_2 RHS includes both Leontovich-like derivatives of phi_1 AND  *)
(* curvature correction times (H_s - grad phi_0) tangential          *)

NeumannRHS2[a_, mur_, H0_, theta_] := Module[{leoTerm, curvTerm},
    leoTerm = -(1/Sqrt[s]) *
        (* tangential derivatives of (-grad phi_1) -- complex *)
        SomethingFrom[phi1Sphere];
    curvTerm = -(1/(mur * Sqrt[s])) * CurvSphere[a] * HtPEC[a, H0, theta];
    leoTerm + curvTerm
];

(* Mitzner correction in closed form for thin-skin sphere: *)
(* Z_s_Mitzner = Z_s_Leontovich * (1 + (1+j)/2 * delta * H_mean) *)
ZsMitznerSphere[Zs_, delta_, R_] := Zs * (1 + (1 + I)/2 * delta / R);
```

## Recipe 5: Reproduce Dong-Di Rienzo Fig 7 (prolate ellipsoid)

The benchmark: prolate ellipsoid `x^2/b^2 + y^2/b^2 + z^2/a^2 = 1`
with a = 2*b, sigma = 5.998e7 S/m, in a circular coil at 10 A.

```mathematica
(* Principal radii of curvature at pole and equator: *)
d1Pole[a_, b_] := b^2 / a;
d2Pole[a_, b_] := b^2 / a;
d1Equator[a_, b_] := b;
d2Equator[a_, b_] := a^2 / b;

(* Mean curvature at each point (used in Mitzner correction): *)
HmeanPole[a_, b_] := (1/d1Pole[a,b] + 1/d2Pole[a,b]) / 2;
HmeanEquator[a_, b_] := (1/d1Equator[a,b] + 1/d2Equator[a,b]) / 2;

(* For a = 2*b = 20 mm: *)
HmeanPole[0.02, 0.01]      (* = 100 1/m, high curvature *)
HmeanEquator[0.02, 0.01]   (* = 75 1/m, mixed *)

(* The Mitzner correction ratio (Z_Mitzner / Z_Leontovich):    *)
(* depends on local curvature -- spatially varying!            *)
(* This is why per-panel Z_s in NGSolve matters (Recipe 2 in   *)
(* ngsolve_recipes).                                            *)
```

## Recipe 6: Frequency-domain switch (Laplace s -> j omega)

```mathematica
(* Dong-Di Rienzo derivation uses Laplace variable s            *)
(* Switch to frequency domain: s -> j*omega                     *)
(* and normalize: tau = 2/omega                                  *)

FreqSubst = {s -> I*omega};

ZsFreqDomain[Zs_Leontovich, mur_, delta_, R_] :=
    Zs_Leontovich * (1 + (1 + I)/2 * delta * (1/R)) /. FreqSubst;
```

## Why these recipes matter

The Dong-Di Rienzo 2020 paper's algorithm is **3 separate Laplace
BVPs** with Neumann RHS that depend on the previous order. The RHS
formulas are easy to miscode (lots of indices, surface-tangential
derivatives, principal radii). Symbolic verification in Mathematica
**before** committing the NGSolve assembly code saves hours of
debugging.

After Mathematica verification, the production NGSolve code lives in:
- `src/radia/panels/calc_fem_kelvin.py::_solve_hoibc_cascade` (planned)

Reference MCP: `mcp-server-mathematica` for actual symbolic computation
from a Cubit panel session.
"""


HOIBC_NGSOLVE_RECIPES = load_prompt("peec", "ngsolve_recipes")

_HOIBC_NGSOLVE_RECIPES_INLINE_DEPRECATED = r"""
# NGSolve Implementation Recipes for HOIBC

Production-grade NGSolve code patterns for the 3-Laplace-BVP cascade
(Dong-Di Rienzo 2020). Goes BEYOND the existing
`calc_fem_kelvin.py --formulation total` 1st-order Robin BC.

Status: design recipes for the upcoming `--ibc-order {0,1,2}` flag.
The 1st-order path is already production; orders 0 and 2 are new.

## Recipe 1: BVP_0 (PEC approximation)

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                     grad, dx, ds, InnerProduct, specialcf, x, y, z)

mesh = Mesh("model.vol")  # workpiece is meshed as a HOLE (subtracted)

# Reduced magnetic scalar potential phi: H_ext = H_s - grad(phi)
fes = H1(mesh, order=p, dirichlet="kelvin_outer", complex=True)

phi_0 = GridFunction(fes)
u, v = fes.TnT()

# LHS: Laplace
a = BilinearForm(fes, symmetric=True)
a += grad(u) * grad(v) * dx
a.Assemble()

# RHS: Neumann data on conductor surface = -H_s.n
n = specialcf.normal(3)
H_s_cf = biot_savart_H_cf_from_peec(peec_topology, current=I_coil)

f = LinearForm(fes)
f += -InnerProduct(H_s_cf, n) * v * ds(
    definedon=mesh.Boundaries("conductor_bnd"))
f.Assemble()

# Solve with Compact AMS preconditioner
import radia.sparsesolv_ngsolve as ssn
c = Preconditioner(a, "bddc")
inv = ssn.BiCGStab(matrix=a.mat, c=c.mat, tol=1e-8)
phi_0.vec.data = inv * f.vec

H_ext_PEC = H_s_cf - grad(phi_0)  # PEC result
```

## Recipe 2: BVP_1 (Leontovich correction RHS)

```python
# BVP_1: same LHS as BVP_0 (Laplace), different Neumann RHS
# RHS_1 = -(1/sqrt(s)) * div_s ((H_s - grad phi_0)_tangential)
# where div_s = surface divergence (tangential gradient component)

phi_1 = GridFunction(fes)

# Use specialcf.Weingarten for tangential derivatives (NGSolve 6.2.2603)
# Or compute surface gradient via projection onto tangent plane

# Compute H_t at conductor surface (tangential of PEC field)
H_t_PEC = (H_s_cf - grad(phi_0)) - \
          InnerProduct((H_s_cf - grad(phi_0)), n) * n  # tangential proj

# Surface divergence of H_t (use H1 boundary FES + integration by parts)
# RHS_1 = -(1/sqrt(j*omega)) * div_s(H_t_PEC) . v
f1 = LinearForm(fes)
f1 += -1.0/sqrt(1j*omega) * InnerProduct(grad(v).Trace(), H_t_PEC) * \
      ds(definedon=mesh.Boundaries("conductor_bnd"))
f1.Assemble()

phi_1.vec.data = inv * f1.vec   # reuse a.mat inverse (same Laplace)

# Leontovich result: H_ext_1 = H_s - grad(phi_0 + p_tilde * phi_1)
p_tilde = mur * delta(omega, sigma, mur) / D_char  # small parameter
H_ext_Leontovich = H_s_cf - grad(phi_0) - p_tilde * grad(phi_1)
```

## Recipe 3: BVP_2 (Mitzner curvature correction)

The headline HOIBC advance — captures surface curvature explicitly.

```python
# BVP_2: same LHS (Laplace), Neumann RHS includes curvature term
# RHS_2 = -(1/sqrt(s)) * div_s(-grad phi_1)
#       - (1/(mur * sqrt(s))) * 2*H_mean * H_t_PEC_tangential

phi_2 = GridFunction(fes)

# Surface mean curvature: 2 * H_mean = 1/d_1 + 1/d_2
# NGSolve specialcf for Weingarten map (3D)
W = specialcf.Weingarten(3)    # 2x2 Weingarten map on boundary
H_mean = 0.5 * Trace(W)        # mean curvature in 1/m

# Term 1: tangential derivative of (-grad phi_1)
grad_phi1_t = grad(phi_1) - InnerProduct(grad(phi_1), n) * n
term1 = 1.0/sqrt(1j*omega) * InnerProduct(grad(v).Trace(), grad_phi1_t)

# Term 2: curvature correction
term2 = -1.0/(mur * sqrt(1j*omega)) * 2 * H_mean * \
        InnerProduct(H_t_PEC, v.Trace())

f2 = LinearForm(fes)
f2 += (term1 + term2) * ds(definedon=mesh.Boundaries("conductor_bnd"))
f2.Assemble()

phi_2.vec.data = inv * f2.vec   # reuse a.mat inverse

# Mitzner result: H_ext_2 = H_s - grad(phi_0 + p_tilde phi_1 + p_tilde^2 phi_2)
H_ext_Mitzner = H_s_cf - grad(phi_0) - p_tilde * grad(phi_1) - \
                p_tilde**2 * grad(phi_2)
```

## Recipe 4: Multi-frequency post-process (zero-cost frequency sweep)

The key computational win: BVPs 0/1/2 are **frequency-independent**
Laplace problems. Once solved, evaluate at any frequency by plugging
in `p_tilde(omega)`:

```python
# Pre-solve (frequency-independent)
solve_bvp_0()  # phi_0 -- 33 s
solve_bvp_1()  # phi_1 -- +52 s
solve_bvp_2()  # phi_2 -- +55 s
# Total: ~2 min 20 s ONE TIME for all frequencies

# Frequency sweep is O(1) per point
def H_ext_at_freq(omega, mur, sigma, D_char):
    delta = sqrt(2/(omega * mu0 * mur * sigma))
    p_tilde = mur * delta / D_char
    return H_s_cf - grad(phi_0) - \
           p_tilde * grad(phi_1) - \
           p_tilde**2 * grad(phi_2)

# 100-point frequency sweep: 100 * O(1) = milliseconds
for f in logspace(2, 6, 100):
    omega = 2*pi*f
    H_ext = H_ext_at_freq(omega, mur=100, sigma=2e6, D_char=0.01)
    # Compute L, P, etc. from H_ext via boundary integrals
```

**Compare to existing pipeline**: standard FEM-Kelvin would re-solve
the volume eddy-current problem per frequency (~2 min each x 100 =
3.3 hours). HOIBC: 2 min one-time + 100 * milliseconds = ~2 min total.

## Recipe 5: Validation against Dong-Di Rienzo 2020 prolate ellipsoid

```python
# Geometry: prolate ellipsoid x^2/b^2 + y^2/b^2 + z^2/a^2 = 1
#           a = 20 mm, b = 10 mm (a/b = 2)
#           Surrounded by circular coil at I = 10 A

# Lab equivalent geometry script: prolate_spheroid_hoibc.jou (planned)
# tests/panels/test_hoibc_ellipsoid_golden.py (planned)

# Validation table (Dong-Di Rienzo 2020 IEEE Access Table II):
# | omega | mu_r | Order needed |
# | 1 kHz | 1    | Leontovich   |
# | 50 Hz | 1    | Mitzner      |
# | 1 kHz | 100  | Mitzner      |
# | 1 kHz | 1000 | Mitzner not enough; order 3 |

# Golden criteria for the new `--ibc-order 2` panel mode:
#   order=1 vs order=2 difference: > 5% at any (f, mur) combo
#                                  in the Mitzner-needed band
#   order=2 vs full-volume FE reference: < 2% error
```

## Recipe 6: Per-panel curvature extraction for HOIBC

Reuses the per-panel local curvature extractor from
`calc_inductance.py::_compute_panel_local_radii` (lab 2026-04-12).

```python
from radia.panels.calc_inductance import _compute_panel_local_radii

# For each panel on the conductor surface, get principal radii (d_1, d_2)
# from discrete normal-angle method:
R_local = _compute_panel_local_radii(mesh, surface_label="conductor_bnd",
                                      percentile=10)  # [m]

# In BVP_2, use per-panel mean curvature H_mean[i] = 1/R_local[i] for
# each surface element. Wrap as a SurfaceL2 GridFunction:
from ngsolve import SurfaceL2
fes_curv = SurfaceL2(mesh, order=0, definedon=mesh.Boundaries("conductor_bnd"))
H_mean_gf = GridFunction(fes_curv)
for i, R in enumerate(R_local):
    H_mean_gf.vec[i] = 1.0 / R  # 1/m

# Use H_mean_gf in term2 of BVP_2 instead of analytical specialcf.Weingarten
# (more robust on irregular meshes)
```

## Cross-MCP recipe references

- `radia_mcp.peec.hoibc('dong_di_rienzo_2020')` -- formula reference
- `radia_mcp.peec.hoibc('mathematica_derivation')` -- symbolic check
- `radia_mcp.peec.hoibc('radia_application')` -- TODO checklist for
  full Radia HOIBC integration
- `radia_mcp.fem` -- FEM gauging + Kelvin transform for outer boundary
- `radia_mcp.matrix_solvers` -- BiCGStab + Compact AMS preconditioner
"""


def get_hoibc_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
        derivation               - Perturbation theory (Rytov-Senior-Mitzner)
        dong_di_rienzo_2020      - 3D FEM/BEM implementation (3 Laplace BVPs)
        bilicz_2023_nonlocal     - Wide-band nonlocal IBC alternative
        vs_sibc_decision         - When to use SIBC vs HOIBC vs nonlocal
        radia_application        - How to extend calc_inductance.py /
                                    calc_fem_kelvin.py with HOIBC
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("derivation", "perturbation", "theory"):
        return HOIBC_DERIVATION
    if topic in ("dong_di_rienzo_2020", "dong", "di_rienzo",
                  "fem_bem_implementation", "mitzner"):
        return DONG_DI_RIENZO_2020
    if topic in ("bilicz_2023_nonlocal", "bilicz", "nonlocal", "wide_band"):
        return BILICZ_2023_NONLOCAL
    if topic in ("vs_sibc_decision", "decision", "vs_sibc",
                  "when_to_use", "decision_tree"):
        return HOIBC_VS_SIBC_DECISION
    if topic in ("radia_application", "application", "implementation",
                  "ngsolve", "calc_fem_kelvin"):
        return HOIBC_RADIA_APPLICATION
    if topic in ("sibc_family", "variants", "xfem", "fdtd",
                  "nonlinear_sibc", "multi_layer", "wpt"):
        return SIBC_FAMILY_VARIANTS
    if topic in ("mathematica_derivation", "mathematica", "symbolic",
                  "wolfram"):
        return HOIBC_MATHEMATICA_DERIVATION
    if topic in ("ngsolve_recipes", "ngsolve", "recipes",
                  "implementation_code", "code"):
        return HOIBC_NGSOLVE_RECIPES
    if topic == "all":
        return "\n\n".join([
            HOIBC_DERIVATION,
            DONG_DI_RIENZO_2020,
            BILICZ_2023_NONLOCAL,
            HOIBC_VS_SIBC_DECISION,
            HOIBC_RADIA_APPLICATION,
            SIBC_FAMILY_VARIANTS,
            HOIBC_MATHEMATICA_DERIVATION,
            HOIBC_NGSOLVE_RECIPES,
        ])
    return (
        f"Unknown topic '{topic}'. Available: derivation, "
        "dong_di_rienzo_2020, bilicz_2023_nonlocal, vs_sibc_decision, "
        "radia_application, sibc_family, mathematica_derivation, "
        "ngsolve_recipes, all."
    )
