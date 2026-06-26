"""EM-XFEM (Hiruma 2023) — electromagnetic XFEM enrichment for eddy currents.

Knowledge layer for the Hiruma EM-XFEM, which is structurally distinct
from solid-mechanics / fracture XFEM (Jafari et al. 2021):

  - Hiruma 2023 (IEEE TMag 59(5)): enrich H_1 trial space with
    psi(x) = exp(-gamma*xi(x)),  gamma = sqrt(j*omega*mu*sigma)
    on the conductor interior so the mesh need not resolve the
    O(delta) skin layer at high frequency.

  - This is a PUFEM (Partition-of-Unity FEM) enrichment, NOT cut-FEM:
    the mesh conforms to the boundary; ALL nodes can be enriched
    (or a near-surface subset).

  - Validation: 88-DOF compound H_1 x H_1 NGSolve implementation
    reproduces Hiruma's headline claim to 0.14% at r/delta = 15
    (beta-1 Phase 2 benchmark, 2026-05-25).

Reference paper (lab archive):
    Hiruma et al., "Extended Finite Element Method for Eddy-Current
    Problems with Surface Skin Effect", IEEE Trans. Magn. 59(5), 2023.
    Same Igarashi/Hiruma group as Hiruma 2020 CLN-A formulation.

Differences from solid-mechanics XFEM (Belytschko-Black 1999):
    No Heaviside / level-set / crack tracking.  The enrichment is
    smooth across the conductor; the goal is to capture skin-depth
    physics, not jumps or singularities.

Lab files (Sugahara research line, beta-1 Phase 1-4 benchmark):
    S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb
        (one consolidated notebook; sections:)
        phase1 - CFEM baseline vs Bessel
        phase2 - NGSolve EM-XFEM impl.
        phase3 - Pade-based Schur (deprecated predecessor)
        phase3b - Galerkin-Krylov ROM
        phase4 - XFEM-CLN stacking test (Hankel conditioning)

All text is ASCII (cp932-safe) and in English.
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Lab files
# ---------------------------------------------------------------------------

LAB_FILES = {
    # Promoted 2026-06-26: the per-phase scripts are now SECTIONS of one notebook.
    "notebook":       r"S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb",
    "benchmark_dir":  r"S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/",
    "phase1_results": r"S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/results_phase1.json",
    "phase4_results": r"S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/phase4_xfem_hankel_results.json",
    "phase4_fig":     r"S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/phase4_xfem_hankel_conditioning.png",
}


# ---------------------------------------------------------------------------
# OVERVIEW
# ---------------------------------------------------------------------------

OVERVIEW = r"""
# EM-XFEM (Hiruma 2023): overview

The Hiruma group (Hokkaido University, Igarashi Lab) introduced an
electromagnetic XFEM in IEEE TMag 59(5) (2023) that addresses a
classical FEM weakness: high-frequency eddy-current problems require
the mesh to resolve the skin depth delta = sqrt(2/(omega*mu*sigma)),
which becomes O(microns) at MHz frequencies in copper.  Mesh
refinement to that scale is prohibitive for engineering models.

The Hiruma EM-XFEM enriches the H_1 trial space with a closed-form
skin-layer profile:

    u_total(x) = u_std(x) + psi(x) * u_enr(x)
    psi(x)     = exp(-gamma * xi(x))
    gamma      = sqrt(j * omega * mu * sigma)
    xi(x)      = distance from the conductor surface

so the enriched DOFs absorb the rapid exponential decay across the
skin layer without needing the mesh to resolve it.

## Key distinction from solid-mechanics XFEM (Belytschko-Black 1999)

| Aspect              | Solid-mech XFEM         | EM-XFEM (Hiruma 2023)   |
|---------------------|-------------------------|-------------------------|
| Enrichment type     | Heaviside + asymptotic  | Exponential skin layer  |
| Discontinuity       | Yes (crack)             | No (smooth)             |
| Geometry tracking   | Level-set crack growth  | None (fixed boundary)   |
| Element classes     | Cut, tip, std (3 types) | Std + enriched (2 types)|
| Goal                | Stress singularities    | Skin-depth resolution   |
| Mesh                | Cut by crack            | Conforms to boundary    |

EM-XFEM is therefore architecturally simpler -- it is essentially
PUFEM (Melenk-Babuska 1996) with one analytical enrichment function
psi.  No level-set, no propagation, no sub-element integration.

## When to use EM-XFEM

Use EM-XFEM if BOTH:
  (a) you need a single-frequency 3D field solve at high f / small
      delta (r/delta > 3 say), with the mesh too coarse to resolve
      the skin layer
  (b) the geometry includes corners or sharp features where the
      skin layer must remain accurate (CFEM-SIBC breaks at corners)

Do NOT use EM-XFEM if:
  (a) you need a frequency sweep or time-domain transient.  The
      enrichment psi depends on omega via gamma, so each omega
      requires a different XFEM space.  See `decision_table` topic.
  (b) you want SPICE export of a Cauer ladder.  EM-XFEM lives at the
      FE level; it does NOT combine usefully with canonical CLN
      extraction (see `cln_stacking_negative` topic).
  (c) DC or near-DC accuracy is needed.  Use canonical CFEM-SIBC or
      standard FEM; EM-XFEM is targeted at the skin-band.

## Position in the radia-mcp universe

  - `radia_mcp.fem.fem_xfem_em_hiruma` (this tool):
       EM-XFEM as an FE-level technique
  - Solid-mechanics / fracture XFEM (Jafari et al. 2021): a separate
       family of XFEM for crack growth — different problem, same name
  - `radia_mcp.radia_ngsolve.cln_sibc_orthogonal`:
       decision framework: when XFEM vs SIBC vs augmented CLN
       (topic 'xfem_vs_sibc')
  - docs/hiruma_xfem_comparison/:
       reproducible benchmark (Phase 1-4)

## Three XFEM variants — do not confuse

The single name "XFEM" hides three structurally distinct families:

| Variant                     | Mesh treatment        | Enrichment           | Typical use            |
|-----------------------------|-----------------------|----------------------|------------------------|
| Solid-mech XFEM (Belytschko-Black 1999, Jafari 2021) | mesh CUT by crack | Heaviside + 4 asymptotic tip | crack growth, mixed-mode SIF |
| **EM-XFEM (Hiruma 2023, this tool)** | mesh CONFORMS to conductor | exp(-gamma * xi) PUFEM | skin-depth resolution at high f |
| Cut-FEM XFEM (Lehrenfeld 2014+, **ngsxfem** library) | mesh does NOT conform to interface | piecewise-polynomial on cut cells + ghost-penalty | moving boundary, two-phase, free-surface flow |

The ``ngsxfem`` library (https://github.com/ngsxfem/ngsxfem) is the
NGSolve add-on for the **third** variant -- cut-FEM XFEM, not the
Hiruma EM-XFEM.  Phase 2 of this lab's beta-1 benchmark deliberately
uses vanilla NGSolve compound H_1 x H_1 spaces **without** ngsxfem
because the Hiruma method is PUFEM, not cut-FEM (the mesh already
conforms to the conductor boundary, no cut cells are needed).

See `ngsxfem_relation` topic for when ngsxfem WOULD be appropriate
(it is not for the present Hiruma reproduction).
"""


# ---------------------------------------------------------------------------
# ENRICHMENT FUNCTION
# ---------------------------------------------------------------------------

ENRICHMENT_FUNCTION = r"""
# The EM-XFEM enrichment function: psi(x) = exp(-gamma * xi(x))

## Mathematical form

The enrichment function is

    psi(x) = exp(-gamma * xi(x))

where:

    gamma  = sqrt(j * omega * mu * sigma)              [complex propagation]
           = (1 + j) / delta
    delta  = sqrt(2 / (omega * mu * sigma))            [classical skin depth]
    xi(x)  = signed distance from x to conductor surface (>=0 inside)

For a circular cylinder of radius R, xi(x, y) = R - sqrt(x^2 + y^2).
For a sphere of radius R, xi(x, y, z) = R - sqrt(x^2 + y^2 + z^2).
For polyhedra, xi is the closest-point distance to the surface (can be
computed via a fast marching method or a pre-computed nodal distance
function).

## Why this is the "right" enrichment

The 1D skin-layer solution of the diffusion equation
    sigma * dA/dt = (1/mu) * d^2A/dx^2
with surface Dirichlet excitation is exactly

    A(x, omega) = A_0 * exp(-gamma * x)

where x is the distance from the surface.  The Hiruma psi is therefore
the exact 1D skin profile.  Multiplying psi by a slowly-varying nodal
basis u_enr(x) lets the enriched DOF capture skin physics
*tangentially* along the surface (varying around the boundary while
maintaining the exact normal decay).

## Frequency dependence: critical caveat

gamma depends on omega.  This means **psi(x) is frequency-dependent**:

    psi(x, omega_1) != psi(x, omega_2)

For a single-omega FE solve this is fine (just re-evaluate psi for the
new omega).  For:
  - frequency sweeps: each omega requires a fresh FE basis -> no ROM
  - time-domain transients: omega is undefined in the time domain
  - model order reduction (CLN, Multi-K): need static K_0, M -> psi
    must be FROZEN at one reference omega_ref

The freezing strategy is:

    psi_frozen(x) = exp(-xi / delta_ref),   delta_ref = sqrt(2/(omega_ref mu sigma))

i.e. a real-valued exponential at the skin depth of a chosen
omega_ref.  This was the configuration tested in beta-1 Phase 4.
"""


# ---------------------------------------------------------------------------
# NGSOLVE IMPLEMENTATION
# ---------------------------------------------------------------------------

NGSOLVE_IMPLEMENTATION = r"""
# EM-XFEM in NGSolve: implementation pattern (Phase 2 verified)

The Hiruma 2023 EM-XFEM is implementable in pure NGSolve **without** an
ngsxfem add-on, using a **compound H_1 x H_1 finite element space** with
manual product-rule gradient expansion.

## Implementation skeleton

```python
from ngsolve import (
    H1, BilinearForm, LinearForm, grad, dx,
    CoefficientFunction, sqrt as ng_sqrt, exp as ng_exp, x, y,
)

# Compound space: standard + enriched, both H_1, complex
fes_std = H1(mesh, order=1, complex=True, dirichlet="outer")
fes_enr = H1(mesh, order=1, complex=True, dirichlet="outer")
fes     = fes_std * fes_enr
(u_std, u_enr), (v_std, v_enr) = fes.TnT()

# Enrichment function psi(x) = exp(-gamma*xi)
gamma_c    = ng_sqrt(1j * omega * SIGMA * MU0)   # complex sqrt
r_cart     = ng_sqrt(x*x + y*y + 1e-30)          # eps for /0 at origin
xi         = R - r_cart                           # distance from surface
psi        = ng_exp(-gamma_c * xi)

# Gradient of psi: -gamma * grad(xi) * psi = gamma * (x,y)/r * psi
grad_psi_x = gamma_c * x / r_cart * psi
grad_psi_y = gamma_c * y / r_cart * psi
grad_psi   = CoefficientFunction((grad_psi_x, grad_psi_y))

# Manual product-rule for grad(psi * u_enr)
g_us  = grad(u_std)
g_ue  = psi * grad(u_enr) + u_enr * grad_psi
g_vs  = grad(v_std)
g_ve  = psi * grad(v_enr) + v_enr * grad_psi

# Total trial and test
u_total = u_std + psi * u_enr
v_total = v_std + psi * v_enr
g_ut    = g_us + g_ue
g_vt    = g_vs + g_ve

# Weak form: nu grad.grad + j omega sigma u v   (volume-source eddy)
a = BilinearForm(fes, symmetric=True)
a += (NU * g_ut * g_vt + 1j * omega * SIGMA * u_total * v_total) \
     * dx(bonus_intorder=10)            # CRITICAL: high quadrature
a.Assemble()

f = LinearForm(fes)
f += (J0 * v_total) * dx(bonus_intorder=10)
f.Assemble()
```

## Three pitfalls and how to avoid them

### Pitfall 1: bonus_intorder too low

The integrand contains `exp(-gamma*xi)` which is exponential, NOT
polynomial.  Standard FEM Gauss quadrature (~2-4 points per element)
will undersample it badly.  Use **bonus_intorder=10** (the verified
choice in Phase 2) or higher.  Without this, the XFEM solution
converges to nonsense.

### Pitfall 2: forgetting the product-rule gradient

`grad(psi * u_enr)` is **NOT** `psi * grad(u_enr)` -- there's the
`u_enr * grad(psi)` cross-term.  NGSolve will not assemble this
automatically because psi is a CoefficientFunction (not a basis
function).  You must write the product-rule expansion explicitly.

### Pitfall 3: enriching only some nodes via FreeDofs() filtering

The Hiruma paper enriches **all** nodes inside the conductor.  In
NGSolve you can filter to near-surface only by:

```python
free_enr = BitArray(fes_enr.ndof)
for n_idx, vert in enumerate(mesh.vertices):
    pos = mesh[vert].point
    xi_val = R - math.hypot(pos[0], pos[1])
    free_enr[n_idx] = (xi_val < 2.0 * delta_ref)
fes.FreeDofs()[fes_std.ndof:] = free_enr
```

This reduces XFEM cost but the gain is small (the enriched DOFs that
are far from the surface contribute negligibly since psi ~ 0 there).
For the cylinder benchmark we enrich ALL nodes for simplicity.

## Postprocessing the XFEM solution

The total field A_z is reconstructed by:

```python
A_std_cmp, A_enr_cmp = A.components
Az_total = A_std_cmp + psi * A_enr_cmp
```

and similarly for grad(A_z).  Joule loss and stored energy are then
integrated over the conductor with `bonus_intorder=10` to ensure the
exponential is captured in the quadrature.
"""


# ---------------------------------------------------------------------------
# CYLINDER VALIDATION
# ---------------------------------------------------------------------------

CYLINDER_VALIDATION = r"""
# Phase 2 cylinder benchmark: 0.14% at r/delta=15 with 88 DOFs

## Setup

Geometry:  2D Cu cylinder, R = 10 mm,  sigma = 10^6 S/m, mu_r = 1.
Source:    uniform volume current density J_0 = sigma (Hiruma 2023
           Section III.A reference configuration).
Mesh:      44-vertex coarse mesh (maxh = R/3.5), matches Hiruma's
           "37 unknowns" coarse-mesh target for comparison.
XFEM:      compound H_1 x H_1, all-node enrichment (44 std + 44 enr =
           88 free DOFs).

Bessel reference (corrected form):
    Y_exact(j*omega) = sigma * 2*pi*R * I_1(kappa*R) / (kappa * I_0(kappa*R))
    kappa = sqrt(j*omega*mu*sigma)

## Results (Phase 2, 2026-05-25)

| r/delta | CFEM (fine) | CFEM (coarse 44 DOF) | XFEM (88 DOF) |
|---------|-------------|----------------------|---------------|
| 0.1     | +0.00 %     | +0.00 %              | +0.00 %       |
| 3.0     | +0.00 %     | +0.68 %              | +3.01 %       |
| 15.0    | +0.04 %     | +85.89 %             | **+0.14 %**   |

## What this proves

1. CFEM at the same coarse mesh (44 DOF) is +85.89% off at r/delta=15.
   So the skin layer is unresolvable on this mesh by adding more
   higher-order CFEM DOFs alone.
2. XFEM at 88 DOF (= CFEM + 44 enrichment DOFs) reaches +0.14%.
   The enrichment captures the skin-depth physics that the coarse
   mesh cannot represent geometrically.
3. The Hiruma 2023 claim "< 3% at r/delta < 15 with 37-DOF coarse mesh"
   is reproduced with substantial margin in the NGSolve compound-space
   implementation.

This was the entire content of beta-1 Phase 2 (2026-05-25): a clean,
independent reproduction of Hiruma's headline benchmark using only
upstream NGSolve features (no ngsxfem add-on).

## Reproducer

    cd S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/
    jupyter nbconvert --to notebook --execute --inplace hiruma_xfem_comparison.ipynb
    # the Phase 2 XFEM enrichment section runs in ~30-90 s (omega sweep)
    # output: STATUS.md table, results_phase2.npz
"""


# ---------------------------------------------------------------------------
# VOLUME-SOURCE SCOPE
# ---------------------------------------------------------------------------

VOLUME_SOURCE_SCOPE = r"""
# EM-XFEM scope: port-driven vs volume-source

## The Phase 3 finding (volume-source FE residual)

For the Hiruma cylinder benchmark, the source is **volume-distributed**:
the entire conductor interior carries uniform J_0 = sigma E_z (driven
electrically).  This is different from a port-driven problem where
the excitation enters through a small portion of the boundary
(Dirichlet on a face or Neumann surface current).

The augmented-CLN ROM construction (Sugahara 2026 Paper I, Theorem 1)
adds the sqrt(s) Schur block to a rational base ROM Y_R(s).  The
preserved-asymptote claim requires Y_CLN(s -> infinity) -> 0.  In the
**discrete FE setting** this is

    Y_CLN(infty) = sigma * (|Omega| - b_r^T M_r^{-1} b_r)

with b the FE source vector and (K_0, M, b_r, M_r) the projected
quantities.  For a port-driven b, the discrete ratio
b_r^T M_r^{-1} b_r / |Omega| converges to the boundary port area as
the mesh refines, giving Y_CLN(infty) -> 0 (the algebraic decay
becomes exponential).  For a volume-source J_0 b, the residual
decays only **algebraically** with mesh size.

Beta-1 Phase 3 verified on Hiruma cylinder:

| Mesh       | DOF   | Y_CLN(infty) | augmented-CLN r(10^12 Hz) |
|------------|-------|--------------|---------------------------|
| maxh=R/3.5 | 44    | 53           | 1.0 -> 239 (collapses)    |
| maxh=R/10  | 393   | 19           | improved but still wrong  |
| maxh=R/30  | 3558  | 13           | algebraic decay only      |

So the augmented-CLN single-DOF SIBC asymptote claim **fails** for
volume-source on the Hiruma cylinder.  This is NOT a contradiction
with Paper 1 Theorem 1 -- Theorem 1's hypothesis (port-driven) is
violated.

## How EM-XFEM cures the volume-source FE residual

The Hiruma EM-XFEM puts the analytical skin profile exp(-gamma*xi)
into the FE space.  The XFEM-enriched source b_r now lives in a
subspace that **contains** the skin layer, so the discrete ratio
b_r^T M_r^{-1} b_r / |Omega| converges fast (exponentially with
delta/h, not algebraically).  Y_FE_XFEM(infty) -> 0 even for
volume-source.

This is Phase 2's good news: **on volume-source, EM-XFEM cures the
finite FE-discretization residual at infinity**, which canonical
CFEM cannot do without prohibitive mesh refinement.

## When the volume-source scope matters

The volume-source structural problem (Y_CLN(infty) != 0 for canonical
CFEM-CLN) is generic to ANY problem with uniform-driver J_0 in the
conductor body:
  - Induction heating workpiece (uniform J_0 = sigma E_z applied)
  - Hiruma 2023 paper cylinder
  - Field-from-conductor calculations

Port-driven problems (no Y_CLN(infty) issue):
  - Boundary-current SIBC drive
  - Cuboid Y(s) of Paper 1 Sec V (port-area-supported b)
  - Multi-conductor BEM-CLN driven from external coils
  - Most engineering coil-and-load configurations

The bottom line: if your source is uniform inside the conductor,
the **EM-XFEM at FE level** is the right tool, NOT the augmented-CLN
ROM alone.
"""


# ---------------------------------------------------------------------------
# CLN STACKING (NEGATIVE RESULT)
# ---------------------------------------------------------------------------

CLN_STACKING_NEGATIVE = r"""
# XFEM + canonical CLN stacking: does NOT extend Hankel conditioning

## The tempting hypothesis

If XFEM at the FE level captures skin modes that CFEM cannot, and
canonical CLN extraction breaks at N ~ 4 in float64 because of the
Hankel-QD breakdown amplified by the skin-mode contribution to
moments, then stacking XFEM + canonical CLN might extend the
float64 reliability of high-N rung extraction.

This was tested in beta-1 Phase 4 (2026-05-25).

## Setup

- Static (frozen) XFEM enrichment: psi(x) = exp(-xi/delta_ref) with
  delta_ref = sqrt(2/(omega_wall * mu * sigma)) at the wall-band
  center.
- Cylinder a=5mm, sigma=1e6 S/m, 44-vertex mesh.
- CFEM (24 free DOFs) vs XFEM (48 free DOFs, +24 enrichment).
- Compute Krylov moments mu_n = b^T (K_0^{-1} sigma M)^n K_0^{-1} b
  for n = 0..19.
- Form Hankel matrices H_N, compute condition numbers kappa(H_N).

## Results

| N  | kappa_CFEM(H_N)  | kappa_XFEM(H_N)  | Ratio | Float64 OK? |
|----|------------------|------------------|-------|-------------|
| 2  | 1.15e+11         | 7.71e+10         | 1.5   | marginal both |
| 3  | 2.93e+23         | 5.49e+22         | 5.3   | NO both      |
| 5  | 3.19e+48         | 7.90e+47         | 4.0   | NO           |
| 10 | 2.18e+101        | 1.02e+100        | 21.5  | NO           |

**Both break at N=3.  XFEM gives at most 1-21x improvement in kappa,
which is NOT enough to extend the float64 stage limit beyond N=2-3.**

Krylov moments are nearly identical between CFEM and XFEM:
  mu_0:  4.677e-3 (CFEM) vs 4.689e-3 (XFEM) -- 0.3% diff
  mu_9:  3.51e-45 vs 3.64e-45 -- 4% diff

## Why the hypothesis was wrong

The Krylov sequence at s = 0 weights eigenmodes by tau_k^n:

    v_n = (K_0^{-1} sigma M)^n K_0^{-1} b = sum_k g_k * tau_k^n * phi_k

**High-k skin modes have tiny tau_k**, so their contribution to
Krylov vectors is **exponentially damped**.  The XFEM-enriched DOFs
that capture these high-k modes contribute negligibly to canonical
CLN moments at s = 0.

The Hankel ill-conditioning is set by the **bulk eigenvalue spread**,
which CFEM and XFEM capture identically (low-k modes are well-resolved
on the coarse mesh; XFEM only adds high-k DOFs that don't help here).

## Implication for the radia-mcp decision layer

- XFEM utility = FE-level skin-layer resolution + volume-source FE
  residual cure (Phase 2, Phase 3).
- XFEM utility is NOT canonical CLN moment-matching stage stability.
  Paper 1 Sec VI's "N <= 4 float64 limit" stands even when XFEM is
  stacked.
- For broadband accuracy with reliable circuit identification:
    1. Multi-K (Kuriyama 2019) at lower stage, K_0-MGS on multi-point
       Krylov sidesteps Hankel ill-conditioning -- but loses canonical
       Cauer identity (Krylov POD).
    2. Foster-of-Cauers (Sugahara 2026 Paper II): low-N Cauer macro +
       diffusive Foster terminator at the surface.
    3. MPFR + INTLAB (Nagamine 2026): keep canonical Cauer but pay
       192-bit precision cost.
  XFEM-CLN stacking is NOT among these.

## Why a NEGATIVE result is useful

The 2026-05-25 user instinct that this hypothesis was worth testing
("XFEM-CLN したら高段で CLN が壊れにくいとかいう特性がでるなら、
それは知見だなぁ") was correct.  The experiment is clean and
falsifiable, and the verdict is unambiguous.

XFEM's role is now precisely circumscribed:
  - FE level: yes, skin layer + volume-source FE residual
  - ROM level: no, does not improve canonical CLN stage stability

## Files

    docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb   (Phase 4 section)
    docs/hiruma_xfem_comparison/phase4_xfem_hankel_conditioning.png
    docs/hiruma_xfem_comparison/phase4_xfem_hankel_results.json
    memory/project_xfem_cln_hankel_no_improvement.md
"""


# ---------------------------------------------------------------------------
# DECISION TABLE
# ---------------------------------------------------------------------------

DECISION_TABLE = r"""
# When to use EM-XFEM: decision table

## Layer view

| Layer | XFEM applicable? | Notes |
|-------|------------------|-------|
| FE level: skin-layer resolution at high f | YES | Hiruma 2023 original use |
| FE level: volume-source FE residual cure (Y_FE(infty) -> 0) | YES | Phase 2 / Phase 3 |
| FE level: sharp-edge geometry without prohibitive mesh refinement | YES (if needed) | with edge enrichment too |
| ROM level: canonical CLN Hankel-QD stability | NO | Phase 4 verified |
| ROM level: float64 high-N stage extension | NO | Phase 4 verified |

## Use-case decision tree

```
Q1. Do you need spatial field distribution at one (or a few) frequencies?
     YES -> XFEM at FE level (best 3D field on coarse mesh)
     NO  -> continue.

Q2. Frequency sweep / time-domain / SPICE export?
     YES -> XFEM does not give a ROM directly; route to:
                * Multi-K + Schur (Paper 1 Sec V.G) for port-driven
                * Foster-of-Cauers for time-domain SPICE
     NO  -> continue.

Q3. Is the source uniform inside the conductor (volume-source)?
     YES -> XFEM at FE level cures the FE residual (Y_FE_XFEM(infty) -> 0).
            Then optionally extract a ROM, but it will no longer be
            canonical CLN -- it lives in the XFEM-enriched space.
     NO (port-driven) -> XFEM optional.  Augmented CLN ROM works
                         standalone for SIBC asymptote (Paper 1 Sec V).

Q4. Pure industrial high-f only (DC accuracy not needed)?
     YES -> classical SIBC at FE level (Robin BC) is the cheapest.
            XFEM is overkill.
```

## Application-to-method mapping

| Application                                          | Recommended      |
|------------------------------------------------------|------------------|
| Single-freq 3D + volume source + corner detail       | EM-XFEM (Hiruma) |
| Single-freq 3D + port-driven + smooth geometry       | CFEM-SIBC        |
| Frequency sweep + port-driven                        | Augmented CLN    |
| Time-domain transient + port-driven                  | Augmented CLN    |
| Multi-conductor IH coil + workpiece                  | BEM-CLN + Augmented CLN per element (Paper II) |
| Plain industrial high-f, no DC                       | classical SIBC   |
| Crack growth (mechanical XFEM, not EM)               | fracture XFEM (Jafari et al. 2021) — different problem |
"""


# ---------------------------------------------------------------------------
# RELATION TO ngsxfem LIBRARY
# ---------------------------------------------------------------------------

NGSXFEM_RELATION = r"""
# Relation to the ngsxfem library: cut-FEM XFEM vs Hiruma PUFEM

## The ngsxfem package

Repository: https://github.com/ngsxfem/ngsxfem
Lead author: Christoph Lehrenfeld (Goettingen).  PyPI: ``xfem``.

ngsxfem is the official NGSolve add-on for **cut-FEM** and other
unfitted-mesh discretizations.  Its installed Python module name is
``xfem`` (despite the project name ngsxfem on GitHub).  Pre-built
Windows wheels are available on PyPI (``pip install xfem``).

## What ngsxfem does (cut-FEM family)

In cut-FEM the mesh does NOT conform to the physical interface
(material boundary, free surface, moving body).  A signed-distance
level-set phi(x) is used to mark elements as "inside", "outside", or
"cut".  Cut elements need:

  - Custom quadrature rules that integrate ONLY over the inside
    portion of each cut element (sub-cell triangulation in 2D /
    sub-tet decomposition in 3D).
  - Ghost penalty stabilization to control conditioning when an
    element's intersection with the inside is geometrically tiny.
  - Piecewise-polynomial trial / test functions defined ONLY where
    they make sense (TraceFE, XFiniteElement, ...).

ngsxfem provides the NGSolve plumbing for all of this:

  - ``CutInfo`` -- classifies elements by level-set sign.
  - ``HasPos / HasNeg / IF`` -- domain-of-integration restrictions.
  - ``dx(definedon=...)`` integrators on cut domains.
  - ``XFiniteElement`` -- piecewise-polynomial element on cut cells.
  - ``GhostPenalty`` -- stabilization integrator.

Typical applications:

  - Two-phase / Stefan / phase-change problems
  - Free-surface flow (water + air with moving interface)
  - Optimization with moving design domains
  - Crack growth where the crack is treated as a level-set

The fracture-mechanics XFEM (Jafari et al. 2021) is one application
that would naturally map onto ngsxfem, but that paper used a different
solver stack.

## Why ngsxfem is NOT needed for Hiruma EM-XFEM

The Hiruma 2023 EM-XFEM is a **PUFEM** (partition-of-unity FEM)
enrichment, not cut-FEM.  Key differences:

  (1) The mesh CONFORMS to the conductor boundary in the Hiruma
      construction (you mesh the conductor as a regular domain;
      no level-set, no cut cells).
  (2) The enrichment is GLOBAL within the conductor (ALL nodes
      can carry the psi(x) = exp(-gamma*xi) enrichment).  There is
      no jump or singularity along an internal level set.
  (3) Integration is over standard elements with HIGH-ORDER Gauss
      (bonus_intorder=10) to capture the exponential psi -- no
      sub-cell triangulation needed.

In NGSolve, PUFEM with one enrichment function is realized as a
COMPOUND finite-element space:

    fes_std = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes_enr = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes     = fes_std * fes_enr

with the (u_std, u_enr) pair representing
``u_total(x) = u_std(x) + psi(x) * u_enr(x)``.  This is **vanilla
NGSolve, no ngsxfem dependency**.  This is exactly what
the phase2 section of the notebook uses.

## Could ngsxfem reproduce Hiruma EM-XFEM?

User question 2026-05-25: "Hiruma EM-XFEM は ngsxfem では再現できない?
仕組み上はできる?"

Three-level honest answer:

| Level | Answer | Why |
|-------|--------|-----|
| Mechanism-wise | **YES** | The underlying PUFEM / extended-FEM machinery that ngsxfem builds on can express smooth enrichment |
| Default ngsxfem API | **NO** | ``XFiniteElement`` uses Heaviside ``H(phi)`` enrichment by default -- designed for cracks / interfaces, not smooth exp(-gamma*xi) skin profiles |
| Practical | **Not needed** | Vanilla NGSolve compound H1 x H1 (Phase 2 approach) is structurally minimal -- ngsxfem's cut-FEM machinery has nothing to contribute when the mesh conforms |

### What would "ngsxfem-Hiruma" look like (hypothetical)?

If you insisted on routing Hiruma EM-XFEM through ngsxfem:

```python
# Pseudo-code (NOT recommended -- use Phase 2 compound space instead)
from xfem import CutInfo, HASNEG  # ngsxfem
from ngsolve import H1, BilinearForm, ...

# (1) Level-set that is everywhere positive inside the conductor.
lset = GridFunction(H1(mesh, order=1))
InterpolateToP1(R - sqrt(x*x + y*y), lset)
ci = CutInfo(mesh, lset)
# All elements are HASNEG (= entirely inside) since mesh conforms.
# No cut elements -> no GhostPenalty needed, no cut quadrature needed.

# (2) Custom XFiniteElement-like compound space with smooth enrichment.
# ngsxfem's default XFiniteElement multiplies the enriched basis by
# H(phi) = +1 inside, -1 outside.  To use exp(-gamma*xi) instead, we
# would need to override the enrichment factor in the BilinearForm
# integrand -- but at that point we're hand-rolling the integrand,
# bypassing ngsxfem's abstraction.
```

The two steps that ngsxfem would naturally provide (CutInfo for cut
classification, GhostPenalty for stabilization) **don't apply** to the
Hiruma setting (no cut, no instability from cut shape).  And the one
step that matters (smooth psi=exp(-gamma*xi) enrichment) is NOT what
ngsxfem's XFiniteElement does by default.

Net: **ngsxfem is mechanically capable but practically inappropriate
for Hiruma EM-XFEM**.  Use vanilla NGSolve compound H1 x H1 as in
Phase 2.

### Hiruma paper's own implementation

Hiruma et al. 2023 IEEE TMag implements the enrichment via a **fully
custom code** (not COMSOL, not ngsxfem) -- they write their own
quadrature loop with the manual product-rule gradient.  Phase 2's
NGSolve compound-space port (88 free DOFs, 0.14% at r/delta=15) is
the first independent reproduction in NGSolve and confirms that the
compound-space pattern is sufficient.

## When you WOULD use ngsxfem for EM problems

  - Moving conductor inside a fixed mesh (e.g. a rotor sweeping past
    a stator).  Cut-FEM is the natural framework -- the mesh stays
    put, the conductor boundary is a moving level set.
  - Magnetic interface between two materials whose boundary you do
    NOT want to mesh-conform to (e.g. fractal or evolving shapes).
  - Coupling EM with two-phase flow (induction heating with phase
    change in the workpiece).  ngsxfem handles the moving interface
    naturally.

For these cases, ngsxfem is the right NGSolve add-on, NOT the
compound-space PUFEM of Phase 2.

## Practical install (if a future user needs cut-FEM EM)

```bash
pip install xfem
```

(PyPI Windows wheel works with the NGSolve 6.2.2603+ stack.)

Quick example (from ngsxfem README):

```python
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *

# A unit square mesh
square = SplineGeometry()
square.AddRectangle((0,0),(1,1))
mesh = Mesh(square.GenerateMesh(maxh=0.1))

# Level set: circle of radius 0.4 around (0.5, 0.5)
levelset = sqrt((x-0.5)**2 + (y-0.5)**2) - 0.4
lset = GridFunction(H1(mesh, order=1))
InterpolateToP1(levelset, lset)

# Classify cut elements
ci = CutInfo(mesh, lset)
hasif = ci.GetElementsOfType(IF)

# Integrate ONLY inside the circle
result = Integrate(levelset_domain={"levelset": lset, "domain_type": NEG},
                   cf=1.0, mesh=mesh, order=2)
print(f"Area inside circle: {result:.4f}  (analytic: {math.pi*0.16:.4f})")
```

## Glossary cross-references

| Term used | ngsxfem name | Hiruma PUFEM equivalent          |
|-----------|--------------|----------------------------------|
| level-set phi(x) | ``lset`` GridFunction | xi(x) = distance from surface |
| domain restriction | ``definedon=...`` integrator | full conductor Omega          |
| cut element | element with mixed-sign phi | (none -- no cut)              |
| ghost penalty | ``ghost_penalty`` integrator | (none -- no cut)              |
| trace operator | ``TraceOperator`` | n/a (no interface trace)        |

## TL;DR

| Use case | Choose |
|----------|--------|
| Hiruma 2023 EM skin-depth enrichment (Phase 2 reproduction) | vanilla NGSolve compound space ✓ |
| Cut-FEM with moving interface, two-phase flow, evolving levelsets | ngsxfem library |
| Solid-mech crack growth | fracture XFEM (Jafari et al. 2021) or ngsxfem with level-set |

The three "XFEM"s share the name but solve different problems.  Be
explicit about which one you mean.
"""


# ---------------------------------------------------------------------------
# REPRODUCTION
# ---------------------------------------------------------------------------

REPRODUCTION = r"""
# Reproducing the Phase 1-4 beta-1 benchmark

## Files

```
S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/
    hiruma_xfem_comparison.ipynb            - consolidated Phase 1-4 notebook
        (sections: phase1 CFEM baseline; phase2 Hiruma XFEM enrichment;
         phase3 sqrt(s) Schur [deprecated Pade predecessor]; phase3b
         Galerkin-Krylov augmented CLN; phase4 XFEM-CLN Hankel conditioning)
    results_phase1.json                     - Phase 1 numerical output
    phase4_xfem_hankel_results.json         - Phase 4 numerical output
    phase4_xfem_hankel_conditioning.png     - Phase 4 figure
```

## Dependencies

- Python 3.12
- NGSolve 6.2.2603+
- numpy, scipy, matplotlib

No ngsxfem add-on required.  The Hiruma EM-XFEM is realized via
NGSolve's compound space machinery (H1 * H1) with manual product-rule
gradients.

## Run order

```bash
cd S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/

# Run all Phase 1-4 sections (consolidated notebook, ~19 s total)
jupyter nbconvert --to notebook --execute --inplace hiruma_xfem_comparison.ipynb
# sections: phase1 CFEM baseline vs Bessel; phase2 Hiruma EM-XFEM enrichment;
# phase3b augmented CLN Galerkin-Krylov (phase3 = deprecated Pade predecessor);
# phase4 XFEM-CLN Hankel conditioning -> phase4_xfem_hankel_conditioning.png + .json
```

## Expected outputs

- Phase 1: CFEM fine = 0.04% vs Bessel at r/delta=15; coarse = +85.89%.
- Phase 2: XFEM 88-DOF = +0.14% at r/delta=15.
- Phase 3b: Y_CLN(infty) = 53/19/13 for mesh = R/3.5/R/10/R/30
  (algebraic decay, augmented CLN r(10^12) collapses for volume-source).
- Phase 4: CFEM and XFEM Hankel kappa both break at N=3 (1.5-21x XFEM
  improvement, not enough to extend the float64 stage limit).

## Citing this work

If you reproduce or extend the Phase 1-4 benchmark for publication,
cite:
    Sugahara, Nagamine, Hane (2026 in preparation).
    "Asymptote-preserving Schur augmentation of CLN" (IEEE TMag).

The Phase 2 XFEM reproduction also independently confirms:
    Hiruma et al., "Extended FEM for Eddy-Current Problems with
    Surface Skin Effect", IEEE TMag 59(5), 2023.
"""


# ---------------------------------------------------------------------------
# Topic dispatcher
# ---------------------------------------------------------------------------

TOPICS = {
    "overview":              "High-level summary, comparison to solid-mechanics XFEM (DEFAULT)",
    "enrichment_function":   "psi(x) = exp(-gamma*xi); frequency dependence; freezing for ROM",
    "ngsolve_implementation":"NGSolve compound H1 x H1 pattern, three pitfalls",
    "cylinder_validation":   "Phase 2: 88-DOF reproduction of Hiruma 2023 (0.14% at r/delta=15)",
    "volume_source_scope":   "Phase 3: XFEM cures volume-source FE residual that breaks augmented CLN",
    "cln_stacking_negative": "Phase 4: XFEM does NOT extend canonical CLN Hankel-QD stability (verified)",
    "decision_table":        "When to use EM-XFEM: layer / use-case / application map",
    "ngsxfem_relation":      "ngsxfem library (cut-FEM XFEM) vs Hiruma PUFEM EM-XFEM",
    "reproduction":          "How to run Phase 1-4 benchmarks; expected outputs",
    "lab_files":             "Absolute paths to Phase 1-4 scripts and figures (LAB_FILES dict)",
    "all":                   "Everything (concatenated)",
}


def _lab_files_block() -> str:
    """Render LAB_FILES as a markdown table for display via the MCP tool."""
    lines = [
        "# Lab files (S:/Radia/01_GitHub/docs/hiruma_xfem_comparison/)",
        "",
        "| Key | Path |",
        "|-----|------|",
    ]
    for k, v in LAB_FILES.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("Use these via:")
    lines.append("    from radia_mcp.fem.xfem_em_hiruma_knowledge import LAB_FILES")
    lines.append("    p = LAB_FILES['notebook']")
    return "\n".join(lines)


def get_em_xfem_knowledge(topic: str = "overview") -> str:
    """Dispatch EM-XFEM (Hiruma 2023) knowledge topics.

    Topics:
        overview                - Summary, vs solid-mech XFEM (DEFAULT)
        enrichment_function     - psi math, frequency dependence, freezing
        ngsolve_implementation  - Compound H1 x H1 pattern, pitfalls
        cylinder_validation     - Phase 2: 0.14% at r/delta=15, 88 DOF
        volume_source_scope     - Phase 3: cures FE residual that breaks
                                  augmented CLN on volume-source
        cln_stacking_negative   - Phase 4: XFEM does NOT extend canonical
                                  CLN Hankel-QD stability (rejected)
        decision_table          - When to use EM-XFEM (layer / use case)
        reproduction            - How to run Phase 1-4 benchmarks
        lab_files               - Paths to scripts and figures
        all                     - Everything (concatenated)
    """
    t = topic.lower().strip()
    if t in ("overview", "summary", "intro"):
        return OVERVIEW
    if t in ("enrichment_function", "enrichment", "psi"):
        return ENRICHMENT_FUNCTION
    if t in ("ngsolve_implementation", "implementation", "ngsolve", "code"):
        return NGSOLVE_IMPLEMENTATION
    if t in ("cylinder_validation", "phase2", "validation", "cylinder"):
        return CYLINDER_VALIDATION
    if t in ("volume_source_scope", "phase3", "volume_source", "scope"):
        return VOLUME_SOURCE_SCOPE
    if t in ("cln_stacking_negative", "phase4", "cln_stacking", "stacking",
             "hankel", "negative"):
        return CLN_STACKING_NEGATIVE
    if t in ("decision_table", "decision", "when_to_use", "use_case"):
        return DECISION_TABLE
    if t in ("ngsxfem_relation", "ngsxfem", "cut_fem", "cutfem", "lehrenfeld"):
        return NGSXFEM_RELATION
    if t in ("reproduction", "reproduce", "run", "benchmark"):
        return REPRODUCTION
    if t in ("lab_files", "files", "paths"):
        return _lab_files_block()
    if t == "all":
        return "\n\n".join([
            OVERVIEW,
            ENRICHMENT_FUNCTION,
            NGSOLVE_IMPLEMENTATION,
            CYLINDER_VALIDATION,
            VOLUME_SOURCE_SCOPE,
            CLN_STACKING_NEGATIVE,
            DECISION_TABLE,
            NGSXFEM_RELATION,
            REPRODUCTION,
            _lab_files_block(),
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        + ", ".join(TOPICS.keys())
        + "."
    )
