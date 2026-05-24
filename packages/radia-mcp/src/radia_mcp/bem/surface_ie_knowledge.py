"""Surface Integral Equations: EFIE / MFIE / CFIE / PMCHWT."""

EFIE = r"""
# Electric Field Integral Equation (EFIE)

For PEC scattering: enforce E_tangential = 0 on PEC surface Γ:

```
(1) -jωμ A(r) + ∇φ(r) = E_inc(r) × n̂    on Γ
where:
    A(r)  = ∫_Γ G_k(r, r') J(r') dS'
    φ(r)  = (1/(jωε)) ∫_Γ G_k(r, r') div_Γ J(r') dS'
    G_k   = e^{-jkr} / (4π r)   (Helmholtz Green) or
          = 1 / (4π r)           (Laplace, for ω → 0)
```

## RWG discretization

Expand J = Σ_n J_n f_n with RWG basis f_n.  Galerkin testing with f_m:

```
Z_mn = jωμ ∫_Γ ∫_Γ f_m(r) · f_n(r') G_k(|r-r'|) dS dS'
     - (1/(jωε)) ∫_Γ ∫_Γ div_Γ f_m(r) div_Γ f_n(r') G_k(|r-r'|) dS dS'

b_m = ∫_Γ f_m(r) · E_inc(r) dS
```

System: Z · J = b.

## Pros and cons

| Pros | Cons |
|------|------|
| Works for open surfaces (sheets) | Interior resonances at PEC modes |
| Well-posed at high freq | LOW-FREQ BREAKDOWN: kernel ill-conditioned as ω→0 |
| Standard formulation | Charge-current scale mismatch (1/ω vs ω) |

## Low-frequency breakdown — the fundamental problem

As ω → 0:
- The A-term (jωμ A) → 0
- The φ-term (∇φ via div J) → ∞ (1/ω blowup)

Z matrix becomes ill-conditioned with condition number ~ 1/ω².
Iterative solvers diverge.

**Fix**: Loop-Star decomposition or Calderón preconditioner
(see `low_freq` topic).

## In ngsolve.bem

```python
from ngsolve.bem import (SingleLayerPotentialOperator,
                          DoubleLayerPotentialOperator)
fes = HDivSurface(mesh, order=1)   # RWG
S = SingleLayerPotentialOperator(fes, kernel='helmholtz', k=k)
# ... assemble EFIE system Z = -jωμ·M_curl + (1/jωε)·M_div
```

References in this folder:
- 04_efie_mfie_cfie/IE_Solution_Maxwell_Zero_to_Microwave.pdf
- 05_low_freq_stabilization/* (broad LF stabilization papers)
"""


MFIE = r"""
# Magnetic Field Integral Equation (MFIE)

Enforce H_tangential continuity at PEC: H_inside · n̂ = -n̂ × J:

```
(1/2) J(r) + n̂(r) × p.v. ∫_Γ ∇G(r,r') × J(r') dS' = n̂(r) × H_inc(r)
```

## Properties vs EFIE

| Property | EFIE | MFIE |
|----------|------|------|
| Works for open surfaces (sheets) | yes | NO (needs closed surface) |
| Identity term (well-conditioned) | NO (purely integral) | YES (I/2 + K) |
| Low-frequency behavior | BREAKS DOWN | also breaks down but differently |
| Conditioning at high freq | hard (no I term) | EASY (I + compact) |
| Accuracy near sharp edges | good | POOR |

## MFIE with RWG: a documented problem

[LOCAL] 04_efie_mfie_cfie/mfie/Inaccuracy_MFIE_RWG.pdf

The standard RWG basis is NOT well-suited to MFIE because RWG is
H(div)-conforming but MFIE needs H(curl)-conforming testing.  Result:
MFIE-RWG gives systematically wrong currents on coarse meshes.

**Fixes** (in this folder):
- BC (Buffa-Christiansen) basis functions — duality with RWG
- Curl-conforming basis (Velasquez 2005, [LOCAL] Solution_MFIE_curl_conforming.pdf)
- Linear-linear basis (Davis-Warnick 2003)

## MFIE at low frequency

[LOCAL] 04_efie_mfie_cfie/mfie/MFIE_at_Very_Low_Freq.pdf
[LOCAL] 05_low_freq_stabilization/Overcoming_LowFreq_Breakdown_MFIE.pdf

Vico-Greengard-Gimbutas-Cools (2013, IEEE TAP 61:1285) showed MFIE
breakdown at low-freq stems from the static MFIE kernel being
non-invertible on simply connected regions.  Fix: add Calderón-style
identity correction.
"""


CFIE = r"""
# Combined Field Integral Equation (CFIE) — interior resonance free

CFIE = α · EFIE + (1-α) · η · MFIE,  with α ∈ (0,1), η = √(μ/ε).

## Why CFIE

Both EFIE and MFIE for PEC have **interior resonance** problems: at
PEC interior cavity resonance frequencies, the surface IE has spurious
null space → singular system matrix.

The combination CFIE avoids ALL interior resonances at NO frequency.

## When to use which

| Case | Use |
|------|-----|
| Open surface (sheet, plate) | EFIE only (MFIE/CFIE require closed) |
| Closed PEC body, no resonance risk | EFIE (simpler) |
| Closed PEC body, high freq scan | **CFIE** (resonance-free) |
| Penetrable dielectric | **PMCHWT** (not CFIE) |

## α choice

- α = 0.5 is a common default
- α → 1: more EFIE flavor (better for edges, worse conditioning)
- α → 0: more MFIE flavor (better conditioning, worse for edges)

For lab-relevant problems (low-freq, PEEC, eddy current), CFIE is
**not commonly used** — the low-frequency breakdown of both EFIE and
MFIE dominates.

## In ngsolve.bem

```python
# Galerkin form:
S = SingleLayerPotentialOperator(fes, kernel='helmholtz', k=k)
D = DoubleLayerPotentialOperator(fes, kernel='helmholtz', k=k)
alpha = 0.5
Z_cfie = alpha * (-1j*omega*mu * S + 1/(1j*omega*eps) * S_div) + \
         (1-alpha) * eta * (0.5 * I + K)
```

(Pseudocode — actual ngsolve.bem API differs; check
`radia_ngsolve.ngsbem_inductance` for usage.)
"""


PMCHWT = r"""
# PMCHWT — Penetrable dielectric / conductor scattering

PMCHWT (Poggio-Miller-Chang-Harrington-Wu-Tsai) couples E and H
unknowns on BOTH sides of a dielectric interface.

## The formulation

On interface Γ between regions 1 (exterior) and 2 (penetrable body):
- Unknowns: J = n̂ × H, M = E × n̂  (both surface currents)
- 4 equations (E and H continuity in both regions)
- 2-by-2 block system

## Why PMCHWT > simpler alternatives

| Formulation | When applicable | Issues |
|-------------|-----------------|--------|
| EFIE-only on PEC | PEC bodies | Doesn't work for dielectric |
| Volume IE | Dielectric bodies | N^3 scaling (volume DOFs) |
| Single-source IE on dielectric | Dielectric bodies | Internal resonances |
| **PMCHWT** | Any penetrable body | Low-freq breakdown |

## Low-freq breakdown — PMCHWT version

[LOCAL] 05_low_freq_stabilization/LowFreq_PMCHWT_Dielectric_Conductive.pdf

At ω → 0, PMCHWT becomes ill-conditioned because:
- M (= E × n̂) → E_static · n̂  scales as O(1)
- J (= n̂ × H) → 0  scales as O(ω) (no eddy current at ω=0)

The mismatch in J vs M scale → bad conditioning.

**Fix (lab-relevant)**: For eddy current problems where the body is
a good conductor (σ >> ωε), USE PMCHWT-LF stabilization instead of
trying to solve the original PMCHWT at low frequencies.

In ngsolve.bem, the recommended **single-trace formulation**
(Weggler-Hiptmair) absorbs this directly.

## Use in lab

For Radia IH and PEEC workflows, PMCHWT is **not directly used** —
the lab uses FEM-SIBC instead (FEM inside the conductor + Robin BC at
surface).  PMCHWT is documented here for completeness and for cases
where pure surface BEM is preferred over hybrid FEM-BEM.
"""


def get_surface_ie_knowledge(topic: str = "efie") -> str:
    """Dispatch surface IE topics.

    Topics:
        efie         - Electric Field Integral Equation (DEFAULT)
        mfie         - Magnetic Field Integral Equation
        cfie         - Combined Field IE (interior resonance free)
        pmchwt       - Penetrable dielectric/conductor (Poggio-Miller-...)
        all          - Everything
    """
    topic = topic.lower().strip()
    if topic == "efie":
        return EFIE
    if topic == "mfie":
        return MFIE
    if topic == "cfie":
        return CFIE
    if topic in ("pmchwt", "dielectric"):
        return PMCHWT
    if topic == "all":
        return "\n\n".join([EFIE, MFIE, CFIE, PMCHWT])
    return (f"Unknown topic '{topic}'. Available: efie, mfie, cfie, pmchwt, "
            "all.")
