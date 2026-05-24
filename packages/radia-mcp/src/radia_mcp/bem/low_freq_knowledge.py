"""Low-frequency stabilization for surface BEM."""

PROBLEM = r"""
# The low-frequency breakdown problem

When ω → 0 in surface IE (EFIE, MFIE, PMCHWT), the discrete system
becomes ill-conditioned:

| Issue | Mechanism |
|-------|-----------|
| EFIE breakdown | A-term scales O(ω), φ-term scales O(1/ω) — scale mismatch |
| MFIE breakdown | Static MFIE kernel non-invertible on simply-connected regions |
| PMCHWT breakdown | J → 0 but M → O(1) — scale mismatch |

Why this matters for LAB:
- IH workpiece at 1-100 kHz is "low frequency" relative to RF
- PEEC coupling at kHz frequencies needs LF-stable formulations
- ngsolve.bem at DC needs LF stabilization

## Two main fix strategies

1. **Loop-Star basis decomposition** (Vecchi, Wilton-Cwik)
   — split current J into divergence-free (Loop) and curl-free (Star)
2. **Calderón preconditioner** (Andriulli-Cools-Bagci 2008)
   — pre-multiply by dual operator
3. **Single-trace formulation** (Weggler-Hiptmair)
   — reformulate to avoid both breakdowns

See sections below.
"""


LOOP_STAR = r"""
# Loop-Star Basis Decomposition

For RWG basis on closed surface:
- Total DOFs = N_edges
- Loops: N_v - 1 (one per closed loop, except the "global" loop)
- Stars: N_T (one per triangle face, except the global constant)
- Plus 2g "harmonic" handles for genus-g surface

So:  RWG_dim = N_edges = (N_v - 1) + N_T + 2g - 1   (Euler)

## The transformation

Transform RWG to (Loop, Star, Harmonic) basis:

```
J = J_loop + J_star + J_harm
```

where:
- div_Γ J_loop = 0  (solenoidal)
- curl_Γ J_star = 0 (irrotational)
- harmonics span the cohomology

## Why this fixes LF breakdown

In the EFIE block matrix on (Loop, Star) basis:

```
Z_EFIE = [ Z_LL    Z_LS  ]
         [ Z_SL    Z_SS  ]
```

- Z_LL scales O(ω)  (loop-loop, A-term only)
- Z_SS scales O(1/ω) (star-star, φ-term only)
- Z_LS scales O(ω⁰) (off-diagonal)

Rescale: multiply Loop block by 1/√ω, Star by √ω → all blocks O(1).

After rescaling, EFIE is well-conditioned at all ω.

References (this folder):
- 05_low_freq_stabilization/LoopStar_Basis_Robust_Prec_EFIE.pdf
- 05_low_freq_stabilization/LoopStar_Decomposition_EFIE.pdf
- 05_low_freq_stabilization/A_Study_Recent_MoM_Accurate_Very_Low_Freq.pdf

## In Radia PEEC

The PEEC framework uses Loop-Star **naturally**:
- Loops = filament currents (one per closed current loop)
- Stars = panel charges (one per panel)

No transformation needed — they're independent unknowns by construction.

```python
from peec_matrices import PyPEECBuilder
# Loops: add_connected_segment defines current path (loop)
# Stars: panel basis on faces
```

See `radia_mcp.peec` for code usage.
"""


CALDERON = r"""
# Calderón Preconditioner — Andriulli et al. 2008

Reference: F.P. Andriulli et al., "A Multiplicative Calderón Preconditioner
for the Electric Field Integral Equation", IEEE Trans. Antennas Propag.
56(8):2398-2412, 2008.  DOI: 10.1109/TAP.2008.926788.

## The idea

The EFIE operator T satisfies the Calderón identity:

```
T · T = -(I/4) + K_compact   (compact operator)
```

So applying T as a left-preconditioner to itself gives a
well-conditioned system:

```
T · (T · J) = T · b
```

→ The Calderón-preconditioned matrix has bounded condition number
independent of mesh refinement AND independent of frequency.

## Practical implementation

Need to discretize T on TWO bases simultaneously:
1. RWG for J on the primal mesh
2. **BC (Buffa-Christiansen)** basis on the dual barycentric mesh

The BC basis is the rotated dual of RWG.  Multiplication T_BC · T_RWG
gives the preconditioned system.

## When the lab uses it

- High-frequency scattering with ngsolve.bem (not lab primary)
- NOT for low-freq eddy current (use Weggler ST or PEEC instead)

## Tradeoff

| Approach | Setup cost | Memory | Robust at LF? | Robust at HF? |
|----------|------------|--------|---------------|---------------|
| Plain EFIE | O(N²) | O(N²) | NO | YES (mostly) |
| Loop-Star EFIE | + O(N) re-basis | same | YES | YES |
| Calderón EFIE | 2× system | 2× | YES | YES |
| PEEC + Loop-Star | by construction | sparse + low-rank | YES | LF only |

For the lab's MQS regime (DC to ~100 kHz), PEEC + Loop-Star natural
decomposition is the production choice.  Calderón is documented for
completeness when ngsolve.bem is used at higher frequencies.
"""


SINGLE_TRACE = r"""
# Single-Trace Formulation (Weggler-Hiptmair)

For low-frequency stable surface BEM, a modern alternative to
Loop-Star and Calderón is the **single-trace formulation** used in
NGSolve's bem module.

## The idea

Instead of carrying both J = n̂ × H and M = E × n̂ as unknowns (as in
PMCHWT), define a single "trace" of the field — typically:

```
T = (n̂ × E, n̂ × H)   on Γ
```

Combine the two traces via a compatibility condition that absorbs the
scale mismatch.  The resulting system is:
- One unknown vector per surface DOF
- Well-conditioned at DC (no LF breakdown)
- Suitable for low-conductor problems

## Why ngsolve.bem uses this

The ngsolve.bem module (developed by L. Weggler, R. Hiptmair, and
collaborators at TU Wien/ETH) ships single-trace as the **default**
formulation for low-freq scenarios.

```python
from ngsolve.bem import (SingleLayerPotentialOperator,
                          DoubleLayerPotentialOperator)
# ngsolve.bem builds single-trace automatically when frequency is low
```

## In the lab workflow

The lab uses `ngsolve.bem` (single-trace) for:
- PEEC inductance validation against full BEM
- Surface-impedance Robin BC on workpiece (lab: IH path)
- Open-boundary Laplace problems

See `radia_mcp.radia_ngsolve.ngsbem_inductance` for code recipes.

## References

[LOCAL] 05_low_freq_stabilization/Boundary element methods for magnetostatic field problems.pdf
[LOCAL] 05_low_freq_stabilization/Far_Field_Scattering_Low_Freq.pdf

Plus Ostrowski-Hiptmair 2021 Two-Step Maxwell which extends this idea
to full Maxwell with tree-cotree gauging.

## Decision: which LF stabilization?

| Use case | Recommended |
|----------|-------------|
| Lab PEEC | Loop-Star (natural) |
| Lab IH FEM-SIBC | (none — already FEM, not surface BEM) |
| Lab ngsolve.bem | Single-trace Weggler-Hiptmair (built-in) |
| External research code | Calderón (most rigorous) |
| Hand-coded EFIE | Loop-Star (simpler than Calderón) |
"""


def get_low_freq_knowledge(topic: str = "problem") -> str:
    """Dispatch low-frequency BEM stabilization topics.

    Topics:
        problem         - Why LF breakdown happens (DEFAULT)
        loop_star       - Loop-Star decomposition
        calderon        - Calderón preconditioner
        single_trace    - Weggler-Hiptmair single-trace (ngsolve.bem)
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("problem", "overview", "breakdown"):
        return PROBLEM
    if topic in ("loop_star", "loopstar", "loop-star"):
        return LOOP_STAR
    if topic in ("calderon", "calderón"):
        return CALDERON
    if topic in ("single_trace", "singletrace", "weggler"):
        return SINGLE_TRACE
    if topic == "all":
        return "\n\n".join([PROBLEM, LOOP_STAR, CALDERON, SINGLE_TRACE])
    return (f"Unknown topic '{topic}'. Available: problem, loop_star, "
            "calderon, single_trace, all.")
