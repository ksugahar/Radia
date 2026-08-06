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
from radia.peec_matrices import PEECBuilder
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


WEGGLER = r"""
# Lucy Weggler's Stabilized Maxwell EFIE (ngsolve.bem)

The **canonical low-frequency BEM stabilisation in ngsolve.bem** is
Lucy Weggler's **product-space formulation** `HDivSurface × SurfaceL2`
with `κ² V_κ` re-scaling.  It is the LF tool the lab uses through
`ngsolve.bem`.

Upstream reference: [Weggler/docu-ngsbem `Maxwell_DtN_Stabilized.ipynb`](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb)

Lab notebook (Sugahara 2026-02, applies Weggler to PEEC + coil
inductance + Calderón preconditioner):
`to_developers/ngsolve/low_freq_efie_ngbem_applications.ipynb`

## The pathology Weggler fixes

The classical EFIE block matrix is

```
[ A_κ    Q_κ   ] [ J ]   [ b_J ]
[ Q_κᵀ   V_κ/κ²] [ Q ] = [ b_Q ]
```

The `V_κ / κ²` term blows up as `κ → 0` -- the matrix becomes
ill-conditioned `O(κ⁻²)`, Krylov fails, the inductance / DC limit is
unreachable from the EFIE solve.

## The Weggler fix: product space + κ² scaling

Reformulate as a product space:

| Component         | Function space  | Physics             |
|-------------------|-----------------|---------------------|
| solenoidal (`J`)  | `HDivSurface`   | loop / inductance   |
| irrotational (`Q`)| `SurfaceL2`     | charge / potential  |

and rescale the `V_κ` block by `κ²` so the system becomes

```
[ A_κ        Q_κ      ]
[ Q_κᵀ       κ² V_κ   ]
```

with **`O(1)` condition number for all `κ`** -- including `κ → 0`.

This is mathematically equivalent to a Loop-Star decomposition but
achieves the separation at the **variational** level, no explicit
basis transformation needed.

## Crucial `κ → 0` reduction (the lab uses this directly)

In the static limit the `A_κ` block becomes `LaplaceSL` on
`HDivSurface`.  For a genus-1 surface (e.g. a coil), the harmonic
subspace (div-free + curl-free) of `HDivSurface` is 1-dimensional --
exactly the circulating-current mode.  Project `LaplaceSL` onto it
and you get the air-core self-inductance:

```python
from ngsolve.bem import LaplaceSL
from ngsolve import HDivSurface
fes  = HDivSurface(coil_surface_mesh, order=2)
L_op = LaplaceSL(fes)          # κ = 0 limit of A_κ
# Hodge decomposition: find the 1D harmonic subspace
H, _ = scipy.linalg.null_space(...)   # see notebook §2
L    = H.T @ L_op @ H          # 1x1 inductance matrix
```

Lab notebook §2 verifies this gives 48.7 nH for a flat circular ring
R = 10 mm, w = 1 mm vs analytical 48.8 nH (0.3 % agreement).

## ngsolve.bem operators that wrap Weggler

| Operator                                            | Role                                                   |
|-----------------------------------------------------|--------------------------------------------------------|
| `HelmholtzSL(HDivSurface)`                          | Vector SLP `A_κ` (solenoidal current)                  |
| `HelmholtzSL(SurfaceL2)`                            | Scalar SLP `V_κ` for stabilized block (κ² rescaling)   |
| `HelmholtzSL(div(HDivSurface), SurfaceL2)`          | Coupling `Q_κ`                                         |
| `LaplaceSL(HDivSurface)`                            | `κ = 0` limit of `A_κ`, harmonic-subspace inductance   |

## Where the lab uses Weggler

| Scenario                                          | Use Weggler? |
|---------------------------------------------------|--------------|
| PEEC inductance validation via `ngsolve.bem`      | YES          |
| Coil self-inductance from genus-1 surface mesh    | YES (`LaplaceSL`+Hodge — notebook §2) |
| IH FEM-SIBC                                       | NO (FEM-side, not BEM solve)          |
| Equivalence-source one-shot evaluator (Stratton-Chu integrate, no solve) | NO direct use -- see `equivalence_source_lf` topic |

## Companion: Calderón preconditioner

Schoeberl + Andriulli's rotated-SL Calderón preconditioner
(`(n̂ ×) ∘ S ∘ (n̂ ×)`) is **complementary** to Weggler's product-
space stabilization:

- Weggler restructures the system (`O(1)` condition number)
- Calderón accelerates the iterative solve on the standard system

Both are demonstrated side-by-side in the lab notebook (§1: Weggler,
§3: Calderón).  See `bem_low_freq('calderon')`.

## References

1. L. Weggler, "Stabilized Maxwell DtN Formulation,"
   ngsbem demo, 2026.
2. J. Ostrowski et al., "ngbem -- A BEM Library for NGSolve," 2024.
3. V. Giunzioni et al., "Low-Frequency Stabilizations of the PMCHWT
   Equation for Dielectric and Conductive Media,"
   *IEEE TAP* 73(8):5725--5740, 2025.  arXiv:2408.01321.
4. S. B. Adrian, A. Dely, D. Consoli, A. Merlini, F. P. Andriulli,
   "Electromagnetic Integral Equations: Insights in Conditioning and
   Preconditioning," *IEEE OJAP* 2:1143--1174, 2021.

[LOCAL] `to_developers/ngsolve/low_freq_efie_ngbem_applications.ipynb`
[LOCAL] `to_developers/ngsolve/NGSolve-BEM_ stabilized Maxwell.pdf`

## Decision: which LF stabilization?

| Use case                          | Recommended                                |
|-----------------------------------|--------------------------------------------|
| Lab PEEC                          | Loop-Star (natural, inherent in PEEC)      |
| Lab IH FEM-SIBC                   | (none — already FEM, not surface BEM)      |
| Lab ngsolve.bem (any LF solve)    | **Weggler product space** (the answer)     |
| Equivalence-source one-shot eval  | Laplace-ML routing (NOT Weggler — see `equivalence_source_lf`) |
| Hand-coded EFIE                   | Loop-Star (simpler than Calderón)          |
| Calderón preconditioner           | Complementary to Weggler — see `calderon`  |
"""


EQUIV_SOURCE_LF = r"""
# Low-frequency equivalence-theorem source — Weggler is NOT the answer

For the equivalence-theorem one-shot evaluator
(`radia.equivalence_source.NearFieldSource`) at low frequency, the
LF stabilisation story is **different from** what `ngsolve.bem`
itself needs.  Two LF "breakdowns" exist; only one of them affects
the equivalence-source evaluator.

## Two distinct LF breakdowns

| Pathology         | Stage                | Fix                              |
|-------------------|----------------------|----------------------------------|
| **Operator-level** `V_κ / κ²` blowup of EFIE → matrix ill-conditioned `O(κ⁻²)` | BEM **solve** (BilinearForm + GMRES) | Lucy Weggler product space (`bem_low_freq('weggler')`) |
| **FMM-level** `J_n(kR) → 0` underflow → multipole expansion fails | Kernel **summation** (multilevel expansion) | Greengard-Huang-Rokhlin LF-FMM, or simpler: route through Laplace ML |

## Why Weggler does NOT directly apply

The equivalence-source evaluator does **not solve a system**.  It
takes a known surface (E, H) and evaluates the Stratton-Chu integral
at obs points.  The classical EFIE matrix-conditioning pathology
that Weggler fixes is irrelevant -- there is no matrix to invert.

## What DOES apply: Laplace ML routing

The equivalence-source evaluator runs into the FMM-level breakdown
when `kR ≪ 1` (the typical Radia low-frequency regime, e.g. IH 10
kHz with R ≈ 1 m gives `kR ~ 2×10⁻⁵`).  Standard Helmholtz FMM
loses precision; the right move is:

1. Route the **real part** of the Stratton-Chu integral through the
   Laplace multilevel expansion (`BiotSavartRegularMLCF` /
   `LaplaceSL`) -- Laplace ML has no breakdown at any frequency.
2. Add the **imaginary part** (`O(kR)` smaller than the real part
   at `kR ≪ 1`) via the direct C++ kernel
   (`_EquivalenceSourceHarmonic`).  Cheap when ω is small.

Conservative cutoff: `kR_max < 0.01`.

Verified: the C++ Phase B kernel at `ω = 1 Hz` (kR ~ 10⁻⁸) agrees
with the Phase A static kernel (`evaluate_static_H`) to 5e-15
relative diff -- machine precision -- so the static-limit reduction
is exact, the imaginary correction is genuinely tiny, and the
Laplace ML + direct hybrid is sound.

## Conceptual unity with Weggler

Both Weggler's product space and the equivalence-source Laplace-ML
routing exploit the same underlying observation:

> **At low frequency, the static (Laplace) physics dominates and
> the wave correction is a small perturbation.**

Weggler's `A_κ` block reduces to `LaplaceSL(HDivSurface)` as
`κ → 0` (notebook §2 -- gives exact air-core inductance).  The
equivalence-source LF rule reduces the harmonic Stratton-Chu integral
to its Laplace ML form as `kR → 0` (same physics, evaluation side).

**Weggler is the SOLVE-side analogue of equivalence-source's
EVALUATE-side Laplace-ML routing.**  Same insight, different
machinery, both correct.

## When Weggler DOES matter for equivalence_source

There is one situation: if a future user wants to **compute** the
surface (E, H) from a **BEM solve** (instead of an FEM solve), the
underlying `ngsolve.bem` solve should use Weggler's product-space
formulation at low frequency.  Today the design assumes the surface
data comes from an FEM solve (e.g. `calc_fem_kelvin.py`), so Weggler
lives on the FEM-coupling side and the equivalence-source layer
never touches it.

## See also

- `bem_low_freq('weggler')` — full Weggler product-space details
- `bem_low_freq('problem')` — the underlying `1/κ²` pathology
- `docs/equivalence_source/FMM_DESIGN.md` §3.3-3.5 — the two-LF-
  breakdown distinction and Laplace-ML routing rule
- `docs/equivalence_source/USER_GUIDE.md` — when to use static vs
  harmonic kernel
"""

# Backwards-compat alias -- legacy single_trace topic name kept
# pointing at the (more accurate) Weggler product-space content.
SINGLE_TRACE = WEGGLER


def get_low_freq_knowledge(topic: str = "problem") -> str:
    """Dispatch low-frequency BEM stabilization topics.

    Topics:
        problem               - Why LF breakdown happens (DEFAULT)
        loop_star             - Loop-Star decomposition
        calderon              - Calderón preconditioner
        weggler               - ★ Lucy Weggler's product-space EFIE
                                stabilization (ngsolve.bem canonical
                                LF fix; HDivSurface × SurfaceL2 with
                                κ² V_κ re-scaling).
        single_trace          - Alias for weggler (legacy name).
        equivalence_source_lf - LF for the equivalence-theorem
                                one-shot evaluator: Weggler does NOT
                                apply (no solve); Laplace ML routing
                                does.  Two-breakdown distinction.
        all                   - Everything
    """
    topic = topic.lower().strip()
    if topic in ("problem", "overview", "breakdown"):
        return PROBLEM
    if topic in ("loop_star", "loopstar", "loop-star"):
        return LOOP_STAR
    if topic in ("calderon", "calderón"):
        return CALDERON
    if topic in ("weggler", "single_trace", "singletrace",
                 "product_space", "stabilized_efie"):
        return WEGGLER
    if topic in ("equivalence_source_lf", "equivalence_source",
                 "nfs_lf", "near_field_source_lf", "stratton_chu_lf",
                 "two_breakdowns", "lf_fmm"):
        return EQUIV_SOURCE_LF
    if topic == "all":
        return "\n\n".join([PROBLEM, LOOP_STAR, CALDERON, WEGGLER,
                            EQUIV_SOURCE_LF])
    return (f"Unknown topic '{topic}'. Available: problem, loop_star, "
            "calderon, weggler, single_trace (alias), "
            "equivalence_source_lf, all.")
