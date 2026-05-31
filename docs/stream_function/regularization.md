# Regularisation choices for FE-direct ψ

The SF problem `A ψ = B_target` with `M = 25` constraints and
`ndof ≈ 1773` (H¹ order 2) is **massively underdetermined** — there
are infinitely many ψ vectors that hit `B_target` exactly.  Choosing
*which* such ψ to return is the role of the **regularisation**.

Five options shipped, all in
[`demo_planar_uniform_fem_psi_advanced.py`](../../examples/stream_function/demo_planar_uniform_fem_psi_advanced.py)
via `--regularize {l2, l2_aca, h1, h1_sigma, inductance_diag, linf}`.

## Quick comparison table (planar uniform Bz benchmark)

| Mode               | Math                                       | RMS      | p2p/mean | Contours | Wall  | Physical meaning                |
|--------------------|--------------------------------------------|----------|----------|----------|-------|---------------------------------|
| `l2`               | `min ‖ψ‖²` s.t. `Aψ=B`                     | 1.78 %   | 7.59 %   | 56       | 0.05s | numpy lstsq, default            |
| `l2_aca`           | same via HACApK ACA+TSVD                   | 1.78 %   | 7.59 %   | 56       | 0.10s | validates ACA+ contract         |
| `h1`               | `min ∫\|∇ψ\|²` s.t. `Aψ=B`                  | 2.09 %   | 6.81 %   | 27       | 1.0s  | min surface current L² norm     |
| `h1_sigma`         | `min ∫(1/σ)\|∇ψ\|²` s.t. `Aψ=B`             | 1.17 %\* | 4.01 %\* | 44       | 1.0s  | true ohmic dissipation          |
| `inductance_diag`  | `min ∑ wᵢ \|∇φᵢ\|²` with `wᵢ ~ size(elᵢ)` | 2.46 %   | 7.75 %   | 27       | 1.0s  | self-inductance proxy           |
| `linf`             | min residual s.t. `\|∇ψ\|_∞ ≤ J_max`        | exp.     | exp.     | depends  | slow  | hot-spot suppression            |

\* For `--sigma-cf "exp(-(x*x + y*y) / 0.02)"` (Gaussian concentration).

## When to pick which

### `l2` / `l2_aca` (default)

When you don't care about smoothness, just want the lowest-norm `ψ`
that hits `B_target`.  `l2_aca` routes through HACApK ACA+TSVD instead
of numpy lstsq — at our M=25 scale this is *slower* (1.12s vs 0.92s)
because ACA+ overhead dominates, but the win is when M >> 25 or each
matrix entry is expensive (= material kernels via Radia MMM, where
each entry is a Radia container solve).

### `h1` — minimum surface current density

```python
min ∫_Γ |∇ψ|² dA   s.t.   A ψ = B
```

Since `|K| = |∇ψ|` (with `K = n̂ × ∇ψ`), this minimises
`∫_Γ |K|² dA` — the L² norm of the surface current density.  This
is the **physically motivated default**: minimum surface current
magnitude.

Implementation (closed-form via Lagrangian):

```python
S = H1 stiffness matrix   = ∫ ∇φᵢ · ∇φⱼ dA
ψ = S⁻¹ Aᵀ (A S⁻¹ Aᵀ)⁻¹ B   # smallest H¹ seminorm hitting B exactly
```

The `(A S⁻¹ Aᵀ)` system is `M × M` (small); the heavy lift is
`solve(S, Aᵀ)` which is `ndof × M`.  No α to tune.

### `h1_sigma` — true ohmic dissipation

```python
min ∫_Γ (1/σ(x, y)) |∇ψ|² dA   s.t.   A ψ = B
```

Weights the H¹ seminorm by `1/σ(x, y)`, the local resistivity.  When
`σ` is large somewhere, the optimiser will preferentially place
current there (low dissipation).  When `σ → 0`, current is pushed AWAY
from that region.

Use cases:

  - **Forbidden regions**: `σ(x, y) = 1 / (1 + 100·IfPos(x²+y²-0.04, 1, 0))`
    pushes all current inside the disc of radius 0.2.
  - **Litz wire / mask**: variable resistivity from a 2D conductivity
    map.
  - **Active shielding**: spatially-varying penalty.

The expression is parsed in `eval()` with `IfPos`, `exp`, `sqrt`,
`sin`, `cos`, `log` and the spatial CFs `x`, `y` in scope.

### `inductance_diag` — self-inductance proxy

Diagonal approximation of the **true self-inductance form**

    L = (μ₀ / 4π) ∫∫ ∇ψ(x) · ∇ψ(y) / |x − y| dA(x) dA(y)

by weighting `|∇ψ|²` by a per-element characteristic length
`size(element)` (NGSolve `specialcf.mesh_size`).  Captures the
leading-order self-energy contribution; the full off-diagonal form
needs surface BEM (= ngsolve.bem `MaxwellSL` on a thin 3D shell, TODO).

Use as "physics-flavoured smoothing" rather than a strict inductance
optimum.  Not actually better than `h1` numerically.

**Full inductance implementation path** (open extension): use
ngsolve.bem 6.2.2604+ `MaxwellSingleLayerPotentialOperator` on a thin
3D shell embedding the source plane; assemble the bilinear form
`∫∫ K(x)·K(y)/(4π|x-y|) ds(x) ds(y)`; combine with the constraint
`A ψ = B` in a Lagrangian.  ~1 week.

### `linf` — peak current cap (experimental)

```python
min ‖A ψ − B‖²    s.t.    max_vertex |∇ψ(v)| ≤ J_max
```

Cap the peak surface current density to avoid hot spots.  Solved via
scipy SLSQP with per-vertex quadratic constraints; warm-started from
the `h1` solution.

**EXPERIMENTAL**: SLSQP at ndof = 1773 H¹ order-2 is slow and may not
satisfy the tight cap.  Tune with `--order 1 --maxh 0.05` for smaller
problem.  Practical fix would be IRLS or ADMM (1 day work).

## Combining with Path-A iteration

All regularisation modes compose with `--compensated-iter` from
[`demo_planar_uniform_fem_psi.py`](../../examples/stream_function/demo_planar_uniform_fem_psi.py).
Path-A iteration RE-USES the inner solver's cached factorisation:

```bash
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1 --order 3 \
    --compensated-iter 100 --compensated-step 0.05
```

Best results combine `h1` regularisation + `order=3` + Path-A.

## High-order ψ sweep

| order p | ndof    | RMS     | p2p/mean | Notes                          |
|---------|---------|---------|----------|--------------------------------|
| 1       | ~500    | 1.02 %  | 3.77 %   | piecewise linear, surprisingly OK |
| 2       | ~1773   | 2.09 %  | 6.81 %   | default                        |
| **3**   | ~3700   | **0.51 %** | **1.83 %** | **sweet spot**             |
| 4       | ~6000   | 2.04 %  | 7.65 %   | regression                     |
| 5       | ~9000   | 2.09 %  | 7.27 %   | plateau                        |

**Non-monotone in p**: the FE smoothness × `nlevels=12` contour
extraction × single-stroke chain interaction has a resolution-matched
sweet spot at p=3 for our default mesh.  This is a benchmark-class
observation, not a generic FE convergence statement.

## Design rule: deformation + near-optimal baseline can backfire

Empirical (planar uniform Bz, same target):

  - order=3 H¹ baseline:                       RMS 0.51 %
  - order=3 H¹ + bump deform (20 trials, 27s): RMS 0.82 %    **(worse!)**
  - order=2 H¹ baseline:                       RMS 2.09 %
  - order=2 H¹ + bump deform (20 trials, 22s): RMS 0.77 %    (better)

When baseline is already sub-0.5 %, the CMA-ES outer loop wastes its
budget searching for a sub-percent improvement and gets stuck.  When
baseline is several percent, deformation helps a lot.

  **RULE**: turn deformation OFF when single-shot accuracy is already
  better than what your CMA-ES budget can credibly improve on.

## Regularisation folded into the ACA+TSVD factorisation

All five regularisations above (and any other SPD `S`) admit a
**single closed form** built on top of an ACA+TSVD of `A`:

> **`ψ = S⁻¹ V · W⁻¹ · Σ⁻¹ · Uᵀ B,    W = Vᵀ S⁻¹ V (k × k)`**

This is the formula implemented by
[`radia.stream_function.RegularizedTSVD`](../../src/radia/stream_function.py)
and exercised by
[`demo_regularized_aca.py`](../../examples/stream_function/demo_regularized_aca.py).

### Derivation

The constrained min-norm problem

    min   ψᵀ S ψ      s.t.    A ψ = B

has the Lagrangian solution `ψ = S⁻¹ Aᵀ (A S⁻¹ Aᵀ)⁻¹ B`.  Substitute
the truncated SVD `A = U Σ Vᵀ` (with `V` having `k` orthonormal columns
from ACA+TSVD):

    Aᵀ           = V Σ Uᵀ                        (N × M)
    S⁻¹ Aᵀ       = (S⁻¹ V) Σ Uᵀ                  (N × M)
    A S⁻¹ Aᵀ     = U Σ (Vᵀ S⁻¹ V) Σ Uᵀ           (M × M, rank k)
                 = U Σ W Σ Uᵀ
    (A S⁻¹ Aᵀ)†  = U Σ⁻¹ W⁻¹ Σ⁻¹ Uᵀ              (Moore-Penrose; rank k)
    ψ            = (S⁻¹ V) Σ Uᵀ · U Σ⁻¹ W⁻¹ Σ⁻¹ Uᵀ B
                 = (S⁻¹ V) · W⁻¹ · Σ⁻¹ · Uᵀ B    (because Uᵀ U = I_k)

So everything reduces to two precomputations:

  1. `Sinv_V = S⁻¹ V` — one sparse / dense solve with `S` and `k` RHS.
  2. `W_inv = (Vᵀ Sinv_V)⁻¹` — one `k × k` inverse.

After that, **each solve costs `O(k · (M + N))`** — a few matvecs and
one `k × k` matvec-with-W⁻¹, no further sparse / dense system solve.

### Sanity check: `S = I` reduces to the standard pseudo-inverse

For `S = I`, since the TSVD has orthonormal `V` columns,
`W = Vᵀ V = I_k`, so `W⁻¹ = I_k` and the formula collapses to

    ψ = V · Σ⁻¹ · Uᵀ B

which is exactly `pseudo_inverse_solve(result, B)`.  This is verified
to machine precision (relative diff `≈ 1e-15`) in
`demo_regularized_aca.py` and the unit-test sanity script in the
implementation history.

### Why this matters

The Path-A compensated iteration re-uses the same `A` and the same
regularisation `S` across `k = 50 – 200` outer iterations, with only
the right-hand side `B = B_target − Iₖ · Bz_chain^(k)` changing each
time.  In the direct `solve_h1` formulation, each iteration costs a
fresh `(A S⁻¹ Aᵀ)` solve at the `O(N² + N · ndof²)` scale.  With the
folded form, each iteration is **a few microseconds** of dense linear
algebra on the cached `(Sinv_V, W_inv)` pair.

| Form                                | Path-A inner-solve cost (per iter) |
|-------------------------------------|------------------------------------|
| Direct `solve_h1`                   | `O(N³)` — N = 1773 free DOFs       |
| `RegularizedTSVD.solve(B)` (cached) | `O(k · N) ≈ 0.04 ms` at k = 25     |

Measured on the planar benchmark: 100-iter Path-A with direct H1 is
~50 s (most spent in 100 dense LU solves); with the cached form,
total iteration wall time drops to ~1 s.

Same value across:

  - **Regularisation sweep**: scan multiple `S` choices on the same
    `A`.  Each new `S` costs one `Sinv_V` + one `k × k` inverse, then
    instant `solve(B)`.
  - **Path-A iteration**: factorise once, fire `solve(B)` per outer
    iter at zero marginal cost.
  - **Material kernels** (Radia MMM, shielded coil):  each `A(i, j)`
    is a `rad.Solve()` + `rad.Fld()` call, so the **ACA+ amortisation
    matters most** — and once it has run, all regularisation choices
    plug into the same cache.

### IRLS for L∞ on top of the cache

The bounded-ratio active-set IRLS (`solve_linf_irls` in
[`demo_planar_uniform_fem_psi_advanced.py`](../../examples/stream_function/demo_planar_uniform_fem_psi_advanced.py))
replaces the slow scipy SLSQP path:

  - ACA+TSVD of `A` runs ONCE.
  - Each IRLS iter rebuilds `S^(k)` (weighted H¹ stiffness with
    per-vertex penalty depending on the previous iterate's gradient
    distribution) and a fresh `RegularizedTSVD.from_stiffness(res,
    S^(k))` cache.  The vertex weights are **bounded-ratio**
    (`weight_ratio_max`, default 20) and **damped** (default 0.5)
    to prevent the classical pure-Lawson IRLS oscillation.
  - An optional final scalar rescale guarantees `max_v |∇ψ| ≤ jmax`
    exactly.

On the planar uniform-Bz benchmark, the H¹ baseline already has
peak/mean gradient `≈ 1.47` (nearly L∞-uniform), so the IRLS gives
only modest (10 – 30 %) peak reduction; the value of the cached form
is more visible on geometries where the H¹ baseline has hot spots
(e.g. surface-deformed coils, asymmetric targets).  The framework
scales without re-factorising `A`.

## Optimising the regularisation SHAPE (fixed surface, ACA+ reused)

The deformation outer loop ([deformation.md](deformation.md)) moves
the *surface* — so `A` changes every Optuna trial and the whole
`RegularizedTSVD` is rebuilt.  The dual move is to keep the surface
fixed and optimise *where the current is allowed to flow* by varying
the conductivity field `σ(x, y)` in the `1/σ`-weighted H¹ seminorm.
Now `A` is **constant**, so the ACA+TSVD base `(U, Σ, V)` is computed
ONCE and reused across every trial — only the cheap fold
`S⁻¹V + W = Vᵀ S⁻¹ V` is rebuilt per trial.

[`demo_reg_hyperparam_aca.py`](../../examples/stream_function/demo_reg_hyperparam_aca.py)
runs CMA-ES over a Gaussian conductivity feature

    σ(x, y) = 1 + amp · exp(−((x−cx)² + (y−cy)²) / width)

(parameters `amp, cx, cy, width`), with the inner solve
`min ∫(1/σ)|∇ψ|² s.t. A ψ = B` folded through the cached base.

### Result (planar uniform Bz, 30 trials, 15 s)

| Quantity            | Uniform σ (= plain H¹) | Optimised σ |
|---------------------|------------------------|-------------|
| single-stroke chain RMS | 2.09 %             | **0.73 %**  |
| wire length         | 18.98 m                | 18.91 m     |

The optimiser finds `amp ≈ −0.81` (an ~80 % conductivity DIP at the
centre), which pushes surface current toward the periphery — the
Helmholtz-like distribution that produces a more uniform central
`Bz`.  This drops chain RMS by **2.9×** with NO geometry change, NO
new target, and NO surface manufacturing cost — purely by reshaping
the regularisation penalty.  It is comparable to what the surface
deformation loop achieves (0.77 %), via an orthogonal and physically
cheaper degree of freedom.

### What the cache buys here

```
ACA+TSVD base (U, Sigma, V):   computed ONCE  (~33 ms)
per trial:  build S(sigma)  ->  S^-1 V  ->  W = V^T S^-1 V  ->  solve
            (~43 ms fold; ACA+ base NOT recomputed)
```

At FE-direct scale (cheap `A(i, j)`) the ~33 ms ACA+ base is small, so
the wall-clock win is modest.  The pattern is what matters: when each
`A(i, j)` is a **material-kernel** evaluation (Radia `rad.Solve()` +
`rad.Fld()`, ~100× slower per entry), the ACA+ base dominates and
reusing it across a 30-trial Optuna sweep is the entire performance
argument for the folded form.

### Cross-reference

  - Theory: [theory.md](theory.md)
  - Single-stroke discretisation: [single_stroke.md](single_stroke.md)
  - Surface deformation: [deformation.md](deformation.md)
  - API: [`radia.stream_function.RegularizedTSVD`](../../src/radia/stream_function.py)
  - Demos: [`demo_regularized_aca.py`](../../examples/stream_function/demo_regularized_aca.py)
    (5-mode sweep), [`demo_reg_hyperparam_aca.py`](../../examples/stream_function/demo_reg_hyperparam_aca.py)
    (σ-shape optimisation, ACA+ reused)
  - MCP topic: `aca_tsvd(topic=regularized)`
