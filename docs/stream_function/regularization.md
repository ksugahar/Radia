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

### Tikhonov (soft data fit) is the SAME core with `+ α I`

The closed form above hits `A ψ = B` **exactly** — it is the `α = 0`,
equality-constrained limit.  Generalised **Tikhonov** instead trades
data misfit against the seminorm,

    min_ψ   ‖A ψ − B‖²  +  α · ψᵀ S ψ          (α > 0)

with normal equation `(AᵀA + α S) ψ = Aᵀ B`.  In the demo
([`demo_sf_to_peec_gx.py`](../../examples/stream_function/demo_sf_to_peec_gx.py))
`--regularize tikhonov` currently solves this as a **separate dense**
`(AᵀA + α I)⁻¹ AᵀB`, disjoint from the `tsvd` ACA+ path.  **It does not
need to be a separate path.**  Substitute the same ACA+TSVD
`A = U Σ Vᵀ` (so `AᵀA = V Σ² Vᵀ`) and look for the solution in the
`S⁻¹V` subspace, `ψ = S⁻¹V y`.  Because `Aᵀ B = V Σ UᵀB ∈ range(V)` and
`V` has orthonormal columns, the `N × N` system collapses to a `k × k`
one (left-multiply by `Vᵀ`, use `VᵀV = I`):

    (α I + Σ² W) y = Σ Uᵀ B,        W = Vᵀ S⁻¹ V      (k × k)

> **`ψ(α) = (S⁻¹V) · (α I + Σ² W)⁻¹ · Σ · Uᵀ B`**

— this is the `α = 0` exact-fit formula with a single `α I` added inside
the small core.  Everything expensive (`U, Σ, V` from ACA+, and the fold
`S⁻¹V`, `W`) is **computed once**; sweeping `α` for an L-curve re-solves
only the `k × k` core (sub-millisecond), reusing the one ACA+
factorisation that `tsvd` already built.  `α I + Σ² W` is invertible for
every `α > 0` (it is similar to `α I + Σ W Σ`, SPD), so no extra care is
needed.

#### One filter, four special cases

The core `(α I + Σ² W)⁻¹ Σ` is just a **spectral filter** on the shared
ACA basis:

| Case                | Reduces to                                       | Filter on mode `i`              |
|---------------------|--------------------------------------------------|---------------------------------|
| `α = 0`             | `S⁻¹V · W⁻¹ · Σ⁻¹ · UᵀB` (exact-fit seminorm)     | `1/σᵢ` (invert every mode)      |
| `S = I`, `α > 0`    | `V · diag(σ/(σ²+α)) · UᵀB` (standard Tikhonov)    | `σᵢ/(σᵢ²+α)` (smooth roll-off)  |
| `S = I`, `α = 0`    | `V · Σ⁻¹ · UᵀB` (plain TSVD pseudo-inverse)       | `1/σᵢ`                          |
| `+ k_mode = r`      | hard-truncate to the top `r` modes               | step `1` (kept) / `0` (dropped) |

So the three "different" regularisers are one object: **TSVD** is the
**hard** spectral cut (drop `σ < tol`), **Tikhonov** is the **smooth**
version of the same idea (`σ/(σ²+α)` rolls each mode off gradually), and
the **seminorm** `S` rotates the filter into the `S`-metric via `W`.  Two
knobs (`α` smooth, `k_mode` hard) and one metric (`S`) — not three
solvers.

#### Why fold it (vs the separate dense Tikhonov solve)

  - **One factorisation for the whole L-curve.**  The α-sweep that locates
    the L-curve corner re-solves only the `k × k` core per α; the dense
    path re-factorises an `N × N` matrix per α.
  - **Stability is explicit and free.**  `A` here is badly conditioned —
    on the planar gradient case σ runs `[1.5e-10, 2.3e-6]`.  Exact-fit
    (`α = 0`) inverts the `1.5e-10` mode (≈ 7e9 amplification of any
    target noise); Tikhonov with `α ≫ σ_min²` **damps** it.  The fold
    makes that a one-line knob: raise `α`, the tiny-σ modes roll off, no
    re-factorise.
  - **Composes with everything.**  Hard-truncate to `k_mode` modes AND
    soft-damp the survivors with `α` AND measure the norm in the H¹ (or
    `1/σ`, or inductance) metric — all on the one cached
    `(U, Σ, V, S⁻¹V, W)`.

Verified numerically: machine precision against the explicit filter
`σ/(σ²+α)` for `S = I` (`≈ 1e-16`), the ACA tolerance against the dense
`(AᵀA + α S)⁻¹ AᵀB` for general `S` (`≈ 2e-5`, = the ACA+ compression
residual of `A`, **constant in α**), and the `α → 0` limit recovers
`RegularizedTSVD.solve` once `α ≪ σ_min²` (`≈ 1e-12` at `α = 1e-30`).

#### Code

`RegularizedTSVD.solve(B)` implements the `α = 0` core
(`Sinv_V @ (W_inv @ (UᵀB / Σ))`).  The Tikhonov core is the **same cached
pieces** plus one `α`:

```python
# folded Tikhonov on the SAME cached ACA+TSVD factorisation (alpha > 0)
U, Sigma, V = base.U, base.S, base.V          # from the ACA+TSVD of A
W    = V.T @ reg.Sinv_V                        # k x k  (= inv(reg.W_inv))
core = alpha * np.eye(k) + (Sigma**2)[:, None] * W
y    = np.linalg.solve(core, Sigma * (U.T @ B))
psi  = reg.Sinv_V @ y                          # == (A^T A + alpha S)^-1 A^T B
```

i.e. a ~6-line `solve(B, alpha=0.0)` extension of `RegularizedTSVD` —
`alpha = 0` keeps the present exact-fit branch; `alpha > 0` swaps the
`W⁻¹·Σ⁻¹` core for the `(αI + Σ²W)⁻¹·Σ` core.  No new factorisation.

#### Three dials, not two — the `aca_eps` pitfall

The "two knobs + one metric" framing above quietly assumes a **third**
dial is pinned: the ACA+ tolerance `aca_eps`, which sets the **rank
`k_aca`** at which `A ≈ U Σ Vᵀ` is built.  That rank is an
**approximation** of `A`, *not* a regulariser — but a **loose** `aca_eps`
drops small-σ directions before `α` / `k_mode` ever see them, so the ACA+
truncation silently regularises out of your control:

| dial      | sets        | role                          | guidance                     |
|-----------|-------------|-------------------------------|------------------------------|
| `aca_eps` | `k_aca`     | **approximation** of `A`      | keep **TIGHT** (1e-10…1e-13) |
| `k_mode`  | hard filter | TSVD regularisation rank      | the auditable hard knob      |
| `α`       | soft filter | Tikhonov ridge `σ/(σ²+α)`     | the auditable soft knob      |

Rule: make `aca_eps` **tighter** than the regularisation you intend, then
steer with `k_mode` / `α` / `S`.  The verification above is the evidence —
at `aca_eps = 1e-10` the cached ψ already carries ~1.5e-3 null-space
noise, even though the optimisation *target* (seminorm + constraint)
stays near machine precision; the specific least-norm **representative**
has drifted because ACA+ pre-truncated it.

(General-`S` Tikhonov `min ‖Aψ − B‖² + α‖Lψ‖²`, `S = LᵀL`, is classically
the **GSVD of the pair `(A, L)`**; the folded core is an ACA-friendly
*implicit* GSVD — it forms only `S⁻¹V` and the `k × k` `W`, never the
dense GSVD.)

#### Choosing α / k_mode — L-curve, Morozov, and what "noise" means in SF

  - **L-curve**: plot `(‖Aψ − B‖, ‖ψ‖ or peak J)` over an `α` (or
    `k_mode`) sweep; the **corner** (maximum curvature) is the design
    point.  The fold makes the whole sweep cheap (one factorisation,
    `k × k` re-solves).
  - **Morozov discrepancy** (stop at `‖Aψ − B‖ = noise`): SF has **no
    measurement noise** — the floor is the **field unreachability** on
    this winding surface (the residual *plateau*).  So "Morozov" here
    means **stop at the achievable-on-this-surface floor**; chasing a
    lower residual past the plateau only inflates `‖ψ‖` / peak J.
  - **GCV** (parameter-free): not wired in; use the L-curve corner.

  **Honest caveat — the α-L-curve is NON-MONOTONIC.**  The iso-contour
  **topology** jumps with `α` (a saddle merges / splits a current
  region), so a smooth L-curve is not guaranteed.  Report **both**
  residual **and** peak J and pick on the corner, never residual alone.
  Measured (Gx single-stroke): plain **TSVD mode-truncation beat the best
  ridge `α`** (8.45 % vs 8.97 %) — the best regulariser is
  target-dependent, which is why the menu (`tsvd` / `tikhonov` / `h1` /
  `inductance`) is kept rather than collapsed to one.

#### Conditioning — the constant-ψ null space

Both `A` **and** a gradient seminorm `S = LᵀL` annihilate an additive
**constant** in ψ: `grad(const) = 0`, and a uniform current potential
carries no current (`K = n̂ × grad ψ = 0`).  So `S` is SPD only on the
quotient and `W = Vᵀ S⁻¹ V` can blow up.  Two fixes, both already in the
code (distinct from the small-σ damping above, which is *data*-side
conditioning):

  - a tiny **ridge** `S += (1e-9 · mean_diag) · I` (the `ridge = 1e-9` in
    [`_inductance_seminorm`](../../src/radia/panels/calc_streamfunction.py))
    lifts the harmless null space;
  - `--confine abe / on` **grounds** it via the reduction matrix `R` (one
    free constant per physical boundary component).

Without either, `W⁻¹` amplifies the constant mode — a DC offset in ψ that
the equal-increment contouring happens to ignore, but which wrecks any
`‖ψ‖`-based L-curve selection.  This is the regularisation-side reason to
run the Verify-First FES checks (slaved DOFs, constant-mode grounding)
**before** trusting an L-curve.

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

## Pushing the (homogeneity, peak-J) Pareto front

The folded Tikhonov α-sweep traces the front for a **fixed** design.  Three
**stackable levers** push the whole front toward the origin (lower peak at
the same homogeneity):

| Lever | Mechanism | Effect (planar gradient hot-spot, order-3 ψ) |
|-------|-----------|----------------------------------------------|
| **1. Tikhonov α** | moves ALONG the front (misfit ↔ seminorm) | the front itself (free α-sweep) |
| **2. L∞ seminorm** (IRLS) | redistributes current within a FIXED `A` | **−18 %** peak (median) at matched homogeneity |
| **3. geometry** (former size) | changes `A` itself — more room for the current | **−34 %** exact-homog peak (former 18 → 42 cm) |
| **4. sheet-metal** (板金) deform | FORM the surface `z=f(x,y)` (changes `A`), optimise the shape | **−17 %** exact-homog peak, whole front **−5…−18 %**, at FIXED average standoff |

**Sheet-metal (板金) deformation — honest decomposition.** Forming the
conductor surface lowers the peak two ways, and they must be separated to
make an honest claim:

  - **Standoff** — forming the sheet *closer* to the target (smaller average
    gap).  Large (≈ −53 % here) but achievable by simply repositioning a
    FLAT former, so it is **not** genuine forming.
  - **Bending** — reshaping at **fixed average standoff** (a **zero-mean**
    deformation).  This is the genuine sheet-metal contribution.

[`demo_pareto_deform.py`](../../examples/stream_function/demo_pareto_deform.py)
defaults to the zero-mean (pure-bending) constraint and optimises the shape
**per homogeneity level** (CMA-ES over a Gaussian-localised polynomial basis,
warm-started across α, `|z| ≤ 5 cm`).  The genuine bending pushes the WHOLE
front down ~ **−5 to −18 %** (best −17 % near exact homogeneity).  Note the
planar field is distance-dominated, so the pure-bending benefit is bounded
(richer DOF does not break it without raising the `|z|` cap).

**Cylinder — the lever flips.**  On a cylinder the *out-of-surface* (radial)
forming is WEAK (~ −3 %, the optimiser barely uses the `|dr|` budget); the
*in-surface* axial bending is the DOMINANT lever — opposite of the plane.  A
length-preserving axial reparametrisation `Z(z) = z + Σ_k b_k sin(kπ(z+L/2)/L)`
redistributes the loop spacing ALONG the surface at FIXED radius, so it is
**100 % genuine forming with NO standoff component** (no zero-mean trick
needed).  It pushes the whole cylinder Gx-fingerprint front down
~ **−10 to −25 %** (best −25 %), done correctly (local-spacing peak +
spacing-weighted seminorm `S`, monotone reparametrisation).  See
[`demo_pareto_cylinder_deform.py`](../../examples/stream_function/demo_pareto_cylinder_deform.py)
(`--target {gx,c2,s2,z2} [--azimuthal]`).

The in-surface lever has TWO directions (axial `Z(z)` and azimuthal
`Φ(φ)`), and **which one matters follows the target's azimuthal order `m`**:

| target | `m` | axial-only | azimuthal benefit (2D vs axial) |
|--------|-----|-----------|----------------------------------|
| `Gx = x` | 1 | −16 % | **+0.0 %** (smooth `cos φ`, no azimuthal hot spot) |
| `Z2 = 2z²−r²` | 0 | −20 % | −1.8 % (axisymmetric → axial dominates) |
| `S2 = xy` | 2 | −11 % | −1.2 % |
| `C2 = x²−y²` (ellipse) | 2 | −5 % | **−5.4 %** (2D ≈ DOUBLES it to ~ −10 %) |

So **axial-only** is the right default for `Gx` / `Z2`; turn on
`--azimuthal` for the high-`m` azimuthal shims (`C2` ellipse / `S2`), where
the azimuthal current has genuine hot spots to redistribute.

The optimal forming direction is geometry-dependent (out-of-surface on the
plane, in-surface on the cylinder; and within the cylinder, axial vs
azimuthal by target symmetry) — the same plane-vs-cylinder reversal seen for
the single-stroke distortion lever in [single_stroke.md](single_stroke.md).

Levers 1–2 reuse one ACA factorisation (the `+ α I` core + the cheap
re-fold `S⁻¹V`); lever 3 rebuilds `A` + ACA per geometry — the EXPENSIVE
outer loop wrapping the cheap inner front.  This is the natural nesting for
a **multi-objective optimiser**: NSGA-II over (geometry, α) traces the
3-objective (homogeneity, peak, **former size**) surface, autonomously
finding that low peak at a given homogeneity *requires* a larger former
(the geometry-peak coupling).  See
[`demo_pareto_geometry_nsga.py`](../../examples/stream_function/demo_pareto_geometry_nsga.py).

On the **cylinder** (Gx fingerprint, the standard MRI/shim geometry) the
geometry lever is the cylinder **length**, and it behaves differently:
instead of monotone diminishing returns it has an **optimum** (≈ 50 cm for a
15 cm-radius former over a ±8 cm DSV, **−37 %** vs a 32 cm cylinder) — once
the length covers the DSV, loops beyond ≈ ±2·DSV barely reach the target so
a longer cylinder stops helping.  The optimal lever direction is
geometry-dependent (cf. the planar vs cylinder sheet-metal-distortion lever
in [single_stroke.md](single_stroke.md)).  See
[`demo_pareto_cylinder.py`](../../examples/stream_function/demo_pareto_cylinder.py).

### Cross-reference

  - Theory: [theory.md](theory.md)
  - Single-stroke discretisation: [single_stroke.md](single_stroke.md)
  - Surface deformation: [deformation.md](deformation.md)
  - API: [`radia.stream_function.RegularizedTSVD`](../../src/radia/stream_function.py)
  - Demos: [`demo_regularized_aca.py`](../../examples/stream_function/demo_regularized_aca.py)
    (5-mode sweep), [`demo_reg_hyperparam_aca.py`](../../examples/stream_function/demo_reg_hyperparam_aca.py)
    (σ-shape optimisation, ACA+ reused),
    [`demo_pareto_tikhonov_aca.py`](../../examples/stream_function/demo_pareto_tikhonov_aca.py)
    (**(homogeneity, peak-J) Pareto front** via the folded Tikhonov α-sweep —
    the direct application of the `+ α I` core: one ACA factorisation, the
    whole front swept at ≈ 50 µs/point),
    [`demo_pareto_geometry_nsga.py`](../../examples/stream_function/demo_pareto_geometry_nsga.py)
    (geometry lever + NSGA-II joint front),
    [`demo_pareto_cylinder.py`](../../examples/stream_function/demo_pareto_cylinder.py)
    (cylinder Gx, length lever with an optimum),
    [`demo_pareto_deform.py`](../../examples/stream_function/demo_pareto_deform.py)
    (planar sheet-metal / 板金 bending lever, zero-mean genuine-forming
    decomposition),
    [`demo_pareto_cylinder_deform.py`](../../examples/stream_function/demo_pareto_cylinder_deform.py)
    (cylinder in-surface axial bending — the dominant cylinder lever, radius
    fixed = 100 % genuine, −10…−25 %)
  - MCP topic: `streamfunction(topic=regularized)`
