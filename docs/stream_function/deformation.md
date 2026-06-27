# Surface deformation outer loop (bilevel optimisation)

> Runnable companion: [`deformation.ipynb`](deformation.ipynb) runs
> `run_deformation_search` live (penalty CMA-ES + NSGA-II Pareto front).

The inner SF design (ACA+TSVD or FE-direct + regularisation) finds the
best ψ on a *fixed* source surface.  Letting the source surface SHAPE
also vary is **bilevel optimisation**:

  - **Inner loop** (linear): given surface, find ψ via min-seminorm /
    ACA+TSVD / Lagrangian.  Fast (~1 sec).
  - **Outer loop** (nonlinear): given deformation params, evaluate the
    inner-loop result's cost.  Use CMA-ES via Optuna.

This is the SA-25-020 "(ACA+)+TSVD + CMA-ES" pattern extended with
**geometric** (not just amplitude) design freedom.

## Implementation

[`demo_planar_uniform_fem_psi_advanced.py`](../../examples/stream_function/demo_planar_uniform_fem_psi_advanced.py)
ships the framework:

```bash
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1 --deform \
    --deform-params bump --deform-trials 20
```

Available `--deform-params` (composable comma-separated):

| Param   | DOF | Description                                  |
|---------|-----|----------------------------------------------|
| `zoff`  | 1   | Source-plane z translation (toward/away from target) |
| `bump`  | 3   | Gaussian bump: amplitude + (x, y) centre     |

The CMA-ES outer loop samples deformation params via Optuna's
`CmaEsSampler`, re-meshes the source surface, runs the inner solve,
returns the target-plane RMS as the cost.

## Measured results (planar uniform Bz benchmark)

| Variant | Trials | Wall | Best RMS    | Best params                  |
|---------|--------|------|-------------|------------------------------|
| `zoff`  | 10     | 7 s  | 1.50 %      | zoff = +20 mm                |
| `bump`  | 20     | 23 s | **0.77 %**  | amp = −22 mm, off-centre     |
| `zoff+bump` | 50 | 33 s | 0.85 %     | competing local minimum      |

Best bump found by CMA-ES: a **concave dip** away from the target,
off-centre.  Physically: redirects current toward the target by creating
a "ramp" that focuses the field.  Best `zoff` alone: source plane lifted
20 mm toward the target (less spreading).

## Composability with regularisation choice

The inner solver can be any of the regularisation modes:

```bash
# 1/sigma weighted inner + bump deform outer
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1_sigma \
    --sigma-cf "exp(-(x*x + y*y) / 0.02)" \
    --deform --deform-params bump --deform-trials 30

# l2_aca (HACApK ACA+) inner + zoff+bump outer
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize l2_aca \
    --deform --deform-params "zoff,bump" --deform-trials 50
```

## When to use deformation freedom

**USE deformation** if the inner-solver baseline RMS is *several percent*
or more.  CMA-ES has measurable improvement headroom (e.g., 2 % → 0.8 %).

**DON'T use deformation** if the baseline is already sub-0.5 %.  CMA-ES
will waste its budget searching a near-optimal landscape and may get
stuck (we measured order=3 H¹ baseline 0.51 % → bump 0.82 %, WORSE).

## Extension: parameterise more of the geometry

The current implementation deforms only a flat plane.  Easy extensions
(future work):

  - Multi-bump surfaces (N Gaussians) — increase deformation DOFs
  - Cylindrical surfaces with axial perturbations (= Gx fingerprint
    coils with optimised cylinder geometry)
  - Conformal surfaces from OCC (NGSolve `Mesh` from CAD)
  - Bi-planar primary + shield with relative spacing as DOF

Each is ~1 day of mesh-generation refactoring on top of the existing
CMA-ES framework.

## Extension: optimise σ as a design variable

Currently `--sigma-cf` is user-given.  Letting CMA-ES choose `σ(x, y)`
parameters (= location and strength of forbidden / allowed regions) is
a natural bilevel:

  - Inner: solve `min ∫(1/σ)|∇ψ|²` s.t. `A ψ = B`
  - Outer: choose σ params to minimise target RMS + manufacturing cost

This is a 1-week implementation on top of what we have.

## Constrained reg-aware loop — `min ψᵀSψ s.t. RMS ≤ ε`

The default outer loop minimises the **field error**.  The reg-aware
mode flips this: minimise the **regularisation norm** `ψᵀSψ` (= an
energy / dissipation / inductance proxy depending on which `S`)
subject to a hard tolerance on field accuracy.  This is the standard
MRI-shim / gradient-coil formulation — *"the lowest-energy coil that
meets the field-uniformity spec"*.

```bash
# Minimise H1 seminorm (∫|∇ψ|² ≈ surface current L² norm)
# subject to target-plane RMS <= 2.5 %
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1 --deform --deform-params bump --deform-trials 20 \
    --minimize-reg --eps-rms 0.025

# Physically meaningful: minimise ohmic dissipation
# (1/σ weighted H1) subject to field accuracy
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1_sigma --sigma-cf "1.0 + exp(-(x*x+y*y)/0.02)" \
    --deform --deform-params "zoff,bump" --deform-trials 30 \
    --minimize-reg --eps-rms 0.02
```

### Penalty form

The constraint enters as a quadratic penalty scaled by the **first
trial's reg_norm** (so the penalty and the objective share units):

    cost = reg_norm + reg_scale · penalty_weight · max(0, RMS/ε − 1)²

For trials inside the feasible set (`RMS ≤ ε`) the cost equals
`reg_norm`; outside it grows quadratically.  `penalty_weight=100`
(default) is enough to push CMA-ES toward the feasible boundary
without numerical issues; tighten via `--reg-penalty 1000` if trials
keep landing outside the cap.

The study tracks the **best feasible** trial (lowest `reg_norm` among
those satisfying the constraint) separately as
`study.user_attrs["best_feasible"]` — that is the engineering answer,
distinct from `study.best_value` (which is the penalised cost CMA-ES
optimised).

### Measured Pareto trade-off (planar uniform Bz, 8 trials, ε=2.5 %)

| Mode                          | RMS    | `ψᵀSψ` (H¹ seminorm) | Note                          |
|-------------------------------|--------|------------------------|-------------------------------|
| Flat baseline (no deform)     | 2.09 % | 2.59e+06               | reference                     |
| `min(RMS)` (CMA-ES)           | **0.77 %** | 2.74e+06           | accuracy-optimal — energy **up** |
| `min(reg) s.t. RMS ≤ 2.5 %`   | 1.82 % | **2.15e+06**           | meets spec, **28 % lower energy** vs `min(RMS)` |

The accuracy-optimal point spends MORE energy than the flat baseline
(2.74e+06 > 2.59e+06) because squeezing RMS from 2 % to 0.8 % requires
the optimiser to add localised current concentrations.  The reg-min
constrained point lives lower-left on the Pareto front: accepts a
~2 % RMS (still inside spec) in exchange for ~25-30 % less surface
current L².

### When to use which

| Goal                                            | Recommended mode |
|-------------------------------------------------|------------------|
| Best possible field uniformity, energy secondary | (default) `--deform` |
| Manufacturing efficiency target (current limit / heat / inductance) with field spec | `--minimize-reg --eps-rms ε` |
| Trace the full Pareto front directly             | `--pareto` (NSGA-II multi-objective, see below) |

### Cache reuse

The cached `RegularizedTSVD` (see [regularization.md](regularization.md))
is rebuilt per trial because **the matrix `A` changes when the
surface deforms** — each Optuna trial calls
`build_fem_matrix_deformed` to re-assemble for the new surface.
What IS amortised: the `_compute_reg_norm` helper builds `S` once
per trial (~50 ms), then a single matvec gives `ψᵀSψ`.

For a fixed surface + reg-aware loop over the regularisation
hyperparameters (e.g. `--sigma-cf` parameters), `A` would not
change and the ACA+TSVD base + `Sinv_V` cache survive across trials
— a clean extension worth ~1 day.

## Multi-objective Pareto front — `--pareto` (NSGA-II)

The penalty form above gives ONE point per Optuna run.  To trace the
entire Pareto curve of `(RMS, ψᵀSψ)` in a single shot, switch to
multi-objective NSGA-II:

```bash
python demo_planar_uniform_fem_psi_advanced.py \
    --regularize h1 --deform --deform-params bump \
    --deform-trials 50 --pareto
```

Implementation:

  - `optuna.samplers.NSGAIISampler(population_size = max(8, trials//4, 24))`
  - `optuna.create_study(directions=["minimize", "minimize"])`
  - `objective(trial) -> (rms, reg_norm)` (tuple, not scalar)
  - `study.best_trials` returns the non-dominated set; the demo
    deduplicates (same params evaluated by NSGA-II across
    generations is common) and sorts by RMS.

Outputs:

  - `demo_pareto_results.json` — every trial plus the non-dominated
    set with params, `(rms, reg_norm)` values, and the surface
    deformation chosen.
  - `demo_pareto_plot.png` — scatter of all 50 trials in
    `(RMS [%], ψᵀSψ)` space with the Pareto front overlaid.

### Measured Pareto front (planar uniform Bz, 50 trials, H¹, bump)

| RMS    | `ψᵀSψ` (H¹) | bump (amp, cx, cy)              | Note         |
|--------|--------------|----------------------------------|--------------|
| 0.58 % | 2.80e+06     | (−0.032, −0.059, +0.007)         | accuracy end |
| 0.61 % | 2.46e+06     | (+0.021, −0.144, +0.141)         | **knee**     |
| 1.04 % | 2.29e+06     | (+0.045, +0.140, −0.121)         |              |
| 1.39 % | 2.27e+06     | (+0.045, −0.059, +0.146)         |              |
| 1.67 % | 2.27e+06     | (+0.045, −0.062, +0.146)         |              |
| 2.13 % | 2.20e+06     | (+0.044, +0.099, +0.081)         | energy end   |

Wall time: **33 s** (50 trials of the inner SF solve).  The curve has
a sharp knee at ~0.6 % RMS — dropping reg_norm by ~12 % (from 2.80e+06
to 2.46e+06) costs only 0.03 % RMS accuracy.  Past the knee the curve
flattens: dropping reg_norm by another ~10 % costs ~1.5 % RMS.  This
is the classical "knee in the L-curve" that gives the engineering
operating point.

For context, plain `min(RMS)` on the same problem with CMA-ES (12 s,
8 trials) lands at 0.77 % RMS, `reg_norm = 2.74e+06`.  NSGA-II's
accuracy end (0.58 %, 2.80e+06) is BETTER on field uniformity but
slightly worse on energy — exactly the expected behaviour: scalar
minimisation collapses to one specific point on the front, while
NSGA-II actively explores the full curve.

### When to pick `--pareto` vs `--minimize-reg`

| Situation                                              | Use            |
|--------------------------------------------------------|----------------|
| Spec is known up-front ("RMS must be ≤ 2 %")           | `--minimize-reg --eps-rms 0.02` |
| Trying to **choose** the spec (= operating point unknown) | `--pareto` (lets you see the curve before committing) |
| Need a single number for a report / DoE                | `--minimize-reg` |
| Need a Pareto-front figure for a paper                  | `--pareto` |
| Budget < 15 trials                                     | `--minimize-reg` (CMA-ES converges faster on a scalar) |
| Budget ≥ 30 trials                                     | either; `--pareto` gives more information per trial |

## A different objective: the (homogeneity, peak current density) front

The `--pareto` front above is **(RMS, ψᵀSψ)** — field accuracy vs the L²
surface-current energy.  A SEPARATE Pareto study optimises **(field
homogeneity, PEAK current density `max|∇ψ|`)** via the folded-Tikhonov
α-sweep, and uses surface deformation as one of its levers — there it is the
**sheet-metal (板金) surface-forming** lever (genuine bending at FIXED average
standoff, with the honest standoff-vs-bending decomposition): planar
out-of-surface −17 %, cylinder in-surface −25 %.  Same surface-reshape +
re-solve-ψ mechanism as this page, different objective (peak, not energy).
See [regularization.md § Pushing the front](regularization.md#pushing-the-homogeneity-peak-j-pareto-front)
and `demo_pareto_deform.py` / `demo_pareto_cylinder_deform.py`.

## Cross-reference

  - Inner solvers: [regularization.md](regularization.md)
  - Pushing the (homogeneity, peak-J) front + sheet-metal forming:
    [regularization.md § Pushing the front](regularization.md#pushing-the-homogeneity-peak-j-pareto-front)
  - Math: [theory.md](theory.md)
  - Optuna BBO recipes: MCP `optuna(topic=lab_applications)`
  - SA-25-020 lineage: MCP `streamfunction(topic=cmaes)`
  - Regularisation closed form + Pareto/板金: MCP `streamfunction(topic=regularized)`
