# API reference — `radia.stream_function`

The public API has two layers: the kernel-agnostic ACA+--QR--TSVD factor and
the improved-DUCAS/Abe current-potential inverse-design contract.  The latter
keeps the physical evidence for every eigenmode instead of reducing TSVD to a
single prefix length.

## `aca_tsvd(M, N, entry, modes=None, kmax=None, aca_eps=1e-4, method="aca_qr_tsvd")`

Truncated SVD of an `M × N` matrix `A`.  Two methods (peer review JIAM-2026-36):
`method="aca_qr_tsvd"` (default) = ACA + the standard "SVD of a low-rank product"
(QR each factor, then one small TSVD); `method="dense"` = the direct dense TSVD
(exact reference).  NO backward compatibility (legacy 2/3 and terse "qr"/"aca"
raise).

**Parameters**

| Name     | Type            | Description                                |
|----------|-----------------|--------------------------------------------|
| `M`      | int             | rows (= field / observation points)        |
| `N`      | int             | columns (= basis sources / DOFs)           |
| `entry`  | callable        | `entry(i, j) -> float` for 0-based i ∈ [0,M), j ∈ [0,N).  Called on demand by ACA+ — O(k_aca·(M+N)) calls, not full M·N. |
| `modes`  | int, optional   | TSVD modes to return (clamped to k_aca).  Default = kmax. |
| `kmax`   | int, optional   | maximum ACA+ rank.  Default = min(M, N).   |
| `aca_eps`| float, optional | ACA+ stopping tolerance (pivot threshold). Default 1e-4. |
| `method` | {"aca_qr_tsvd","dense"} | "aca_qr_tsvd" (default) = ACA + QR + TSVD (fast); "dense"/"tsvd" = direct dense TSVD (exact reference, small problems).  NO backward compat: legacy 2/3 and terse "qr"/"aca" raise. |

**Returns** `StreamTSVD` with `U (M, modes)`, `S (modes,)`,
`V (N, modes)` row-major NumPy arrays, `k_aca`, `method`.

## `pseudo_inverse_solve(result, B, k_mode=None)`

Least-norm pseudo-inverse solve via TSVD: `φ = V · diag(1/S) · Uᵀ · B`.

**Parameters**

  - `result` — `StreamTSVD` from `aca_tsvd`.
  - `B` — `(M,)` target field values.
  - `k_mode` — int, modes to use (≤ `result.modes`).  Default = all.

**Returns** `ndarray (N,)` — the basis coefficients φ.

The cached `result` can be REUSED across many calls with different
right-hand sides (e.g., Path-A iteration), since the factorisation is
independent of `B`.

## `solve(M, N, entry, B, modes=None, k_mode=None, kmax=None, aca_eps=1e-4, method="aca_qr_tsvd")`

Convenience: `aca_tsvd` then `pseudo_inverse_solve` in one call.

**Returns** `(φ, result)` — basis coefficients and the cached
factorisation (for reuse).

## Improved DUCAS / Abe current-potential solve

### `solve_abe_current_potential(response, target_field, **options)`

Solves the weighted node-current-potential equations

> `W_B B_TG = W_B A R diag(delta') q`, and
> `T = T0 + R diag(delta') q`.

This is not a blind first-`k` TSVD.  Every magnetic-field eigenmode records
its singular value, target strength `u_i.T @ (W_B B_TG)`, normalized strength,
target correlation, peak potential contribution and rejection reason.  Modes
may be non-contiguous through `allowed_modes`; accumulation stops when the
requested physical peak-to-peak and/or RMS field residual is reached.

The returned `AbeCurrentPotentialSolution` contains the full potential,
reconstructed field, checked `A @ T` residual, selected zero-based modes,
`StreamTSVD` factor and `AbeCurrentPotentialModeDiagnostics`.  Supplying the
same `precomputed_factor` makes repeated right-hand-side correction a cheap
mode sum without another ACA/TSVD factorisation.

Key options are:

| Option | Meaning |
|---|---|
| `reduction=R` | independent-to-full boundary/connection constraint map |
| `field_weights=W_B` | positive field-accuracy multipliers |
| `node_potential_scales` / `independent_potential_scales` | improved-DUCAS geometric current-potential weights |
| `initial_potential=T0` | smooth/manufacturable seed; unselected high-order content is retained |
| `allowed_modes` | arbitrary zero-based mode mask/list, including non-contiguous symmetry modes |
| `minimum_mode_strength` | threshold on `abs(u_i.T W_B B_TG)/sqrt(M)` |
| `residual_peak_to_peak`, `residual_rms` | physical field stopping criteria |
| `maximum_abs_potential` | fail-loud engineering feasibility limit |

### `abe_nearest_field_distance_scales(node_points, field_points, ...)`

Constructs the improved-DUCAS starting weight `delta_i proportional to d_i^2`,
where `d_i` is the distance from current-potential node `i` to its nearest
magnetic-field evaluation point (Abe 2013, eq. 20).  This is deliberately not
an element-area weight: the response assembly already contains area.

### `abe_reduce_node_potential_scales(R, node_potential_scales)`

Averages full-node scales into independent potentials through `R` and
normalizes the largest scale to one (Abe 2013, eqs. 24--25).  Sparse `R` stays
sparse.

### `solve_abe_bounded_current_potential(response, target_field, **options)`

Repeats the same precomputed current-potential solve after enforcing lower and
upper potential/material bounds.  With `lower_potential=0`, negative inferred
material is removed, the actual post-clip error field is recomputed, and that
new error is solved again.  The return value
`AbeBoundedCurrentPotentialSolution` includes clip, potential-change and
physical-residual histories and fails loud on stagnation or capacity limits.
At least one of `residual_peak_to_peak` or `residual_rms` is required so the
bounded iteration has an explicit physical acceptance criterion.
Apply equality constraints first and pass the reduced response: clipping
full-node potentials independently could violate `R`.

The panel exposes the same route with `--inverse-method abe`, one-based
`--abe-allowed-modes`, optional `--abe-initial-potential`, and
`--abe-node-weights nearest-distance2`.  The latter is currently limited to
`--order 1`; Radia does not reconstruct high-order NGSolve DOF coordinates in
Python.  Result JSON separates `abe_solve_converged` from
`abe_residual_target_met`; the latter is null unless a physical residual target
was explicitly supplied.

## `radia_field_kernel(obs_points, sources, component=2, field="b")`

Build a matrix-entry callback from Radia's existing field computation:
`A(i, j) = (component of field) at obs_points[i] from sources[j]`.

**Parameters**

  - `obs_points` — `(M, 3)` array of 3D observation points.
  - `sources` — length-N sequence of Radia object handles.
  - `component` — 0/1/2 for x/y/z (default 2 = z).
  - `field` — Radia field id passed to `radia.Fld` (default "b").

**Returns** `entry(i, j) -> float`.

Works for **ANY source family** — coils (Biot-Savart from `ObjFlmCur`),
permanent magnets (`ObjRecMag`), soft iron (HDiv-VIM solved materials),
shielded systems via Radia container — because
`radia.Fld` is the universal field evaluator.

## `StreamTSVD` dataclass

Returned by `aca_tsvd`.  Has attributes `U`, `S`, `V`, `k_aca`,
`method`, plus properties `M`, `N`, `modes`.

## `RegularizedTSVD` dataclass

Folds an SPD stiffness `S` into the ACA+TSVD pseudo-inverse to give a
**cached** closed-form regularised solver:

> `ψ = (S⁻¹ V) · W⁻¹ · diag(1/Σ) · Uᵀ B,    W = Vᵀ S⁻¹ V (k × k)`

| Attribute  | Type        | Description                              |
|-----------|-------------|------------------------------------------|
| `base`     | StreamTSVD  | The underlying ACA+TSVD of `A`           |
| `Sinv_V`   | ndarray     | Precomputed `S⁻¹ V`  shape `(N, k)`      |
| `W_inv`    | ndarray     | Precomputed `(Vᵀ Sinv_V)⁻¹`  shape `(k, k)` |

**Classmethod** `from_stiffness(base, S)` — precomputes `Sinv_V` (one
dense or sparse solve with `k` RHS) and `W_inv` (one `k × k` inverse).
`S` may be a dense `ndarray` or a `scipy.sparse` matrix.

**Method** `solve(B, k_mode=None)` — apply the cached factorisation:
`O(k · (M + N))` matvec.  If `k_mode < base.modes`, re-inverts the
top-left `k × k` block of `W` (because `W_inv` depends on `k`).

**Special case `S = I`**: `W = I_k` (because TSVD `V` columns are
orthonormal), so the formula reduces exactly to `pseudo_inverse_solve`.
Verified to machine precision (relative diff `≈ 1e-15`).

**Path-A pattern** (the canonical use):

```python
from radia.stream_function import aca_tsvd, RegularizedTSVD

res = aca_tsvd(M, N, entry, modes=M, aca_eps=1e-10)   # ONCE
reg = RegularizedTSVD.from_stiffness(res, S_free)     # ONCE
psi = reg.solve(B_target)

for it in range(n_iter):
    residual = B_target - I_w * Bz_chain(psi)
    delta = reg.solve(residual)                       # 0.04 ms each
    psi   = psi + step * delta
```

See [`demo_regularized_aca.py`](demo_regularized_aca.py)
for a 5-mode sweep on the planar uniform-Bz benchmark, and
[regularization.md](regularization.md) for the derivation.

## `pseudo_inverse_solve_regularized(result, B, S, k_mode=None)`

One-shot helper: builds `RegularizedTSVD.from_stiffness(result, S)`
and applies it once.  Use `RegularizedTSVD` directly when the same
`S` is reused across multiple `B` (= Path-A or regularisation sweep).

**Returns** `ndarray (N,)` — the regularised `ψ` satisfying
`A ψ = B` with minimum `ψᵀ S ψ`.

## Higher-level pipeline functions

These live in the demo files and are not part of the formal API, but
are the entry points users typically want:

### From `demo_planar_uniform_fem_psi.py`

  - `build_fem_matrix(plane_half, maxh, order, dirichlet_bc, targets)`
    Mesh the source plane + build the `M × ndof` Biot-Savart matrix
    via `LinearForm` per target.

  - `solve_tikhonov(A, B, fes, alpha=0.0, regularize="h1")`
    Three regularisation modes (L2, H1, Tikhonov) via a dense
    `(AᵀA + α S)⁻¹ AᵀB` solve.  Mathematically identical to folding the
    same `α` into the ACA+TSVD core — `ψ(α) = S⁻¹V·(αI + Σ²W)⁻¹·Σ·UᵀB` —
    which reuses one factorisation across an α-sweep; see
    [regularization.md](regularization.md) § "Tikhonov is the SAME core
    with `+ α I`".

  - `sample_psi_grid(psi_vec, fes, mesh, plane_half, n_sample)`
    Sample a GridFunction ψ on a regular grid for marching-squares.

### From `demo_planar_uniform_fem_psi_advanced.py`

  - `build_fem_matrix_deformed(..., surface_z_fn=None)`
    Same but with surface deformation `z = surface_z_fn(x, y)`.

  - `solve_h1(A, B, fes, sigma_cf_expr=None)`
    Min H¹ seminorm s.t. `A ψ = B`, with optional `1/σ` weighting.

  - `solve_inductance_diagonal(A, B, fes)`
    Diagonal proxy for self-inductance minimisation.

  - `solve_linf(A, B, fes, jmax, init_psi=None)`
    L∞ peak-current cap via scipy SLSQP (experimental).

  - `solve_l2_aca(A, B, fes, aca_eps=1e-10)`
    L² min-norm via HACApK ACA+TSVD (validates the callback contract).

  - `solve_and_evaluate(plane_half, maxh, order, target_args, B0, ...)`
    Full inner-loop wrapper for the deformation outer loop.

  - `run_deformation_search(args, trials=20)`
    Optuna CMA-ES on surface deformation parameters.

### Chain construction (from `demo_planar_uniform_coil.py`,
`demo_sf_to_peec_gx.py`)

  - `contour_polylines_xy(psi_xy, x_grid, y_grid, n_levels)` —
    marching-squares contour extraction in `(x, y)` space.
  - `contour_polylines_phi_z(psi_zphi, phi_grid, z_grid, n_levels)` —
    same in `(φ, z)` space with periodic-φ wrap (cylinder).
  - `single_stroke_spiral_xy(polylines, n_blend=1)` — Kuijpers
    Method-1 chain in `(x, y)`.
  - `single_stroke_kuijpers_chain_phi_z(polylines, a, n_blend=1)` —
    same on cylinder surface (4-lobe topology).
  - `single_stroke_lobe_chain_phi_z(polylines, a, n_arc=8)` —
    4-quadrant lobe-aware chain (LEGACY, kuijpers preferred).
  - `single_stroke_chain_phi_z(polylines, n_arc=8)` — global greedy
    NN chain (LEGACY, kuijpers preferred).

### Field evaluation

  - `bz_at(path_3d, current, obs)` — Bz at observation points from a
    3D wire path (Biot-Savart via `radia.biot_savart.h_segments_batch`).

## Cross-reference

  - Theory: [theory.md](theory.md)
  - Demos: [demos.md](demos.md)
  - MCP topic: `streamfunction(topic=api)` for the formal SA-25-020 manuscript
    API.
