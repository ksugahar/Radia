# ESIM Implementation: Architecture, Karl Loop, Convergence Analysis

**Audience.** Reviewers / readers who need to evaluate the algorithmic
correctness of the Radia ESIM implementation, the three-solver
dispatch architecture, and the convergence properties of the outer
Karl iteration.

**Companion documents.**
- [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) — formulation and weak forms.
- [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) — numerical benchmarks.
- [`SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md) — why this implementation architecture (scalar BIE on path A) is the right one to publish.
- [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) — focused note on the PEEC vs BEM-A coil R discrepancy.

---

## 1. Three-Solver Dispatch Architecture

Radia ships three independent coupled-solver paths that all consume
the same 1-D cell-problem solver but differ in how the **outer**
electromagnetic field is represented:

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 4 (CLI):   calc_inductance.py / calc_fem_kelvin.py /     │
│                   calc_fem_coilmesh.py                          │
│  Layer 3 (Karl):  outer-loop wrapper, Z_s update                │
│  Layer 2 (Field): BEM-SIBC / HCurl-Kelvin / HCurl-AV            │
│  Layer 1 (Cell):  ESIMFiniteSlabSolver(geometry='cylinder')     │
└─────────────────────────────────────────────────────────────────┘
```

| Path | Coil model | Workpiece model | Outer DOFs | Wall-time |
|---|---|---|---|---|
| `calc_inductance.py` | PEEC filament bundle (Biot-Savart) **or** BEM-A surface RWG | Scalar BIE-SIBC on workpiece surface | ~10²–10³ | ~5 s |
| `calc_fem_kelvin.py` | PEEC filament bundle (line-integral RHS) | HCurl A on volumetric workpiece + Kelvin transformation | ~10⁴ | ~70 s |
| `calc_fem_coilmesh.py` | Volumetric A-V compound FES on coil mesh | HCurl A on workpiece + Robin BC | ~10⁵ | ~100 s |

All three call the **same** `ESIMFiniteSlabSolver` instance for the
cell problem, ensuring consistent saturation behaviour.  The Karl loop
implementations differ in line count but are logically identical (see
§ 3 below).

### 1.1 Why three paths?

Each path is correct for a different physics regime + accuracy goal:

| Path | Sweet spot |
|---|---|
| `calc_inductance.py` (PEEC + scalar BEM-SIBC) | Fast P_wp + ΔL screening at 1-way coupling.  Thin-skin workpiece. |
| `calc_fem_kelvin.py` | Workpiece geometry that doesn't fit the scalar-potential assumption (thick skin, large μ_r contrast, nontrivial topology). |
| `calc_fem_coilmesh.py` | Full 3-D Maxwell reference; coil ohmic losses; geometry-resolved workpiece-coil proximity. |

For the IGTE digest, `calc_inductance.py` is the primary tool
(production speed for the publication's parameter sweep);
`calc_fem_coilmesh.py` serves as the volumetric A-V reference for
each headline number.

---

## 2. The 1-D Cell Problem Solver

Two classes ship in [`src/radia/esim_cell_problem.py`](../../src/radia/esim_cell_problem.py):

| Class | File line | Role |
|---|---|---|
| `ESIMCellProblemSolver` | 816 | Legacy infinite-domain solver |
| `ESIMFiniteSlabSolver` | 339 | Finite-thickness / curvature-aware solver (used by all four production scripts) |

### 2.1 ESIMFiniteSlabSolver — discretisation

**PDE family** (lines 352, 573):

- Slab: $\rho\,\partial_z^2 H + j\omega\mu(|H|)\,H = 0$ on $z \in [0, a]$,
  BCs $H(0) = H_t$ (driven), $H'(a) = 0$ (insulated)
- Cylinder: $(\rho/r)\,\partial_r[r\,\partial_r H] + j\omega\mu(|H|)\,H = 0$
  on $r \in [0, R]$, BCs $H'(0) = 0$ (regularity), $H(R) = H_t$
- Finite-slab (2-sided): as slab + Robin / mirror BC

**Mesh** (line 428):

```python
self.mesh_points = np.linspace(0, half_thickness, n_nodes)
```

Default `n_nodes = 100` (line 366).  Uniform; geometric grading is a
roadmap item (see [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md)
§ 2.1).

**Finite-difference stencil — interior node (cylinder)**:

```python
r_ip = r[i] + h/2                  # r_{i+1/2}
r_im = r[i] - h/2                  # r_{i-1/2}
coef_p = ρ * r_ip / (r[i] * h²)
coef_m = ρ * r_im / (r[i] * h²)
diag_main[i] = -(coef_p + coef_m) + jω * μ_dist[i]
diag_upper[i] = coef_p
diag_lower[i-1] = coef_m
```

(Lines 594-602.) This is the **standard finite-volume conservative
scheme** on a uniform radial mesh, with cell-centred fluxes at
$r_{i\pm 1/2}$ — provably second-order accurate in the smooth interior.

**Finite-difference stencil — axis at r = 0**:

```python
coef_0 = 4 * ρ / h²                # L'Hôpital limit of (1/r) ∂_r[r ∂_r]
diag_main[0] = -coef_0 + jω * μ_dist[0]
diag_upper[0] = coef_0
```

(Lines 609-611.)  Symmetry `H'(0) = 0` is imposed via the ghost-cell
identity `H[-1] = H[1]`, which combined with `(1/r)∂_r[r∂_r H] → 2∂_r^2 H`
as `r → 0` gives the `4ρ/h²` coefficient.  This is essential to
recover the analytical Bessel solution (§ 3.1 of
[`CROSS_VALIDATION.md`](CROSS_VALIDATION.md)).

**Linear solve.**  The tridiagonal system is assembled as a `scipy.sparse`
CSR matrix and solved with `scipy.sparse.linalg.spsolve`.  Memory and
work scale linearly in `n_nodes`; one solve is `~50 μs` at
`n_nodes = 100`.

### 2.2 Inner Picard iteration (cell-internal nonlinearity)

The cell-problem nonlinearity comes from `μ(|H|)`.  Picard linearises
by fixing `μ` from the previous iterate:

```python
for iteration in range(max_iter):
    H_new = self._solve_linear_system(H0, mu_dist)
    mu_new = np.array([self._get_mu(abs(h)) for h in H_new])
    mu_dist = (1 - relaxation) * mu_dist + relaxation * mu_new
    rel_change = max|mu_new - mu_dist| / max|mu_dist|
    if rel_change < tol:
        break
```

(Lines 508-541 of [`esim_cell_problem.py`](../../src/radia/esim_cell_problem.py).)

**Defaults**: `tol = 1e-6`, `max_iter = 50`, `relaxation = 0.5`.

For monotone BH curves (all production materials), Picard converges
linearly with rate ≤ `relaxation × L_local` where `L_local` is the
local Lipschitz of `μ(|H|)`.  Convergence is fast (~5-10 iter) when
`|H|` stays on one side of the BH knee; convergence slows to ~0.5/iter
when `|H|` straddles the knee.

### 2.3 Analytical seed for the cell

To minimise inner Picard iterations, the cell is seeded with the
constant-μ analytical solution:

```python
# Cylinder (line 489-499):
γ = (1 + 1j) / δ
H_seed = H0 * scipy.special.iv(0, γ * z) / scipy.special.iv(0, γ * a)

# Slab (line 504-505):
γ = (1 + 1j) / δ
H_seed = H0 * cosh(γ * (a - z)) / cosh(γ * a)
```

For high-`ξ` regimes (`scipy.special.iv` overflows), the seed falls
back to a thin-skin exponential `H ~ H0 * exp(-(a-z)/δ)` (line 506-511).

### 2.4 Output dict

`ESIMFiniteSlabSolver.solve(H0)` returns (lines 547-559):

| Key | Type | Meaning |
|---|---|---|
| `Z` | complex | Effective surface impedance `Z = 2(P' + jQ')/|H_0|²` [Ω] |
| `P_prime` | float | Active power density `½ Re(Z) |H_0|²` [W/m²] |
| `Q_prime` | float | Reactive power density `½ Im(Z) |H_0|²` [VAR/m²] |
| `P_magnetic` | float | Hysteretic / grain-eddy magnetic loss density [W/m²] (only nonzero when `complex_mu` is supplied) |
| `R_ratio` | float | AC-to-DC resistance ratio `R_ac / R_dc` (diagnostic output from ESIM; the panel itself uses `cylinder_ac_impedance` from `analytical_formulas.conductor_impedance` for the Bessel per-filament `Zs_fil` injection — see [R_MISMATCH_PEEC_VS_BEMA.md](R_MISMATCH_PEEC_VS_BEMA.md)) |
| `H_solution` | complex array | H(r) profile, length `n_nodes` |
| `mu_final` | float | Surface-region mean μ from the converged iteration |
| `converged` | bool | True if `rel_change < tol` reached |
| `iterations` | int | Number of Picard iter to convergence |
| `xi`, `delta` | float | ξ = R/δ and δ for diagnostic logging |

The `R_ratio` is computed by direct numerical integration of the eddy
current dissipation:

```python
# Slab (line 706-716):
R_ratio = a * ∫_0^a |∂_z H|² dz / |H(0) - H(a)|²

# Cylinder (line 680-704):
R_ratio = ∫_0^R |∂_r H|² · r dr / [R · |H(0)|² / 2]
```

The closed-form Dowell reference for circular cross-section is at
[`esim_cell_problem.py:787-813`](../../src/radia/esim_cell_problem.py#L787-L813)
and matches the numerical `R_ratio` to ~10⁻³ in `1 ≤ ξ ≤ 10`.

---

## 3. The Outer Karl Iteration

The Karl iteration wraps any of the three outer solvers with the
fixed-point operator:

$$
Z_s \xrightarrow{\text{outer solve}} H_t \xrightarrow{\text{cell solve}} Z_s'.
$$

**Pseudo-code** (common to all three paths; line numbers vary):

```python
# Seed
Z_s = esim_solver.solve(H0=5.0, max_iter=5)['Z']    # conservative low-H seed

for k in range(max_iter):                            # outer Karl loop
    res = outer_solve(Z_s)                           # BEM / FEM solve
    H_t_rms = res['H_t_rms']                         # mesh-RMS tangential H
    Z_s_old = Z_s
    Z_s_new = esim_solver.solve(H_t_rms)['Z']
    Z_s = relax * Z_s_new + (1 - relax) * Z_s_old    # damped Picard
    dZ = |Z_s - Z_s_old| / |Z_s_old|
    if dZ < tol and k > 0:                           # avoid spurious seed convergence
        break

# Final re-solve at converged Z_s so post-proc reflects it
res_final = outer_solve(Z_s)
```

This pattern is implemented at:
- [`calc_inductance.py:588-712`](../../src/radia/panels/calc_inductance.py#L588-L712)
- [`calc_fem_kelvin.py:550-941`](../../src/radia/panels/calc_fem_kelvin.py#L550-L941)
- [`calc_fem_coilmesh.py:275-422`](../../src/radia/panels/calc_fem_coilmesh.py#L275-L422)

**Defaults** (in panel / CLI args): `max_iter = 15`, `tol = 1e-3`,
`relax = 0.5`, seed at `H_0 = 5 A/m`.

### 3.1 Why "Karl"?

"Karl iteration" is lab shorthand referring to **Karl Hollaus**, the
first author of the canonical ESIM scalar-potential formulation paper:

> K. Hollaus, M. Kaltenbacher, J. Schöberl,
> *A Nonlinear Effective Surface Impedance in a Magnetic Scalar
> Potential Formulation*, IEEE Trans. Magn., 2025.
> DOI: [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932).

External readers / reviewers should see this as "Hollaus-type Picard
relaxation"; the lab nickname is preserved in the code for historical
continuity.

### 3.2 The H_t extraction step

Each outer path extracts a single scalar `H_t_rms` to feed back into
the cell.  The formula varies:

| Path | H_t extraction |
|---|---|
| `calc_inductance.py` (BIE) | `H_t_rms² = (φᵀ K φ) / A_wp` where φ is the scalar potential, K is the stiffness matrix |
| `calc_fem_kelvin.py` (HCurl) | `H_t_rms² = |jω/Z_s|² · ∫ |A_t|² dS / A_wp`, A_t = A − (A·n)n on workpiece BND |
| `calc_fem_coilmesh.py` (HCurl A-V) | Same as `calc_fem_kelvin.py` |

The BIE formula derives from the variational identity for the scalar
potential; the HCurl formulas use the Robin BC inversion
`H_t = (jω/Z_s) A_t`.

### 3.3 Per-DOF Karl ("per-panel ESIM")

For workpieces with strong spatial variation of `|H_t|` (saturation
pattern), a single mesh-RMS `H_t_rms` is too coarse.  The per-DOF
variant ships in all three paths (BIE: v4.47.2+; FEM: v4.55+).

**Per-DOF `|H_t|` extraction (v4.67.0+)**: triangle-wise P1 gradient
`∇_s φ = Σ_j φ_j ∇N_j` per triangle, area-weighted average to
vertices.  See `radia.bem_sibc_solver.extract_H_t_per_dof_grad`:

```python
# calc_inductance.py:849-851 (v4.67.0+):
from radia.bem_sibc_solver import extract_H_t_per_dof_grad
phi_vec = np.asarray(res_bem["phi_vec"])
H_t_per = extract_H_t_per_dof_grad(phi_vec, wp_mesh)

# Per-DOF ESIM call
Z_s_new = np.array([esim_solver.solve(H_t_per[i])['Z']
                    for i in range(bem.ndof)])
Z_s = anderson.step(Z_s_old, Z_s_new)            # ndarray update
```

> **DO NOT** use the Galerkin localization
> `|H_t|_i² ∝ φ_i (Kφ)_i` (which we shipped in v4.47.2 - v4.66.x).
> It samples the **surface Laplacian**, not the gradient norm.  On
> the steel cylinder benchmark (I=100 A, f=50 kHz) it mis-places
> the saturation hot-spot and gave `P_per = 45.4 W` vs the correct
> `P_per = 18.75 W` — i.e. it **flipped the sign of the per-vs-scalar
> disagreement** (we initially reported +48 % "scalar
> under-estimates"; the correct sign is −38.5 %, scalar
> over-estimates).  The inline comment at
> `calc_inductance.py:842-848` warns about this in the source.

`Z_s` becomes an `ndarray[ndof]`; the BIE solver accepts this via row-
scaling: `A[i, :] = 0.5 M - DL + γ[i] · SL · M⁻¹ · K[i, :]` where
`γ[i] = jω/Z_s[i]`.

**Cost.**  Per-DOF Karl multiplies the cell-solve cost by `N_DOF`.  At
166 DOFs (production gapped-torus sample) the overhead is ~7 s per
Karl iteration; at 5k DOFs it dominates and Anderson acceleration
becomes important.

### 3.4 Convergence criterion

Both scalar and per-DOF Karl use:

$$
dZ = \frac{\max_i |Z_s^{(k+1)}[i] - Z_s^{(k)}[i]|}{\max_i |Z_s^{(k)}[i]|}
$$

with iteration `k > 0` required to avoid spurious termination at the
seed (which has `dZ = 1.0` by definition of the relaxation update,
even if `Z_s` happened to match after the first cell call).

**Per-DOF `dZ` is a harsher criterion than scalar `dZ`.**  For
per-element Karl (`--esim-per-panel`) the `max_i` in `dZ_max` is
taken over the **worst** DOF (typically a corner / hot-spot where
the cell-solver Lipschitz constant is closest to 1).  Per-DOF runs
routinely hit `max_iter` with `dZ_max` still in the 0.1-0.5 range
even though scalar-Karl on the same problem converges in <10 iter.
A `max_iter` cap on `dZ_max` is therefore **not by itself** evidence
of divergence — but it is **not by itself** evidence of a usable
result, either.  You must inspect the trajectory.

**Decision rule.**  Plot `esim_history` with
[`examples/ih_esim_benchmark/plot_karl_history.py`](../../examples/ih_esim_benchmark/plot_karl_history.py),
which overlays `dZ` (or `dZ_max`), `Z_s_abs` (with min/max band for
per-panel), and `|H_t|` per iteration.  Read it as:

| Trajectory signature | What it means | Action |
|---|---|---|
| `dZ` decays monotonically, crosses `esim_tol` | Karl converged | use the result |
| `dZ` non-monotone, but mean `Z_s_abs` and `H_t_rms` settle over last 3-5 iter (variation < ~5 %) | per-DOF noise on a useful plateau | accept the `max_iter` cap; report P_wp with the caveat |
| `dZ` non-monotone AND `Z_s_abs`, `H_t_rms` still drifting at iter N | Karl genuinely under-relaxed or BH knee straddled | lower `--esim-relax` to 0.3, raise `--esim-max-iter`; do NOT publish the iter-N number |
| `dZ` oscillates / grows | true divergence (`α L > 1` somewhere) | drop `--esim-relax` to 0.2, check BH curve monotonicity |

**P_wp robustness vs damping / iteration count (pre-2026-05-24, with
Galerkin localization).**  Two v4.47.2-era runs at the headline
benchmark (steel cylinder, 50 kHz, I_port = 100 A) gave:

| Run | `--esim-relax` | `--esim-max-iter` | Iters used | `P_wp` [W] | Last-5-iter drift, `<\|Z_s\|>` | Last-5-iter drift, `<\|H_t\|>` |
|---|---|---|---|---|---|---|
| v4 | 0.5 | 15 | 15 (capped) | 45.143 | 4.37 % | 8.44 % |
| v5 | 0.3 | 30 | 30 (capped) | 45.196 | 0.54 % | 0.88 % |

> **Note (2026-05-24)**: the `P_wp ≈ 45 W` values in this table are
> from the Galerkin localization `|H_t|² ∝ φ_i(Kφ)_i` (which we
> shipped until v4.66.x).  After switching the per-DOF extractor to
> the triangle-wise P1 gradient in v4.67.0 (§ 3.3 above), the same
> benchmark gives **P_per = 18.75 W** (sweep_v2; see
> [`CROSS_VALIDATION.md` § 6b](CROSS_VALIDATION.md#6b-per-element-vs-scalar-z_s-the-headline-contribution)).
> The convergence behaviour also changes: with
> `--esim-anderson-m 5` (production default, v4.68+) per-DOF Karl
> reaches `dZ_max < 1e-3` in **7-30 iter** across the IGTE 32-case
> sweep, no longer hitting the `~5e-2` per-DOF noise floor that was
> typical with plain damped Picard.

The two-run `P_wp` agreement to 0.12 % across damping settings is
still informative — it shows the **integrated quantity is robust
even when per-DOF `dZ` has not converged** — but the absolute value
has shifted with the gradient-extraction fix.

The per-DOF `dZ_max` remained in the 0.06-0.20 range across both
runs under damped Picard alone — this is the per-element ESIM noise
floor on the hot-spot DOFs and **cannot** be driven below `~5e-2` by
damping alone.  Anderson Type-II acceleration (m=5, v4.68+) is what
breaks past it; see § 6b.2 in `CROSS_VALIDATION.md` for the
production convergence numbers.

---

## 4. The Robin BC on the Workpiece Surface

The outer BEM/FEM solves all impose the Leontovich SIBC

$$
\mathbf{n} \times \mathbf{E} = Z_s\,\mathbf{H}_t
$$

on the workpiece surface, parametrised by `Z_s` from the cell solver.
The implementation differs per path:

### 4.1 BIE (calc_inductance.py)

The scalar BIE on the workpiece is `(SL + γ M) φ = φ_inc` where
`γ = jω/Z_s`, SL is the single-layer operator, M the mass on the H1
basis.

In `bem_sibc_solver.ScalarBIESIBCSolver.solve(phi_inc, Z_s, omega)`,
`Z_s` is folded into the system matrix `A_sys = 0.5 M - DL + γ · SL · M⁻¹ · K`
(see [`bem_sibc_solver.py:401-412`](../../src/radia/bem_sibc_solver.py#L401-L412)).
Per-DOF `Z_s` is supported by element-wise multiplication on the diagonal.

### 4.2 HCurl + Robin (calc_fem_kelvin.py)

The curl-curl bilinear form is

$$
a(A, A') = \int_\Omega \nu(r')\,\mathrm{curl}\,A \cdot \mathrm{curl}\,A'\,d\Omega
         + \int_{\Gamma_{\mathrm{sibc}}} \frac{j\omega}{Z_s}\,A_t \cdot A'_t\,dS
$$

(see [`calc_fem_kelvin.py:565`](../../src/radia/panels/calc_fem_kelvin.py#L565)).

The Robin coefficient `jω/Z_s` is built as an NGSolve `CoefficientFunction`;
for per-DOF Karl it becomes a `GridFunction` on the SIBC boundary FES.

Kelvin transformation (lines 240-252) replaces the outer Dirichlet
truncation: the Kelvin sphere `r' = R_kelvin` is the boundary of an
inverted-coordinate domain where `ν → NU_0 · (r'/R)²`, making the
exterior contribution to the integral finite without a far-field
truncation.

### 4.3 HCurl A-V compound (calc_fem_coilmesh.py)

The compound FES `HCurl(A) × H1(φ_coil)` (lines 196-202) supports the
A-V formulation:

$$
\mathbf{J}_{\mathrm{coil}} = -j\omega\,\sigma_{\mathrm{coil}}\,(A + \nabla\varphi)
$$

with Dirichlet `φ = 1` at the source face and `φ = 0` at the sink.
The workpiece Robin BC is identical to the calc_fem_kelvin form,
folded into the same compound bilinear form.

The coil dissipation `P_coil = 0.5 σ⁻¹ ∫|J|²` is computed
post-hoc and adds to `P_wp` for the total dissipation.

---

## 5. Per-DOF Backend Support Matrix

| Path | Scalar Z_s | Per-DOF Z_s | Notes |
|---|---|---|---|
| `calc_inductance.py` `--wp-bem-backend intree-dense` | YES | YES | Production-ready (v4.47+) |
| `calc_inductance.py` `--wp-bem-backend hacapk` | YES | NO | Per-DOF requires dense matrix access |
| `calc_fem_kelvin.py` | YES | YES (v4.55+) | Robin term becomes a CF on the SIBC boundary FES |
| `calc_fem_coilmesh.py` | YES | YES (v4.55+) | Same Robin-CF pattern |

The HACApK backend gap is the only remaining hole; closing it would
require either (a) wrapping `Z_s` into a diagonal preconditioner inside
the HACApK GMRES, or (b) extending the ACA assembler to support
position-dependent kernels.  Tracked as roadmap.

---

## 6. Lipschitz Estimation and Damping Choice

The damped Picard iteration

$$
Z_s^{(k+1)} = (1 - \alpha)\,Z_s^{(k)} + \alpha\,T(Z_s^{(k)}),
\qquad T = \mathcal{E} \circ H_t,
$$

is contractive iff `α·L < 1` where `L` is the local Lipschitz of `T`.
We do not have a closed-form bound on `L`; we estimate it empirically
from the iteration history.

**Empirical estimate**: log-fit the `dZ/Z` sequence:

$$
\log\,dZ^{(k)} = \log\,dZ^{(0)} + k \cdot \log(\alpha L) \quad\Rightarrow\quad
\alpha L \approx \frac{dZ^{(k)}}{dZ^{(k-1)}}.
$$

For the production gapped-torus + steel-cylinder benchmark at 50 kHz
(Karl iteration history in
[`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) § 3):

| k | dZ/Z | Ratio |
|---|---|---|
| 1 | 1.7e-2 | (seed) |
| 2 | 8e-3   | 0.47 |
| 3 | 3e-3   | 0.375 |
| 4 | 1e-3   | 0.33 |
| 5 | 3e-4   | 0.30 |

The contraction factor ratchets down from ~0.5 toward ~0.3 as the
iteration approaches the fixed point — i.e. `αL` decreases.  With
`α = 0.5` this gives empirical `L ∈ [0.6, 1.0]`, comfortably in the
contraction regime.

**Implication for `--esim-relax`:**

- For `L < 1` (most engineering regimes), any `α ≤ 1` converges; `α = 1`
  is theoretically faster but oscillates in iter 1-2 due to the
  conservative seed undershooting the operating point.
- For `L ≈ 1.5-2` (deep saturation, sharp BH knee), `α = 1` diverges;
  use `α = 0.3` for safety.
- For `L ≫ 2` (degenerate BH curves; data error), no damping rescues
  Picard — replace with Newton or Anderson.

---

## 7. Anderson Acceleration (Planned)

Anderson acceleration (AA) replaces the scalar damping with a
quasi-Newton update that maintains a history window:

$$
Z_s^{(k+1)} = Z_s^{(k)} - \sum_{i=0}^{m} \gamma_i^*
              \left(Z_s^{(k-i)} - Z_s^{(k-i-1)}\right)
            + \beta\,F^{(k)},
$$

where `F^{(k)} = T(Z_s^{(k)}) - Z_s^{(k)}` is the residual and
`γ_i^*` solve

$$
\min_\gamma \left|F^{(k)} + \sum_i \gamma_i (F^{(k-i)} - F^{(k)})\right|^2.
$$

For Karl on a single complex scalar, `m = 2-3` is typical and gives a
2-4× iter-count reduction in the saturated regime.  For per-DOF Karl,
AA is applied independently per DOF (no cross-coupling), so memory
scales `O(m · N_DOF)`.

Implementation lives at [`roadmap § 7 of MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md);
not yet wired into production paths.  Expected first-implementation
target: `calc_inductance.py` scalar path, ~3 days effort.

---

## 8. Failure Modes and Diagnostics

| Symptom | Likely cause | Diagnostic |
|---|---|---|
| Karl diverges (dZ ↑) | `α·L > 1` — try `--esim-relax 0.3` | inspect `esim_history.dZ`, look for oscillation |
| Karl plateaus at `dZ ≈ 0.5` | non-monotone BH curve | print `μ(H)` from the loaded BH file at sample H values |
| Karl converges to wrong `Z_s` | wrong sign on Robin BC | compare to scalar SIBC at low-H regime; `Z_s_esim ≈ Z_s_dowell` should hold |
| Karl 50+ iter without convergence | `α` too small OR BH knee straddled | try `α = 0.7` with `max_iter = 30`; check `H_t_rms` vs the BH knee |
| Karl hits `max_iter` but `Z_s_abs` / `H_t_rms` look plateaued | per-DOF noise floor on worst DOF, not divergence (see § 3.4) | plot `esim_history` with [`plot_karl_history.py`](../../examples/ih_esim_benchmark/plot_karl_history.py); if integrated quantities are monotone, accept the cap or raise `--esim-tol` to 5e-3 |
| Cell Picard inner-loop diverges | `tol` too tight for the BH curve smoothness | raise cell `tol` to 1e-4 (currently hard-coded — see [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) § 2.3) |
| Per-DOF Karl gives wildly different result from scalar | mass-lumping `M_lump` close to zero on some DOFs | check `min(M_lump)` — should be `>10⁻⁶ × max(M_lump)`; coarsen mesh if not |

Every Karl iteration is logged to `esim_history` in the JSON output
([`USAGE.md`](USAGE.md) § 4); inspecting this history is the first step
in any post-mortem.

---

## 9. Performance Characterisation

| Operation | Cost | Bottleneck |
|---|---|---|
| Cell-problem `solve(H_t)` | ~50 μs at `n_nodes=100` | tridiagonal `spsolve` |
| BIE outer solve, 166 DOFs | ~0.2 s | dense LU |
| BIE outer solve, 5k DOFs (HACApK) | ~5 s | ACA + GMRES |
| FEM-Kelvin solve, 12k DOFs (pardiso) | ~12 s | sparse direct |
| FEM-coilmesh solve, 38k DOFs (pardiso) | ~25 s | sparse direct |
| Karl loop (typical, 5 iter, BIE) | ~1 s + post-proc | outer-solve cost × iter count |
| Karl loop (typical, 5 iter, FEM-coilmesh) | ~125 s | re-assembly + solve per iter |
| Per-DOF Karl (166 DOFs × 5 iter) | +7 s | cell-solve `N_DOF × iter` overhead |
| Per-DOF Karl (5k DOFs × 5 iter) | +200 s | dominates total time |

All numbers are LAB (Windows, MKL, NGSolve 6.2.2603, radia 4.55.3).
Production deploys (100号機, mdx) match within 10%.

---

## 10. Code Layout Summary

```
src/radia/
  esim_cell_problem.py        # Layer 1: cell problem (FD tridiagonal)
  em_material.py              # Material wrapper + create_esim_solver factory
  bem_sibc_solver.py          # Layer 2 (BIE): ScalarBIESIBCSolver
  peec_bundle.py              # Layer 2 (PEEC coil): build_loop_bundle_impedance
  workpiece_surface.py        # Telegen reciprocity helpers
  panels/
    calc_inductance.py        # Layer 4: BEM-SIBC dispatch
    calc_fem_kelvin.py        # Layer 4: HCurl + Kelvin
    calc_fem_coilmesh.py      # Layer 4: HCurl A-V volumetric coil
  radia_ih.py                 # Layer 3 GUI (PySide6) that drives all four CLIs

examples/
  ih_esim_benchmark/          # benchmark.py + analytical_bessel_baseline.py + results.json
  induction_heating/          # esim_demo.py / esim_induction_heating_demo.py
  # (canonical ESIM lives in src/radia/esim_cell_problem.py; the old
  #  examples/effective_surface_impedance research scripts were removed 2026-06-27)

tests/
  test_esim_integration.py    # cell-problem + ESI table + coupled solver tests
  panels/golden/              # peec_bem_*.json + fem_kelvin_*.json + fem_coilmesh_*.json

packages/radia-mcp/src/radia_mcp/
  radia_ngsolve/knowledge/esim.py    # general ESIM MCP knowledge
  ih/sibc_knowledge.py               # IH-specific SIBC + Karl MCP knowledge
```

The cell solver is a single ~1600-line Python file with no external
dependencies beyond `numpy` and `scipy.sparse` + `scipy.special.iv`.
This makes the ESIM core easy to port, audit, and verify
independently — a deliberate design choice motivated by the academic
publication path.

---

**Document version**: 2026-05-18 (radia v4.55.3+).
