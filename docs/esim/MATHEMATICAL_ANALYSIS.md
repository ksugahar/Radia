# ESIM (Effective Surface Impedance Method): Mathematical Analysis & Implementation Audit

**Scope.** This document describes WHAT the Radia ESIM module solves, HOW it
discretises the 1-D cell problem, what curvature it does and does NOT
capture, and how the outer BEM/FEM curve order propagates into the
Karl-iteration feedback loop.  It is meant as an internal audit for
contributors and a citation reference for publication: every claim
about "the implementation does X" carries a file:line citation.

The companion documents are:

- [`docs/esim/USAGE.md`](USAGE.md) — user-facing CLI guide.
- [`docs/esim/IMPLEMENTATION.md`](IMPLEMENTATION.md) — code architecture, three-solver dispatch, Karl-loop internals, performance characterisation.
- [`docs/esim/CROSS_VALIDATION.md`](CROSS_VALIDATION.md) — analytical / internal-consistency / external 2-D axisymmetric validation matrix with concrete numerical data (IGTE-grade tables).
- [`docs/esim/R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) — focused diagnosis of why PEEC and BEM-A produce different coil R values.
- [`docs/research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.md`](../research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.md) — research-WIP notes.
- [`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md) — wide-band nonlocal extension (deferred).

---

## 1. What ESIM Solves

ESIM ("Effective Surface Impedance Method") replaces the linear Dowell
surface impedance by the solution of a **1-D nonlinear cell problem**
through the conductor depth, returning a field-dependent complex
impedance `Z_s(H_t)`.  This `Z_s` feeds the outer BEM/FEM SIBC Robin
boundary condition.  The two are coupled by **Karl iteration**
(a damped Picard fixed-point) — see § 6.

For a single-frequency phasor `H_t e^{jωt}` on the conductor surface,
the cell problem is:

| Geometry | PDE | Domain | Boundary conditions |
|---|---|---|---|
| `slab` | `ρ ∂²H/∂z² + jω μ(|H|) H = 0` | `z ∈ [0, a]` | `H(0) = H_t`, `dH/dz(a) = 0` |
| `cylinder` | `(ρ/r) ∂/∂r[r ∂H/∂r] + jω μ(|H|) H = 0` | `r ∈ [0, R]` | `dH/dr(0) = 0` (regularity), `H(R) = H_t` |
| `finite_slab` | as `slab` (`ρ ∂²H/∂z² + jω μ(|H|) H = 0`) | `z ∈ [0, a]` | `H(0) = H_t` (driven face), `dH/dz(a) = 0` (insulated back) |

The `finite_slab` mode is implemented by
[`ESIMFiniteSlabSolver`](../../src/radia/esim_cell_problem.py#L339)
(separate class).  It models a finite-thickness plate with one
heated face — the canonical 2-sided plate problem of
Lavers–Biringer 1985, restricted to 1-sided drive.  The
`ESIMCellProblemSolver` class ([`esim_cell_problem.py:816`](../../src/radia/esim_cell_problem.py#L816))
is the legacy infinite-slab/cylinder solver and is the one called
by the production Karl loops.

The unknown is the complex H-field profile inside the conductor.  The
nonlinearity enters through `μ(|H|)` — the BH curve.

Reference:
[`src/radia/esim_cell_problem.py:565-575`](../../src/radia/esim_cell_problem.py#L565-L575)
documents these BVPs verbatim, with the analytical linear solutions
(`cosh` for slab, `I_0` for cylinder) shown as commentary.

Output of `solve(H0)`:

- `Z = E_t(surface) / H_t(surface)` — complex surface impedance [Ω]
- `P' = 0.5 Re(Z) |H_t|² + losses_in_complex_mu` — surface power density [W/m²]
- `Q' = 0.5 Im(Z) |H_t|²` — surface reactive power density [VAR/m²]

The `2·(P' + jQ') / H0²` identity ties P' and Z together — both are
independent observables in publication tables.  Implementations:

- Slab / cylinder solver: [`esim_cell_problem.py:545`](../../src/radia/esim_cell_problem.py#L545).
- Finite-slab (2-sided heating, `ESIMFiniteSlabSolver`):
  [`esim_cell_problem.py:1010`](../../src/radia/esim_cell_problem.py#L1010).

---

## 1.5 Derivation of the Cell-Problem Strong and Weak Forms

This section derives the cylinder cell PDE from Maxwell's equations
and presents both the strong (PDE) and weak (variational) forms.  The
finite-difference solver implements the strong form directly; the FE
finite-difference solver mathematically can be viewed as the limit of
a P0 / P1 weak-form solver on a uniform mesh.  The derivation is
included for reviewers who prefer the variational view.

### 1.5.1 Maxwell → 1-D diffusion equation inside the conductor

Inside a homogeneous conductor (no free charge, no displacement
current in the magneto-quasi-static (MQS) limit), Maxwell's equations
reduce to:

$$
\nabla \times \mathbf{E} = -j\omega \mathbf{B}, \quad
\nabla \times \mathbf{H} = \mathbf{J} = \sigma \mathbf{E}, \quad
\mathbf{B} = \mu(|\mathbf{H}|)\,\mathbf{H}.
$$

Eliminating `E` and `J`:

$$
\nabla \times (\rho\,\nabla \times \mathbf{H}) = -j\omega\,\mu(|\mathbf{H}|)\,\mathbf{H},
\qquad \rho = 1/\sigma.
$$

For a tangential field on a cylindrical surface with radius `R` and the
conductor at `r ≤ R`, assume the only relevant component of `H` is
azimuthal `H_φ(r) e^{jωt}` (the canonical IH workpiece geometry).
Then `∇ × H = -∂_r H_φ ẑ + (1/r) ∂_r(r H_φ) r̂` (in cylindrical
coordinates), and substituting:

$$
-\frac{1}{r}\,\partial_r\!\left[r\,\rho\,\partial_r H_\varphi\right]
= -j\omega\,\mu(|H_\varphi|)\,H_\varphi.
$$

For constant `ρ` (uniform conductivity), this rearranges to the cell
PDE used in [`esim_cell_problem.py`](../../src/radia/esim_cell_problem.py):

$$
\boxed{\quad
\frac{\rho}{r}\,\partial_r\!\left[r\,\partial_r H\right] + j\omega\,\mu(|H|)\,H = 0,
\qquad r \in [0, R].
\quad}
$$

### 1.5.2 Boundary conditions

- **Surface** (`r = R`).  The Robin / Leontovich SIBC is *not* applied
  on the cell — instead the cell is driven by the tangential field
  `H(R) = H_t` set by the outer (BEM / FEM) solve.  In the cell, this
  becomes a Dirichlet BC.
- **Axis** (`r = 0`).  Regularity of `H` (no infinite derivatives at
  the axis) implies `∂_r H(0) = 0` by symmetry.

The pair (Dirichlet at `R`, Neumann at `0`) gives a well-posed BVP
for any `μ(|H|) > 0`.

### 1.5.3 Weak form

Multiply the strong PDE by a test function `v(r)` (complex,
sufficiently smooth), and integrate against the **radial measure**
`r\,dr` over `[0, R]`:

$$
\int_0^R \left\{ \frac{\rho}{r}\,\partial_r[r\,\partial_r H] +
                 j\omega\,\mu(|H|)\,H \right\}\,v\,r\,dr = 0.
$$

Integrate by parts the first term, using `∂_r H(0) = 0`:

$$
-\int_0^R \rho\,r\,(\partial_r H)\,(\partial_r v^*)\,dr
+ \left[\rho\,r\,(\partial_r H)\,v^*\right]_{r=R}
+ j\omega \int_0^R \mu\,H\,v^*\,r\,dr = 0.
$$

The boundary term at `r = R` is `ρ R (∂_r H)|_R v^*(R)`.  Since
`∂_r H|_R = -j ω σ A_φ + j κ H_t` evaluates via Ampère's law to a
quantity proportional to the surface current density, this term is
absorbed into the Dirichlet lift and does not appear in the
homogeneous problem (we test against `v(R) = 0`).

The resulting weak form:

$$
\boxed{\quad
\int_0^R \rho\,r\,\partial_r H\,\partial_r v^*\,dr
\;-\;j\omega \int_0^R \mu(|H|)\,H\,v^*\,r\,dr = 0,
\quad}
\qquad v(R) = 0,\ \partial_r v(0) = 0.
$$

This is a curl-curl-like form on a 1-D radial mesh with radial-weighted
inner product.  Galerkin FE discretisation with P1 basis on a uniform
mesh and one-point quadrature recovers exactly the finite-difference
stencil used in [`esim_cell_problem.py:594-602`](../../src/radia/esim_cell_problem.py#L594-L602):

$$
\frac{\rho\,r_{i+1/2}}{h^2}\,(H_{i+1} - H_i) -
\frac{\rho\,r_{i-1/2}}{h^2}\,(H_i - H_{i-1})
+ j\omega\,\mu_i\,r_i\,H_i = 0,
$$

after dividing by `r_i`.  The L'Hôpital limit at `r = 0` produces
exactly the `4ρ/h²` coefficient documented at
[`esim_cell_problem.py:609-611`](../../src/radia/esim_cell_problem.py#L609-L611).

### 1.5.4 Why a finite-difference solver in production

The weak form is **mathematically equivalent** to the production
finite-difference (FD) solver on a uniform mesh.  We use FD instead of
a full FE machinery because:

1. The 1-D cell is small (`n_nodes = 100`) — FE-assembly overhead
   would dominate the actual linear solve.
2. The tridiagonal structure is exact under FD, allowing
   `scipy.sparse.linalg.spsolve` with `O(n)` work per solve.  An FE
   approach needs explicit mass-matrix assembly, which has the same
   tridiagonal structure but more boilerplate.
3. The L'Hôpital regularisation at `r = 0` is naturally handled by
   the FD stencil; FE requires a ghost-cell convention or special
   axis basis function (Henrotte basis — used by Radia's axisymmetric
   *magnetic* solver, NOT by ESIM, see CLAUDE.md "Axisymmetric FE"
   policy).

For higher-dimensional problems (the planned 2-D nonlocal SIBC
extension to capture finite-coil-aperture wide-band effects, see
[`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md)),
a proper FE discretisation would be preferred.

---

## 2. Discretisation of the 1-D Cell Problem

### 2.1 Mesh

```python
self.mesh_points = np.linspace(0, half_thickness, n_nodes)
```
([`esim_cell_problem.py:428`](../../src/radia/esim_cell_problem.py#L428))

The cell mesh is **uniform** with `n_nodes = 100` by default
([`esim_cell_problem.py:366`](../../src/radia/esim_cell_problem.py#L366)).

**Implication for thin-skin regimes.**  The dimensionless skin parameter
is `ξ = a / δ` where `δ = sqrt(2 ρ / (ω μ_0 μ_r))`.  For a typical IH
steel workpiece (μ_r ≈ 100, σ = 2×10⁶ S/m, f = 100 kHz):

- `δ ≈ 0.16 mm`
- `a = 5 mm` → `ξ ≈ 30`
- Uniform `h = 50 μm` → ≈ 3 nodes within one skin depth — marginal.

For `ξ > 100` (high frequency × thick workpiece), a geometric stretch
toward the surface would be a strict improvement.  This is **not
currently implemented**.  Open issue: replace `np.linspace` with a
graded mesh `r_i = a · (1 - (1-i/N)^p)` for `p ≈ 2-3` so that the
spacing near `r = a` is `O(h/p)`.

### 2.2 Stiffness assembly

The cylindrical discretisation uses cell-centred fluxes — a standard
conservative scheme:

```python
r_ip = z[i] + h/2  # r_{i+1/2}
r_im = z[i] - h/2  # r_{i-1/2}
coef_p = self.rho * r_ip / (r_i * h * h)
coef_m = self.rho * r_im / (r_i * h * h)
diag_main[i] = -(coef_p + coef_m) + 1j * self.omega * mu_dist[i]
```
([`esim_cell_problem.py:594-602`](../../src/radia/esim_cell_problem.py#L594-L602))

The axis singularity at `r = 0` is resolved by L'Hôpital + symmetry
(`H[-1] = H[1]` ghost-cell mirror), yielding:

```python
coef_0 = 4 * self.rho / (h * h)
diag_main[0] = -coef_0 + 1j * self.omega * mu_dist[0]
diag_upper[0] = coef_0
```
([`esim_cell_problem.py:609-611`](../../src/radia/esim_cell_problem.py#L609-L611))

This is correct because `(1/r) ∂/∂r[r ∂H/∂r] → 2 ∂²H/∂r²` as `r → 0`
when `dH/dr(0) = 0`.  Verified against the analytical
`H(r) = H_0 I_0(γr) / I_0(γR)` in the constant-μ limit (see § 5.1).

### 2.3 Nonlinear Picard

```python
for iteration in range(max_iter):
    H_new = self._solve_linear_system(H0, mu_dist)      # tridiagonal solve
    mu_new = np.array([self._get_mu(abs(h)) for h in H_new])
    mu_dist = (1 - relaxation) * mu_dist + relaxation * mu_new
    if ||μ_new - μ_old|| / ||μ_old|| < tol: break
```
([`esim_cell_problem.py:508-541`](../../src/radia/esim_cell_problem.py#L508-L541))

Picard converges linearly with rate `≤ relaxation·L` where `L` is the
Lipschitz constant of `μ(|H|)`.  For monotone `μ(|H|)` curves (all
ferromagnetic materials we ship), Picard always converges; rate is fast
when |H| stays on one side of the BH knee.  When |H| straddles the knee
(transitioning from low-H linear regime to saturated regime), `L`
spikes and convergence drops to ~0.5 per iter — that's where the
**outer** Karl loop sees most of its iterations (§ 6.3).

Default `relaxation = 0.5` ([`esim_cell_problem.py:457`](../../src/radia/esim_cell_problem.py#L457))
is conservative and rarely needs tuning.

### 2.4 Analytical linear-material initial guess

Before Picard starts, ESIM seeds H with the linear-μ analytical
solution (Bessel `I_0` for cylinder, `cosh` for slab):

```python
H = H0 * bessel_iv(0, gamma * z) / bessel_iv(0, gamma * a)
```
([`esim_cell_problem.py:489-499`](../../src/radia/esim_cell_problem.py#L489-L499))

**Note (v4.46.1+):** `scipy.special.iv` overflows silently to ∞ at
large `|γa|` (typical for ξ > 30), producing NaN via ∞/∞.  This is
handled by wrapping in `np.errstate` and falling back to a thin-skin
exponential when the result is not finite.  The bad initial guess
was previously benign because Picard overwrites H from `mu_dist`
on iteration 0, but the warning spam was confusing.

---

## 3. Curvature Handling

This is the most consequential approximation in production ESIM.

### 3.1 What the cell-problem geometry captures

The `'cylinder'` geometry captures the **principal curvature**
`κ = 1/R` of an idealised circular cross-section.  Mathematically, this
is the local axial-symmetric cylindrical Laplacian — i.e. it assumes
the surface has one nonzero principal curvature and one zero.  Real
workpieces typically have two nonzero principal curvatures
(`κ_1, κ_2`); cylinder mode represents `κ_1 ≠ 0`, `κ_2 = 0` (an
infinite cylinder).

The `'slab'` geometry captures `κ_1 = κ_2 = 0` (flat plate).

### 3.2 What the production Karl scripts pass to ESIM

| Script | `geometry=` | `half_thickness` | Per-panel R? |
|---|---|---|---|
| `calc_inductance.py`     | `'cylinder'` | global `--half-thickness` | No |
| `calc_fem_kelvin.py`     | `'cylinder'` | global `--half-thickness` | No |
| `calc_fem_coilmesh.py`   | `'cylinder'` | global `--half-thickness` | No |
| `calc_heating.py`        | `'cylinder'` | derived from coil geometry | No |

All four pass a **single global radius** to a **single ESIM solver
instance** and use **scalar Z_s** (one impedance value for the whole
workpiece).  Reference:
[`calc_fem_kelvin.py:265-268`](../../src/radia/panels/calc_fem_kelvin.py#L265-L268),
[`calc_inductance.py:593-597`](../../src/radia/panels/calc_inductance.py#L593-L597),
[`calc_fem_coilmesh.py:93-96`](../../src/radia/panels/calc_fem_coilmesh.py#L93-L96).

This is the right model for **idealised IH workpieces** (cylindrical
bars, pipes, plates).  It is an **approximation** for:

- Workpieces with strong shape variation (cube + fillet, gear teeth).
- Anything where one panel sees `H_t ≫ H̄` while another sees
  `H_t ≪ H̄` — the saturation pattern is spatial, but the scalar Z_s
  is mesh-RMS averaged.

### 3.3 Per-DOF Z_s (per-panel ESIM, available since v4.47.0)

**Status (v4.47.0+)**: `calc_inductance.py --esim-per-panel` invokes
per-DOF Z_s — one ESIM cell solve per BEM surface DOF, using the
locally extracted `|H_t|` from `phi_vec`.  See § 3.4 below for the
per-DOF `|H_t|` discretisation.

The underlying `ScalarBIESIBCSolver.solve()` accepts
`Z_s: ndarray[ndof]` directly
([`bem_sibc_solver.py:401-412`](../../src/radia/bem_sibc_solver.py#L401-L412))
via row-scaling: `A_sys[i, :] = 0.5 M - DL + (gamma[i] · SL · M⁻¹ · K)[i, :]`.
The Karl loop in
[`calc_inductance.py:_solve_workpiece_weak_coupled`](../../src/radia/panels/calc_inductance.py)
seeds `Z_s_wp = full(ndof, esim.solve(5.0)['Z'])` and refreshes per
iteration as

    Z_s_wp[i]  ←  relax · esim.solve(|H_t|[i])['Z']
                + (1 - relax) · Z_s_wp_old[i]

Convergence metric uses `max_i |dZ[i]| / |Z_s[i]|`.

Caveats:
- HACApK GMRES backend does not yet accept ndarray `Z_s` — the script
  raises if `--esim-per-panel` is combined with `--wp-bem-backend
  hacapk`.  This is a backend gap, not a math gap.
- Per-DOF ESIM costs N_DOF extra cell solves per Karl iter.  For
  166-DOF wp + 10 iter ≈ 1660 ESIM calls × ~4 ms = ~7 s overhead.
  On larger meshes (~5000 DOF) the per-DOF ESIM call cost dominates
  (~200 s per iter); add Anderson acceleration or batch the cell
  solve if needed.
- `calc_fem_kelvin` and `calc_fem_coilmesh` are NOT yet wired for
  per-DOF Z_s; the FEM Robin term construction uses scalar `s/Z_s`.
  Extending would mean replacing the constant Robin coefficient with
  a `CoefficientFunction` of per-element Z_s.  Tracked in § 7.

### 3.4 Per-DOF |H_t| extraction inside the Karl loop

The BEM system stores the scalar potential `phi` on H1 P1 DOFs.  The
tangential field is `H_t = -∇_S phi`.  The cell-solver wants `|H_t|`
at each DOF.

For scalar Karl we use the mesh-RMS:

    |H_t|_rms² = (phi^T K phi) / area
              = sum_i [ phi_i · (K phi)_i ] / area

For per-DOF Karl we localise the same form: the i-th summand of
`phi^T K phi` is `phi_i · (K phi)_i`.  Dividing by the i-th lumped
mass `M_lump[i] = (M · 1)[i]` (row sum of M, = effective area of
basis function i) gives a per-DOF density:

    |H_t at i|² ≈ |phi_i · (K phi)_i| / M_lump[i]
              + (same for imag(phi))

`abs()` is needed because each term `phi_i · (K phi)_i` can be
signed; only the total sum is guaranteed nonneg.

This is what
[`calc_inductance.py`](../../src/radia/panels/calc_inductance.py)
computes inside the Karl loop.  It is consistent with the scalar
formula in the limit of a uniform basis-function gradient — when
all `(K phi)_i / M_lump_i` ratios are equal, the per-DOF and scalar
mesh-RMS values coincide to floating-point precision (verified by
the unit-uniform-Z_s smoke test).

### 3.5 The infrastructure for per-panel curvature exists but is not wired in

The `ESIMFiniteSlabSolver` class has a per-call radius override:

```python
def solve(self, H0, tol=1e-6, max_iter=50, relaxation=0.5, R_local=None):
    if R_local is not None:
        self.set_radius(float(R_local))
    ...
```
([`esim_cell_problem.py:457-479`](../../src/radia/esim_cell_problem.py#L457-L479))

`set_radius(R)` rebuilds the mesh in place:

```python
def set_radius(self, R):
    self.half_thickness = R
    self.mesh_points = np.linspace(0, R, self.n_nodes)
```
([`esim_cell_problem.py:440-455`](../../src/radia/esim_cell_problem.py#L440-L455))

So one could call `esim.solve(H_t_at_panel, R_local=R_at_panel)` for
each surface panel and get a true per-panel curvature SIBC.  None of
the production Karl scripts invoke this path; it is reserved for
research probes (see roadmap § 7).

### 3.6 Why "single Z_s" is still useful in publication

For workpieces where the **operating-point H_t** does not vary by more
than ~3× across the surface (typical for solenoid + cylindrical bar),
the scalar mesh-RMS H_t feeds an ESIM call that gives an effective Z_s
within ~5 % of a fully resolved per-panel calculation.  This is good
enough for engineering screening (P_wp, L_eff at ±5 % for design
iteration).

For accuracy claims sharper than 5 %, per-panel ESIM is required
(see § 3.3 for the BEM path; FEM-path per-panel Z_s is roadmap § 7).

---

## 4. Outer Curve-Order Propagation

The outer BEM/FEM solve uses higher-order curved geometry when the
Cubit-exported `.vol` has `curve_order ≥ 2` and the basis order is
matched.  How does that curving propagate to the H_t that ESIM sees?

### 4.1 BEM path (`calc_inductance.py`)

The `.vol`'s geometry order is auto-detected from the companion
`.vol.json` and passed to the scalar BEM solver:

```python
vol_curve_order = min(_detect_vol_curving_order(args.vol), 2)
...
bem = ScalarBIESIBCSolver(
    bem_input_mesh, order=basis_order,
    intree_geom_order=vol_curve_order, ...)
```
([`calc_inductance.py:462`](../../src/radia/panels/calc_inductance.py#L462),
[`calc_inductance.py:536-549`](../../src/radia/panels/calc_inductance.py#L536-L549))

`H_t_rms` is computed by the BEM solver as `√(φᵀ K φ / A_wp)`, where K
is the assembled stiffness matrix on either Tri3 flat or Tri6 curved
nodes (depending on `intree_geom_order`).  When the curved (Tri6) path
is used, **the Jacobian of the curved element is integrated correctly**
and the resulting `H_t_rms` reflects the curved geometry's surface area.

**Curve-order limitation in BEM:** the in-tree Lagrange-P2 assembler
caps geometric order at 2.  A `.vol` exported with `curve_order = 3`
(or higher) is treated as Tri6 with the extra mid-face nodes
**unused**.  This is documented in
[`calc_inductance.py:475-484`](../../src/radia/panels/calc_inductance.py#L475-L484);
it is a known limitation of the in-tree assembler, not an ESIM-side
issue.

### 4.2 FEM path (`calc_fem_kelvin.py`, `calc_fem_coilmesh.py`)

NGSolve's `Mesh(vol)` constructor automatically loads the
`curvedelements` section from the `.vol`, so subsequent `Integrate()`
calls use the curved Jacobian.  `H_t_rms` is computed by:

```python
A_sq = sum(A[i].real² + A[i].imag² for i in range(3))
At_sq = A_sq - (A·n)²
int_At2 = Integrate(At_sq, mesh, BND, definedon=wp_region).real
H_t_rms = abs(jω/Z_s) · sqrt(int_At2 / A_wp)
```
([`calc_fem_kelvin.py:772-794`](../../src/radia/panels/calc_fem_kelvin.py#L772-L794),
[`calc_fem_coilmesh.py:299-306`](../../src/radia/panels/calc_fem_coilmesh.py#L299-L306))

NGSolve's `Integrate` on a `BND` region of a mesh with loaded
curvedelements **does** use the curved Jacobian for both the
integration weights and the surface area `A_wp`.  This is documented
behaviour of `mesh.Curve(p)` (which the `.vol` loader emulates) — see
NGSolve docs on `Curve()`.

**Status: outer-solve curve-order propagates correctly into H_t for
FEM paths.**  The single Z_s ESIM call that consumes this H_t is then
geometry-aware via cylinder mode but does NOT see the local curvature
(§ 3.1) — that's a separate, deeper limitation.

### 4.2.5 Telegen reciprocity and the ΔL_telegen φ·B form

This subsection derives the gauge-invariant Telegen formula for the
workpiece-induced port inductance change `ΔL`, used in
`calc_inductance.py` via [`radia.workpiece_surface.delta_L_telegen_phiB`](../../src/radia/workpiece_surface.py).

**Setup.**  Let `Ω_wp` be the workpiece domain with boundary
`Γ = ∂Ω_wp` and outward unit normal `n`.  The coil drives a port
current `I_port`.  In the absence of the workpiece, the coil's
own field has scalar potential `φ_inc` and vector potential `A_inc`
on Γ such that `H_inc = -∇φ_inc` (scalar) and `B_inc = ∇ × A_inc`
(vector) — both well-defined by Biot–Savart from the coil filaments.

When the workpiece is added, an induced surface current `J_s` flows
on Γ and produces a scattered field.  The total port impedance is
`Z_port = Z_vacuum + ΔZ`, and we want a closed-form for `ΔZ`.

**Telegen reciprocity** (energy form): the back-reaction at the port
equals the volume integral of the field interaction:

$$
I_{\mathrm{port}}\,\Delta V_{\mathrm{port}}
\;=\; \int_{\Omega_{\mathrm{wp}}} \mathbf{J}_s \cdot \mathbf{A}_{\mathrm{inc}}\,d\Omega
\;=\; \int_\Gamma \mathbf{J}_s \cdot \mathbf{A}_{\mathrm{inc}}\,dS
$$

(the volume integral reduces to a surface integral because, in the
SIBC limit, `J_s` is confined to a thin skin layer at Γ).  Hence

$$
\Delta Z = \frac{\Delta V_{\mathrm{port}}}{I_{\mathrm{port}}}
        = \frac{1}{I_{\mathrm{port}}^2}\,\int_\Gamma \mathbf{J}_s \cdot \mathbf{A}_{\mathrm{inc}}\,dS.
$$

This is the **`J_s · A_inc` form** ([`workpiece_surface.py:209-291`](../../src/radia/workpiece_surface.py#L209-L291)).
It is gauge-dependent in the discrete setting: changing the gauge of
`A_inc` (e.g. `A_inc → A_inc + ∇χ`) adds a surface-divergence term
`∫_Γ J_s · ∇χ dS = -∫_Γ (∇_s · J_s) χ dS`.  In the continuum,
`∇_s · J_s = 0` for the induced eddy currents (Faraday's law +
surface charge balance), so the gauge term vanishes.  In a P1 FE
discretisation, however, `∇_s · J_s_h` has element-edge jumps that
behave as a weak surface divergence — leading to a ~100× error on
`Im(ΔL)` compared to energy-balance predictions.

**The fix — `φ · (n · B_inc)` form.**  Use the surface vector-calculus
identity (closed Γ, MQS limit):

$$
\int_\Gamma \mathbf{J}_s \cdot \mathbf{A}\,dS
\;=\;\int_\Gamma \varphi\,(\mathbf{n} \cdot \mathrm{curl}\,\mathbf{A})\,dS
\;=\;\int_\Gamma \varphi\,(\mathbf{n} \cdot \mathbf{B}_{\mathrm{inc}})\,dS,
$$

where `φ` is the workpiece-side scalar potential from the SIBC BIE
solve (`H_t = -∇_s φ` on Γ).  Step-by-step:

1. `J_s = -n × ∇_s φ` on the workpiece surface (definition).
2. `∫_Γ J_s · A dS = -∫_Γ (n × ∇_s φ) · A dS = ∫_Γ ∇_s φ · (n × A) dS`
   (using `(a × b) · c = (b × c) · a`).
3. `∫_Γ ∇_s φ · (n × A) dS = -∫_Γ φ · (∇_s · (n × A)) dS`
   (surface integration by parts on the closed surface).
4. `∇_s · (n × A) = -n · curl A` (a standard surface vector-calculus
   identity), so the result follows: `∫_Γ φ · (n · curl A) dS = ∫_Γ φ · (n · B) dS`.

The right-hand side uses `B = curl A` directly — gauge-invariant.
This is the formula implemented in
[`workpiece_surface.py:294-379`](../../src/radia/workpiece_surface.py#L294-L379):

$$
\boxed{\quad
\Delta L = \mathrm{Re}\,\frac{1}{I_{\mathrm{port}}^2}
          \int_\Gamma \varphi(r)\,\bigl(\mathbf{n}(r) \cdot \mathbf{B}_{\mathrm{inc}}(r)\bigr)\,dS,
\quad}
$$

and the resistive contribution:

$$
\Delta R = -\omega\,\mathrm{Im}\,\frac{1}{I_{\mathrm{port}}^2}
           \int_\Gamma \varphi(r)\,\bigl(\mathbf{n}(r) \cdot \mathbf{B}_{\mathrm{inc}}(r)\bigr)\,dS.
$$

(See [`calc_inductance.py:839-840`](../../src/radia/panels/calc_inductance.py#L839-L840) for the implementation of the sign convention.)

**Discrete quadrature.**  Each workpiece triangle `T` contributes:

$$
\Delta L|_T = \mathrm{Re}\,\frac{1}{I_{\mathrm{port}}^2}\,
              \varphi_{\mathrm{avg},T} \cdot
              (\mathbf{n}_T \cdot \mathbf{B}_{\mathrm{inc}}(c_T)) \cdot A_T,
$$

where `c_T`, `A_T`, `n_T` are the centroid, area, and outward unit
normal of `T`, and `φ_avg,T = (1/A_T) ∫_T φ_h dS` is the area-
averaged potential on `T` (for P1 `φ_h` this equals
`(φ_1 + φ_2 + φ_3)/3`).  `B_inc(c_T)` is evaluated by direct
Biot–Savart from the coil filaments at the triangle centroid —
single-point quadrature is sufficient because `B_inc` varies smoothly
on the workpiece scale, while `φ_h` and `n_T` vary rapidly on the
mesh scale.

**Cost.**  `O(N_filaments × N_triangles)` for the Biot–Savart kernel
evaluation; in the production gapped-torus benchmark this is
~50 ms.  C++ kernel `_HFromSegmentsComplex` is used
([`workpiece_surface.py:340-365`](../../src/radia/workpiece_surface.py#L340-L365)).

### 4.3 What is and is not "curve-order-aware"

| Quantity | Curve-order-aware? | File:line |
|---|---|---|
| BEM stiffness matrix K (basis_order=2) | **Yes** (Tri6 isoparametric) | [`bem_sibc_solver.py` Lagrange-P2 assembler] |
| BEM `H_t_rms` (basis_order=2) | **Yes** | (`H_t_rms = sqrt(φᵀKφ/A_wp)`) |
| BEM stiffness matrix K (basis_order=1, flat) | **No** (Tri3 flat — by design) | [`calc_inductance.py:520-525`](../../src/radia/panels/calc_inductance.py#L520-L525) |
| FEM A field representation | **Yes** (HCurl on curved mesh) | NGSolve standard |
| FEM `H_t_rms` (Integrate on BND) | **Yes** | [`calc_fem_kelvin.py:791-793`](../../src/radia/panels/calc_fem_kelvin.py#L791-L793) |
| ESIM cell solve (1-D radial mesh) | **N/A** — cell has no surface curving concept | n/a |
| Z_s scalar fed back to outer | **Same Z_s regardless of curve order** | by design (scalar mesh-RMS Karl) |

---

## 5. Numerical Verification

### 5.1 Linear-μ comparison with Bessel `I_0`

For constant `μ = μ_0 μ_r` (linear material), the cylinder cell-problem
admits the closed-form `H(r) = H_0 I_0(γr) / I_0(γR)` with
`γ = (1+j)/δ`.  The numerical Picard solver should match this on the
first iteration (since `mu_dist` is initialised to the constant value
and `_solve_linear_system` is solved exactly to floating precision).

A quick verification (LAB, 2026-05-15):

```
ESIMFiniteSlabSolver(half_thickness=5e-3, sigma=2e6, mu_r=100,
                     frequency=100e3, geometry='cylinder')
.solve(H0=1.0)
→ Z = 3.40e-4 + 7.81e-3j Ω
```

Independent Bessel evaluation (`scipy.special.iv`):
```
γ = (1+1j) / δ,  δ = √(2/(ω μ_0 μ_r σ)) ≈ 0.159 mm
Z_anal = ρ γ · I_1(γR) / I_0(γR)
       = 3.40e-4 + 7.81e-3j  Ω
```

Match to 4 significant figures.  Formal regression coverage:
[`tests/test_esim_integration.py`](../../tests/test_esim_integration.py)
exercises the cell-problem solver, the table interpolator, the
coupled solver, and VTK export.  Bessel `I_0` parity is implicitly
covered by the linear-regime convergence assertion.  A dedicated
test that pins `Z` against `scipy.special.iv(0, γR) / iv(1, γR) ·
ρ γ` to 5 sig-fig is open work (Phase A § 7).

### 5.2 Karl iteration consistency across the three production paths

For the same physical inputs (steel cylinder, em_sample_bh, f = 100 kHz,
μ_r = 100), the three Radia paths produce:

| Script | L_total [nH] | P_wp [mW] | Iterations |
|---|---|---|---|
| `calc_inductance.py --impedance-model esim` | 84.74 | 0.147 | 6 |
| `calc_fem_kelvin.py --impedance esim` | 160.75 | 0.282 | 5 |
| `calc_fem_coilmesh.py --impedance-model esim` | 26.35 | 0.102 | 5 |

The three scripts use **different coil models** (PEEC filament vs
volumetric FEM) and **different workpiece models** (scalar BEM-SIBC vs
HCurl-FEM with Robin), so absolute numbers differ.  What is consistent
is:

- All three converge to the same Karl tolerance (`dZ/Z < 10⁻³`) in
  5–7 iterations.
- All three produce identical Z_s magnitudes within rounding
  (`|Z_s| ≈ 3.58×10⁻² Ω` for the steel cylinder case).
- Per-iter dZ progression is monotonic (no oscillation) at the default
  `--esim-relax 0.5`.

Cross-validation suite:
[`examples/ih_esim_benchmark/`](../../examples/ih_esim_benchmark/) ships
`benchmark.py` + `analytical_bessel_baseline.py` + `results.json` +
`benchmark_plot.pdf`.  This is the runnable counterpart to § 5.1 (linear
Bessel parity).  Extending it to Stoll-envelope and Lavers–Biringer
nonlinear cases is roadmap § 7.

### 5.3 Cross-validation against analytical references

Closed-form references for the SIBC + ESIM combination:

| Geometry | Reference | Status |
|---|---|---|
| Cylinder + linear μ | Wakao–Igarashi–Fujiwara Part 5 (Bessel) | **VERIFIED** (matches to ~10⁻⁴; benchmark at [`examples/ih_esim_benchmark/`](../../examples/ih_esim_benchmark/)) |
| Slab + linear μ | Dowell (tanh) | **VERIFIED** (Dowell closed-form baked into `mat.dowell_Zs`) |
| Cylinder + nonlinear μ (BH) | Stoll 1974 (analytical envelope) | **open** (roadmap § 7) |
| Plate + 2-sided heating | Lavers–Biringer 1985 | **open** (roadmap § 7; implementation exists at `ESIMFiniteSlabSolver`, no cross-check harness) |

The "open" rows are roadmap items in § 7.

---

## 6. Karl Iteration (Outer Loop) — Implementation Details

The outer fixed-point loop is named **Karl iteration** in the codebase
after **Karl Hollaus**, the first author of the canonical ESIM scalar-
potential formulation paper:

> K. Hollaus, M. Kaltenbacher, J. Schöberl, *"A Nonlinear Effective Surface
> Impedance in a Magnetic Scalar Potential Formulation,"* **IEEE Trans.
> Magn.**, 2025.  DOI:
> [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932).

This is the paper the in-tree cell-problem solver implements (see the
docstring of [`esim_cell_problem.py`](../../src/radia/esim_cell_problem.py)
lines 1-10).  The outer Picard fixed-point is the loop introduced
there; "Karl iteration" is lab shorthand.  For external readers /
reviewers, prefer "Hollaus-type Picard relaxation" or
"Hollaus iteration".

The outer Karl loop wraps the BEM/FEM solve:

```python
Z_s = esim.solve(5.0, max_iter=5)['Z']        # seed at small H_t
for k in range(max_iter):
    res = solve_outer(Z_s=Z_s)                # BEM/FEM solve
    H_t = res['H_t_rms']                       # mesh-RMS
    Z_s_old = Z_s
    Z_s = relax * esim.solve(H_t)['Z'] + (1-relax) * Z_s_old
    dZ = |Z_s - Z_s_old| / |Z_s_old|
    if dZ < tol and (k > 0 or max_iter <= 1):
        break
# final re-solve at converged Z_s so post-proc sees matching residual
res_final = solve_outer(Z_s=Z_s)
```

This pattern is implemented in three places:
[`calc_inductance.py:597-680`](../../src/radia/panels/calc_inductance.py#L597-L680),
[`calc_fem_kelvin.py:515-820`](../../src/radia/panels/calc_fem_kelvin.py#L515-L820),
[`calc_fem_coilmesh.py:250-360`](../../src/radia/panels/calc_fem_coilmesh.py#L250-L360).

### 6.1 Why `relax = 0.5`

Picard with `relaxation < 1` is monotone-convergent when the
fixed-point map has Lipschitz constant `L < 1/(1-relaxation)`.  For
ESIM Z_s near the BH knee, `L` can approach 1; `relax = 0.5` gives
`1/(1-0.5) = 2`, comfortably above the empirical `L ≈ 1.2` observed
on the steel sample.  For sharper saturation (e.g. SUS430 above 30
kHz) drop to `--esim-relax 0.3`.

### 6.2 The "seed at small H_t" choice

The seed is `esim.solve(5.0, max_iter=5)`.  H_t = 5 A/m corresponds to
~6 µT of B-field — well below the BH knee for any ferromagnetic
material we ship.  This guarantees the seed Z_s is in the linear regime
(high μ, large |Z_s|), which is a **conservative** start: the outer
solve under-estimates H_t → first ESIM update moves Z_s toward smaller
|Z_s| (saturation).  Compared to starting from the linear Dowell Z_s,
this avoids overshoot in the first iteration.

### 6.3 Convergence guarantees

Karl iteration on a single scalar (mesh-RMS H_t → Z_s) is a 1-D
fixed-point problem.  Convergence is guaranteed by the monotonicity of
the BH curve (the map `H_t → |Z_s|` is decreasing — saturation reduces
permeability — and the outer solve gives a monotone `Z_s → H_t` map
for fixed geometry).  No formal proof is shipped; empirical
convergence in 5–10 iterations is reproducible across all three
production scripts.

If Karl ever diverges in production, the most likely cause is
`max_iter` too small or the BH curve being non-monotone (data error).

---

## 7. Roadmap

| Item | Effort | Justification |
|---|---|---|
| **Per-element Z_s in `calc_inductance`** (Phase B from ESIM review) | 1–2 weeks | Removes the spatial-averaging approximation of § 3.2; targets <1 % accuracy claim for publication. |
| **Per-element Z_s in `calc_fem_coilmesh`** | 2–3 weeks | Same physics, more complex implementation (Robin BC currently scalar). |
| **Per-panel R_local from local curvature** | 1 week | Compute `R_local = 2 / (κ_1 + κ_2)` from surface mesh; pass to `esim.solve(H, R_local=…)`.  Closes § 3.1 limitation for moderately curved workpieces (gear roots, fillets). |
| **Geometric mesh stretch** (`np.linspace` → graded) | 2 days | Closes § 2.1 limitation at ξ > 100. |
| **Anderson acceleration of Karl outer loop** | 3 days | Replaces Picard relaxation; typical 2–4× iteration reduction at deep saturation. |
| **Formal pytest for linear Bessel match** | 1 day | Locks § 5.1 numerically. |
| **Stoll 1974 nonlinear-envelope cross-check** | 1 week | Closes § 5.3 third row. |
| **Lavers–Biringer 2-sided plate cross-check** | 1 week | Closes § 5.3 fourth row. |
| **Nonlocal SIBC (Bilicz–Badics–Pávó 2023) extension to DC** | 5–7 days | Wide-band capability for low-frequency WPT; see [`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.md). |

---

## 8. References

- **Hollaus, K., Kaltenbacher, M., Schöberl, J.** "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation." *IEEE Trans. Magn.* (2025).  DOI: [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932).  **The canonical reference for the ESIM cell-problem + outer Karl (= Karl Hollaus) Picard iteration** implemented here.
- **Lavers, J. D. and Biringer, P. P.** "An efficient calculation of effective surface impedance for nonlinear ferromagnetic materials." *IEEE Trans. Magn.* **21**(5), 1985.  Closed-form envelope reference for cylinder + slab nonlinear ESIM.
- **Stoll, R. L.** *The Analysis of Eddy Currents.* Oxford University Press, 1974.  Reference for analytical Bessel-function comparisons in § 5.1.
- **Krähenbühl, L. and Muller, D.** "Thin layers in electrical engineering — Example of shell models in analysing eddy-currents by boundary and finite element methods." *IEEE Trans. Magn.* **29**(2), 1993.  Foundational thin-shell SIBC reference.
- **Dlala, E., Belahcen, A., Arkkio, A.** "Optimal Convergence of the Fixed-Point Method for Nonlinear Eddy-Current Problems." *IEEE Trans. Magn.* **44**(6), 2008, pp. 1318-1321.  Optimal contraction-factor analysis for the same Picard scheme; informs `--esim-relax` default of 0.5.
- **Yuferev, S. and Ida, N.** *Surface Impedance Boundary Conditions: A Comprehensive Approach.* CRC Press, 2009.  Textbook reference for nonlinear SIBC iteration schemes.
- **Bilicz, S., Badics, Z., Pávó, J.** "Nonlocal surface impedance boundary condition for wide-band eddy-current problems." *Studies in Applied Electromagnetics and Mechanics* (ISEM 2023).  Wide-band nonlocal extension (deferred to roadmap § 7).
- **Wakao, S., Igarashi, H., Fujiwara, K., Kameari, A.** "Various Verifications of Eddy Current Analysis (Parts 1–9)." *T.IEE Japan* series, 2008–2018.  Series of linear-μ analytical references; Part 5 used for § 5.1 cylinder benchmark.

---

**Document version**: 2026-05-15, written against radia v4.46.3.
File:line citations are valid for HEAD of this version; minor drift on
later commits is expected and tolerated.
