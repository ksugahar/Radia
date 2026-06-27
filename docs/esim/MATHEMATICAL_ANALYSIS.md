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
- [`docs/esim/SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md) — why scalar BIE + curved Tri6 + per-element ESIM is the right combination (vs vector BEM-A / FEM-Kelvin / FEM-coilmesh), with the error-order match argument and the per-iteration Karl cost comparison.
- [`docs/esim/R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) — focused diagnosis of why PEEC and BEM-A produce different coil R values.
- [`docs/research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.ipynb`](../research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.ipynb) — research-WIP notes.
- [`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb) — wide-band nonlocal extension (deferred).

---

## 0. Nomenclature

For self-contained reading.  All quantities are SI; phasor convention
`e^{+jωt}` (engineering / Hollaus 2025 / NGSolve convention).

### Fields and potentials

| Symbol | Meaning | Units |
|---|---|---|
| **H** | Magnetic field intensity (vector) | A/m |
| **B** | Magnetic flux density (vector) | T |
| **E** | Electric field (vector) | V/m |
| **J** | Volume current density (vector) | A/m² |
| **J_s** | Surface current density on workpiece | A/m |
| **A** | Magnetic vector potential | T·m |
| φ | Scalar magnetic potential (on workpiece surface for BIE) | A |
| **H_t** | Tangential H on workpiece surface | A/m |
| H_0, H_t_rms | Driven surface H magnitude / mesh-RMS value | A/m |
| **A_t** | Tangential A on workpiece surface (= A − (A·n)n) | T·m |

### Material parameters

| Symbol | Meaning | Units |
|---|---|---|
| σ | Electrical conductivity | S/m |
| ρ | Electrical resistivity (= 1/σ) | Ω·m |
| μ_0 | Vacuum permeability (4π·10⁻⁷) | H/m |
| μ_r | Relative permeability (real, possibly H-dependent) | 1 |
| μ | Total permeability (= μ_0 μ_r) | H/m |
| μ' − jμ" | Complex permeability (real + hysteretic loss) | H/m |
| BH curve | Tabulated H[A/m] → B[T] for ferromagnetic materials | (table) |

### Geometric / numerical parameters

| Symbol | Meaning | Units |
|---|---|---|
| R | Cylinder radius (cell-problem half-thickness) | m |
| a | Slab half-thickness | m |
| δ | Linear skin depth `√(2ρ / (ω μ_0 μ_r))` | m |
| ξ | Dimensionless skin parameter `R/δ` or `a/δ` | 1 |
| ω | Angular frequency (= 2πf) | rad/s |
| f | Frequency | Hz |
| n_peri | PEEC perimeter filament count per coil cross-section | 1 |
| n_nodes | Cell-problem radial mesh resolution (default 100) | 1 |

### Solver outputs

| Symbol | Meaning | Units | Source |
|---|---|---|---|
| Z_s | Complex effective surface impedance | Ω | cell solver |
| P' | Active surface power density `½ Re(Z_s)|H_t|²` | W/m² | cell solver |
| Q' | Reactive surface power density `½ Im(Z_s)|H_t|²` | VAR/m² | cell solver |
| P_wp | Workpiece total dissipation `∫_Γ P' dS` | W | outer solver |
| L_coil | Self-inductance of coil (vacuum) | H | outer solver |
| ΔL | Workpiece-induced port inductance change | H | Telegen reciprocity |
| L_total | L_coil + ΔL (with workpiece) | H | outer solver |
| R_coil | AC resistance of coil | Ω | outer solver |
| ΔR | Workpiece-induced port resistance change | Ω | Telegen reciprocity |
| R_total | R_coil + ΔR | Ω | outer solver |
| R_ac/R_dc | AC-to-DC resistance ratio (per Dowell) | 1 | cell solver |

### Iteration parameters

| Symbol | Meaning | Default |
|---|---|---|
| α | Karl-iteration under-relaxation parameter `--esim-relax` | 0.5 |
| tol_Karl | Karl convergence threshold `--esim-tol` | 1e-3 |
| max_iter_Karl | Karl max iterations `--esim-max-iter` | 15 |
| tol_cell | Cell-problem inner Picard tolerance | 1e-6 (hard-coded) |
| max_iter_cell | Cell-problem inner Picard max iter | 50 |
| L (Lipschitz) | Local contraction-factor of `T = E ∘ H_t` | empirical ≈ 1 |

### Acronyms

| | |
|---|---|
| ESIM | Effective Surface Impedance Method |
| SIBC | Surface Impedance Boundary Condition |
| BIE | Boundary Integral Equation |
| BEM | Boundary Element Method |
| PEEC | Partial Element Equivalent Circuit |
| FEM | Finite Element Method |
| HCurl | H(curl) Nédélec edge-element finite-element space |
| MQS | Magneto-Quasi-Static (no displacement current) |
| RWG | Rao–Wilton–Glisson (HDivSurface RT₀ basis) |
| RT₀ | Lowest-order Raviart–Thomas (= RWG) |
| Telegen | Telegen's reciprocity theorem |
| Karl iter. | Hollaus-type outer Picard iteration (see § 6) |

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
is the original infinite-slab/cylinder solver and is the one called
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

For an infinite axially-uniform cylindrical conductor (no `z` or `φ`
dependence) driven by an azimuthal surface field, symmetry requires
`H = H_φ(r) φ̂`.  The cylindrical curl with this ansatz is

$$
\nabla \times \mathbf{H} = \frac{1}{r}\,\partial_r\!\bigl[r\,H_\varphi\bigr]\,\hat{z}
$$

(no `r̂` or `ẑ` components from a purely-`φ̂` field that depends only
on `r`).  Then `J = σ E = (σ/r) ∂_r[r H_φ] ẑ`, and applying
`∇ × (∇ × H) = -∇²H` plus `∇ × J = jωμ H` gives, after one more
cylindrical-curl operation:

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

For the homogeneous Dirichlet problem we test against `v ∈ H¹_0([0,R])`
with `v(R) = 0`, so the boundary term at `r = R` vanishes
identically.

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

*Caveat on the FE↔FD equivalence statement above.*  Strictly, the
equivalence holds when the stiffness term uses midpoint quadrature
on each element (which gives the `r_{i±1/2}` weights) AND the mass
term uses vertex / lumped quadrature (which gives `μ_i r_i`).  Mixing
two quadrature rules in this way is non-standard but well-defined; a
pure-Galerkin P1 implementation with consistent (Gauss) quadrature
gives a structurally identical tridiagonal stencil with slightly
different coefficients (differ by `O(h)`) and converges to the same
continuum limit.

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
[`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb)),
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

All three pass a **single global radius** to a **single ESIM solver
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
              = ∫_Γ |∇_s phi|² dS / area

For per-DOF Karl we need a **physically-correct local gradient norm**
at each vertex.  The current implementation
([`extract_H_t_per_dof_grad`](../../src/radia/bem_sibc_solver.py)
since v4.69.x+) uses the manual triangle-wise P1 gradient:

For each surface triangle with vertices `(p_0, p_1, p_2)`, the
linear interpolant of `phi` has constant in-plane gradient:

    ∇_S phi |_tri = Σ_{j=0,1,2} phi_j · ∇N_j

with the P1 basis-function gradients

    ∇N_0 = (p_2 - p_1) × n̂ / (2 area),
    ∇N_1 = (p_0 - p_2) × n̂ / (2 area),
    ∇N_2 = (p_1 - p_0) × n̂ / (2 area),

(`n̂` = outward unit normal).  Each `∇N_i` is in-plane perpendicular
to the opposite edge with magnitude `1/h_i`, the altitude from
vertex `i`.  Then `|∇phi|²` per triangle is summed area-weighted
into the three incident vertices, giving the per-vertex value
`|H_t|_i²`.

**Bug history (v4.66.0 → v4.69.x)**: an earlier implementation used
a Galerkin localization `|H_t at i|² ≈ |phi_i · (K phi)_i| / M_lump_i`,
treating each diagonal entry of `phi^T K phi` as the per-DOF density.
Although this integrates to the correct total `phi^T K phi`, the
per-DOF distribution **samples the surface Laplacian of `phi`, not
its gradient norm**.  The Laplacian peaks at edge / corner DOFs
where the curvature of `phi` is high, not at the saturation
hot-spots where the field `|H_t|` is large.  Feeding the Laplacian-
sample into the cell solver returned wrong local `Z_s` values, and
inflated the per-element vs scalar `P_wp` gap by ~2-3× in the IH
regime (sign-reversed in some cases).  The fix (commit `630527d4`)
replaced the Galerkin sample with the triangle-gradient formula
above and matches the v4.67.0 fix already applied to the
spatial `q_surf` output.  See
[`docs/esim/CROSS_VALIDATION.md`](CROSS_VALIDATION.md) § 6d for the
gap-reframing analysis.

In the limit of uniform `|∇phi|` over the mesh (e.g. linear-`mu_r`
regime with uniform external `H_t`), the per-DOF formula reduces to
the scalar mesh-RMS to floating-point precision.

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

### 4.2.5 Lorentz reciprocity / reaction integral for ΔL (the φ·B form)

This subsection derives the gauge-invariant formula for the
workpiece-induced port impedance change `ΔZ`, implemented in
`calc_inductance.py` via [`radia.workpiece_surface.delta_L_telegen_phiB`](../../src/radia/workpiece_surface.py).

**Naming note.**  The codebase function name uses "telegen" as lab
shorthand.  The correct literature term for the integral identity
below is the **Lorentz reciprocity theorem** in its
**reaction-integral form** (Rumsey 1954; Harrington 1961, *Time-
Harmonic Electromagnetic Fields*, §3.8).  Telegen's reciprocity is a
distinct, network-theoretic result.  IGTE-paper readers should cite
Harrington / Rumsey.

**Setup.**  Let `Ω_wp` be the workpiece domain with boundary
`Γ = ∂Ω_wp` and outward unit normal `n`.  The coil drives a port
current `I_port`.  In the absence of the workpiece, the coil's
own field has scalar potential `φ_inc` and vector potential `A_inc`
on Γ such that `H_inc = -∇φ_inc` (scalar) and `B_inc = ∇ × A_inc`
(vector) — both well-defined by Biot–Savart from the coil filaments.

When the workpiece is added, an induced surface current `J_s` flows
on Γ and produces a scattered field.  The total port impedance is
`Z_port = Z_vacuum + ΔZ`, and we want a closed-form for `ΔZ`.

**Reaction-integral identity** (Lorentz reciprocity; in the MQS
limit the field interaction reduces to the magnetic-flux coupling
form): the back-reaction at the port is

$$
\Delta V_{\mathrm{port}}\,I_{\mathrm{port}}^{-1}
\;=\; \frac{j\omega}{I_{\mathrm{port}}^2}\,
\int_\Gamma \mathbf{J}_s \cdot \mathbf{A}_{\mathrm{inc}}\,dS
$$

(the volume integral collapses to a surface integral because, in the
SIBC limit, `J_s` is confined to a thin skin layer at Γ).  Defining
the auxiliary complex quantity

$$
\Lambda \;\equiv\; \frac{1}{I_{\mathrm{port}}^2}\,
\int_\Gamma \mathbf{J}_s \cdot \mathbf{A}_{\mathrm{inc}}\,dS
\qquad [\,\mathrm{H}\,]
$$

(units of henries), we have `ΔZ = jω Λ`, hence:

$$
\Delta L \;=\; \mathrm{Re}\,\Lambda, \qquad
\Delta R \;=\; -\omega\,\mathrm{Im}\,\Lambda.
$$

(Verification: `Λ = ΔL_real + j·(ΔL_imag)`; `ΔZ = jω Λ = -ω·Im(Λ) + jω·Re(Λ)`;
matching with `ΔZ = ΔR + jωΔL` gives `ΔR = -ω·Im(Λ)` and
`ΔL = Re(Λ)`.  The code variable named `delta_L_complex` is `Λ`, see
[`workpiece_surface.py:374-379`](../../src/radia/workpiece_surface.py#L374-L379)
and the sign convention at
[`calc_inductance.py:834-835`](../../src/radia/panels/calc_inductance.py#L834-L835).)

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
identity (purely geometric — holds for any divergence-free surface
current `J_s = n × H_t` with `H_t = -∇_s φ`, regardless of frequency
or material):

$$
\int_\Gamma \mathbf{J}_s \cdot \mathbf{A}\,dS
\;=\;\int_\Gamma \varphi\,(\mathbf{n} \cdot \mathrm{curl}\,\mathbf{A})\,dS
\;=\;\int_\Gamma \varphi\,(\mathbf{n} \cdot \mathbf{B}_{\mathrm{inc}})\,dS,
$$

where `φ` is the workpiece-side scalar potential from the SIBC BIE
solve (`H_t = -∇_s φ` on Γ).  Step-by-step (all quantities applied to
the smooth Biot–Savart `A_inc`, NOT to a P1-discretised `A_inc`; the
gauge-failure of the J·A form arises precisely because P1 H_h has
edge-jumps that violate Step 3's smoothness premise):

1. `J_s = n × H_t = -n × ∇_s φ` on the workpiece surface (definition).
2. `∫_Γ J_s · A_t dS = -∫_Γ (n × ∇_s φ) · A_t dS = ∫_Γ ∇_s φ · (n × A_t) dS`
   (using `(a × b) · c = (b × c) · a`).  The normal component
   `(A · n) n` integrates to zero against `n × ∇_s φ` since
   `n × (A·n)n = 0`, so we may replace `A` by `A_t` without loss.
3. `∫_Γ ∇_s φ · (n × A_t) dS = -∫_Γ φ · (∇_s · (n × A_t)) dS`
   (surface integration by parts on the closed C¹ surface Γ; valid
   for the smooth Biot–Savart `A_inc` field used here).
4. `∇_s · (n × A_t) = -n · curl A` (standard surface vector-calculus
   identity on C¹ closed surfaces), so the result follows:
   `∫_Γ φ · (n · curl A) dS = ∫_Γ φ · (n · B) dS`.

The right-hand side uses `B = curl A` directly — gauge-invariant.
This is the formula implemented in
[`workpiece_surface.py:294-379`](../../src/radia/workpiece_surface.py#L294-L379)
as the auxiliary quantity `Λ` introduced earlier in this section:

$$
\boxed{\quad
\Lambda \;=\; \frac{1}{I_{\mathrm{port}}^2}\int_\Gamma \varphi(r)\,
              \bigl(\mathbf{n}(r) \cdot \mathbf{B}_{\mathrm{inc}}(r)\bigr)\,dS,
\quad
\Delta L = \mathrm{Re}\,\Lambda,
\quad
\Delta R = -\omega\,\mathrm{Im}\,\Lambda.
\quad}
$$

(The sign-tracking is implemented at
[`calc_inductance.py:834-835`](../../src/radia/panels/calc_inductance.py#L834-L835):
`delta_L_nH = delta_L_complex.real * 1e9` and
`delta_R_mOhm = -delta_L_complex.imag * omega * 1e3`, where
`delta_L_complex = Λ`.)

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

## 4.4 HCurl Weak Form (FEM-Kelvin and FEM-coilmesh Paths)

This section formally derives the HCurl weak form used by
`calc_fem_kelvin.py` and `calc_fem_coilmesh.py`, from Maxwell's
equations to the Galerkin assembly.  The BIE weak form has already
been described under § 1.5 (the workpiece scalar BIE is dual to the
cell-problem weak form on the surface).

### 4.4.1 Strong form (MQS A formulation with SIBC)

In the MQS limit (no displacement current), the magnetic vector
potential `A` satisfies:

$$
\nabla \times (\nu\,\nabla \times \mathbf{A}) + j\omega\,\sigma\,\mathbf{A}
= \mathbf{J}_{\mathrm{coil}}, \qquad \nu = 1/\mu.
$$

In the air domain `Ω_air` (σ = 0): `∇ × (ν ∇ × A) = J_coil`.
In the workpiece volume `Ω_wp`: replaced by the SIBC on `Γ_wp = ∂Ω_wp`,
which couples `A_t` to the workpiece scalar surface impedance via the
Leontovich relation `n × E = Z_s H_t`.

**Derivation of the Robin coefficient `jω/Z_s`.**  In the MQS Coulomb
gauge, `E = -jω A`.  Substitute into Leontovich:

$$
\mathbf{n} \times \mathbf{E} = Z_s\,\mathbf{H}_t
\;\Longrightarrow\;
-j\omega\,(\mathbf{n} \times \mathbf{A}) = Z_s\,\mathbf{H}_t.
$$

The natural BC for the curl-curl form is `n × (ν ∇×A)`.  Using the
identity `n × (∇×A) = μ_0 H × n × n` (Ampère on the surface, with
`n` the outward normal of `Γ_wp`), and the standard surface identity
`n × A = -n × n × A_t = A_t` (after sign-tracking the orientation),
we get:

$$
\mathbf{H}_t = -\frac{j\omega}{Z_s}\,(\mathbf{n} \times \mathbf{A})
            = +\frac{j\omega}{Z_s}\,\mathbf{A}_t.
$$

Equivalently, in terms of the curl on the conductor side:

$$
\boxed{\quad
\mathbf{n} \times (\nu\,\nabla \times \mathbf{A}) =
\frac{j\omega}{Z_s}\,\mathbf{A}_t \quad\text{on } \Gamma_{\mathrm{wp}}.
\quad}
$$

This is the Robin coefficient that appears in the weak form below.
The factor `jω/Z_s` is implemented at [`calc_fem_kelvin.py:559`](../../src/radia/panels/calc_fem_kelvin.py#L559)
(scalar case) and [`calc_fem_kelvin.py:590`](../../src/radia/panels/calc_fem_kelvin.py#L590)
(per-DOF case).

### 4.4.2 Weak form (Galerkin H(curl))

Multiply by a test function `v ∈ H(curl; Ω)`, integrate over Ω, use
integration by parts:

$$
\int_\Omega \nu\,(\nabla \times \mathbf{A}) \cdot (\nabla \times \mathbf{v})\,d\Omega
\;+\; \int_{\Gamma_{\mathrm{wp}}} \frac{j\omega}{Z_s}\,\mathbf{A}_t \cdot \mathbf{v}_t\,dS
\;=\; \int_\Omega \mathbf{J}_{\mathrm{coil}} \cdot \mathbf{v}\,d\Omega.
$$

The boundary term on the outer Dirichlet `Γ_∞` is dropped by choice
of trial space (`A_t = 0` on `Γ_∞` or its Kelvin-transformed image).

This is exactly the bilinear form assembled at
[`calc_fem_kelvin.py:565`](../../src/radia/panels/calc_fem_kelvin.py#L565):

```python
a_bf = BilinearForm(fes)
a_bf += nu_cf * curl(u) * curl(v) * dx                    # curl-curl
a_bf += (1j*omega / Z_s_cf) * u.Trace() * v.Trace() * ds(sibc_bnd)  # Robin
```

with the Kelvin transformation realised through `nu_cf` taking the
scaled value `NU_0 (r'/R)²` inside the Kelvin-inverted sphere
([`calc_fem_kelvin.py:240-252`](../../src/radia/panels/calc_fem_kelvin.py#L240-L252)).

### 4.4.3 Per-DOF Robin coefficient

In the per-DOF Karl variant, `Z_s` becomes a function of position on
`Γ_wp` — represented as an H1(`Γ_wp`)-valued GridFunction:

```python
gf_Zs = GridFunction(fes_Zs_on_sibc_bnd)
# ... per-DOF ESIM call to populate gf_Zs ...
robin_coeff = 1j * omega / gf_Zs    # CoefficientFunction broadcasting
a_bf += robin_coeff * u.Trace() * v.Trace() * ds(sibc_bnd)
```

NGSolve broadcasts the per-DOF coefficient onto each quadrature point
at assembly time.  Convergence of the per-DOF Karl loop is on
`max_i |Z_s_new[i] − Z_s_old[i]| / |Z_s_old[i]|` (see [`IMPLEMENTATION.md`](IMPLEMENTATION.md) § 3.4).

### 4.4.4 Compound A-V form (calc_fem_coilmesh)

For volumetric coil, the bilinear form gains a coil-side conductivity
term, *active only on the coil sub-domain* `Ω_coil ⊂ Ω`
(i.e. `σ(r) = σ_coil` on `Ω_coil`, `σ(r) = 0` elsewhere):

$$
\int_{\Omega_{\mathrm{coil}}} j\omega\,\sigma_{\mathrm{coil}}\,
(\mathbf{A} + \nabla\varphi) \cdot (\mathbf{v} + \nabla\psi)\,d\Omega.
$$

The trial / test space pair `(φ, ψ) ∈ H¹(Ω_coil)` carries source-port
Dirichlet BCs `φ = 1` (source face) / `φ = 0` (sink face).  These
Dirichlet conditions also serve as the gauge fix for the otherwise
non-unique scalar component (a pure constant added to `φ` is killed
by the two-face Dirichlet pinning).  The compound FES is
`HCurl(A) × H¹(φ_coil)` and the linear solve is a single
block-direct factorisation (see
[`calc_fem_coilmesh.py:216-222`](../../src/radia/panels/calc_fem_coilmesh.py#L216-L222)).

---

## 4.5 Complex-μ Permeability (Hysteretic Loss)

The ESIM cell solver supports two BH-relation modes (the design is in
[`ComplexPermeabilityInterpolator`](../../src/radia/esim_cell_problem.py#L49)
and `ESIMCellProblemSolver.__init__`, lines 866-875):

### 4.5.1 Real BH curve (saturation only)

A two-column table `H[A/m] → B[T]` is loaded via
[`em_material.load_bh_file`](../../src/radia/em_material.py#L81).
The cell solver interpolates `μ(|H|) = B(|H|) / |H|` (cubic spline if
≥4 points, linear otherwise — see [`esim_cell_problem.py:75`](../../src/radia/esim_cell_problem.py#L75)).
This branch models **saturation only**; the loss is purely Joule
(`P' = ½ Re(Z_s) |H_t|²`).

Use case: structural / electrical steel, lamination cores, IH workpieces.

### 4.5.2 Complex permeability `μ' − jμ"` (hysteretic + grain eddy loss)

For lossy ferrite cores, the constitutive relation is

$$
B = \mu(H) \cdot H, \qquad \mu = \mu' - j\mu",
$$

where the imaginary part `μ"` captures both hysteretic loss and
grain-level eddy-current loss inside the conductor.  The cell PDE is
unchanged in form but `μ` is now complex:

$$
\frac{\rho}{r}\,\partial_r[r\,\partial_r H] + j\omega\,\mu(|H|)\,H = 0,
\qquad \mu = \mu'(|H|) - j\mu"(|H|).
$$

The dissipation formula gains an extra "magnetic" term:

$$
P' = \tfrac{1}{2}\,\rho\,\int_0^R |\partial_r H|^2 r\,dr
   + \underbrace{\tfrac{1}{2}\omega\,\int_0^R \mu"(|H|)\,|H|^2 r\,dr}_{P_{\text{magnetic}}}.
$$

(See [`esim_cell_problem.py:1100-1158`](../../src/radia/esim_cell_problem.py#L1100-L1158)
where `P_magnetic` is computed separately and returned in the result dict.)

**Input format** (cell-solver `complex_mu` argument):

1. Constant tuple `(μ'_r, μ"_r)` — constant complex permeability.
2. Table `[[H, μ'_r, μ"_r], ...]` — H-dependent.
3. Dict `{'H': [...], 'mu_prime': [...], 'mu_double_prime': [...]}` —
   numpy-friendly form.

Use case: Mn-Zn / Ni-Zn ferrites, ferromagnetic insulators where
the bulk loss is non-negligible vs the eddy loss.

**Not yet exposed via CLI**: the production CLIs accept only `--bh-file`
(real BH curve) at present.  Complex-μ requires direct Python use of
the cell solver.  This is a roadmap item.

---

## 4.6 R_ac/R_dc Output of the Finite-Slab / Cylinder Solver

`ESIMFiniteSlabSolver.solve(H0)` also returns an **R_ac/R_dc ratio**:

$$
\frac{R_{\mathrm{ac}}}{R_{\mathrm{dc}}} =
\frac{P_{\mathrm{ohmic, ac}}}{P_{\mathrm{ohmic, dc}}}
$$

where `P_ohmic,ac` is the actual eddy-current dissipation at the
given `H_t` and `P_ohmic,dc` is the equivalent DC dissipation
`I² R_dc = I² ρ L / A`.  This ratio is used inside `peec_bundle`
to scale per-filament Dowell resistances (when `Zs_fil` is provided
to `solve_loop_bundle`).

Implementation:
[`esim_cell_problem.py:706-744`](../../src/radia/esim_cell_problem.py#L706-L744)
(numerical integration of `|∂_r H|² r dr`).

Closed-form Dowell reference for circular cross-section is at
[`esim_cell_problem.py:787-813`](../../src/radia/esim_cell_problem.py#L787-L813):

$$
\frac{R_{\mathrm{ac}}}{R_{\mathrm{dc}}}\bigg|_{\mathrm{Dowell}} =
\xi \cdot \frac{\sinh(2\xi) + \sin(2\xi)}{\cosh(2\xi) - \cos(2\xi)},
\qquad \xi = a/\delta.
$$

The numerical R_ac/R_dc from the BVP solver matches the Dowell
closed form to ~10⁻³ in the regime `1 ≤ ξ ≤ 10`; outside that
range either the Dowell limit (`ξ → 0`: R_ratio → 1; `ξ → ∞`:
R_ratio → ξ) or the BVP solver alone is preferred.

---

## 4.7 Limitations and Regime of Validity

The Radia ESIM implementation is validated and expected to perform
well **inside** a specific physics regime.  Outside it, the user
should not rely on the results.

### 4.7.1 Frequency regime

| Frequency | Validity | Reason |
|---|---|---|
| DC – 1 Hz | **Outside** (use FEM A-V with `--impedance-model sibc` and `mu_r=1`) | δ ≫ workpiece size; no skin layer to apply SIBC to |
| 10 Hz – 10 kHz | OK if workpiece is conductive (δ ~ workpiece size); marginal otherwise | Transition regime — verify with `calc_fem_coilmesh` |
| 10 kHz – 1 MHz | **In regime** | Production IH frequencies; thin-skin approx valid |
| 1 MHz – 100 MHz | OK if δ ≪ workpiece size, but verify | Displacement current still negligible for most engineering scales |
| > 100 MHz | **Outside** | Wave effects in the air (radiation, ε₀ matters); use full-wave EFIE/MFIE (`ngsolve.bem`) |

The MQS approximation (no displacement current) is hard-coded
throughout Radia (see CLAUDE.md "Green's Function: Laplace Kernel
Only").  For high-frequency wave problems, use ngsolve.bem with the
Helmholtz kernel.

### 4.7.2 Material regime

| Material class | Validity | Caveat |
|---|---|---|
| Non-magnetic conductor (Cu, Al, brass, μ_r = 1) | **In regime** | Use `--impedance-model sibc`; ESIM is over-engineered here |
| Linear ferromagnetic (μ_r constant > 1) | **In regime** | Use `--impedance-model sibc --mu-r 100` |
| Soft-magnetic with single-valued BH curve (silicon steel, electrical steel, soft iron) | **In regime** | Use `--impedance-model esim --bh-file ...` |
| Hard-magnetic / hysteretic materials (rare-earth magnets) | **Outside** | ESIM assumes single-valued μ(|H|); does not model branching BH loops |
| Lossy ferrite (Mn-Zn, Ni-Zn) | **In regime** but CLI not yet exposed | Use the Python cell-solver directly with `complex_mu` argument |
| Anisotropic materials (laminated cores) | **Outside** | ESIM cell solver is scalar / isotropic; for laminated cores use the homogenisation method of Hollaus/Hannukainen/Hiptmair (separate solver) |

### 4.7.3 Geometric regime

| Geometric feature | Validity | Caveat |
|---|---|---|
| Idealised cylindrical / planar workpieces | **In regime** | Single global radius / thickness OK |
| Complex shapes with locally varying curvature | OK with caveat | Single scalar Z_s assumes a single global radius; per-panel R is roadmap (see § 3.4) |
| Sharp corners (gear roots, fillets) | **Marginal** | Cell solver assumes locally-1D field; corners violate this. Mesh refinement in the corner region + FEM A-V comparison is recommended |
| Very thin workpieces (`a < δ`) | **Outside** the thin-skin approximation that SIBC requires | Use `--impedance-model sibc` with thin-shell formulation, OR use FEM A-V volumetric workpiece mesh |
| Very thick workpieces (`a > 100 δ`) | **In regime** but cell mesh too coarse | Bump `n_nodes` to 2000+ (see § 5.1 cross-check) |

### 4.7.4 Coupling-mode regime

| Coupling assumption | Validity | When it fails |
|---|---|---|
| Weak coupling (coil current fixed, workpiece reacts) | **In regime** for most IH | Strong-coupling cases (deep saturation + strong workpiece back-reaction onto coil current distribution) need `calc_fem_coilmesh.py` full A-V |
| Linear coil (Dowell skin formula in filament bundle) | **In regime** | Strong coil saturation (NOT a coil mode in Radia — coil materials are linear by design, see CLAUDE.md) |

### 4.7.5 Known accuracy caveats

- Scalar mesh-RMS H_t → single Z_s gives ±5 % accuracy for typical IH
  designs.  For < 5 % accuracy use `--esim-per-panel`.
- Telegen ΔL is gauge-invariant in continuum but ~1 % error from
  centroid quadrature on coarse workpiece meshes.
- The PEEC ↔ BEM-A R_coil discrepancy (1.3 – 3 ×) is **expected** —
  see [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md).
- For publication accuracy, always anchor against
  `calc_fem_coilmesh.py` (volumetric A-V); it is the highest-fidelity
  reference inside Radia.

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
[`docs/ih_esim_benchmark/`](../ih_esim_benchmark/) ships
`benchmark.py` + `analytical_bessel_baseline.py` + `results.json` +
`benchmark_plot.pdf`.  This is the runnable counterpart to § 5.1 (linear
Bessel parity).  Extending it to Stoll-envelope and Lavers–Biringer
nonlinear cases is roadmap § 7.

### 5.3 Cross-validation against analytical references

Closed-form references for the SIBC + ESIM combination:

| Geometry | Reference | Status |
|---|---|---|
| Cylinder + linear μ | Wakao–Igarashi–Fujiwara Part 5 (Bessel) | **VERIFIED** (matches to ~10⁻⁴; benchmark at [`docs/ih_esim_benchmark/`](../ih_esim_benchmark/)) |
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

### 6.4 Per-DOF Lipschitz and the noise-floor failure mode

When the outer loop is **per-element** (`--esim-per-panel`), the
fixed-point map becomes vector-valued:

$$
\mathbf{G} : \mathbb{C}^{n_{\mathrm{DOF}}} \to \mathbb{C}^{n_{\mathrm{DOF}}},
\qquad
\mathbf{G}(\mathbf{Z}_s)[i] = \mathcal{E}\bigl(|H_t|_i(\mathbf{Z}_s)\bigr),
$$

where `|H_t|_i(Z_s)` is the per-DOF extraction (§3.4) applied to
`phi_vec(Z_s)`.  The damped-Picard contraction condition becomes a
**componentwise** statement: at iterate `k`, the LOCAL Lipschitz
constant at DOF `i` is

$$
L_i^{(k)} \;\approx\; \left| \frac{\partial \mathcal{E}}{\partial |H_t|}
\right|_{|H_t|_i^{(k)}}
\cdot
\left\| \frac{\partial |H_t|_i}{\partial \mathbf{Z}_s} \right\|_{\!\infty}.
$$

The first factor is the slope of the cell-solver `Z_s` vs `|H_t|`
envelope; it is **bounded** away from 1 except at the BH knee, where
`mu_r` drops sharply.  The second factor is the BIE sensitivity at
DOF `i`; for hot-spot DOFs (where `|H_t|_i` is near `|H_t|_{max}`)
this can be `2`–`3` in our IH benchmarks.

**Consequence**: the damped-Picard update with a SINGLE `alpha` cannot
simultaneously satisfy `alpha · L_i < 1` at every DOF.  Picking
`alpha = 0.3` puts most DOFs deep in contraction but lets the
hot-spot DOFs sit at `alpha · L_i ≈ 0.6–0.9`.  In the current dense
IGTE sweep, Anderson acceleration with memory 5 and relaxation 0.5
converges 49 of 54 per-DOF cases.  Five high-current cases reach the
30-iteration cap and are excluded from quantitative headline maxima.

### 6.5 Anderson acceleration

The natural remedy is **Anderson-type-II acceleration** (Anderson
1965; Walker-Ni 2011), which augments damped Picard with a
least-squares combination of the last `m` iterates to suppress the
slowly-decaying modes that constant-alpha damping cannot reach.

Setup.  Let `f_k = G(x_k) - x_k` be the residual at iteration `k`,
and `m_k = min(m, k)` the available memory depth.  Define

$$
\Delta X_k = [\, x_{k-m_k+1} - x_{k-m_k},\; \ldots,\; x_k - x_{k-1}\,], \quad
\Delta F_k = [\, f_{k-m_k+1} - f_{k-m_k},\; \ldots,\; f_k - f_{k-1}\,].
$$

Both have shape `(n_DOF, m_k)`.  Solve the unconstrained real
least-squares

$$
\gamma_k \;=\; \mathrm{arg\,min}_{\gamma \in \mathbb{R}^{m_k}}
   \;\bigl\| \Delta F_k \,\gamma - f_k \bigr\|_2,
$$

with the complex columns of `ΔF_k` and `f_k` stacked as
`[\Re\,;\, \Im]` to keep `gamma` real (the convention for complex
fixed-point iterations, Walker-Ni 2011 § 2.2 remark).  The
damped-Anderson update is

$$
x_{k+1} \;=\; x_k \;+\; \alpha\, f_k
   \;-\; \bigl(\Delta X_k + \alpha\, \Delta F_k\bigr)\,\gamma_k.
$$

`m = 0` recovers plain damped Picard.  `alpha = 1` (no damping) gives
undamped Anderson.

**Memory depth `m`.**  Empirically `m = 3–5` is sufficient for the IH
benchmark; `m = 1` is too short to capture the per-DOF "ringing"
mode; `m > 7` increases LSQ cost without further iteration savings
on this problem.

**Damping `alpha`.**  Keep `alpha = 0.5–0.7` even with Anderson on
— the first few iterations have `m_k < m` so the early steps still
need damping to avoid overshoot.

**Why Anderson succeeds where damped Picard fails.**  Damped Picard
applies the SAME contraction factor `alpha` to every mode of the
linearised map `G'(x*)`; Anderson adapts its effective factor
per-mode via the least-squares combination.  For the IH per-element
problem the slow modes are localised at hot-spot DOFs; Anderson
captures these in the `ΔF_k` history and zeros them out in the LSQ.

**Implementation.**  The :class:`AndersonAccelerator` in
[`src/radia/esim_anderson.py`](../../src/radia/esim_anderson.py) is
a drop-in replacement for the `alpha * G(x) + (1-alpha) * x` update
in the Karl loop.  Wired into `calc_inductance.py` behind the CLI
flag `--esim-anderson-m N` (default `N = 0` = plain damped Picard).

### 6.6 Safeguarding (Walker-Ni 2011 § 3.4)

Naive Anderson on per-element ESIM hits a known pathology: the LSQ
matrix `dF_k` can become ill-conditioned at intermediate plateaus
(early `dZ_max` valleys that are NOT true fixed points), causing
the next step to over-extrapolate and the trajectory to drift to a
different basin of attraction.  Two safeguards are implemented in
:class:`AndersonAccelerator`:

**Step-clipping.**  After computing the Anderson correction
`corr = -(dX + alpha * dF) @ gamma`, if `||corr|| > step_clip *
||alpha * f_k||` (default `step_clip = 2.0`), scale `corr` down so
the inequality holds.  Prevents pathological over-extrapolation
when `dF_k` is rank-deficient.  Reports `n_clips` per run.

**Relative-residual restart.**  Track the previous iteration's
relative inf-norm residual `rel_resid = max_i |f_i| / max_i |x_i|`.
If the current iteration's `rel_resid` exceeds the previous by
`restart_growth` (default `2.0`), the most recent Anderson step
diverged — clear the history (keep the current iterate) so the next
step is plain damped Picard, then Anderson rebuilds from scratch.
Reports `n_restarts` per run.

**Why "vs previous iter" not "vs min-ever".**  Per-element ESIM has
an intrinsic per-DOF noise floor (§ 6.4): the dZ_max never drops
below ~5e-3 on the IH benchmark.  Comparing to min-ever triggers a
restart on every iteration after the first deep valley, killing the
Anderson history and effectively disabling acceleration (verified
empirically: 18/30 iter restarts with min-ever vs 2/30 with prev-iter).

**Validation on IH per-element benchmark** (steel cylinder, dense
108-case sweep): the production setting is `--esim-anderson-m 5
--esim-relax 0.5 --esim-max-iter 30`.  At `I_port = 100 A` and
`f = 50 kHz`, the per-DOF model converges in 7 iterations and gives
`P_wp = 18.75 W`; the corresponding scalar uniform model gives
`P_wp = 30.51 W`.

Safeguarded Anderson with prev-iter restart criterion brings the
per-DOF `dZ_max` from the ~0.1 noise floor down to **5e-3**, within
5× of the formal `tol = 1e-3` cutoff — and crucially the trajectory
is **monotone** in `dZ_max` at the end (last 5 iter values
strictly decreasing).  Extending `--esim-max-iter` to 50 is
expected to cross `1e-3`.

For typical IH benchmarks set `--esim-anderson-m 5` with default
safeguards; expect convergence behavior matching the bottom row of
the table above.

---

## 6.7 Axisymmetric Volumetric FEM Truth Reference

To evaluate ESIM's *absolute* accuracy (as opposed to scalar-vs-
per-element internal disagreement) we need a reference solver that
resolves the volumetric eddy current inside the workpiece — not a
SIBC Robin BC.  The 2D axisymmetric A_phi formulation
([`src/radia/panels/calc_axisym_volumetric.py`](../../src/radia/panels/calc_axisym_volumetric.py))
provides this.

**Formulation.**  Time-harmonic Maxwell in axisymmetric `(r, z)`,
single component `A_phi(r, z)` complex-valued.  The system is

    (1/mu) [ d_r(r A) d_r(r v) / r + r d_z A d_z v ] + j w sigma r A v
        = J_phi_src r v

(NGSolve weak form, complex H1, `dirichlet="axis|outer|top|bot"`).
The workpiece-side material is `sigma > 0` (eddy current resolved
volumetrically) and `mu_r` from a BH curve (linear or nonlinear).
The coil is a current-density source `J_phi_src = I/A_coil` in a
small ring region.  Workpiece power is

    P_wp = (1/2) Re int_workpiece sigma (omega A_phi)^2 dV_axisym

where `dV_axisym = 2 pi r dr dz`.

**Linear-mu Bessel validation** (see CROSS_VALIDATION.md § 5b for the
detailed table).  Long cylinder R=5mm H=200mm, mu_r=100, sigma=2e6,
f=50kHz:

| Method           | P_wp   |
|------------------|--------|
| FEM volumetric   | 0.182 W |
| Bessel I_0/I_1   | 0.193 W |
| Agreement        | -5.7 % (end-effect) |

**IGTE-geometry linear-mu validation** (R_wp=5mm, H_wp=10mm,
R_coil=20mm, mu_r=100, sigma=2e6, f=50kHz, I=100A):

| Method           | P_wp    | |H_t|_side  |
|------------------|---------|------------|
| FEM volumetric   | 1.01 W  | 1421 A/m   |
| Bessel-from-FEM  | 0.997 W | (input)    |
| Agreement        | +1 %    |            |

For the same geometry under nonlinear BH (the IGTE benchmark
operating point), the comparison requires a nonlinear FEM
extension (Picard on local mu(|B|)) — open item.

**Why this matters.**  Once nonlinear FEM is operational, the
scalar-vs-per-element ESIM gap reported in
[`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) § 6d can be split into:

  - which side of the gap is closer to the FEM truth, AND
  - what the absolute error of EACH method is.

This is the only way to determine whether scalar SIBC's over-prediction
relative to the per-DOF sweep is the dominant absolute error, or whether
both reduced models deviate comparably from the volumetric solution.

---

## 7. Roadmap

| Item | Effort | Justification |
|---|---|---|
| **Per-element Z_s in `calc_inductance`** (Phase B from ESIM review) | 1–2 weeks | Removes the spatial-averaging approximation of § 3.2; targets <1 % accuracy claim for publication. |
| **Per-element Z_s in `calc_fem_coilmesh`** | 2–3 weeks | Same physics, more complex implementation (Robin BC currently scalar). |
| **Per-panel R_local from local curvature** | 1 week | Compute `R_local = 2 / (κ_1 + κ_2)` from surface mesh; pass to `esim.solve(H, R_local=…)`.  Closes § 3.1 limitation for moderately curved workpieces (gear roots, fillets). |
| **Geometric mesh stretch** (`np.linspace` → graded) | 2 days | Closes § 2.1 limitation at ξ > 100. |
| **Anderson acceleration of Karl outer loop** | ✅ done (radia ≥ 4.55.9) | Implemented in [`src/radia/esim_anderson.py`](../../src/radia/esim_anderson.py); enable with `--esim-anderson-m 5`.  Closes the per-DOF noise floor on the IH benchmark (§ 6.5). |
| **Formal pytest for linear Bessel match** | 1 day | Locks § 5.1 numerically. |
| **Stoll 1974 nonlinear-envelope cross-check** | 1 week | Closes § 5.3 third row. |
| **Lavers–Biringer 2-sided plate cross-check** | 1 week | Closes § 5.3 fourth row. |
| **Nonlocal SIBC (Bilicz–Badics–Pávó 2023) extension to DC** | 5–7 days | Wide-band capability for low-frequency WPT; see [`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb). |

---

## 8. References

- **Hollaus, K., Kaltenbacher, M., Schöberl, J.** "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation." *IEEE Trans. Magn.* (2025).  DOI: [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932).  **The canonical reference for the ESIM cell-problem + outer Karl (= Karl Hollaus) Picard iteration** implemented here.
- **Lavers, J. D. and Biringer, P. P.** "An efficient calculation of effective surface impedance for nonlinear ferromagnetic materials." *IEEE Trans. Magn.* **21**(5), 1985.  Closed-form envelope reference for cylinder + slab nonlinear ESIM.
- **Stoll, R. L.** *The Analysis of Eddy Currents.* Oxford University Press, 1974.  Reference for analytical Bessel-function comparisons in § 5.1.
- **Krähenbühl, L. and Muller, D.** "Thin layers in electrical engineering — Example of shell models in analysing eddy-currents by boundary and finite element methods." *IEEE Trans. Magn.* **29**(2), 1993.  Foundational thin-shell SIBC reference.
- **Dlala, E., Belahcen, A., Arkkio, A.** "Optimal Convergence of the Fixed-Point Method for Nonlinear Eddy-Current Problems." *IEEE Trans. Magn.* **44**(6), 2008, pp. 1318-1321.  Optimal contraction-factor analysis for the same Picard scheme; informs `--esim-relax` default of 0.5.
- **Anderson, D. G.** "Iterative Procedures for Nonlinear Integral Equations." *J. ACM* **12**(4), 1965, pp. 547-560.  Original Anderson-acceleration paper.
- **Walker, H. F. and Ni, P.** "Anderson Acceleration for Fixed-Point Iterations." *SIAM J. Numer. Anal.* **49**(4), 2011, pp. 1715-1735.  The modern derivation (Type II, ΔX-ΔF formulation) used in [`src/radia/esim_anderson.py`](../../src/radia/esim_anderson.py).
- **Yuferev, S. and Ida, N.** *Surface Impedance Boundary Conditions: A Comprehensive Approach.* CRC Press, 2009.  Textbook reference for nonlinear SIBC iteration schemes.
- **Bilicz, S., Badics, Z., Pávó, J.** "Nonlocal surface impedance boundary condition for wide-band eddy-current problems." *Studies in Applied Electromagnetics and Mechanics* (ISEM 2023).  Wide-band nonlocal extension (deferred to roadmap § 7).
- **Wakao, S., Igarashi, H., Fujiwara, K., Kameari, A.** "Various Verifications of Eddy Current Analysis (Parts 1–9)." *T.IEE Japan* series, 2008–2018.  Series of linear-μ analytical references; Part 5 used for § 5.1 cylinder benchmark.

---

**Document version**: 2026-05-15, written against radia v4.46.3.
File:line citations are valid for HEAD of this version; minor drift on
later commits is expected and tolerated.
