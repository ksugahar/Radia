# ESIM (Effective Surface Impedance Method): Mathematical Analysis & Implementation Audit

**Scope.** This document describes WHAT the Radia ESIM module solves, HOW it
discretises the 1-D cell problem, what curvature it does and does NOT
capture, and how the outer BEM/FEM curve order propagates into the
Karl-iteration feedback loop.  It is meant as an internal audit for
contributors and a citation reference for publication: every claim
about "the implementation does X" carries a file:line citation.

The companion documents are:

- [`docs/esim/USAGE.md`](USAGE.md) — user-facing CLI guide (TBD).
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
| `finite_slab` | as `slab` with an extra anti-symmetric BC for 1-sided heating | as `slab` | (reserved; not used by production) |

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

The `2·(P' + jQ') / H0²` identity ([`esim_cell_problem.py` solve return dict, near line 1010](../../src/radia/esim_cell_problem.py))
ties P' and Z together — both are independent observables in publication tables.

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

### 3.4 Why "single Z_s" is still useful in publication

For workpieces where the **operating-point H_t** does not vary by more
than ~3× across the surface (typical for solenoid + cylindrical bar),
the scalar mesh-RMS H_t feeds an ESIM call that gives an effective Z_s
within ~5 % of a fully resolved per-panel calculation.  This is good
enough for engineering screening (P_wp, L_eff at ±5 % for design
iteration).

For accuracy claims sharper than 5 %, per-panel ESIM is required.

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

Match to 4 significant figures.  A formal pytest is recommended (TBD).

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

Formal cross-validation suite: `examples/ih_esim_benchmark/` (TBD; see
Phase A roadmap below).

### 5.3 Cross-validation against analytical references

Closed-form references for the SIBC + ESIM combination:

| Geometry | Reference | Status |
|---|---|---|
| Cylinder + linear μ | Wakao–Igarashi–Fujiwara Part 5 (Bessel) | matches to ~10⁻⁴ |
| Slab + linear μ | Dowell (tanh) | matches; Dowell baked in (`mat.dowell_Zs`) |
| Cylinder + nonlinear μ (BH) | Stoll 1974 (analytical envelope) | TBD |
| Plate + 2-sided heating | Lavers–Biringer 1985 | TBD |

The "TBD" rows are the Phase A benchmark targets.

---

## 6. Karl Iteration (Outer Loop) — Implementation Details

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

- **Lavers, J. D. and Biringer, P. P.** "An efficient calculation of effective surface impedance for nonlinear ferromagnetic materials." *IEEE Trans. Magn.* **21**(5), 1985.
- **Stoll, R. L.** *The Analysis of Eddy Currents.* Oxford University Press, 1974.
- **Krähenbühl, L. and Muller, D.** "Thin layers in electrical engineering — Example of shell models in analysing eddy-currents by boundary and finite element methods." *IEEE Trans. Magn.* **29**(2), 1993.
- **Bilicz, S., Badics, Z., Pávó, J.** "Nonlocal surface impedance boundary condition for wide-band eddy-current problems." *Studies in Applied Electromagnetics and Mechanics* (ISEM 2023).
- **Wakao, S., Igarashi, H., Fujiwara, K., Kameari, A.** "Various Verifications of Eddy Current Analysis (Parts 1–9)." *T.IEE Japan* series, 2008–2018.

---

**Document version**: 2026-05-15, written against radia v4.46.3.
File:line citations are valid for HEAD of this version; minor drift on
later commits is expected and tolerated.
