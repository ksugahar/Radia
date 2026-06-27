# Kelvin Transformation for Open Boundary Magnetostatics

The Kelvin transformation maps an unbounded exterior domain to a bounded
computational domain via sphere inversion, enabling standard FEM on a
finite mesh while exactly representing far-field decay. This document
covers the theory, API, and practical workflow as implemented in Radia.

---

## 1. Overview

In magnetostatics, the field extends to infinity. Truncation (Dirichlet/ABC)
introduces artificial error. The Kelvin transformation avoids this by
inverting the exterior `|r| > R` bijectively onto a bounded interior
`|r'| < R`, then modulating the material properties so that the FEM
energy in the transformed domain equals the physical energy exactly.

**Key properties**:
- No truncation error (exact open boundary)
- No PML, no infinite elements, no special basis functions
- Standard FEM discretization on the transformed domain
- Works with any formulation (A, Omega, H) and any element type

**Implementation status**: complete 5-layer Python API on NGSolve, with
Cubit and OCC geometry paths, 3D/axisymmetric/2D support, and PEEC
coil source integration.

---

## 2. Theory

### 2.1 Kelvin Map

3D sphere inversion centered at `c` with radius `R`:

```
phi: r -> r' = c + R^2 / |r - c|^2 * (r - c)
```

Properties:
- **Involutive**: `phi(phi(r)) = r`
- Maps exterior `|r - c| > R` to interior `|r' - c| < R`
- Interface `|r - c| = R` is fixed pointwise
- Physical infinity maps to `r' = c` (the GND point)

### 2.2 Jacobian

Let `rho' = |r' - c|` and `n = (r' - c) / rho'`. The Jacobian of the
map in computational coordinates:

```
J = (rho'/R)^2 * H        where H = I - 2 n n^T  (Householder reflection)
```

The inverse Jacobian (by involution):

```
J^{-1} = (R/rho')^2 * H
```

Key properties: `det(H) = -1` (orientation-reversing), `H^2 = I`,
`|det(J^{-1})| = (R/rho')^6`.

### 2.3 Pullback of Differential Forms

| Form degree | Physical quantity | Pullback formula | Exponent |
|-------------|-------------------|------------------|----------|
| 0-form (scalar) | Potential Omega | `Omega_comp(r') = Omega_phys(phi(r'))` | 0 |
| 1-form (vector) | Vector potential A | `A_comp = (R/rho')^2 H A_phys` | 2 |
| 2-form (pseudo) | Flux density B | `B_comp = -(R/rho')^4 H B_phys` | 4 |

The Householder `H` flips the radial component and preserves tangential
components. For tangential A (e.g. azimuthal from a z-directed B),
`H A = A` and only the scalar factor `(R/rho')^2` survives.

The minus sign in the B pullback arises from `det(H) = -1` via the
Levi-Civita pseudovector identity.

### 2.4 Material Modulation

The bilinear energy equivalence condition
`W_phys = W_comp` yields modified material properties in the
transformed domain. This is the central result of Nagamine CEFC 2026.

#### 3D Spherical (conformal, orientation-reversing)

```
nu' = (rho'/R)^2 * nu_0          [A-formulation, HCurl]
mu' = (R/rho')^2 * mu_0          [Omega / H-formulation, H1]
```

These are pointwise reciprocals: `mu' * nu' = mu_0 * nu_0 = 1`.

**Derivation outline** (Nagamine eqs. 5-9):

1. Pullback of orthonormal 1-form basis: `k*(e^{r'}) = -(R/r)^2 e^r`
2. Hodge operator: `*(k*(B')) = -(R/r)^2 k*(*(B'))`
3. Inner product: `g(k*(w'), k*(w')) = (R/r)^4 g'(w', w')`
4. Bilinear integrand: `nu <dw, dA> = nu <dw', dA'>' (R/r)^8`
5. Volume pullback: `dOmega = k*(-R^6/r'^6 dOmega')` with `sgn(k) = -1`
6. Equating energies: `nu' = nu (r'/R)^2`

#### Axisymmetric (r, z)

Same as 3D spherical, with `rho' = sqrt(r^2 + (z - z_offset)^2)` being
the 3D radial distance from the Kelvin sphere center. The axisymmetric
`2*pi*r` weight is NOT absorbed into the Kelvin factor:

```python
a += nu_cf * grad(u) * grad(v) * r_coord * dx    # correct
a += nu_cf * grad(u) * grad(v)           * dx    # WRONG: off by O(r)
```

#### 2D Cylindrical (non-conformal, anisotropic)

Cylindrical inversion `k(rho, phi, z) = (R^2/rho, phi, z)` yields an
anisotropic material tensor (Nagamine eq. 12):

```
nu' = diag(1, 1, (rho'/R)^4) * nu
mu' = diag(1, 1, (R/rho')^4) * mu
```

Only the axial (z) slot is modulated. Which slot enters the bilinear
form depends on the field components:

| Case | B/H components | Kelvin factor | Frequency |
|------|----------------|---------------|-----------|
| **In-plane** | Hx, Hy or Bx, By | **1 (identity)** | Common |
| **Axial** | Hz or Bz | (rho'/R)^4 for nu, (R/rho')^4 for mu | Rare |

### 2.5 Material Modulation vs Solution Pullback

Two DIFFERENT Kelvin factors arise. Do not confuse them:

| Concept | Factor | Used in |
|---------|--------|---------|
| Material nu (bilinear form) | `(rho'/R)^2 * nu_0` | FEM assembly |
| Material mu (reciprocal) | `(R/rho')^2 * mu_0` | Omega/H-formulation |
| Solution A (1-form pullback) | `(R/rho')^2 * H * A` | Source field evaluation |
| Solution B (2-form pullback) | `-(R/rho')^4 * H * B` | B-field recovery |

### 2.6 Energy Invariance Derivation

The bilinear energy equivalence `W_phys = W_comp` determines `nu_kelvin`.
Substituting `|B_comp|^2 = (R/rho')^8 |B_phys|^2` (from 2-form pullback)
and `dV_phys = (R/rho')^6 dV_comp`:

```
nu_0 (R/rho')^6  =  nu_kelvin (R/rho')^8     on each computational cell
=> nu_kelvin = nu_0 (rho'/R)^2
```

This is the canonical Nagamine CEFC 2026 result.

### 2.7 Curl Verification (uniform B background)

For `B_phys = B_0 z_hat`, `A_phys = (B_0/2)(-y, x, 0)` (azimuthal,
tangential). At `r' = (rho', 0, 0)`, Householder `H` is identity for
tangential A_phys, so `A_comp = (R/rho')^2 (B_0/2)(0, R^2/rho', 0)
= (B_0/2) R^4 / rho'^3 y_hat`. Generalizing to arbitrary r':

```
A_comp(r') = (B_0/2) R^4 / rho'^4 (-r'_y, r'_x, 0)
```

Direct computation of `curl_z(A_comp)`:

```
∂(A_comp_y)/∂(r'_x) at (rho', 0, 0)
   = (B_0/2) R^4 [1/rho'^4 - 4 r'_x^2/rho'^6] = -(3 B_0/2) R^4/rho'^4
-∂(A_comp_x)/∂(r'_y)
   = (B_0/2) R^4 / rho'^4

curl_z(A_comp) = (-3/2 + 1/2) (B_0 R^4 / rho'^4) = -(R/rho')^4 B_0
```

This matches the 2-form pullback `B_comp = -(R/rho')^4 H B_phys` —
internal consistency between the 1-form A pullback and the 2-form
B pullback, confirming the differential-form framework.

### 2.8 Numerical Validation

Toroidal current loop (Nagamine CEFC 2026 Section III):
- Major radius a = 0.1 m, wire radius b = 0.01 m
- Analytical exterior dipole energy: 3.333e-8 J
- FEM on transformed domain with `nu' = (rho'/R)^2 nu_0`: 3.344e-8 J
- Error: **+0.33%**

Reproducible test: `tests/test_pullback_dipole_exterior_energy.py`
(2026-05-04, p=3 HCurl + bonus_intorder=8 → +0.02%).

Stone-bridge test (mu_r=100 iron cylinder + ring coil, 2026-04-17):
- 3D Kelvin vs 2D axisymmetric reference: **1.15% max error** (h=4mm, p=3)
- Iron interior: 0.3-0.5%, near-surface: 1.0-1.15%, exterior: <0.1%

---

## 3. Geometry: Sugahara Two-Sphere Convention

The implementation uses two identical, offset spheres — NOT a concentric
shell:

```
Inner sphere (physical):  center = origin,  radius = R_K
Outer sphere (Kelvin):    center = offset,  radius = R_K  (same R)
```

Periodic boundary conditions identify the two sphere surfaces, creating
a 1:1 node correspondence. The center of the outer sphere (image of
physical infinity) is the GND vertex.

### 3.1 Why Two Separate Spheres?

- Concentric shell (R_inner to R_outer) is NOT Kelvin — the transformed
  material modulation `(rho'/R)^2` would be applied to the wrong domain
- Two offset spheres allow standard FEM meshing of each sphere independently
- Periodic BC handles the surface identification automatically
- The offset must be large enough that the spheres do not overlap

### 3.2 GND Vertex

| Formulation | GND Required? | Purpose |
|-------------|---------------|---------|
| H1 (Omega, scalar) | **Essential** | Uniqueness (Dirichlet at infinity) |
| HCurl (A, vector) | Optional | Gauge regularization suffices |

GND is placed at the outer sphere center (`offset` point), which
corresponds to physical infinity under Kelvin inversion.

---

## 4. API Reference

The Python API is organized in 5 layers:

### Layer 0: Math Primitives (`radia.kelvin_source`)

| Function | Description |
|----------|-------------|
| `kelvin_map_3d(points, center, R)` | 3D sphere inversion `r' = c + R^2(r-c)/|r-c|^2` |
| `is_in_kelvin_exterior_domain(points, center, R)` | Boolean test for Kelvin exterior membership |
| `kelvin_pullback_vector(A, points, center, R)` | 1-form pullback: `(R/rho')^2 H A` |
| `kelvin_pullback_B_pseudovector(B, points, center, R)` | 2-form pullback: `-(R/rho')^4 H B` |
| `kelvin_factor_scalar(points, center, R)` | Scalar magnitude `(R/rho')^2` |

#### Material Modulation CoefficientFunctions

| Function | Formula | Use case |
|----------|---------|----------|
| `kelvin_nu_factor_3d_cf(center, R)` | `(rho'/R)^2` | A-formulation, HCurl, 3D |
| `kelvin_mu_factor_3d_cf(center, R)` | `(R/rho')^2` | Omega/H-form, H1, 3D |
| `kelvin_nu_factor_axisym_cf(z_offset, R)` | `(rho'/R)^2` | Axisymmetric A-form |
| `kelvin_mu_factor_axisym_cf(z_offset, R)` | `(R/rho')^2` | Axisymmetric Omega/H-form |
| `kelvin_factor_2d_inplane_cf()` | **1** | 2D in-plane (common case) |
| `kelvin_nu_factor_2d_axial_cf(offset, R)` | `(rho'/R)^4` | 2D axial nu_zz (rare) |
| `kelvin_mu_factor_2d_axial_cf(offset, R)` | `(R/rho')^4` | 2D axial mu_zz (rare) |
| `build_material_cf(mesh, default, kelvin_cf, overrides)` | Material-indexed CF builder | All formulations |

#### Source Evaluation Helpers

| Function | Description |
|----------|-------------|
| `eval_Omega_physical_from_gf(gf, points, center, R)` | Omega GridFunction at physical points |
| `eval_A_physical_from_gf(gf, points, center, R)` | A GridFunction with 1-form pullback |
| `eval_B_physical_from_gf(gf, points, center, R)` | B = curl(A) with 2-form pullback |
| `eval_H_from_radia_in_kelvin(radia_obj, points, center, R)` | Radia H-field with Kelvin pullback |
| `eval_B_from_radia_in_kelvin(radia_obj, points, center, R)` | Radia B-field with Kelvin pullback |
| `biot_savart_A_at_points(filaments, obs, n_gauss)` | Vector potential from filament bundles |
| `A_s_at_obs_with_kelvin(filaments, obs, center, R, mode)` | Biot-Savart with automatic Kelvin |

### Layer 1: Geometry Helpers (`radia.kelvin_geometry`)

| Function | Description |
|----------|-------------|
| `add_kelvin_exterior_domain(inner_shape, offset, R_K)` | Build two-sphere OCC geometry with periodic identification |

Returns `(geometry, info)` where info contains `center`, `R_K`,
`outer_shape`, `gnd_vertex`, face names for periodic BC.

### Layer 2: Mesh-Aware CFs (`radia.kelvin_material`)

| Function | Description |
|----------|-------------|
| `make_kelvin_nu_cf(mesh, R_K, offset, nu_0, kelvin_mats)` | NGSolve CF for Kelvin-modulated nu |
| `make_kelvin_aware_A_s_cf(mesh, A_factory, R_K, offset)` | Wraps external A_s with Kelvin pullback |

### Layer 3: FEM Drivers (`radia.kelvin_solver`)

| Function | Description |
|----------|-------------|
| `solve_full_A_kelvin()` | Full-A HCurl FEM with volume J source |
| `solve_reduced_A_kelvin()` | Reduced-A with external A_s (PEEC filaments) |

### Layer 4: Validation (`radia.kelvin_validate`)

| Function | Description |
|----------|-------------|
| `compare_against_radia_self_inductance()` | End-to-end inductance comparison vs Radia |

### Panel Integration (`radia.panels`)

| Module | Description |
|--------|-------------|
| `calc_fem_kelvin.py` | 3D FEM-SIBC solver for Radia-NGSolve Cubit panel |
| `add_kelvin.py` | Automate Kelvin exterior for Cubit 3D, OCC 3D, 2D axisym |

---

## 5. Cubit Workflow

### 5.1 Steps

1. **Create inner sphere** (physical domain: coil + air + iron)
   ```
   create sphere radius 0.06
   ```

2. **Create outer sphere** (Kelvin domain, same radius, offset in space)
   ```
   create sphere radius 0.06
   move volume 2 x 0 y 0 z 0.15
   ```

3. **Webcut both** along a matching plane for mesh surface copy
   ```
   webcut volume 1 with plane zplane offset 0
   webcut volume 2 with plane zplane offset 0.15
   ```

4. **Copy mesh surface** (MANDATORY for 1:1 periodic node correspondence)
   ```
   mesh surface <inner_kelvin_face>
   copy mesh surface <inner_face> onto surface <outer_face> \
     source curve ... target curve ... source vertex ... target vertex ...
   ```

5. **Mesh volumes**
   ```
   volume all scheme tetmesh
   mesh volume all
   ```

6. **Assign blocks and sidesets**
   ```
   block 1 add volume 1; block 1 name "air"
   block 2 add volume 2; block 2 name "kelvin"
   ```

7. **Export** — `export netgen` automatically detects `air`+`kelvin`
   blocks and creates periodic identification + `kelvin_int`/`kelvin_ext`
   boundary labels
   ```
   export netgen "model.vol" order 3 overwrite
   ```

### 5.2 Auto-Detection in export netgen

When blocks named `air` and `kelvin` are present, the C++ plugin:
- Labels the shared air|kelvin interface as `kelvin_int`
- Labels the outer kelvin boundary as `kelvin_ext`
- Computes translation offset from vertex centroids
- Creates periodic identification pairs via `ident.Add()`

No manual sideset setup is required for Kelvin boundaries.

### 5.3 Copy Mesh Surface: Why It Is Mandatory

Without `copy mesh surface`, Cubit generates independent tet meshes on
each sphere surface, producing different triangulations. Periodic BC
requires 1:1 vertex correspondence, so mismatched meshes cause a crash
or silent inaccuracy.

### 5.4 PEEC Coil Source (peec block)

When a block named `peec` is present, `export netgen` automatically:
- Exports the peec block's volumes as a STEP file (`model_peec.step`)
- Records the path in the companion JSON (`"peec_step": "model_peec.step"`)
- The peec volume is **not meshed** — it is excluded from the .vol

The STEP file is used by `calc_fem_kelvin.py` to extract PEEC filaments
via `filaments_from_step()`, which are then injected as a Biot-Savart
source in the reduced-A formulation.

#### Cubit .jou Example

```
reset

# Air sphere
create sphere radius 0.06

# Kelvin sphere (same radius, offset)
create sphere radius 0.06
move volume 2 z 0.15

# Coil torus (PEEC — NOT meshed)
create torus major radius 0.03 minor radius 0.003
subtract volume 3 from volume 1 keep

# Mesh air + kelvin only
volume 1 scheme tetmesh
volume 1 size 0.015
mesh volume 1
volume 2 scheme tetmesh
volume 2 size 0.020
mesh volume 2

# Blocks
block 1 add volume 1
block 1 name "air"
block 2 add volume 2
block 2 name "kelvin"
block 3 add volume 3
block 3 name "peec"

# Export → model.vol + model_peec.step
export netgen "model.vol" order 2 overwrite
```

#### Python Solve

```python
from calc_fem_kelvin import solve_fem

result = solve_fem(
    vol_file="model.vol",
    peec_step="model_peec.step",  # auto-detected from .vol.json if omitted
    frequency=0,
    I_total=1.0,
)
print(f"L = {result['L']*1e9:.2f} nH")
```

---

## 6. Usage Examples

### 6.1 3D A-Formulation with Kelvin

```python
from ngsolve import *
from radia.kelvin_source import (
    kelvin_nu_factor_3d_cf,
    build_material_cf,
)

mesh = Mesh("model.vol")

nu0 = 1 / (4e-7 * pi)
nu_cf = build_material_cf(
    mesh, nu0,
    kelvin_nu_factor_3d_cf(center=(0, 0, 0.15), R=0.06),
    overrides={"iron": nu0 / 1000},   # mu_r = 1000
)

fes = HCurl(mesh, order=2)
u, v = fes.TnT()
a = BilinearForm(fes)
a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
```

### 6.2 Axisymmetric with Z-Offset

```python
from radia.kelvin_source import kelvin_nu_factor_axisym_cf

nu_cf = build_material_cf(
    mesh, nu0,
    kelvin_nu_factor_axisym_cf(z_offset=0.15, R=0.06),
)

# Remember the r-weight!
a += nu_cf * grad(u) * grad(v) * x * dx
```

### 6.3 Reduced-A with PEEC Coil Source

```python
from radia.kelvin_source import A_s_at_obs_with_kelvin

# Define filament coil
filaments = [(path_points, current)]

# Evaluate A_s at observation points (auto Kelvin in exterior)
A_s = A_s_at_obs_with_kelvin(
    filaments, obs_points,
    center=(0, 0, 0.15), R=0.06,
    factor_mode='pullback'   # exact 1-form pullback (not scalar approx)
)
```

---

## 7. Reduced Potential Formulations + Kelvin

This section consolidates the reduced potential formulation theory with
Kelvin transformation: H-formulation (scalar potential), A-formulation
(vector potential), and the (ν - ν₀) form pitfall. Source: internal
note "Sugahara, Nagamine, Kameari 2026" (paper not submitted), extended
in this repository to A-formulation.

### 7.1 Reduced Potential Decomposition

In a reduced potential formulation, the total field is decomposed into
a known background contribution and an unknown perturbation:

| Formulation | Decomposition | Background `H_s / A_s` typically known |
|---|---|---|
| Ω (H-form) | `H = H_s - ∇Ω` | uniform applied B, Biot-Savart |
| A-form | `A = A_s + A_r` | Biot-Savart of PEEC coil |
| T-Ω | `T = T_s + T_r`, `Ω = Ω_s + Ω_r` | applied + induced |

The PERTURBATION (Ω, A_r, T_r, Ω_r) is solved by FEM; the BACKGROUND
is provided analytically.

### 7.2 H-formulation Weak Form

From `∇·B = 0` and `H = H_s - ∇m`:

```
-∇·(μ ∇m) = -∇·(μ H_s)
```

After integration by parts, the weak form is:

```
∫_Ω ∇v · (μ ∇m) dΩ = ∫_Ω ∇v · (μ H_s) dΩ - ∫_∂Ω v (n · μ H_s) dΓ    (1)
```

The boundary integral is essential for finite-domain (non-Kelvin)
analysis. With Kelvin + periodic BC, it vanishes automatically (see
§3.2). The volume term reduces in the kext to:

```
a(m, v) = ∫_Ω' (R/r')^2 ∇'v · (μ ∇'m) dΩ'    (5)
f(v)    = -∫_Ω' (r'/R)^2 ∇'v · (μ H_s) dΩ'   (6)
```

These can be absorbed by using modulated `μ' = (R/r')^2 μ_0` and
transformed background `H_s' = -(r'/R)^2 H_s` (see §7.5).

### 7.3 A-formulation Weak Form

From `∇×(ν ∇×A) = J`, decomposing `A = A_s + A_r`:

```
∇×(ν ∇×A_r) = J - ∇×(ν ∇×A_s)
```

Weak form (no Kelvin):

```
∫_Ω ν (∇×A_r)·(∇×v) dΩ = ∫_Ω J·v dΩ - ∫_Ω ν (∇×A_s)·(∇×v) dΩ     (7)
```

If `A_s` satisfies vacuum Ampere `∇×(ν₀ ∇×A_s) = J` everywhere, the
RHS simplifies via algebra to the popular `(ν - ν₀)` form:

```
∫_Ω ν (∇×A_r)·(∇×v) dΩ = -∫_Ω (ν - ν₀) (∇×A_s)·(∇×v) dΩ           (7')
```

This is standard in eddy-current FEM with iron yoke / magnetic body.

### 7.4 Background Field Transformation Rules (3D)

The metric-tensor approach (Sugahara 2022 / Nagamine CEFC 2026) gives
a 1-form transformation rule for any 1-form background field
(`H_s` for H-formulation, `A_s` for A-formulation):

```
H_s'(r') = -(r'/R)^2 H_s(r' - offset)    (3D)
A_s'(r') = -(r'/R)^2 A_s(r' - offset)    (3D)
H_s'(r') = -H_s(r' - offset)             (2D — sign flip only)
```

**Key features**:
- Evaluated at LOCAL (offset-relative) coordinates, NOT at the
  Kelvin-mapped physical point r_phys = T(r').
- NO Householder reflection (this is NOT the proper 1-form geometric
  pullback; that's Convention A used for PEEC source pullback).
- Vanishes at offset (rho' → 0 = physical infinity), no singularity.
- Sign flip ensures matching at the periodic Kelvin boundary
  rho' = R where inner and exterior normals are opposite.

#### Why local (offset-relative) evaluation?

For position-dependent backgrounds (e.g., `A_s = (B₀/2)(-y, x, 0)`):
- Global coords would inject a spurious offset-dependent term
  (e.g., the `o_x = 0.15 m` constant dominates A_s in kext).
- Local coords keep the natural symmetric gauge centered at offset
  — the Kelvin sphere represents physical infinity centered there.
- For uniform fields (`H_s = (0,0,1)`), local vs global is identical.

API: `radia.kelvin_material.make_reduced_potential_background_cf`.

### 7.4.1 Kameari Canonical Pattern (boundary-integral source, NOT bulk)

⭐ **Reference**: Kameari (2025/10/14 slides; private slide deck).

The historical "(ν-ν₀) bulk source" reduced-A form (next subsection §7.5)
is one route. Kameari's canonical NGSolve+Kelvin recipe takes a
DIFFERENT route: source coupling via **boundary integral on the
inner-Kelvin interface ∂Ω**, not via volume integral with (ν-ν₀).

**Three reductions** (Kameari slide 4):

| Method | Inner | Outer (Kelvin) | Use case |
|---|---|---|---|
| Ω-Ω_r | H = -∇Ω_t | H = -∇Ω_r + H_s | Linear magnetic, scalar |
| **A-Ω_r** | B = ∇×A_t | H = -∇Ω_r + H_s | **Recommended for sphere/cuboid** |
| A-A_r | B = ∇×A_t | B = ∇×(A_r + A_s) | Vector source |
| A-φ-A_r | A_t,φ in cond | A_r in air/Kelvin | Eddy current (TEAM 7) |

**Weak forms (Kameari slides 7-9)**:

```python
# A-Ω_r weak form (slide 8) — for vector A in conductor, Ω_r outside
a += 1/mu * curl(N) * curl(A) * dx                       # bulk
f += N.Trace() * Cross(H_s, n) * ds(inner_kelvin_bdry)   # ∂Ω boundary integral

# Ω-Ω_r weak form (slide 7) — for scalar Ω
a += mu * grad(omega) * grad(Omega) * dx
f += -omega.Trace() * (B_s * n) * ds(inner_kelvin_bdry)
```

The applied field (`H_s` or `B_s`) enters via a **boundary integral on
the interface** between inner physical region and Kelvin region (or
between inner and outer truncation). NOT in the bulk.

**Performance** (Kameari slide 17, magnetic sphere a=1m, μ_r=1000,
B_0=1T, theory Bz0=2.9940 T, **coarse mesh**):

| Method | Order | Bz0 | Error | Wm |
|---|---|---|---|---|
| Ω-Ω_r | 2 | 3.3813 | 12.9% | — |
| Ω-Ω_r | 3 | 2.9995 | 0.184% | 1883.8 |
| Ω-Ω_r | 4 | 2.9985 | 0.029% | 1868.4 |
| A-A_r | 3 | 2.9928 | 0.040% | 1867.8 |
| **A-Ω_r** | **3** | **2.9940** | **0.001%** ⭐ | **1867.6** |

**Independence from Kelvin radius rk** (slide 18): rk = 100 (very large)
gives the same accuracy as small rk with adaptive refinement. The
Kelvin region size is a free parameter.

**Recommendation**: for "isolated conductor + applied uniform B" problems,
use Kameari's A-Ω_r with Order ≥ 3 and boundary-integral source. The
(ν-ν₀) bulk-source approach (next §7.5) requires extra care and was
the source of the v11-v14 cuboid CLN issues.

### 7.5 The (ν − ν₀) Form Pitfall with Kelvin (CRITICAL)

The popular `(ν - ν₀)` reduced-A simplification (eq. 7' above) is
**INVALID when combined with Kelvin pullback** for the source A_s in
the kext.

**Why?** The simplification requires `A_s` to satisfy `∇×(ν₀ ∇×A_s) = J`
in the WHOLE domain (Ampere in vacuum). When A_s in kext is the
Kelvin pullback, it satisfies `∇×(ν' ∇×A_s) = 0` (the pulled-back
vacuum Maxwell equation), NOT the `ν₀` version. The algebraic
substitution is no longer valid.

**Symptom**: applying the wrong form gives wildly inflated FEM
inductance — empirically +43% on a torus + Kelvin benchmark.

**Correct form** (returning to eq. 7 directly):

```
a(A_r, v) = (J, v)_inner - ∫_full ν · (∇×A_s) · (∇×v) dV
         = - ∫_kext ν' · (∇×A_s_pullback) · (∇×v) dV
```

The inner contribution cancels via Ampere's law (where A_s = Biot-Savart
satisfies ν₀ Maxwell as expected).

**Validation** (archived classic source
`examples/.../Coil_3D_A_HCurl_PEEC_source.py`, preserved in
`docs/kelvin/kelvin_classic_demos_results.json`, 2026-05-04):

| Linear form | L (FEM) | Analytical | Error |
|---|---|---|---|
| Old: `-(ν - ν₀) (∇×A_s) (∇×v) dx` | 127.06 nH | 88.55 nH | **+43.48%** |
| New: `-ν' (∇×A_s) (∇×v) dx("kelvin")` | 93.97 nH | 88.55 nH | **+6.12%** |

Reference baseline (J-source, no PEEC): +4.80%.

### 7.6 Two Conventions for "Background Field Transformation"

These are EASY to confuse. Both arise from Kelvin, but for different
purposes:

| Convention | Formula (3D, uniform H_s = ẑ) | Use case | API |
|---|---|---|---|
| **A. Proper 1-form pullback** | `H_comp(r') = (R/ρ')² × H × H_phys(r_phys)` | Evaluate physical field at computational frame point (PEEC source pullback into kext for inductance / energy) | `kelvin_pullback_vector`, `make_kelvin_aware_A_s_cf` |
| **B. Reduced-potential background** | `H_s'(r') = -(ρ'/R)² × H_phys(r' - offset)` | Define background field for reduced-potential weak form | `make_reduced_potential_background_cf` |

**Convention A** has 1/ρ'^3 singularity at offset for unbounded
H_phys (e.g., uniform field). Used only for sources that DECAY at
infinity (PEEC coils).

**Convention B** is bounded everywhere (vanishes at offset). Used for
globally-defined backgrounds (uniform B, dipole, quadrupole at
infinity).

### 7.7 Axisymmetric A-formulation (special case)

For axisymmetric problems, a scalar variable transformation simplifies
the FE system. Using cylindrical coords `(r, θ, z)` with θ-symmetry:

```
A = A_θ(r, z) e_θ        (only azimuthal component)
u = r A_θ                (scalar variable)
B_r = -(1/r) ∂u/∂z
B_z =  (1/r) ∂u/∂r
```

**Strong form** for u:
```
-∇·((ν/r) ∇u) = J_θ
```

**Reduced form** with `u = u_s + u_r`:
```
-∇·((ν/r) ∇u_r) = J_θ + ∇·((ν/r) ∇u_s)
```

If `u_s` is the vacuum source (e.g., `u_s = B_0 r^2/2` for uniform B_z),
the (ν - ν₀) form simplification yields:

```
-∇·((ν/r) ∇u_r) = ∇·(((ν_0 - ν)/r) ∇u_s)
```

**Weak form** with axisymmetric `r dr dz` weight:
```
∫ (ν/r) ∇u_r · ∇v · r dr dz = -∫ ((ν_0 - ν)/r) ∇u_s · ∇v · r dr dz + boundary
```

The 1/r and r factors cancel, leaving a clean H1-like form (this is a
special feature of the axisymmetric `u = rA_θ` formulation).

**Z-offset Kelvin variant**: for axisym problems the Kelvin "sphere" is
implemented by offsetting in z (not the standard 3D sphere offset).
This preserves r-coordinate values across inner / kext, retaining
axisymmetry after the transformation:

```
inner: (r, z)         →   kext: (r, z + z_offset)
ρ' = sqrt(r^2 + (z' - z_offset)^2) = sqrt(r^2 + z^2)
ν_outer = ν_0 (ρ'/a)^2     (a = Kelvin radius)
```

API: `radia.kelvin_source.kelvin_nu_factor_axisym_cf(z_offset, R)`.

### 7.8 Anisotropic Permeability Support

The H-formulation bilinear form `∫ μ ∇u · ∇v dΩ` supports tensor μ
directly (no special handling needed beyond using a tensor-valued
CoefficientFunction):

```python
mu_tensor = CoefficientFunction((
    (mu_rr, 0),
    (0, mu_zz)
), dims=(2,2))
a += InnerProduct(mu_tensor * grad(u), grad(v)) * dx
```

This is essential for **2D cylindrical Kelvin** which produces an
anisotropic `μ' = diag(1, 1, (R/ρ')^4) μ` (only z-component modulated;
see §2.4).

### 7.9 Numerical Validation (H-formulation)

Magnetic sphere µ_r = 100 in uniform applied field, R_K = 1 m:

| Method | Dipole | Quadrupole |
|---|---|---|
| Finite domain (no Kelvin) | boundary-dependent error | boundary-dependent error |
| **Kelvin transformation** | matches analytical | matches analytical |

Toroidal current loop (a = 0.1 m, b = 0.01 m, R_K = 1 m):
- Dipole exterior energy analytical: 3.333e-8 J
- FEM with `ν' = (ρ'/R)² ν₀`: 3.344e-8 J (**+0.33%** vs analytical)

The Kelvin formulation reproduces analytical results accurately for
both background field types.

#### Analytical solutions for benchmark cases

For a magnetic sphere of radius `a` and relative permeability `μ_r`:

**Dipole background** (`H_s = (0, 0, 1)`):
- Interior: `H_z,pert = -1 + 3/(μ_r + 2)`
- Exterior: `H_pert ∝ a³/r³` (decay slow)

**Quadrupole background** (`H_s = -∇(xz) = (-z, 0, -x)`):
- Coefficient `B = -2(μ_r - 1)/(2μ_r + 5)`
- Interior: `H_pert = B(-z, 0, -x)`
- Exterior: `H_pert ∝ a⁵/r⁵` (decay faster than dipole)

For µ_r = 100, a = 0.5 m: B = -0.96618. The faster decay of
quadrupole means smaller external boundary errors in the non-Kelvin
case, but Kelvin transformation handles both cases uniformly.

---

## 8. Known Limitations

1. **Source in Kelvin exterior**: PEEC filaments must lie in the physical
   interior domain. Sources extending into the Kelvin exterior are not
   supported.

2. **Mesh refinement near rho'=0**: Material nu vanishes (mu diverges)
   at the Kelvin center. Use `bonus_intorder >= 4` and consider adaptive
   refinement near rho'=0.

3. **Gapped coils with Omega formulation**: Omega-reduced-Omega fails
   when `div(J) != 0` at coil terminals (~12% error). Use full
   A-formulation for gapped coil geometries.

4. **Radia MMM as reference**: Radia's Method of Magnetic Moments may
   show 4-8% variation across mesh densities for Kelvin validation.
   Use 2D axisymmetric FEM as ground truth instead.

5. **Order 3+ Kelvin instability**: Higher polynomial orders (p >= 3)
   on the Kelvin domain may require regularization or adapted meshes
   near rho'=0 to avoid ill-conditioning.

---

## 9. References

### Primary

1. **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based
   Formulation of Kelvin Transformation in Electromagnetic Field
   Analysis," CEFC 2026 (Thessaloniki), id 350.
   - Derives material transformation via pullbacks of differential forms
     and bilinear energy functionals.
   - Canonical result: `nu' = (rho'/R)^2 nu_0` for 3D spherical;
     `nu' = diag(1, 1, (rho'/R)^4) nu` for 2D cylindrical.
   - Validated: toroidal dipole energy +0.33%.

2. **K. Sugahara**, "Electromagnetic analysis of eddy current testing
   with Kelvin transformation," IEEE Trans. Magn. 58(9), 1-6, Sept. 2022.
   - Original A-formulation derivation (cited as ref [3] in Nagamine).

3. **K. Sugahara, H. Nagamine, A. Kameari** (internal note, 2026):
   "Kelvin Transformation for Open Boundary Problems in Reduced
   Potential Formulation". Reflected in
   [Reduced_Potential_Kelvin.md](Reduced_Potential_Kelvin.md).
   - Paper not submitted; content lives in this repository.
   - H-formulation rule: `H'_s = -(rho'/R)^2 H_s` (3D).
   - Reduced-A extension and (nu - nu_0) form pitfall: derived in
     this repository's docs/kelvin/, not in the original digest.

### Classical / Background

4. **A. Bossavit**, *Computational Electromagnetism*, Academic Press, 1998.
   - Differential geometric framework for EM (ref [4] in Nagamine).

5. **E.M. Freeman, D.A. Lowther**, "An open boundary technique using
   a modified Kelvin transformation," IEEE Trans. Magn. 24(6), 1988.
   - Early Kelvin FEM implementation (ref [2] in Nagamine).

6. **S.K.M. Wong, I.R. Ciric**, "Method of conformal transformation
   for open-boundary electromagnetic field problems," COMPEL 4(3), 1985.
   - Open boundary via conformal map (ref [1] in Nagamine).

7. **O. Kuwahara, T. Takeda**, "Kelvin transformation for unbounded
   FEM," 1990.
   - Early treatment referenced in EDDY_CURRENT_METHODS.md.

8. **A. Nabizadeh, R. Ramamoorthi, A. Chern**, "Kelvin transformations
   for simulations on infinite domains," ACM Trans. Graphics 40(4), 2021.
   - General k-form pullback framework (computer graphics perspective).

---

## 10. File Index

### Source Code

| File | Layer | Description |
|------|-------|-------------|
| `src/radia/kelvin_source.py` | L0 | Math primitives, pullbacks, material CFs, Biot-Savart |
| `src/radia/kelvin_geometry.py` | L1 | Two-sphere OCC geometry builder |
| `src/radia/kelvin_material.py` | L2 | Mesh-aware material CoefficientFunctions |
| `src/radia/kelvin_solver.py` | L3 | Full-A and reduced-A FEM drivers |
| `src/radia/kelvin_validate.py` | L4 | End-to-end validation harness |
| `src/radia/panels/calc_fem_kelvin.py` | Panel | 3D FEM-SIBC Cubit panel solver |
| `src/radia/panels/add_kelvin.py` | Panel | Kelvin exterior automation |

### Documentation

| File | Description |
|------|-------------|
| **`docs/kelvin/KELVIN_TRANSFORMATION.md`** (this file) | Comprehensive Kelvin theory + API + workflow + reduced potential (consolidated 2026-05-04 from 8 prior docs) |
| `docs/kelvin/POLICY.md` | NGSolve version + environment policy |
| `docs/kelvin/api_plan.md` | API design plan and milestones |
| `docs/kelvin/Supplement/ErrorEstimator.md` | Equilibrated error estimator theory (adaptive mesh) |
| `docs/kelvin/Supplement/CG-smoother.md` | CG-smoother acceleration for error estimation |
| `docs/kelvin/Supplement/cg_smoother_demo.ipynb` | Runnable CG-smoother showcase (executed, embedded outputs) |
| `docs/kelvin/kelvin_examples_migration.ipynb` | Result-bearing initial migration ledger for the former Kelvin examples tree: 226 Python files classified into docs / validation_test / src-api / memory lanes, synchronized with JSON |
| `docs/kelvin/kelvin_classic_demos.ipynb` | Result-bearing source map for the 37 classic A/H/Omega/Radia-IEM Kelvin demos; the standalone example scripts were pruned after full-source JSON preservation |
| `docs/kelvin/kelvin_adaptive_mesh_archive.ipynb` | Result-bearing full-source archive for 59 AdaptiveMesh scripts; the first 45 repetitive `order=*` runners were pruned after source-hash preservation |
| `docs/kelvin/kelvin_dtn_spectrum_archive.ipynb` | Result-bearing full-source archive for 122 DtN-spectrum scripts; the standalone act scripts were pruned after source-hash preservation |
| `docs/kelvin/kelvin_remaining_examples_archive.ipynb` | Result-bearing final archive for the last 22 Kelvin example scripts; Cubit p-convergence moved to `validation_test`, the rest preserved as full-source records |
| `docs/kelvin/Supplement/experiment_cg_smoother_equilibration.py` | Test script for CG-smoother |
| `docs/kelvin/legacy_assets/kelvin_transformation/CONVENTION.md` | Canonical convention declaration (one-page) |
| `docs/solver/FEM_KELVIN_PLAN.md` | Implementation plan and status |

### Tests

| File | Description |
|------|-------------|
| `tests/test_kelvin_source.py` | 11 unit tests (involution, pullback, Biot-Savart) |
| `docs/kelvin/kelvin_classic_demos.ipynb` | Rendered docs layer for the 37 pruned classic formulation demos (A-form, H-form, Omega-ReducedOmega, Radia IEM vs FEM sphere), synchronized with full-source JSON |
| `validation_test/cubit/kelvin_1_4_p_convergence/` | Promoted Cubit p-convergence fixture covered by `validation_test/cubit/test_kelvin_1_4_p_convergence.py` |
| `docs/kelvin/kelvin_dtn_spectrum_archive.ipynb` | Open-boundary / DtN research prototypes archived with full source JSON; stable behavior is promoted into `src/radia/open_boundary` and `validation_test` |
| `docs/kelvin/kelvin_adaptive_mesh_archive.ipynb` | Archived source for AdaptiveMesh studies; repetitive per-order runners are no longer kept as standalone example scripts |
| `docs/kelvin/kelvin_remaining_examples_archive.ipynb` | Final examples cleanup ledger; Kelvin examples now contain no standalone Python scripts |

### MCP Knowledge

| File | Description |
|------|-------------|
| `src/radia/mcp_server/radia_ngsolve/kelvin_knowledge.py` | MCP server knowledge base |
