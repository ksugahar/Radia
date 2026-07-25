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

#### Where every sign comes from: straight vs twisted forms

The Kelvin inversion is **orientation-reversing** (`det J < 0`), and that
is the single origin of every minus sign in this document. In premetric
electromagnetism the field quantities split by **orientation type**:

| type | forms | behaviour under an orientation-reversing map |
|---|---|---|
| **straight** (inner-oriented) | `V` (0-form), `a`, `e` (1-forms), `b` (2-form) | pull back with the Jacobian alone |
| **twisted** (outer-oriented, densities) | **`phi_m` (0-form, the magnetic scalar potential)**, `h` (1-form), `d`, `j` (2-forms), `rho`, `U_m` (3-forms) | pick up an **extra** `sgn(det J) = -1` |

Form **degree** and **orientation parity** are independent axes: the
magnetic scalar potential is a 0-form *and* twisted, which is why its
pullback carries a minus but **no** metric factor (exponent 0 in the
table above). For the Kelvin map `det Dk = -R^6/rho'^6 < 0`, so
`s_k = -1` and

```
phi_m = -k* phi_m'  ,   H = -k* H'  ,   J = -k* J'  ,   rho = -k* rho'
V     =  k* V'      ,   A =  k* A'  ,   B =  k* B'
```

So the minus in the reduced-potential rule `H_s' = -(rho'/R)^2 H_s`
(§7.4) is *not* a fitting constant: `H` is a **twisted** 1-form, and the
inversion reverses orientation. `A` is straight, which is why its own
pullback (§2.3 table) carries no such factor, and `B`'s minus enters
instead through the pseudovector (Hodge-dual) representation.

The **energy is orientation-blind** — the pairing `<B, H>` contains one
straight and one twisted factor, so the two sign flips cancel. That is
why the material modulation `nu' = (rho'/R)^2 nu_0` (§2.4) carries **no**
sign at all, while the field rules do.

Cross-reference: `differential_forms_maxwell('twisted')` on
mcp-server-radia-differential-forms (Bossavit / Kameari diagram).

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

#### Convention B has no valid 0-form flavour — use the twisted pullback

It is tempting to apply the same `-(rho'/R)^2` factor to a background
**potential** so a T-Omega weak form can consume `Omega_s`. **That is
wrong.** A 0-form pullback carries **no** metric factor at all (§2.3
table, exponent 0); the `(rho'/R)^2` is a 1-form weight. Measured cost of
making that substitution: a factor **4/3** (§7.9, route B-0).

The rule that *does* work is the genuine **twisted 0-form pullback**:

```
Omega_s'(r') = - Omega_s( k(r') ) ,   k(r') = (R/rho')^2 (r' - offset)
```

— no metric factor, sign only, evaluated at the **mapped** point. The
minus is `s_k = sgn(det Dk) = -1`, because the magnetic scalar potential
is a **twisted** 0-form (§2.3).

Its partner for the field is the twisted 1-form pullback

```
H_s'(r') = -(R/rho')^2 (I - 2 n n^T) H_s( k(r') )
```

i.e. the **radial** component keeps its sign and the **tangential**
components flip. Because pullback commutes with the exterior derivative
(`g*(dw) = d(g*w)`) and both quantities carry the *same* twist factor,

```
H_s'  ==  -grad'( Omega_s' )      EXACTLY
```

— verified to `1.4e-10` (finite-difference limited) in
`tests/test_reduced_potential_background.py`. **This is the property
Convention B lacks**, and it is what makes the potential route usable.

API: `radia.kelvin_material.make_kelvin_aware_Omega_s_cf` and
`make_kelvin_aware_H_s_cf`. (`A` and `B` are *straight* forms and take no
extra minus — that is the existing `make_kelvin_aware_A_s_cf`.)

**Regularity, and the one real limitation.** For a source that **decays**
at infinity the pullback is regular at the offset: a dipole
`|H_s| ~ 1/r^3` maps to `|H_s'| ~ rho'/R^4`, and `Omega_s' = O(rho'^2)`
(measured: `7.0e-3 -> 1.6e-5` as `rho'` goes `0.4 -> 0.02`). "Vanishing at
infinity" becomes "vanishing at the Kelvin centre". For a background that
does **not** decay — a uniform field applied at infinity — `Omega_s'`
diverges like `R^2/rho'^2` (measured: `2.5 -> 20` as `rho'` goes
`0.4 -> 0.05`), because the uniform-field potential is genuinely unbounded
at infinity. There is **no** bounded 0-form representative in that case;
drive it through the 1-form Convention B route instead.

So: **real coil / dipole sources -> the potential route works**; a uniform
background applied at infinity is the one case that must stay on the field
route.

Reference: <https://www.ele.kindai.ac.jp/laboratory/sugahara/elemag/geometry09.php>
(twisted-form sign table).

`make_reduced_potential_scalar_cf` is retained only because the T-Omega
design note proposes it and its behaviour is contract-locked; it must not
be used to drive a weak form.

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
infinity). It comes in a 1-form and a 0-form flavour:

| Convention B flavour | Formula (3D) | Weak form consumes | API |
|---|---|---|---|
| **B-1 (1-form)** | `H_s'(r') = -(rho'/R)^2 H_s(r'-offset)` | a background **field**, `int mu H_s . grad v` | `make_reduced_potential_background_cf` |
| **B-0 (0-form)** | `Omega_s'(r') = -(rho'/R)^2 Omega_s(r'-offset)` | a background **potential** (T-Omega `Omega_s`) | `make_reduced_potential_scalar_cf` |

B-1 and B-0 are **not** two views of one field: B-1 is not curl-free, so
it has no potential, and differentiating B-0 does not give B-1 (§7.4).
Mixing them costs a factor 4/3 — measured in §7.9.

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

#### What the exterior source is worth (3-route golden, 2026-07-23)

Magnetic sphere `mu_r = 100`, `a = 0.5 m`, uniform applied `H_0 z_hat`,
`R_K = 1 m`, reduced scalar potential `H = H_s - grad(Omega)`, analytical
`H_in = 3 H_0/(mu_r + 2)`. One mesh, one bilinear form; the routes differ
**only** in what `H_s` is inside the Kelvin exterior:

| route | `H_s'` in kext | err (p=2, maxh .14) | err (p=3, maxh .11) | ratio to B-1 |
|---|---|---|---|---|
| **Z** | `0` (source dropped) | −32.83% | −33.33% | **2/3** |
| **B-1** | `-(rho'/R)^2 H_0 z_hat` | +0.74% | **+0.002%** | 1 |
| **B-0** | `-grad(Omega_s')` | +34.33% | +33.34% | **4/3** |

Readings:

1. **B-1 is exact** (+0.002% at p=3): the 1-form Convention B rule is the
   correct way to carry a background field into the Kelvin exterior.
2. **Dropping the exterior source loses exactly one third.** The 2/3 and
   4/3 factors are mesh-independent (0.666669 / 1.333338 at p=3), so they
   are structural, not discretisation error.
3. **Differentiating the 0-form rule into a field overshoots by 4/3** —
   the numerical face of `curl(H_s' 1-form B) != 0` from §7.4.

Golden: `validation_test/kelvin_source/test_kelvin_exterior_source_routes.py`.
Contract locks (symbolic identities, curl, boundedness):
`tests/test_reduced_potential_background.py`.

#### A-Phi (A-V) eddy formulation at p = 2 (measured 2026-07-25)

Conducting sphere `a/delta = 2` in a uniform AC field vs the analytic Smythe
moment, two lanes on one mesh (every verified-recipe element shared):

| | p=1 | p=2 | p=3 | h-sweep at p=2 (ne 10k → 28k → 68k) |
|---|---|---|---|---|
| **A\*** (plain A-method) | 2.889% | 0.473% | 0.442% (saturated) | 0.473 → 0.243 → 0.141% (~O(h)) |
| **A-Phi** (mixed) | 2.889% | **0.053%** | 0.005% | 0.053 → 0.011 → 0.004% |

`nograds=True` (required for Periodic Kelvin) removes the gradient test
functions that enforce discrete charge conservation; the explicit Phi block
restores them. SIBC-driven problems stay plain-A; a volume conductor at
p >= 2 takes the A-Phi block. Golden:
`validation_test/kelvin_source/test_aphi_kelvin_eddy.py`.

#### Rendered walkthrough

`docs/kelvin/kelvin_exterior_source_and_aphi.ipynb` (executed, with the
synchronized `kelvin_exterior_source_and_aphi_results.json`) re-runs all
three goldens of this section — the exterior-source routes, the twisted
0-form pullback contracts, and the A\*/A-Phi p-sweep — and embeds the
outputs.

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

4. **Reference solver choice**: Old moment-path Radia comparisons showed
   mesh-density variation for Kelvin validation.  Use 2D axisymmetric FEM or
   HDiv-VIM / reduced-FEM cross-checks as the reference path instead.

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
| `docs/kelvin/Supplement/experiment_cg_smoother_equilibration.py` | Test script for CG-smoother |
| `docs/kelvin/CONVENTION.md` | Canonical convention declaration (one-page) |
| `docs/solver/FEM_KELVIN_PLAN.md` | Implementation plan and status |

### Tests

| File | Description |
|------|-------------|
| `tests/test_kelvin_source.py` | 11 unit tests (involution, pullback, Biot-Savart) |
| `validation_test/cubit/kelvin_1_4_p_convergence/` | Promoted Cubit p-convergence fixture covered by `validation_test/cubit/test_kelvin_1_4_p_convergence.py` |
| `src/radia/open_boundary/` | Productionized DtN / CLN open-boundary APIs |

### MCP Knowledge

| File | Description |
|------|-------------|
| `src/radia/mcp_server/radia_ngsolve/kelvin_knowledge.py` | MCP server knowledge base |
