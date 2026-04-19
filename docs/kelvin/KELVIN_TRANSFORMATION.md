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

### 2.6 Numerical Validation

Toroidal current loop (Nagamine CEFC 2026 Section III):
- Major radius a = 0.1 m, wire radius b = 0.01 m
- Analytical exterior dipole energy: 3.333e-8 J
- FEM on transformed domain with `nu' = (rho'/R)^2 nu_0`: 3.344e-8 J
- Error: **+0.33%**

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

7. **Export** — `radia_export netgen` automatically detects `air`+`kelvin`
   blocks and creates periodic identification + `kelvin_int`/`kelvin_ext`
   boundary labels
   ```
   radia_export netgen "model.vol" order 3 overwrite
   ```

### 5.2 Auto-Detection in radia_export netgen

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

When a block named `peec` is present, `radia_export netgen` automatically:
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
radia_export netgen "model.vol" order 2 overwrite
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

## 7. Known Limitations

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

## 8. References

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

3. **K. Sugahara, H. Nagamine, A. Kameari**, "Kelvin Transformation
   for Open Boundary Problems in Reduced Potential Formulation"
   (digest, 2026).
   - Companion tutorial: H-formulation rule `H'_s = -(rho'/R)^2 H_s`.

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

## 9. File Index

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
| `examples/kelvin_transformation/CONVENTION.md` | Canonical convention declaration |
| `examples/kelvin_transformation/docs/pullback_derivation_3D.md` | Full pullback derivation |
| `examples/kelvin_transformation/docs/Kelvin_2D.md` | 2D cylindrical Kelvin |
| `examples/kelvin_transformation/docs/Kelvin_3D.md` | 3D H-formulation |
| `docs/solver/FEM_KELVIN_PLAN.md` | Implementation plan and status |

### Tests

| File | Description |
|------|-------------|
| `tests/test_kelvin_source.py` | 11 unit tests (involution, pullback, Biot-Savart) |
| `examples/kelvin_transformation/A-formulation/` | 12 A-formulation examples |
| `examples/kelvin_transformation/H-formulation/` | 14 H-formulation examples |
| `examples/kelvin_transformation/Omega_ReducedOmega/` | Omega + coil source examples |
| `examples/kelvin_transformation/AdaptiveMesh/` | h/p adaptive refinement |

### MCP Knowledge

| File | Description |
|------|-------------|
| `src/radia/mcp_server/radia_ngsolve/kelvin_knowledge.py` | MCP server knowledge base |
