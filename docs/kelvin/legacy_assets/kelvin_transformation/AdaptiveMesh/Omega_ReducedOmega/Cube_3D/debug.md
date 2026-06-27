# Bug Report: External Domain (air_outer) ZZ Error Estimator

## Issue
The estimated error for the exterior domain (air_outer) in the 2x2 plot is always blue, indicating near-zero or incorrect error values.

## Root Causes

There are **TWO bugs** that need to be fixed:

### Bug 1: Wrong r_prime calculation for Kelvin transformation

The `r_prime` (distance from center of exterior domain) was calculated incorrectly:

**Original Code (BUGGY)**:
```python
r_prime_sq = x**2 + y**2 + (z - offset_z)**2  # WRONG!
```

**Fixed Code**:
```python
# For exterior domain: center at (2*kelvin_radius, 0, 0)
center_x_ext = 2 * kelvin_radius
r_prime_sq = (x - center_x_ext)**2 + y**2 + z**2
```

The exterior domain is centered at `(2*kelvin_radius, 0, 0) = (6, 0, 0)`, not offset in the z-direction.

### Bug 2: Wrong flux in compute_error_estimator function

The `compute_error_estimator` function was computing flux incorrectly for air_outer region.

### Original Code (BUGGY)
```python
def compute_error_estimator(mesh, fes, H_pert_cf, mu0, mu_r):
    """Compute ZZ-type error estimator using H(div) flux recovery."""
    B_bounded_dict = {
        "magnetic": (mu_r * mu0) * H_pert_cf,
        "air_total": mu0 * H_pert_cf,
        "air_inner": mu0 * H_pert_cf,
        "air_outer": mu0 * H_pert_cf  # BUG: Should use actual B_pert for exterior
    }
    flux = CoefficientFunction([B_bounded_dict[mat] for mat in mesh.GetMaterials()])
    # ... rest of the function
```

### Problems

1. **Wrong permeability for air_outer**: The exterior domain uses Kelvin-transformed permeability `mu_kelvin = (R/r')^2 * mu0`, but the error estimator uses plain `mu0`. This leads to vastly underestimated flux values.

2. **Using generic H_pert_cf**: The `H_pert_cf` coefficient function is defined with region-specific values, but when used directly with `mu0` for air_outer, it gives incorrect B values because:
   - For air_outer, `H_pert = grad(Omega)` (no Hs subtraction)
   - The actual B field is `B = mu_kelvin * H_pert`

3. **Why all blue**: Since `r_prime` in the Kelvin region is calculated as `sqrt(x^2 + y^2 + (z - offset_z)^2)` where `offset_z = 3.0`, and the exterior domain is centered at `x = 2*kelvin_radius = 6.0`, the coordinates in the exterior domain give values that don't match the actual field distribution. Using constant `mu0` instead of `mu_kelvin` leads to nearly zero error estimates.

## Fix

The function signature needs to accept `B_pert_cf` (which is already correctly computed per-region) instead of `H_pert_cf` and `mu0/mu_r`:

### Fixed Code
```python
def compute_error_estimator(mesh, fes, B_pert_cf):
    """Compute ZZ-type error estimator using H(div) flux recovery."""
    flux = B_pert_cf  # Use the correctly computed B_pert per region

    recovery_order = max(1, fes.globalorder - 1)
    fes_flux = HDiv(mesh, order=recovery_order)
    gf_flux = GridFunction(fes_flux)

    sigma = fes_flux.TrialFunction()
    tau = fes_flux.TestFunction()
    a_flux = BilinearForm(fes_flux)
    a_flux += InnerProduct(sigma, tau) * dx
    a_flux.Assemble()

    f_flux = LinearForm(fes_flux)
    f_flux += InnerProduct(flux, tau) * dx
    f_flux.Assemble()

    gf_flux.vec.data = a_flux.mat.Inverse(fes_flux.FreeDofs(), inverse="sparsecholesky") * f_flux.vec

    err = InnerProduct(flux - gf_flux, flux - gf_flux)
    element_errors = Integrate(err, mesh, element_wise=True)
    return element_errors
```

### Call Site Change
```python
# Before (BUGGY):
element_errors = compute_error_estimator(mesh, fes, H_pert_cf, mu0, mu_r)

# After (FIXED):
element_errors = compute_error_estimator(mesh, fes, B_pert_cf)
```

## Why This Fix Works

The `B_pert_cf` coefficient function is already correctly defined per-region in `solve_omega_formulation`:

```python
B_pert_dict = {
    "magnetic": (mu_r * mu0) * (grad(gfOmega) - Hs),
    "air_total": mu0 * (grad(gfOmega) - Hs),
    "air_inner": mu0 * grad(gfOmega),
    "air_outer": mu_kelvin * grad(gfOmega)  # Correctly uses mu_kelvin!
}
B_pert_cf = CoefficientFunction([B_pert_dict[mat] for mat in mesh.GetMaterials()])
```

This properly accounts for:
- The Kelvin-transformed permeability in `air_outer`
- The different H_pert definitions (with/without Hs) for total vs reduced regions
- Consistent flux recovery across all domains

### Bug 3: Periodic boundary identification AFTER Glue (Root Cause)

**Symptom**: `air_outer: 0.0000e+00 J` - exterior domain has zero energy even after fixing Bugs 1 and 2.

**Root Cause**: The periodic boundary between `air_inner` (Kelvin boundary) and `air_outer` was being identified AFTER the `Glue()` operation. In NGSolve/OCC, face references change after `Glue()`, so `face.Identify()` must be called BEFORE `Glue()`.

**Original Code (BUGGY)**:
```python
# Glue first
geo = Glue([mag_cube, air_total, air_inner, outer_sphere, vertex])
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh))

# Try to find faces AFTER Glue - FAILS!
kelvin_int_face = None
kelvin_ext_face = None
for face in geo.faces:
    if face.name == "kelvin_int":
        kelvin_int_face = face
    elif face.name == "kelvin_ext":
        kelvin_ext_face = face

if kelvin_int_face is not None and kelvin_ext_face is not None:
    kelvin_ext_face.Identify(kelvin_int_face, "periodic", IdentificationType.PERIODIC)
else:
    print("WARNING: Could not find periodic faces!")  # Always prints this!
```

**Fixed Code**:
```python
# Create shapes
air_inner = ...
outer_sphere = ...

# Identify periodic faces BEFORE Glue
print("Identifying periodic boundaries BEFORE Glue...")
kelvin_int_face = None
kelvin_ext_face = None

for face in air_inner.faces:
    fc = face.center
    dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
    if dist > kelvin_radius * 0.8:
        kelvin_int_face = face
        face.name = "kelvin_int"
        break

for face in outer_sphere.faces:
    fc = face.center
    dist = sqrt((fc.x - center_x)**2 + fc.y**2 + fc.z**2)
    if dist > kelvin_radius * 0.8:
        kelvin_ext_face = face
        face.name = "kelvin_ext"
        break

if kelvin_int_face is not None and kelvin_ext_face is not None:
    kelvin_ext_face.Identify(kelvin_int_face, "periodic", IdentificationType.PERIODIC)
    print("  Periodic identification applied BEFORE Glue!")

# THEN Glue
geo = Glue([mag_cube, air_total, air_inner, outer_sphere, vertex])
```

**Reference**: the historical `CubeMesh.py` source is preserved in
`docs/kelvin/kelvin_remaining_examples_archive_results.json`.  Its correct
pattern was:
```python
external_domain.faces[0].Identify(Omega_domain.faces[0], "ud0", IdentificationType.PERIODIC)
geo = Glue([iron, A_domain, Omega_domain, external_domain])  # Glue AFTER Identify
```

## Design Policy: Kelvin Transform Offset Direction

Kelvin変換で外部領域を有界領域にマッピングする際のオフセット方向のポリシー:

| モデルタイプ | オフセット方向 | 外部領域中心 |
|-------------|---------------|-------------|
| 3D全体/1/8モデル | z方向 | `(0, 0, 2*R)` |
| 軸対称 (2D) | y方向 | `(0, 2*R)` |

**理由**:
- 3Dモデル: z軸が通常「上方向」または対称軸となるため、z方向にオフセット
- 軸対称モデル: r-z平面で計算し、y座標がz軸に対応するため、y方向にオフセット

**重要**: オフセット方向は座標系の定義に依存する。1/8モデル（第1象限、x>=0, y>=0, z>=0）の場合、x方向にオフセットすることも可能だが、対称性の観点からz方向（または軸対称ではy方向）を推奨する。

## Files Modified

All 9 Python files in Cube_3D have been fixed:
1. `order=2/Refine_with_zz_estimator/Cube_3D_adaptive_with_Kelvin.py`
2. `order=2/Refine_all_elements/Cube_3D_adaptive_with_Kelvin.py`
3. `order=2/maxh/Cube_3D_maxh_with_Kelvin.py`
4. `order=3/Refine_with_zz_estimator/Cube_3D_adaptive_with_Kelvin.py`
5. `order=3/Refine_all_elements/Cube_3D_adaptive_with_Kelvin.py`
6. `order=3/maxh/Cube_3D_maxh_with_Kelvin.py`
7. `order=4/Refine_with_zz_estimator/Cube_3D_adaptive_with_Kelvin.py`
8. `order=4/Refine_all_elements/Cube_3D_adaptive_with_Kelvin.py`
9. `order=4/maxh/Cube_3D_maxh_with_Kelvin.py`
