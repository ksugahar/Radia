# Omega-Reduced Omega Method with Kelvin Transformation

## Overview

This document describes the correct implementation of field calculations and energy evaluation for the Omega-Reduced Omega magnetostatic formulation with Kelvin transformation.

---

## 1. Formulation Summary

### 1.1 Domain Structure

| Region | Material | Formulation | Description |
|--------|----------|-------------|-------------|
| `magnetic` | Magnetic cylinder (mu_r) | Total | H = grad(Omega) |
| `air_inner` | Air (mu0) | Reduced | H = grad(Omega) + Hs |
| `air_outer` | Kelvin-transformed exterior | Reduced | H = grad(Omega) + Hs_kelvin |

### 1.2 Boundary Conditions

```python
# Dirichlet BC: Omega = 0 on symmetry planes and GND
fes_before = H1(mesh, order=fe_order, dirichlet="GND|sym_z|sym_z_ext")
fes = Periodic(fes_before)  # Periodic BC on Kelvin boundary

# Neumann BC: cylinder boundary
# Set Omega_s on cylinder boundary, apply normal flux condition
gfOmega.Set(Omega_s, BND, mesh.Boundaries("cylinder"))
f += (normal * Bs) * psi * ds("cylinder")
```

### 1.3 Source Field Definitions

```python
# Source potential and field
Omega_s = H0 * z
Hs = CoefficientFunction((0.0, 0.0, H0))
Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

# Kelvin-transformed source field (continuity at r' = R)
# Uses (R/r')² which equals 1 at boundary
Hs_kelvin = CoefficientFunction((0.0, 0.0, (kelvin_radius / r_prime)**2 * H0))
Bs_kelvin = mu0 * Hs_kelvin
```

---

## 2. Field Calculation (solve_omega_formulation)

### 2.1 Perturbation Fields (H_pert, B_pert)

| Region | H_pert | B_pert |
|--------|--------|--------|
| `magnetic` | `grad(Omega) - Hs` | `mu_r * mu0 * (grad(Omega) - Hs)` |
| `air_inner` | `grad(Omega)` | `mu0 * grad(Omega)` |
| `air_outer` | `grad(Omega)` | `mu_kelvin * grad(Omega)` |

**Key point**: In magnetic region (Total formulation), subtract Hs. In air regions (Reduced formulation), H_pert = grad(Omega) directly.

### 2.2 Total Fields (H_total, B_total)

| Region | H_total | B_total |
|--------|---------|---------|
| `magnetic` | `grad(Omega)` | `mu_r * mu0 * grad(Omega)` |
| `air_inner` | `grad(Omega) + Hs` | `mu0 * (grad(Omega) + Hs)` |
| `air_outer` | `grad(Omega) + Hs_kelvin` | `mu_kelvin * (grad(Omega) + Hs_kelvin)` |

### 2.3 Implementation Code

```python
def solve_omega_formulation(mesh, fe_order):
    # ... (FEM setup and solve) ...

    # Kelvin permeability
    r_prime = sqrt(x**2 + y**2 + (z - offset_z)**2 + 1e-20)
    mu_kelvin = (kelvin_radius / r_prime)**2 * mu0

    # Source fields by region
    Hs_kelvin = CoefficientFunction((0.0, 0.0, (kelvin_radius / r_prime)**2 * H0))

    # H_pert by region
    H_pert_dict = {
        "magnetic": grad(gfOmega) - Hs,    # Total: subtract Hs
        "air_inner": grad(gfOmega),         # Reduced: no subtraction
        "air_outer": grad(gfOmega)          # Kelvin: no subtraction
    }
    H_pert_cf = CoefficientFunction([H_pert_dict[mat] for mat in mesh.GetMaterials()])

    # B_pert by region
    B_pert_dict = {
        "magnetic": (mu_r * mu0) * (grad(gfOmega) - Hs),
        "air_inner": mu0 * grad(gfOmega),
        "air_outer": mu_kelvin * grad(gfOmega)
    }
    B_pert_cf = CoefficientFunction([B_pert_dict[mat] for mat in mesh.GetMaterials()])

    # Return fields dictionary
    fields = {
        'H_pert_cf': H_pert_cf,
        'B_pert_cf': B_pert_cf,
        'H_total_cf': H_total_cf,
        'B_total_cf': B_total_cf,
        'Hs_cf': Hs_cf,
        'Bs_cf': Bs_cf,
        'mu_kelvin': mu_kelvin,
    }
    return fes, gfOmega, Mu, fields
```

---

## 3. Perturbation Energy Calculation

### 3.1 Energy Formula by Region

The perturbation energy uses `W = 0.5 * mu * |H_pert|²`:

| Region | Formula | Notes |
|--------|---------|-------|
| `magnetic` | `0.5 * mu_r * mu0 * \|grad(Omega) - Hs\|²` | Use separate GridFunction |
| `air_inner` | `0.5 * mu0 * \|grad(Omega) - grad(Omega_s)\|²` | Subtract boundary contribution |
| `air_outer` | `0.5 * mu_kelvin * \|grad(Omega)\|²` | Use mu_kelvin directly |

### 3.2 Important: air_inner Region

For `air_inner`, the perturbation field must account for the boundary condition contribution:

```python
# Create separate GridFunctions for air regions
fesOr = H1(mesh, order=order, definedon="air_inner|air_outer")
Orr = GridFunction(fesOr)
Oxr = GridFunction(fesOr)
Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr.Set(Omega_s, BND, mesh.Boundaries("cylinder"))

# H_pert in air_inner = grad(Omega) - grad(Omega_s from boundary)
H_pert_reduced = grad(Orr) - grad(Oxr)
energy_air_inner = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)
```

### 3.3 Implementation Code

```python
# ===== Perturbation Energy Calculation =====

# Source field and potential
Hs = CoefficientFunction((0.0, 0.0, H0))
Omega_s = H0 * z

# --- Total region (magnetic cylinder) ---
H_pert_total = grad(gfu) - Hs
energy_magnetic = Integrate(0.5 * (mu_r * mu0) * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)

# --- Reduced region (air_inner) ---
# Use separate GridFunctions to compute H_pert = grad(Omega) - grad(Omega_s)
fesOr = H1(mesh, order=order, definedon="air_inner|air_outer")
Orr = GridFunction(fesOr)
Oxr = GridFunction(fesOr)
Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr.Set(Omega_s, BND, mesh.Boundaries("cylinder"))
H_pert_reduced = grad(Orr) - grad(Oxr)
energy_air_inner = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)

# --- Kelvin region (air_outer) ---
H_pert_kelvin = grad(Orr)
energy_air_outer = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)

# Total perturbation energy (1/8 model -> full = 8x)
energy_1_8 = energy_magnetic + energy_air_inner + energy_air_outer
energy_full = 8 * energy_1_8

# Per-region energy (1/8 -> full = 8x)
energy_magnetic_full = 8 * energy_magnetic
energy_air_inner_full = 8 * energy_air_inner
energy_air_outer_full = 8 * energy_air_outer

# Interior/Exterior energy split (for convergence analysis)
# Interior = physical domain (magnetic + air_inner)
# Exterior = Kelvin-transformed domain (air_outer only)
energy_interior_full = energy_magnetic_full + energy_air_inner_full
energy_exterior_full = energy_air_outer_full
```

### 3.4 History Tracking

```python
# History tracking with per-region energy
history = {
    'ndof': [],
    'elements': [],
    'elements_interior': [],   # magnetic + air_inner elements
    'elements_exterior': [],   # air_outer elements (Kelvin region)
    'elements_magnetic': [],   # magnetic elements only
    'elements_air_inner': [],  # air_inner elements only
    'elements_air_outer': [],  # air_outer elements only
    'error': [],
    'Hz_error_percent': [],
    'energy': [],
    'energy_interior': [],     # magnetic + air_inner energy
    'energy_exterior': [],     # air_outer energy
    'energy_magnetic': [],     # magnetic region only
    'energy_air_inner': [],    # air_inner region only
    'energy_air_outer': []     # air_outer region only
}

# Record history (inside iteration loop)
history['energy_magnetic'].append(energy_magnetic_full)
history['energy_air_inner'].append(energy_air_inner_full)
history['energy_air_outer'].append(energy_air_outer_full)
```

### 3.5 Energy Output with Error Percentage

```python
# Print per-region energy with analytical comparison
print(f"  Energy (8×1/8): {energy_full:.6e} J (analytical: {W_analytical:.6e} J)")
print(f"    magnetic: {energy_magnetic_full:.6e} J (ana: {W_magnetic_analytical:.6e} J, err: {(energy_magnetic_full/W_magnetic_analytical-1)*100:+.2f}%)")
print(f"    air_inner: {energy_air_inner_full:.6e} J (ana: {W_air_inner_analytical:.6e} J, err: {(energy_air_inner_full/W_air_inner_analytical-1)*100:+.2f}%)")
print(f"    air_outer: {energy_air_outer_full:.6e} J (ana: {W_air_outer_analytical:.6e} J, err: {(energy_air_outer_full/W_air_outer_analytical-1)*100:+.2f}%)")
```

---

## 4. Output Data Structure

### 4.1 Shared output_data Dictionary

Both VTK and PNG outputs use a common dictionary:

```python
output_data = {
    'mesh': mesh,
    'gfu': gfu,
    # Perturbation fields (from solve_omega_formulation)
    'B_pert_cf': B_pert_cf,
    'H_pert_cf': H_pert_cf,
    # Total fields
    'B_total_cf': B_total_cf,
    'H_total_cf': H_total_cf,
    # Source fields
    'Bs_cf': Bs_cf,
    'Hs_cf': Hs_cf,
    # Error estimator
    'element_errors': element_errors,
}
```

### 4.2 VTK Output

```python
def output_vtk(mesh, iteration, output_data):
    coefs = [output_data['gfu']]
    names = ["Omega"]

    # Add perturbation, total, and source fields
    for key, name in [('B_pert_cf', 'B_pert'), ('H_pert_cf', 'H_pert'),
                      ('B_total_cf', 'B_total'), ('H_total_cf', 'H_total'),
                      ('Bs_cf', 'Bs'), ('Hs_cf', 'Hs')]:
        if key in output_data and output_data[key] is not None:
            coefs.append(output_data[key])
            names.append(name)

    vtk = VTKOutput(mesh, coefs=coefs, names=names, filename=vtk_path)
    vtk.Do()
```

### 4.3 PNG Output (Convergence Plot)

```python
def generate_convergence_plot(iter_num, history, output_data):
    mesh = output_data['mesh']
    element_errors = output_data['element_errors']
    # ... generate 2x2 plot ...
```

### 4.4 MAT File Output (with per-region data)

```python
mat_data = {
    'ndof': array(history['ndof']),
    'elements': array(history['elements']),
    'error': array(history['error']),
    'energy': array(history['energy']),
    'energy_interior': array(history['energy_interior']),
    'energy_exterior': array(history['energy_exterior']),
    # Per-region energy tracking
    'energy_magnetic': array(history['energy_magnetic']),
    'energy_air_inner': array(history['energy_air_inner']),
    'energy_air_outer': array(history['energy_air_outer']),
    # Analytical values for comparison
    'W_analytical': W_analytical,
    'W_in_analytical': W_in_analytical,
    'W_out_analytical': W_out_analytical,
    'W_magnetic_analytical': W_magnetic_analytical,
    'W_air_inner_analytical': W_air_inner_analytical,
    'W_air_outer_analytical': W_air_outer_analytical
}
sio.savemat(mat_filename, mat_data)
```

---

## 5. Kelvin Transformation Details

### 5.1 Coordinate Transformation

- Physical space: `r > R` (exterior domain)
- Kelvin space: `r' < R` (computational domain)
- Mapping: `r = R²/r'` (inversion)

### 5.2 Permeability Transformation (3D)

```python
mu_kelvin = (kelvin_radius / r_prime)**2 * mu0  # = (R/r')² * mu0
```

At boundary (r' = R): `mu_kelvin = mu0` (continuous)
At center (r' → 0): `mu_kelvin → ∞` (singularity)

### 5.3 Source Field Transformation

```python
# Correct: (R/r')² ensures continuity at boundary
Hs_kelvin = (kelvin_radius / r_prime)**2 * H0

# WRONG: -(r'/R)² causes discontinuity!
# Hs_kelvin = -(r_prime / kelvin_radius)**2 * H0  # DO NOT USE
```

---

## 6. Simulation Parameters

### 6.1 DOF Limit

```python
# Stop simulation when DOF reaches 100,000
if prev_ndof >= 1e5:
    print(f"\n  DOF limit reached ({prev_ndof} >= 1e5), stopping without computing.")
    break
```

### 6.2 Convergence Plot Settings

```python
# X-axis limit for DOF plots
ax.set_xlim(1e2, 1e5)
```

### 6.3 ZZ Error Estimator

**Important**: Use bounded B field for ZZ error estimator to avoid `mu_kelvin` singularity at Kelvin center (r' → 0).

**Algorithm**: Percentage-based marking (no normalization)
- Compute unnormalized element errors
- Mark elements with `error >= theta * max_error` for refinement
- This ensures all regions are refined proportionally to their actual error

**Key Implementation Details**:
1. **L2 projection (not interpolation)**: Use Galerkin projection to H(div) space
2. **Recovery order**: `recovery_order = fes.globalorder - 1` (one order lower than solution space)
3. **Bounded flux**: Use `mu0` instead of `mu_kelvin` in air_outer to avoid singularity

```python
def compute_error_estimator(mesh, fes, H_pert_cf, mu0, mu_r):
    # Create bounded B field (use mu0 for air_outer, not mu_kelvin)
    B_bounded_dict = {
        "magnetic": (mu_r * mu0) * H_pert_cf,
        "air_inner": mu0 * H_pert_cf,
        "air_outer": mu0 * H_pert_cf  # Bounded: use mu0 instead of mu_kelvin
    }
    flux = CoefficientFunction([B_bounded_dict[mat] for mat in mesh.GetMaterials()])

    # H(div) recovery using L2 projection (order - 1 for ZZ estimator)
    recovery_order = max(1, fes.globalorder - 1)
    fes_flux = HDiv(mesh, order=recovery_order)
    gf_flux = GridFunction(fes_flux)

    # L2 projection: solve (sigma, tau) = (flux, tau) for all tau in H(div)
    # This is the correct ZZ-type recovery (NOT interpolation with gf_flux.Set())
    sigma = fes_flux.TrialFunction()
    tau = fes_flux.TestFunction()
    a_flux = BilinearForm(fes_flux)
    a_flux += InnerProduct(sigma, tau) * dx
    a_flux.Assemble()

    f_flux = LinearForm(fes_flux)
    f_flux += InnerProduct(flux, tau) * dx
    f_flux.Assemble()

    gf_flux.vec.data = a_flux.mat.Inverse(fes_flux.FreeDofs(), inverse="sparsecholesky") * f_flux.vec

    # Unnormalized error (no normalization!)
    err = InnerProduct(flux - gf_flux, flux - gf_flux)
    element_errors = Integrate(err, mesh, element_wise=True)
    return element_errors

def mark_elements_by_threshold(element_errors, theta):
    """Mark elements with error >= theta * max_error."""
    max_error = max(element_errors)
    cutoff = theta * max_error
    return [i for i, err in enumerate(element_errors) if err >= cutoff]
```

**Parameters**:
- `theta = 0.3`: Refine elements with error >= 30% of max error (fixed threshold)
- `recovery_order = fes.globalorder - 1`: Use one order lower for ZZ recovery
- `max_elements_for_plot = 50000`: Skip mesh plot if elements exceed this limit

**Adaptive Theta Marking** (for ~2x DOF growth per iteration):

This algorithm dynamically adjusts theta to achieve approximately 2x DOF increase per iteration.

**Algorithm Overview**:
1. Calculate target number of elements to mark based on `expansion_factor`
2. Use binary search to find theta that marks approximately that many elements
3. Mark elements with error >= theta * max_error

**Key Formula**:
```
target_marked = current_ne * (target_ratio - 1) / expansion_factor
```

Where `expansion_factor = 28` is an empirical constant for 3D tetrahedral meshes that accounts for:
- Each marked element splits into 8 children (bisection refinement)
- Neighbor elements must also be refined for mesh conformity
- Total effect: marking n elements results in ~28n additional elements

```python
def mark_elements_adaptive_theta(element_errors, current_ne, target_ratio=2.0):
    """Mark elements with dynamically adjusted theta to achieve target DOF ratio.

    Uses binary search to find theta that marks approximately the right number
    of elements to achieve target_ratio increase in DOF.

    Args:
        element_errors: Array of element-wise error estimates
        current_ne: Current number of elements
        target_ratio: Target ratio for DOF increase (default 2.0)

    Returns:
        Tuple of (marked element indices, theta used)
    """
    max_error = max(element_errors)
    if max_error <= 0:
        return [], 1.0

    # Empirical model for 3D tetrahedral mesh refinement:
    # Due to mesh conformity (green/red refinement), marking n elements
    # typically results in many additional elements (including neighbors).
    # expansion_factor = 28 achieves ~2x DOF growth per iteration
    expansion_factor = 28.0
    target_marked = int(current_ne * (target_ratio - 1) / expansion_factor)
    target_marked = max(1, min(target_marked, current_ne))

    # Binary search for theta
    theta_low, theta_high = 0.0, 1.0
    best_theta = 0.5
    best_marked = []

    for _ in range(20):  # 20 iterations is enough for convergence
        theta = (theta_low + theta_high) / 2
        cutoff = theta * max_error
        marked = [i for i, err in enumerate(element_errors) if err >= cutoff]
        n_marked = len(marked)

        if n_marked == target_marked:
            return marked, theta
        elif n_marked < target_marked:
            theta_high = theta  # Need lower theta to mark more
            if n_marked > len(best_marked):
                best_marked = marked
                best_theta = theta
        else:
            theta_low = theta  # Need higher theta to mark fewer
            if abs(n_marked - target_marked) < abs(len(best_marked) - target_marked):
                best_marked = marked
                best_theta = theta

    # Return closest result
    cutoff = best_theta * max_error
    marked = [i for i, err in enumerate(element_errors) if err >= cutoff]
    return marked, best_theta
```

**Usage in Main Loop**:
```python
# Adaptive theta marking (dynamically adjust theta to achieve ~2x DOF increase)
target_dof_ratio = 2.0
marked, theta_used = mark_elements_adaptive_theta(element_errors, mesh.ne, target_dof_ratio)
max_err = max(element_errors)
print(f"  Marked {len(marked)}/{mesh.ne} elements (theta={theta_used:.3f})")

for el in mesh.Elements():
    mesh.SetRefinementFlag(el, False)
for el_nr in marked:
    mesh.SetRefinementFlag(ElementId(VOL, el_nr), True)

mesh.Refine()
mesh.Curve(order)
```

**Expected Results** (order=2, Sphere_3D):
| Iter | DOFs | Ratio | theta |
|------|------|-------|-------|
| 1 | 186 | - | - |
| 2 | 349 | 1.9x | 0.79 |
| 3 | 705 | 2.0x | 0.38 |
| 4 | 1618 | 2.3x | 0.46 |
| 5 | 3099 | 1.9x | 0.46 |
| 6 | 6930 | 2.2x | 0.43 |
| 7 | 13490 | 1.9x | 0.38 |
| 8 | 27961 | 2.1x | 0.41 |
| 9 | 55348 | 2.0x | 0.28 |
| 10 | 129728 | 2.3x | 0.22 |

---

## 7. Analytical Energy Formulas (3D Sphere)

### 7.1 Energy Partition by Region

For a magnetic sphere of radius `a` with Kelvin boundary at `R`:

| Region | Formula | Description |
|--------|---------|-------------|
| `magnetic` | `W_mag = 0.5 * mu_r * mu0 * H_pert² * V_sphere` | Energy inside sphere |
| `air_inner` | `W_air_in = mu0 * m² / (12*pi) * (1/a³ - 1/R³)` | Energy from a to R |
| `air_outer` | `W_air_out = mu0 * m² / (12*pi*R³)` | Energy from R to infinity |

Where:
- `H_pert = -(mu_r - 1)/(mu_r + 2) * H0` (perturbation field inside sphere)
- `V_sphere = (4/3)*pi*a³` (sphere volume)
- `m = 4*pi*a³ * (mu_r - 1)/(mu_r + 2) * H0` (dipole moment)

### 7.2 Implementation Code

```python
# Perturbation field inside magnetic sphere
H_pert_analytical = -(mu_r - 1) / (mu_r + 2) * H0
V_sphere = (4.0/3.0) * pi * sphere_radius**3

# Magnetic region energy
W_magnetic_analytical = 0.5 * mu_r * mu0 * H_pert_analytical**2 * V_sphere

# Dipole moment
m_dipole = 4 * pi * sphere_radius**3 * (mu_r - 1) / (mu_r + 2) * H0

# Split exterior energy into air_inner (a to R) and air_outer (R to infinity)
a = sphere_radius
R = kelvin_radius
W_air_inner_analytical = mu0 * m_dipole**2 / (12 * pi) * (1/a**3 - 1/R**3)
W_air_outer_analytical = mu0 * m_dipole**2 / (12 * pi * R**3)

# Total energy
W_analytical = W_magnetic_analytical + W_air_inner_analytical + W_air_outer_analytical
```

### 7.3 Energy Comparison Note

**IMPORTANT**: The simulation regions map to analytical regions as follows:
- `energy_magnetic` ↔ `W_magnetic_analytical`
- `energy_air_inner` ↔ `W_air_inner_analytical`
- `energy_air_outer` ↔ `W_air_outer_analytical`

Do NOT compare `energy_interior` (magnetic + air_inner) with `W_in_analytical` (magnetic only)!

---

## 8. Summary of Key Points

1. **Magnetic region**: `H_pert = grad(Omega) - Hs` (subtract source)
2. **Air regions**: `H_pert = grad(Omega)` in Reduced formulation (no subtraction)
3. **air_inner energy**: Use `grad(Orr) - grad(Oxr)` to account for boundary condition
4. **air_outer energy**: Use `mu_kelvin * |grad(Orr)|²` with mu_kelvin = (R/r')² * mu0
5. **Kelvin source field**: Use `(R/r')²` (not `-(r'/R)²`) for continuity at boundary
6. **Output**: VTK and PNG share common `output_data` dictionary
7. **Per-region energy tracking**: Track `energy_magnetic`, `energy_air_inner`, `energy_air_outer` separately
8. **Analytical comparison**: Compare each region's energy with corresponding analytical formula

---

## 9. Plot Structure (2x2 Convergence Plot)

### 9.1 Overview

The convergence plot is a 2×2 figure showing mesh, error distribution, and convergence history:

| Position | Content | Description |
|----------|---------|-------------|
| Top-left (ax1) | Interior domain mesh | All physical regions (magnetic + air_inner, or magnetic + air_total + air_inner for 4-region) |
| Top-right (ax2) | Exterior domain mesh | Kelvin-transformed region (air_outer only) |
| Bottom-left (ax3) | DOF vs Error | Log-log convergence plot with theoretical slope |
| Bottom-right (ax4) | DOF vs Energy | Per-region energy convergence |

### 9.2 Top-left: Interior Domain (ax1)

Shows all physical regions inside the Kelvin boundary:

**3-region structure** (magnetic + air_inner + air_outer):
```python
# Interior = magnetic + air_inner
if mat_name in ["magnetic", "air_inner"]:
```

**4-region structure** (magnetic + air_total + air_inner + air_outer):
```python
# Interior = magnetic + air_total + air_inner
if mat_name in ["magnetic", "air_total", "air_inner"]:
```

**Plot elements:**
- Mesh triangles colored by ZZ error (log10 scale)
- Cylinder/sphere boundary (red solid line)
- Total-Reduced interface (magenta dashed, for 4-region)
- Kelvin boundary (green dashed)

**Axis limits:**
```python
ax1.set_xlim(-0.05, kelvin_radius + 0.05)
ax1.set_ylim(-0.05, kelvin_radius + 0.05)
```

### 9.3 Top-right: Exterior Domain (ax2)

Shows Kelvin-transformed exterior region only:

```python
# Exterior = air_outer only
if mat_name == "air_outer":
```

**Important:** The exterior domain is offset in z-direction by `offset_z`:
```python
ax2.set_ylim(offset_z - 0.05, offset_z + kelvin_radius + 0.05)
ax2.plot(r_kelvin_plot, z_kelvin_plot + offset_z, 'g--', ...)  # Kelvin boundary at offset_z
```

### 9.4 Bottom-left: DOF vs Error (ax3)

Log-log plot of convergence:
```python
ax3.loglog(history['ndof'], history['error'], 'ko-', ...)
ax3.set_xlim(1e2, 1e5)
# Theoretical slope O(N^{-p/3}) for 3D
err_line = err_ref * (N_ref / ndof_line) ** (order / 3)
ax3.loglog(ndof_line, err_line, 'r--', label=f'$O(N^{{-{order}/3}})$')
```

### 9.5 Bottom-right: DOF vs Energy (ax4)

Per-region energy plot:

**3-region:**
```python
ax4.semilogx(history['ndof'], history['energy_magnetic'], 'rs-', label='magnetic')
ax4.semilogx(history['ndof'], history['energy_air_inner'], 'go-', label='air\\_inner')
ax4.semilogx(history['ndof'], history['energy_air_outer'], 'bo-', label='air\\_outer')
ax4.semilogx(history['ndof'], history['energy'], 'k^-', label='Total')
```

**4-region:**
```python
ax4.semilogx(history['ndof'], history['energy_magnetic'], 'rs-', label='magnetic')
ax4.semilogx(history['ndof'], history['energy_air_total'], 'mo-', label='air\\_total')
ax4.semilogx(history['ndof'], history['energy_air_inner'], 'go-', label='air\\_inner')
ax4.semilogx(history['ndof'], history['energy_air_outer'], 'bo-', label='air\\_outer')
ax4.semilogx(history['ndof'], history['energy'], 'k^-', label='Total')
```

### 9.6 Implementation Template

```python
def generate_convergence_plot(iter_num, history, output_data):
    mesh = output_data['mesh']
    element_errors = output_data['element_errors']

    fig = plt.figure(figsize=(14, 12), dpi=150)

    # Kelvin boundary circle for plots
    theta_circle = linspace(0, pi/2, 50)
    r_kelvin_plot = kelvin_radius * sin(theta_circle)
    z_kelvin_plot = kelvin_radius * cos(theta_circle)

    # Mesh size tolerance for y=0 plane detection
    y_tol = maxh_initial * 0.1
    max_elements_for_plot = 50000
    skip_mesh_plot = mesh.ne > max_elements_for_plot

    # ===== Top-left: Interior domain =====
    ax1 = plt.subplot(2, 2, 1)
    if not skip_mesh_plot:
        triangles_interior = []
        error_interior = []
        for el_idx, el in enumerate(mesh.Elements(VOL)):
            mat_name = el.mat
            if mat_name in ["magnetic", "air_inner"]:  # or ["magnetic", "air_total", "air_inner"]
                verts = [mesh[v].point for v in el.vertices]
                y_coords = [v[1] for v in verts]
                if min(y_coords) < y_tol:
                    face_verts = [v for v in verts if abs(v[1]) < y_tol]
                    if len(face_verts) >= 3:
                        xz_coords = [(v[0], v[2]) for v in face_verts[:3]]
                        triangles_interior.append(xz_coords)
                        error_interior.append(element_errors[el_idx])
        # Plot with PolyCollection...

    ax1.set_xlim(-0.05, kelvin_radius + 0.05)
    ax1.set_ylim(-0.05, kelvin_radius + 0.05)
    ax1.set_title('Interior domain (magnetic + air\\_inner) on $y=0$')

    # ===== Top-right: Exterior domain =====
    ax2 = plt.subplot(2, 2, 2)
    if not skip_mesh_plot:
        triangles_exterior = []
        error_exterior = []
        for el_idx, el in enumerate(mesh.Elements(VOL)):
            if el.mat == "air_outer":
                # Same extraction logic...
        # Plot with PolyCollection...

    ax2.plot(r_kelvin_plot, z_kelvin_plot + offset_z, 'g--', label='Kelvin boundary')
    ax2.set_xlim(-0.05, kelvin_radius + 0.05)
    ax2.set_ylim(offset_z - 0.05, offset_z + kelvin_radius + 0.05)
    ax2.set_title('Exterior domain (air\\_outer, Kelvin) on $y=0$')

    # ===== Bottom-left: DOF vs Error =====
    ax3 = plt.subplot(2, 2, 3)
    ax3.loglog(history['ndof'], history['error'], 'ko-', ...)
    ax3.set_xlim(1e2, 1e5)

    # ===== Bottom-right: DOF vs Energy =====
    ax4 = plt.subplot(2, 2, 4)
    ax4.semilogx(history['ndof'], history['energy_magnetic'], 'rs-', ...)
    ax4.semilogx(history['ndof'], history['energy_air_inner'], 'go-', ...)
    ax4.semilogx(history['ndof'], history['energy_air_outer'], 'bo-', ...)
    ax4.semilogx(history['ndof'], history['energy'], 'k^-', ...)
    ax4.set_xlim(1e2, 1e5)
```

---

## 10. Y=0 Cross-Section Detection for 3D Mesh Plotting

### 10.1 Problem

The original approach uses vertex positions to detect elements on the y=0 plane:
```python
y_tol = maxh_initial * 0.1  # Fixed tolerance
for el_idx, el in enumerate(mesh.Elements(VOL)):
    verts = [mesh[v].point for v in el.vertices]
    y_coords = [v[1] for v in verts]
    if min(y_coords) < y_tol:
        face_verts = [v for v in verts if abs(v[1]) < y_tol]
        if len(face_verts) >= 3:  # Need 3+ vertices on y=0
            ...
```

**Problem**: As the mesh is refined, elements get smaller and fewer have 3+ vertices exactly on y=0.

### 10.2 Solution: Plane Intersection Algorithm

The correct approach computes the **actual intersection** of each tetrahedron with the y=0 plane.

**Algorithm**:
1. For each tetrahedron, check if it intersects y=0 (some vertices above, some below or on y=0)
2. Find intersection points by computing where each edge crosses y=0
3. Sort intersection points by angle to form a proper polygon
4. The result is a triangle or quadrilateral representing the exact cross-section

```python
def compute_y0_cross_section(verts):
    """
    Compute the cross-section of a tetrahedron with the y=0 plane.
    Returns a list of (x, z) coordinates forming the cross-section polygon,
    or None if the tetrahedron doesn't intersect y=0.
    """
    y_coords = [v[1] for v in verts]
    y_min, y_max = min(y_coords), max(y_coords)

    # Check if tetrahedron intersects y=0 plane
    if y_min > 0 or y_max < 0:
        return None  # No intersection

    # Collect intersection points
    cross_points = []

    # Check each edge for intersection with y=0
    edges = [
        (0, 1), (0, 2), (0, 3),
        (1, 2), (1, 3), (2, 3)
    ]

    for i, j in edges:
        y0, y1 = verts[i][1], verts[j][1]

        # Check if edge crosses y=0
        if (y0 <= 0 <= y1) or (y1 <= 0 <= y0):
            if abs(y1 - y0) < 1e-12:
                # Edge is on y=0 plane, add both endpoints
                cross_points.append((verts[i][0], verts[i][2]))
                cross_points.append((verts[j][0], verts[j][2]))
            else:
                # Interpolate to find intersection point
                t = -y0 / (y1 - y0)
                x = verts[i][0] + t * (verts[j][0] - verts[i][0])
                z = verts[i][2] + t * (verts[j][2] - verts[i][2])
                cross_points.append((x, z))

    # Remove duplicate points
    unique_points = []
    for p in cross_points:
        is_dup = False
        for q in unique_points:
            if abs(p[0] - q[0]) < 1e-10 and abs(p[1] - q[1]) < 1e-10:
                is_dup = True
                break
        if not is_dup:
            unique_points.append(p)

    if len(unique_points) < 3:
        return None

    # Sort points by angle around centroid for proper polygon ordering
    cx = sum(p[0] for p in unique_points) / len(unique_points)
    cz = sum(p[1] for p in unique_points) / len(unique_points)

    def angle_key(p):
        from math import atan2
        return atan2(p[1] - cz, p[0] - cx)

    unique_points.sort(key=angle_key)

    return unique_points

# Usage in mesh plotting
polygons = []
errors = []

for el_idx, el in enumerate(mesh.Elements(VOL)):
    if el.mat not in target_materials:
        continue

    verts = [mesh[v].point for v in el.vertices]
    cross_section = compute_y0_cross_section(verts)
    if cross_section is not None:
        polygons.append(cross_section)
        errors.append(element_errors[el_idx])
```

### 10.3 Advantages

1. **Exact intersection**: Computes true y=0 cross-section, not approximation
2. **No gaps**: Every element crossing y=0 is captured
3. **Handles all cases**: Works for triangular and quadrilateral cross-sections
4. **Mesh refinement independent**: Works at any mesh resolution

### 10.4 Implementation Notes

- Edge intersection uses linear interpolation: `t = -y0 / (y1 - y0)`
- Points are sorted by angle to form a proper polygon for PolyCollection
- Cross-section can be 3 or 4 points (triangle or quadrilateral)

---

## Files Using This Implementation

- `Cylinder_3D/order=2/Refine_with_zz_estimator/Cylinder_3D_adaptive_with_Kelvin.py`
- `Cylinder_3D/order=2/Refine_all_elements/Cylinder_3D_adaptive_with_Kelvin.py`
- `Cylinder_3D/order=2/metric_based/Cylinder_3D_metric_with_Kelvin.py`
- `Cylinder_3D/order=3/*/Cylinder_3D_*.py`
- `Cylinder_3D/order=4/*/Cylinder_3D_*.py`
- `Sphere_3D/order=2/*/Sphere_3D_*.py`
- `Sphere_3D/order=3/*/Sphere_3D_*.py`
- `Sphere_3D/order=4/*/Sphere_3D_*.py`
