"""
Axisymmetric Omega-Reduced Omega Method for Magnetostatics with Kelvin Transformation
Problem: Magnetic sphere (mu_r=100) in uniform z-directed background field

Based on Omega_ReducedOmega.py implementation:
- Sign convention: B = mu * grad(Omega), H = grad(Omega)
- Source potential: Omega_s = H0 * z (so grad(Omega_s) = H_s)
- Total region (magnetic sphere): No source term in weak form
- Reduced region (air): Source term from Omega_s

Kelvin transformation:
- Maps infinite exterior domain to finite half-circle
- Permeability transformation: mu'(rho') = (R/rho')^2 * mu0
- Periodic BC couples interior (rho=R) with exterior (rho'=R)

Analytical solution:
- Inside sphere: Hz = -3/(mu_r + 2) * H0 (negative, demagnetization effect)
- Outside sphere: dipole + uniform field

Author: Claude Code
Date: 2025-12-23
Updated: Following Omega_ReducedOmega.py implementation exactly
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Axisymmetric Omega-Reduced Omega Method with Kelvin Transform")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
sphere_radius = 0.5      # Magnetic sphere radius [m]
kelvin_radius = 1.0      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1.0  # [A/m]

# Mesh parameters
maxh = 0.03              # Mesh size
fe_order = 3             # Finite element order

# Offset for exterior domain (z-direction for axisymmetric)
offset_z = 3.0

print(f"\nProblem parameters:")
print(f"  Sphere radius: {sphere_radius} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")

# ============================================================
# Geometry Definition using HALF-CIRCLES (axisymmetric, r >= 0)
# x = r (radial), y = z (axial)
# ============================================================
print("\nCreating half-circle geometry with periodic boundary conditions...")

print(f"Using Kelvin transformation with periodic BC:")
print(f"  - Inner domain: 0 <= r < R = {kelvin_radius} m at (0, 0)")
print(f"  - Outer domain: 0 <= r' < R = {kelvin_radius} m at (0, {offset_z})")
print(f"  - Transformation radius R = {kelvin_radius} m")

# ===== INTERIOR DOMAIN (half-circle, r >= 0) =====
outer_half_int = Circle((0, 0), kelvin_radius).Face()
cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_int = outer_half_int - cutter_int  # Keep x >= 0 part

inner_half_int = Circle((0, 0), sphere_radius).Face()
inner_half_int = inner_half_int - cutter_int  # Keep x >= 0 part

# Air region (half-annulus)
air_inner = outer_half_int - inner_half_int

# Name boundaries
for edge in air_inner.edges:
    x_center = edge.center.x
    dist = sqrt(edge.center.x**2 + edge.center.y**2)
    if x_center < 1e-6:  # On z-axis (r = 0)
        edge.name = "axis_int"
    elif abs(dist - kelvin_radius) < kelvin_radius * 0.2:
        edge.name = "kelvin_int"
    else:
        edge.name = "sphere"
air_inner.faces.name = "air_inner"

# Magnetic material (half-circle)
for edge in inner_half_int.edges:
    if edge.center.x < 1e-6:
        edge.name = "axis_int"
    else:
        edge.name = "sphere"
inner_half_int.faces.name = "magnetic"

# ===== EXTERIOR DOMAIN (half-circle, r' >= 0, offset in z-direction) =====
outer_half_ext = Circle((0, offset_z), kelvin_radius).Face()
cutter_ext = MoveTo(-kelvin_radius - 1, offset_z - kelvin_radius - 1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_ext = outer_half_ext - cutter_ext  # Keep x' >= 0 part (r' >= 0)

for edge in outer_half_ext.edges:
    r_center = edge.center.x  # r' coordinate
    if r_center < 1e-6:  # On axis (r' = 0)
        edge.name = "axis_ext"
    else:
        edge.name = "kelvin_ext"
outer_half_ext.faces.name = "air_outer"

# ===== GND VERTEX (center of exterior domain - at infinity) =====
vertex = Vertex(Pnt(0, offset_z, 0))
vertex.name = "GND"

# Glue all domains
shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

# ===== IDENTIFY ALL PERIODIC BOUNDARY EDGE PAIRS =====
print("\nIdentifying periodic boundaries...")

# Find all kelvin edges
kelvin_int_edges = []
kelvin_ext_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_int_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_ext_edges.append(edge)

print(f"  Found {len(kelvin_int_edges)} interior kelvin edges")
print(f"  Found {len(kelvin_ext_edges)} exterior kelvin edges")

# Print edge details
for i, edge in enumerate(kelvin_int_edges):
    print(f"    kelvin_int[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")
for i, edge in enumerate(kelvin_ext_edges):
    print(f"    kelvin_ext[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# Match edges for z-offset Kelvin
# Inner z>0 <-> Outer z>offset_z, Inner z<0 <-> Outer z<offset_z
if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    matched_pairs = 0
    for int_edge in kelvin_int_edges:
        int_z = int_edge.center.y  # z coordinate
        # Find matching exterior edge
        for ext_edge in kelvin_ext_edges:
            ext_z = ext_edge.center.y - offset_z  # Relative z from center
            # Match: int_z > 0 with ext_z > 0, int_z < 0 with ext_z < 0
            if (int_z > 0 and ext_z > 0) or (int_z < 0 and ext_z < 0):
                int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                print(f"  Identified: int(z={int_z:.3f}) <-> ext(z-offset={ext_z:.3f})")
                matched_pairs += 1
                break
    print(f"  Total matched pairs: {matched_pairs}")

# Create geometry
geo = OCCGeometry(shape, dim=2)

print(f"\nGeometry created:")
print(f"  Magnetic sphere radius: {sphere_radius} m at (0, 0)")
print(f"  Inner air domain: half-circle of radius {kelvin_radius} m")
print(f"  Outer air domain: half-circle of radius {kelvin_radius} m at (0, {offset_z})")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.7))

# Apply curved elements for high-order accuracy
mesh.Curve(fe_order)

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")
print(f"  Curved elements: order {fe_order}")

# ============================================================
# Finite Element Space
# ============================================================
print("\nSetting up finite element space...")

# H1 space with Dirichlet BC at GND (infinity)
fes_before = H1(mesh, order=fe_order, dirichlet_bbnd="GND")
fes = Periodic(fes_before)  # Apply periodic BC

# Check if Periodic BC is working (FreeDofs should decrease)
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  Finite element order: {fe_order}")
print(f"  Number of DOFs: {fes.ndof}")
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")

if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

# Trial and test functions
Omega = fes.TrialFunction()
psi = fes.TestFunction()

# Coordinate functions (x = r, y = z)
r_coord = x  # Radial coordinate
z_coord = y  # Axial coordinate

# ============================================================
# Material Properties
# ============================================================
print("\nSetting up material properties...")

# Kelvin-transformed permeability (Nagamine CEFC 2026 canonical):
#   mu_ext = mu_0 * (R/rho')^2  for 3D spherical (conformal) Kelvin
# rho' = sqrt(r^2 + (z - offset_z)^2) in axisym (r,z) meridional plane.
# See examples/kelvin_transformation/CONVENTION.md.
from radia.kelvin_source import kelvin_mu_factor_axisym_cf, build_material_cf

mu_kelvin_factor = kelvin_mu_factor_axisym_cf(
    z_offset=offset_z, R=kelvin_radius,
    r_coord=r_coord, z_coord=z_coord,
)
Mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)

print(f"  air_inner: mu = mu0")
print(f"  air_outer: mu = (R/rho')^2 * mu0 [Nagamine CEFC 2026]")
print(f"  magnetic: mu = {mu_r} * mu0")

# ============================================================
# Following Omega_ReducedOmega.py implementation
# ============================================================
print("\nSetting up weak form (following Omega_ReducedOmega.py)...")

# Detect which domain we're in (offset in z-direction)
is_exterior = IfPos(z_coord - offset_z/2, 1.0, 0.0)

# r-weight for axisymmetric formulation
# For z-offset Kelvin: r' = r (same for both domains)
r_weight_inner = IfPos(r_coord - 1e-10, r_coord, 1e-10)
r_weight_outer = r_weight_inner  # r' = r for z-offset Kelvin
r_weight = r_weight_inner  # Same r weight for both domains

# ===== Source Field and Potential (Omega_ReducedOmega.py convention) =====
# Convention: B = mu * grad(Omega), H = grad(Omega)
# Source potential: Omega_s = H0 * z (so grad(Omega_s) = (0, H0) = H_s)

# Source potential Omega_s (line 50 of Omega_ReducedOmega.py: Ov = Ofield(coil))
# For uniform field H0 in z-direction: Omega_s = H0 * z
Omega_s = H0 * z_coord

# Source magnetic field H_s = (0, H0) and B_s = mu0 * H_s
Hs = CoefficientFunction((0.0, H0))
Bs = CoefficientFunction((0.0, mu0 * H0))

# For Kelvin-transformed exterior domain (z-offset):
# H_s' = -(rho'/R)^2 * H0 (field reversal under Kelvin transform)
rho_prime = sqrt(r_coord**2 + (z_coord - offset_z)**2 + 1e-20)
Hz_exterior = -(rho_prime / kelvin_radius)**2 * H0
Hs_exterior = CoefficientFunction((0.0, Hz_exterior))
Bs_exterior = CoefficientFunction((0.0, mu0 * Hz_exterior))

# ===== Bilinear Form (lines 78-88) =====
# a += Mu*(grad(omega)*grad(psi))*dx(total_region)
# a += Mu*(grad(omega)*grad(psi))*dx(reduced_region)
# For Kelvin: a += Mu*fac*(grad(omega)*grad(psi))*dx("Kelvin")

a = BilinearForm(fes)
# Total region (magnetic sphere)
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("magnetic")
# Reduced region (air_inner)
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_inner")
# Kelvin region (air_outer) with Kelvin factor already in Mu
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_outer")

a.Assemble()

# Dirichlet BC on Total/Reduced interface: Omega = Omega_s on sphere boundary
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

# Linear form
f = LinearForm(fes)
# Source term in Reduced region only
f += Mu * grad(gfOmega) * grad(psi) * r_weight * dx("air_inner")
f.Assemble()

# Neumann boundary condition on Total/Reduced interface
normal = specialcf.normal(mesh.dim)
f += (normal * Bs) * psi * r_weight_inner * ds("sphere")
f.Assemble()

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
gfu = gfOmega

print("  Solution converged")

# ============================================================
# Post-processing (following Omega_ReducedOmega.py lines 118-143)
# ============================================================
print("\nPost-processing (following Omega_ReducedOmega.py)...")

# Following lines 118-143:
# fesOt = H1(mesh, order=feOrder, definedon=total_region)
# fesOr = H1(mesh, order=feOrder, definedon=reduced_region+"|Kelvin")
# Ot = GridFunction(fesOt)
# Orr = GridFunction(fesOr)
# Oxr = GridFunction(fesOr)
#
# Ot.Set(gfOmega, VOL, definedon=total_region)
# Orr.Set(gfOmega, VOL, definedon=reduced_region+"|Kelvin")
# Oxr.Set(Ov, BND, mesh.Boundaries(total_boundary))
#
# Bt = grad(Ot) * Mu
# Br = (grad(Orr) - grad(Oxr)) * mu0
# BField = Bt + Br + Bs

fesOt = H1(mesh, order=fe_order, definedon="magnetic")
fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")

Ot = GridFunction(fesOt)
Orr = GridFunction(fesOr)
Oxr = GridFunction(fesOr)

Ot.Set(gfu, VOL, definedon="magnetic")
Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))

# B field computation (lines 137, 142, 143)
# Bt = grad(Ot) * Mu  (Total region - NO minus sign!)
# Br = (grad(Orr) - grad(Oxr)) * mu0  (Reduced region perturbation)
# BField = Bt + Br + Bs

Bt = grad(Ot) * Mu  # Total region
Br = (grad(Orr) - grad(Oxr)) * mu0  # Reduced region perturbation

# Source field Bs is only in Reduced region
Bs_dict = {
    "air_inner": Bs,
    "air_outer": Bs_exterior,
    "magnetic": CoefficientFunction((0.0, 0.0))
}
Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

# Total B field
BField = Bt + Br + Bs_cf

# H field: H = B/mu
HField = BField / Mu

# ============================================================
# Analytical Solution
# ============================================================
# For a magnetic sphere (mu_r) in uniform external field H0:
# Inside sphere: H_in = 3/(mu_r + 2) * H0
#
# Note: The field INSIDE the sphere is in the SAME direction as H0,
# but reduced in magnitude. This is because the sphere concentrates
# flux lines, increasing B but reducing H inside.
#
# The "demagnetization" interpretation:
# H_in = H0 - N*M where N=1/3 for sphere
# For high mu_r: H_in ~= 3*H0/mu_r (very small)

Hz_analytical_interior = 3.0 / (mu_r + 2) * H0  # = 0.029412 A/m for mu_r=100

print("\n" + "=" * 60)
print("VALIDATION RESULTS")
print("=" * 60)
print(f"\nAnalytical solution:")
print(f"  Interior Hz = 3/(mu_r+2) * H0 = {Hz_analytical_interior:.6f} A/m")
print(f"  (Same direction as H0, but reduced magnitude)")

# Evaluate at interior points
print(f"\nInterior (magnetic sphere - Total potential region):")
print(f"  In Total region: H = grad(Omega) (Omega_ReducedOmega.py convention)")
print()


for r_val in [0.1, 0.2, 0.3, 0.4]:
    try:
        mip = mesh(r_val, 0)
        # In Total region: H = grad(Omega)
        Hz_numerical = grad(gfu)[1](mip)
        error = abs(Hz_numerical - Hz_analytical_interior) / abs(Hz_analytical_interior) * 100
        print(f"  r={r_val}: Hz={Hz_numerical:.6f}, analytical={Hz_analytical_interior:.6f}, error={error:.2f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")

# Evaluate at exterior points
print(f"\nExterior (air - Reduced potential region):")
print(f"  Following Omega_ReducedOmega.py: H = H_s + grad(Omega_r) = H_s + (grad(Omega) - grad(Omega_s))")
print()

for r_val in [0.6, 0.7, 0.8, 0.9]:
    # Analytical: Hz = H0 * (1 - (mu_r-1)/(mu_r+2) * (a/r)^3)
    Hz_analytical_ext = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r_val)**3)
    try:
        mip = mesh(r_val, 0)
        # For Reduced region: H = H_s + (grad(Omega) - grad(Omega_s))
        # = H0 + grad(Omega)[1] - H0 = grad(Omega)[1]
        # Wait, that's not right. Let's use B/mu instead
        B_val = BField(mip)
        Hz_numerical = B_val[1] / mu0
        error = abs(Hz_numerical - Hz_analytical_ext) / abs(Hz_analytical_ext) * 100
        print(f"  r={r_val}: Hz={Hz_numerical:.6f}, analytical={Hz_analytical_ext:.6f}, error={error:.2f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")

# ============================================================
# Profile along r-axis
# ============================================================
print("\nComputing r-axis profile...")

r_profile = linspace(0.02, kelvin_radius - 0.02, 100)
Hz_numerical = zeros(len(r_profile))
Hz_analytical_profile = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    try:
        mip = mesh(r_val, 0)
        if r_val < sphere_radius:
            # Total region: H = grad(Omega)
            Hz_numerical[i] = grad(gfu)[1](mip)
        else:
            # Reduced region: use BField / mu0
            B_val = BField(mip)
            Hz_numerical[i] = B_val[1] / mu0
    except:
        Hz_numerical[i] = nan

    # Analytical solution
    if r_val < sphere_radius:
        Hz_analytical_profile[i] = Hz_analytical_interior
    else:
        Hz_analytical_profile[i] = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r_val)**3)

# Compute error statistics
valid_idx = ~isnan(Hz_numerical)
interior_idx = valid_idx & (r_profile < sphere_radius)
exterior_idx = valid_idx & (r_profile >= sphere_radius)

print(f"\n  Validation results (r-axis, total field Hz):")
print(f"  " + "-" * 40)

if sum(interior_idx) > 0:
    interior_error = Hz_numerical[interior_idx] - Hz_analytical_profile[interior_idx]
    max_err_int = abs(interior_error).max()
    rms_err_int = sqrt(mean(interior_error**2))
    rel_err_int = rms_err_int / abs(Hz_analytical_profile[interior_idx][0]) * 100
    print(f"  Interior (r < {sphere_radius} m):")
    print(f"    Max error: {max_err_int:.6e} A/m")
    print(f"    RMS error: {rms_err_int:.6e} A/m ({rel_err_int:.3f}%)")

if sum(exterior_idx) > 0:
    exterior_error = Hz_numerical[exterior_idx] - Hz_analytical_profile[exterior_idx]
    max_err_ext = abs(exterior_error).max()
    rms_err_ext = sqrt(mean(exterior_error**2))
    print(f"  Exterior (r >= {sphere_radius} m):")
    print(f"    Max error: {max_err_ext:.6e} A/m")
    print(f"    RMS error: {rms_err_ext:.6e} A/m")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                             'bf': 'serif:bold', 'fontset': 'cm'})

fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

# Plot 1: R-axis profile (total field)
ax1 = axes[0, 0]
ax1.plot(r_profile, Hz_numerical, 'b-', linewidth=2, label='NGSolve (Omega method)')
ax1.plot(r_profile, Hz_analytical_profile, 'r--', linewidth=1.5, label='Analytical')
ax1.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax1.axhline(0, color='gray', linestyle='-', alpha=0.3)
ax1.set_xlabel('$r$ (m)', fontsize=11)
ax1.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax1.set_title('R-axis Profile (Total Field)', fontsize=12)
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: Error along r-axis
ax2 = axes[0, 1]
if sum(valid_idx) > 0:
    with errstate(divide='ignore', invalid='ignore'):
        error_profile = abs(Hz_numerical[valid_idx] - Hz_analytical_profile[valid_idx]) / abs(Hz_analytical_profile[valid_idx]) * 100
    error_profile = where(isfinite(error_profile), error_profile, 0)
    ax2.semilogy(r_profile[valid_idx], error_profile + 1e-10, 'b-', linewidth=2)
ax2.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax2.set_xlabel('$r$ (m)', fontsize=11)
ax2.set_ylabel('Relative Error (%)', fontsize=11)
ax2.set_title('Error Distribution (R-axis)', fontsize=12)
ax2.grid(True, alpha=0.3)

# Plot 3: 2D total H field in interior domain
ax3 = axes[1, 0]

# Create 2D grid for visualization
r_grid = linspace(0.02, kelvin_radius - 0.02, 50)
z_grid = linspace(-kelvin_radius + 0.02, kelvin_radius - 0.02, 50)
[rr, zz] = meshgrid(r_grid, z_grid)

Hr_field = zeros(rr.shape)
Hz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        r_dist = sqrt(r_grid[ir]**2 + z_grid[iz]**2)
        if r_dist < kelvin_radius - 0.02:
            try:
                mip = mesh(r_grid[ir], z_grid[iz])
                if r_dist < sphere_radius:
                    # Total region: H = grad(Omega)
                    Hr_field[iz, ir] = grad(gfu)[0](mip)
                    Hz_field[iz, ir] = grad(gfu)[1](mip)
                else:
                    # Reduced region: use BField / mu0
                    B_val = BField(mip)
                    Hr_field[iz, ir] = B_val[0] / mu0
                    Hz_field[iz, ir] = B_val[1] / mu0
            except:
                Hr_field[iz, ir] = nan
                Hz_field[iz, ir] = nan
        else:
            Hr_field[iz, ir] = nan
            Hz_field[iz, ir] = nan

strm = ax3.streamplot(rr, zz, Hr_field, Hz_field,
                      color='black', linewidth=1.0, density=1.5,
                      arrowsize=0.8, arrowstyle='->')

# Draw sphere and Kelvin boundaries
theta_circle = linspace(0, pi, 100)
r_sphere = sphere_radius * sin(theta_circle)
z_sphere = sphere_radius * cos(theta_circle)
ax3.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
ax3.plot(r_sphere, z_sphere, 'r-', linewidth=2, label='Magnetic sphere')

r_kelvin_plot = kelvin_radius * sin(theta_circle)
z_kelvin_plot = kelvin_radius * cos(theta_circle)
ax3.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')

ax3.set_xlabel('$r$ (m)', fontsize=11)
ax3.set_ylabel('$z$ (m)', fontsize=11)
ax3.set_title('Total Field $\\mathbf{H}$ (Omega-Reduced Omega)', fontsize=12)
ax3.legend(loc='upper right', fontsize=8)
ax3.set_aspect('equal')
ax3.set_xlim(0, kelvin_radius + 0.1)
ax3.set_ylim(-kelvin_radius - 0.1, kelvin_radius + 0.1)

# Plot 4: Comparison table
ax4 = axes[1, 1]
ax4.axis('off')

# Create comparison text
comparison_text = f"""
Omega-Reduced Omega Method (Axisymmetric) + Kelvin
==================================================

Formulation (following Omega_ReducedOmega.py):
  B = mu * grad(Omega)  (no minus sign)
  Total region (magnetic): no source term
  Reduced region (air): source from Omega_s

Problem: Magnetic sphere in uniform z-field

Parameters:
  Sphere radius: {sphere_radius} m
  Kelvin radius: {kelvin_radius} m
  mu_r = {mu_r}
  H0 = {H0} A/m

Results at r=0.1:
  Analytical Hz: {Hz_analytical_interior:.6f} A/m
  Numerical Hz: {Hz_numerical[5]:.6f} A/m
  (Interior field is negative - demagnetization)

Periodic BC verification:
  FreeDofs before: {freedof_before}
  FreeDofs after: {freedof_after}
  Reduction: {freedof_before - freedof_after}
"""
ax4.text(0.05, 0.95, comparison_text, transform=ax4.transAxes,
         fontsize=10, fontfamily='monospace', verticalalignment='top')

plt.tight_layout()

# Save figure
png_file = os.path.splitext(__file__)[0] + ".png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved: {png_file}")

# Open the plot
try:
    os.startfile(png_file)
except:
    pass

print("\n" + "=" * 60)
print("Computation completed")
print("=" * 60)
