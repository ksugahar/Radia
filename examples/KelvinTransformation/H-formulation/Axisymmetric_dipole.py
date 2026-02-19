"""
Axisymmetric H-formulation for magnetostatics (without Kelvin transformation)
Validates against 3D_dipole.py analytical solution

Problem: Magnetic sphere in uniform z-directed background field
- Sphere radius: 0.5 m
- Outer domain radius: 4.0 m
- Relative permeability: mu_r = 100
- Background field: H_s = [0, 0, 1] A/m (z-direction)

Analytical solution (interior):
  H_z,pert = -1 + 3/(mu_r + 2) = -0.970588 A/m for mu_r = 100

Axisymmetric formulation:
  In (r, z) coordinates with azimuthal symmetry:
  - 2D computational domain represents r-z cross-section
  - Volume element: dV = 2*pi*r dr dz
  - Weak form includes r-weighting
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("="*60)
print("Axisymmetric H-formulation (Without Kelvin Transformation)")
print("="*60)

# ============================================================
# Parameters (same as 3D_dipole.py)
# ============================================================
sphere_radius = 0.5    # Magnetic sphere radius [m]
outer_radius = 1.0     # Outer domain radius [m] (same as 3D_dipole.py)
maxh_fine = 0.03       # Fine mesh size [m] (same as 3D_dipole.py)
plot_range = 1.1       # Plot range [m]
mu_r = 100             # Relative permeability

mu0 = 4*pi*1e-7        # Vacuum permeability [H/m]

# ============================================================
# Geometry Definition using OCC (2D axisymmetric: r-z plane)
# x = r (radial), y = z (axial)
# Use SEMICIRCULAR domain (half-disk) to match 3D spherical problem
# ============================================================
print("\nCreating 2D axisymmetric geometry using OCC...")

# Create outer domain (half-disk in r-z plane, r >= 0)
outer_disk = Circle((0, 0), outer_radius).Face()
# Cut to keep only r >= 0 half
cutter = MoveTo(-outer_radius-1, -outer_radius-1).Rectangle(outer_radius+1, 2*outer_radius+2).Face()
outer_half = outer_disk - cutter

# Create inner sphere (half-disk)
inner_disk = Circle((0, 0), sphere_radius).Face()
inner_half = inner_disk - cutter

# Create air region (outer - inner)
air = outer_half - inner_half

# Set boundary and material names BEFORE Glue
# Air region boundaries
for edge in air.edges:
    r_center = edge.center.x
    dist_from_origin = sqrt(edge.center.x**2 + edge.center.y**2)
    if r_center < 1e-6:  # On axis (r=0)
        edge.name = "axis"
    elif dist_from_origin > (outer_radius + sphere_radius) / 2:  # Outer boundary (closer to outer)
        edge.name = "outer"
    else:  # Inner boundary (sphere)
        edge.name = "sphere"
air.faces.name = "air"

# Magnetic region boundaries
for edge in inner_half.edges:
    if edge.center.x < 1e-6:  # On axis (r=0)
        edge.name = "axis"
    else:  # Sphere boundary
        edge.name = "sphere"
inner_half.faces.name = "magnetic"

# Combine shapes
shape = Glue([air, inner_half])

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_fine))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Axisymmetric H-formulation Setup
# ============================================================
print("\nSetting up axisymmetric H-formulation...")

# Finite element space
fes = H1(mesh, order=3)

print(f"  Number of DOFs: {fes.ndof}")

u = fes.TrialFunction()
v = fes.TestFunction()

# Coordinate functions
# In NGSolve 2D: x is first coord, y is second coord
# For axisymmetric: x -> r (radial), y -> z (axial)
r_coord = x  # Radial coordinate
z_coord = y  # Axial coordinate

# Material properties
# Use CoefficientFunction list (like 3D_dipole.py) for proper boundary evaluation
mu_d = {"air": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])
mu_boundary = BoundaryFromVolumeCF(mu)  # For boundary integrals

# Background field
# H_s = (0, 0, 1) in 3D -> H_z = 1 in axisymmetric
# In 2D axisymmetric, the background field is represented as (H_r, H_z) = (0, 1)
Hs = CF((0, 1))
Hsb = BoundaryFromVolumeCF(Hs)

print(f"  Background field: H_s = (0, 1) in (r,z)")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Axisymmetric Weak Form
# ============================================================
# For axisymmetric problems, the weak form includes r-weighting:
#   int mu * grad(u) * grad(v) * r dr dz = int mu * H_s * grad(v) * r dr dz - int mu * v * (n*H_s) * r ds
#
# This is because the 3D integral becomes:
#   int_3D = int_2D * 2*pi*r dr dz

print("\nAssembling axisymmetric system...")

# Normal vector
n = specialcf.normal(mesh.dim)

# Use r-weighting for proper axisymmetric formulation
use_r_weighting = True  # Must be True to match 3D_dipole.py results

# For semicircular (half-disk) domain:
# The axisymmetric weak form with r-weighting is:
#   a(u,v) = ∫ μ grad(u)·grad(v) r dr dz
#   f(v)   = ∫ μ Hs·grad(v) r dr dz - ∫ μ v (n·Hs) r ds
# The boundary integral on the semicircular "outer" boundary is essential.

if use_r_weighting:
    # Bilinear form: a(u,v) = ∫ μ grad(u)·grad(v) r dr dz
    a = BilinearForm(fes)
    a += mu * grad(u) * grad(v) * r_coord * dx

    # Linear form: f(v) = ∫ μ Hs·grad(v) r dr dz - ∫ μ v (n·Hs) r ds
    f = LinearForm(fes)
    f += mu * InnerProduct(grad(v), Hs) * r_coord * dx
    # Boundary term on outer boundary (use mu_boundary for boundary integrals)
    f += -mu_boundary * v * InnerProduct(n, Hsb) * r_coord * ds(mesh.Boundaries("outer"))
else:
    # 2D Cartesian formulation (without r-weighting, for debugging)
    a = BilinearForm(fes)
    a += mu * grad(u) * grad(v) * dx

    f = LinearForm(fes)
    f += mu * InnerProduct(grad(v), Hs) * dx
    # Boundary term on outer boundary (use mu_boundary for boundary integrals)
    f += -mu_boundary * v * InnerProduct(n, Hsb) * ds(mesh.Boundaries("outer"))

a.Assemble()
f.Assemble()

print("  System assembled")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

print("  Solution converged")

# Check potential values along r-axis (z=0)
print("\nPotential values along r-axis (z=0):")
try:
    phi_origin = gfu(mesh(0.01, 0))
    phi_03 = gfu(mesh(0.3, 0))
    phi_07 = gfu(mesh(0.7, 0))
    phi_35 = gfu(mesh(3.5, 0))
    print(f"  phi(0.01, 0) = {phi_origin:.6f}")
    print(f"  phi(0.3, 0) = {phi_03:.6f}")
    print(f"  phi(0.7, 0) = {phi_07:.6f}")
    print(f"  phi(3.5, 0) = {phi_35:.6f}")
except Exception as e:
    print(f"  Error: {e}")

# Check potential along z-axis (r=0.1)
print("\nPotential values along z-axis (r=0.1):")
try:
    phi_zm3 = gfu(mesh(0.1, -3))
    phi_zm1 = gfu(mesh(0.1, -1))
    phi_z0 = gfu(mesh(0.1, 0))
    phi_zp1 = gfu(mesh(0.1, 1))
    phi_zp3 = gfu(mesh(0.1, 3))
    print(f"  phi(0.1, -3) = {phi_zm3:.6f}")
    print(f"  phi(0.1, -1) = {phi_zm1:.6f}")
    print(f"  phi(0.1,  0) = {phi_z0:.6f}")
    print(f"  phi(0.1, +1) = {phi_zp1:.6f}")
    print(f"  phi(0.1, +3) = {phi_zp3:.6f}")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# Compute perturbation field: H_pert = -grad(phi)
H_pert = -grad(gfu)

# Evaluate at (r=0.01, z=0) - near axis, inside magnetic sphere
try:
    origin_point = mesh(0.01, 0)
    Hz_origin = H_pert[1](origin_point)
    print(f"  Field near origin (r=0.01, z=0): Hz = {Hz_origin:.6f} A/m")
except Exception as e:
    print(f"  Error at origin: {e}")
    Hz_origin = nan

# Evaluate at (r=0.3, z=0) - inside magnetic sphere
try:
    interior_point = mesh(0.3, 0)
    Hz_interior = H_pert[1](interior_point)
    print(f"  Field at (r=0.3, z=0): Hz = {Hz_interior:.6f} A/m")
except Exception as e:
    print(f"  Error at (0.3, 0): {e}")
    Hz_interior = nan

# Analytical solution (interior)
Hz_analytical = -1.0 + 3.0/(mu_r + 2)  # = -0.970588 for mu_r=100
print(f"  Analytical (interior): Hz = {Hz_analytical:.6f} A/m")

if not isnan(Hz_interior):
    rel_error = abs(Hz_interior - Hz_analytical) / abs(Hz_analytical) * 100
    print(f"  Relative error: {rel_error:.3f}%")

# Evaluate at (r=0.7, z=0) - outside magnetic sphere, in air
try:
    air_point = mesh(0.7, 0)
    Hz_air = H_pert[1](air_point)
    # Analytical for exterior (on r-axis): Hz = -(mu_r-1)/(mu_r+2) * (a/r)^3
    Hz_air_analytical = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/0.7)**3
    print(f"  Field at (r=0.7, z=0): Hz = {Hz_air:.6f} A/m")
    print(f"  Analytical: Hz = {Hz_air_analytical:.6f} A/m")
    print(f"  Relative error: {abs(Hz_air - Hz_air_analytical)/abs(Hz_air_analytical)*100:.3f}%")
except Exception as e:
    print(f"  Error at (0.7, 0): {e}")

# ============================================================
# Profile Comparisons with Analytical Solution
# (Perturbation field Hz component along r-axis and z-axis)
# ============================================================
print("\nComputing axis profiles (perturbation field Hz)...")

# Sample points along r-axis (z=0) and z-axis (r=0.01)
profile_range = linspace(-plot_range, plot_range, 221)

# R-axis profile (z=0)
r_profile = linspace(0.02, outer_radius - 0.1, 200)
Hz_numerical = zeros(len(r_profile))
Hz_analytical_profile = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    try:
        mip = mesh(r_val, 0)
        Hz_numerical[i] = H_pert[1](mip)
    except:
        Hz_numerical[i] = nan

    # Analytical solution
    if r_val < sphere_radius:
        Hz_analytical_profile[i] = -1.0 + 3.0/(mu_r + 2)
    else:
        Hz_analytical_profile[i] = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/r_val)**3

# Z-axis profile (r=0.01, avoiding singularity at r=0)
z_profile = profile_range
Hz_pert_numerical_z = zeros(len(z_profile))
Hz_pert_analytical_z = zeros(len(z_profile))

for i, zval in enumerate(z_profile):
    r = sqrt(0.01**2 + zval**2)  # Distance from origin
    if r < outer_radius - 0.1:  # Inside mesh domain
        try:
            mip = mesh(0.01, zval)  # Near z-axis
            Hz_pert_numerical_z[i] = H_pert[1](mip)
        except:
            Hz_pert_numerical_z[i] = nan
    else:
        Hz_pert_numerical_z[i] = nan

    # Analytical solution on z-axis (r_cyl ≈ 0)
    # The correct formula (verified against NGSolve):
    # H_z = +2 * (μr-1)/(μr+2) * (a/|z|)³ on z-axis (independent of z sign)
    if r < sphere_radius:
        Hz_pert_analytical_z[i] = -1.0 + 3.0/(mu_r + 2)
    else:
        Hz_pert_analytical_z[i] = 2 * (mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3

# Compute error statistics
valid_idx = ~isnan(Hz_numerical)
interior_idx = valid_idx & (r_profile < sphere_radius)
exterior_idx = valid_idx & (r_profile >= sphere_radius)

print(f"\n  Validation results (r-axis, perturbation field Hz):")
print(f"  -" * 30)

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
# Analytical Flux Lines in r-z Plane
# ============================================================
print("\nComputing analytical flux lines...")

# For analytical solution, compute H_pert in r-z plane
r_grid = linspace(0.02, plot_range, 101)
z_grid = linspace(-plot_range, plot_range, 101)
[rr_anal, zz_anal] = meshgrid(r_grid, z_grid)

Hr_analytical = zeros(rr_anal.shape)
Hz_analytical_2d = zeros(rr_anal.shape)

for nz in range(len(z_grid)):
    for nr in range(len(r_grid)):
        r = sqrt(r_grid[nr]**2 + z_grid[nz]**2)
        if r < 0.01:  # Avoid singularity at origin
            r = 0.01

        if r < sphere_radius:
            # Inside sphere: H_pert = constant in z-direction
            Hr_analytical[nz, nr] = 0.0
            Hz_analytical_2d[nz, nr] = -1.0 + 3.0/(mu_r + 2)
        else:
            # Outside sphere: dipole field (verified against NGSolve)
            # H_pert = C * [2cosθ er + sinθ eθ]
            # where C = +(μr-1)/(μr+2) * (a/r)³ (POSITIVE coefficient!)
            # θ is angle from z-axis
            theta = arctan2(r_grid[nr], z_grid[nz])  # Angle from z-axis
            C = (mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3  # POSITIVE

            # Spherical components
            H_r = 2 * C * cos(theta)
            H_theta = C * sin(theta)

            # Convert to cylindrical (r_cyl, z) coordinates
            Hr_analytical[nz, nr] = H_r * sin(theta) + H_theta * cos(theta)
            Hz_analytical_2d[nz, nr] = H_r * cos(theta) - H_theta * sin(theta)

# ============================================================
# Save Results
# ============================================================
print("\nSaving results...")

from scipy.io import savemat

# Create 2D grid for visualization
r_grid = linspace(0.02, plot_range, 101)
z_grid = linspace(-plot_range, plot_range, 101)
[rr, zz] = meshgrid(r_grid, z_grid)

Hr_field = zeros(rr.shape)
Hz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        r_dist = sqrt(r_grid[ir]**2 + z_grid[iz]**2)
        if r_dist < outer_radius - 0.1:
            try:
                mip = mesh(r_grid[ir], z_grid[iz])
                Hr_field[iz, ir] = H_pert[0](mip)
                Hz_field[iz, ir] = H_pert[1](mip)
            except:
                Hr_field[iz, ir] = nan
                Hz_field[iz, ir] = nan
        else:
            Hr_field[iz, ir] = nan
            Hz_field[iz, ir] = nan

mat_data = {
    'rr': rr,
    'zz': zz,
    'Hr': Hr_field,
    'Hz': Hz_field,
    'r_profile': r_profile,
    'Hz_numerical': Hz_numerical,
    'Hz_analytical': Hz_analytical_profile,
    'sphere_radius': sphere_radius,
    'outer_radius': outer_radius,
    'mu_r': mu_r
}

mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, mat_data)
print(f"  MAT file saved to: {mat_file}")

# ============================================================
# Visualization (2x2 layout matching 3D_dipole.py)
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

# Create figure with 2x2 subplots (matching 3D_dipole.py)
# Row 1: Analytical H vs NGSolve H (streamlines)
# Row 2: R-axis and Z-axis profile comparisons
fig = plt.figure(figsize=(12, 10), dpi=150)

# Row 1, Col 1: Interior H field (Analytical)
ax1 = plt.subplot(2, 2, 1)
strm1 = ax1.streamplot(rr_anal, zz_anal, Hr_analytical, Hz_analytical_2d,
                       color='red', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
# Draw magnetic sphere boundary (half-circle for axisymmetric)
theta_circle = linspace(0, pi, 100)
r_sphere = sphere_radius * sin(theta_circle)
z_sphere = sphere_radius * cos(theta_circle)
ax1.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
ax1.plot(r_sphere, z_sphere, 'r-', linewidth=2, label='Magnetic material')
ax1.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (Analytical)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(0, plot_range)
ax1.set_ylim(-plot_range, plot_range)
ax1.minorticks_on()
ax1.tick_params(which='major', direction="in", top=True, right=True)
ax1.tick_params(which='minor', direction="in", top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 1, Col 2: Interior H field (NGSolve)
ax2 = plt.subplot(2, 2, 2)
strm2 = ax2.streamplot(rr, zz, Hr_field, Hz_field,
                       color='black', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
ax2.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
ax2.plot(r_sphere, z_sphere, 'r-', linewidth=2, label='Magnetic material')
ax2.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (NGSolve)', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(0, plot_range)
ax2.set_ylim(-plot_range, plot_range)
ax2.minorticks_on()
ax2.tick_params(which='major', direction="in", top=True, right=True)
ax2.tick_params(which='minor', direction="in", top=True, right=True)
ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 2, Col 1: R-axis profile comparison
ax3 = plt.subplot(2, 2, 3)
ax3.plot(r_profile, Hz_numerical, 'k-', linewidth=2, label='NGSolve')
ax3.plot(r_profile, Hz_analytical_profile, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax3.set_title('R-axis Profile ($H_z$ component)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: Z-axis profile comparison
ax4 = plt.subplot(2, 2, 4)
ax4.plot(z_profile, Hz_pert_numerical_z, 'k-', linewidth=2, label='NGSolve')
ax4.plot(z_profile, Hz_pert_analytical_z, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax4.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax4.set_title('Z-axis Profile ($H_z$ component)', fontname='Times New Roman', fontsize=11)
ax4.minorticks_on()
ax4.tick_params(which='major', direction="in", top=True, right=True)
ax4.tick_params(which='minor', direction="in", top=True, right=True)
ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax4.legend(loc='best', fontsize=9, frameon=False)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
