"""
Axisymmetric A-formulation for magnetic sphere in uniform field
Simple version WITHOUT Kelvin transformation (for validation)

A-formulation for axisymmetric problems:
  Variable: u = r * A_theta (modified vector potential)
  Equation: -div(nu/r * grad(u)) = 0 (no current source)

Relation to B field:
  Br = -(1/r) * du/dz
  Bz = (1/r) * du/dr

Problem: Magnetic sphere (mu_r) in uniform external field H0 (z-direction)

Boundary condition approach:
  u = u_s on outer boundary, where u_s = B0 * r^2 / 2 gives uniform Bz = B0
  u = 0 on axis (r = 0)

Analytical solution for sphere:
  Interior: Hz = 3/(mu_r + 2) * H0  (uniform)
  Exterior: perturbation field decays as 1/r^3

Note: Uses half-circle (r >= 0) geometry for proper axisymmetric treatment.
Boundary naming is done AFTER mesh generation using coordinate-based identification.
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("="*60)
print("Axisymmetric A-formulation")
print("Magnetic Sphere in Uniform Field (Half-circle version)")
print("="*60)

# ============================================================
# Parameters
# ============================================================
R_sphere = 0.5          # Sphere radius [m]
R_outer = 5.0           # Outer boundary radius [m]
maxh = 0.05             # Mesh size [m]

mu0 = 4*pi*1e-7         # Vacuum permeability [H/m]
nu0 = 1/mu0             # Vacuum reluctivity [m/H]
mu_r = 100              # Relative permeability
H0 = 1.0                # Applied field [A/m]
B0 = mu0 * H0           # Applied flux density [T]

print(f"\nParameters:")
print(f"  Sphere radius: {R_sphere} m")
print(f"  Outer boundary: {R_outer} m")
print(f"  Relative permeability: {mu_r}")
print(f"  Applied field: H0 = {H0} A/m")
print(f"  Applied flux density: B0 = {B0:.6e} T")

# ============================================================
# Analytical Solution
# ============================================================
Hz_int_analytical = 3.0 / (mu_r + 2) * H0
Hz_pert_int = Hz_int_analytical - H0  # Perturbation field inside

print(f"\nAnalytical solution (3D sphere):")
print(f"  Interior Hz = {Hz_int_analytical:.6f} A/m")
print(f"  Interior Hz_pert = {Hz_pert_int:.6f} A/m")

# ============================================================
# Geometry - Half-circle using WorkPlane
# ============================================================
print("\nCreating half-circle geometry...")

# Use WorkPlane to create half-circle (right half, x >= 0)
# For axisymmetric: x = r (radial), y = z (axial)
wp = WorkPlane()

# Outer half-circle: arc from (0, -R_outer) to (0, R_outer) via (R_outer, 0)
outer_half = wp.MoveTo(0, -R_outer).Arc(R_outer, 180).LineTo(0, -R_outer).Close().Face()
outer_half.name = 'air'

# Sphere half-circle
wp = WorkPlane()
sphere_half = wp.MoveTo(0, -R_sphere).Arc(R_sphere, 180).LineTo(0, -R_sphere).Close().Face()
sphere_half.name = 'magnetic'
sphere_half.maxh = maxh/2

# Subtract sphere from outer to get air region
air_half = outer_half - sphere_half
air_half.name = 'air'

# Glue together
shape = Glue([air_half, sphere_half])
geo = OCCGeometry(shape, dim=2)

# ============================================================
# Mesh Generation
# ============================================================
print("Generating mesh...")
ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)

# Name boundaries based on coordinates BEFORE creating NGSolve mesh
# This requires direct access to netgen mesh
for seg in ngmesh.Elements1D():
    # Get segment vertices
    pts = [ngmesh[v] for v in seg.vertices]
    p0, p1 = pts[0], pts[1]

    # Center of segment
    cx = (p0[0] + p1[0]) / 2
    cy = (p0[1] + p1[1]) / 2
    dist = sqrt(cx**2 + cy**2)

    # Check if on axis (x ≈ 0)
    if abs(cx) < 0.001 and abs(p0[0]) < 0.001 and abs(p1[0]) < 0.001:
        ngmesh.SetBCName(seg.index - 1, 'axis')
    # Check if on outer arc (dist ≈ R_outer)
    elif abs(dist - R_outer) < 0.1:
        ngmesh.SetBCName(seg.index - 1, 'outer')
    # Otherwise interface
    else:
        ngmesh.SetBCName(seg.index - 1, 'interface')

mesh = Mesh(ngmesh)

print(f"  Number of elements: {mesh.ne}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# Count boundary types
bnds = mesh.GetBoundaries()
axis_count = sum(1 for b in bnds if 'axis' in b)
outer_count = sum(1 for b in bnds if 'outer' in b)
print(f"  'axis' boundaries: {axis_count}")
print(f"  'outer' boundaries: {outer_count}")

# ============================================================
# Finite Element Space
# ============================================================
print("\nSetting up A-formulation...")

# Dirichlet BC: u = 0 on axis, u = u_s on outer
fes = H1(mesh, order=3, dirichlet="axis|outer")
print(f"  Number of DOFs: {fes.ndof}")

u = fes.TrialFunction()
v = fes.TestFunction()

# ============================================================
# Material Properties
# ============================================================
# r-weight for axisymmetric formulation (r = x in r-z plane)
r_weight = IfPos(x - 1e-10, x, 1e-10)

# Build reluctivity by region
materials = mesh.GetMaterials()
print(f"  Materials: {materials}")

is_magnetic = CoefficientFunction([1 if 'magnetic' in mat else 0 for mat in materials])
nu_cf = IfPos(is_magnetic - 0.5, nu0 / mu_r, nu0)

# ============================================================
# Weak Form
# ============================================================
# Bilinear form: integral( nu/r * grad(u) . grad(v) dx )
a_form = BilinearForm(fes)
a_form += nu_cf / r_weight * grad(u) * grad(v) * dx

# Linear form: zero (no current source)
f = LinearForm(fes)

print("\nAssembling system...")
a_form.Assemble()
f.Assemble()
print("  System assembled")

# ============================================================
# Boundary Condition
# ============================================================
# u_s = B0 * r^2 / 2 gives uniform Bz = B0
# On outer boundary (arc at r = R_outer): u = B0 * r^2 / 2 = B0 * x^2 / 2
# On axis (r = 0): u = 0

print("\nSetting boundary condition...")

gfu = GridFunction(fes)
# Set u = B0 * x^2 / 2 on outer boundary (u = 0 on axis automatically via x^2)
gfu.Set(B0 * x**2 / 2, BND)

# Solve using residual method
res = gfu.vec.CreateVector()
res.data = f.vec - a_form.mat * gfu.vec
gfu.vec.data += a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
print("  Solution converged")

# ============================================================
# Validation
# ============================================================
print("\n" + "="*60)
print("VALIDATION: Field comparison with analytical solution")
print("="*60)

# In axisymmetric r-z coordinates: x = r, y = z
# Bz = (1/r) * du/dr = (1/x) * du/dx
# Br = -(1/r) * du/dz = -(1/x) * du/dy

r_sample = 0.05  # Slightly off axis
print(f"\nOn z-axis (r = {r_sample} m):")
print(f"{'z (m)':<10} {'Hz_num (A/m)':<15} {'Hz_ana (A/m)':<15} {'Error (%)':<10}")
print("-"*50)

for z_val in [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0]:
    try:
        mip = mesh(r_sample, z_val)
        grad_u = grad(gfu)
        Bz_num = grad_u[0](mip) / r_sample
        nu_val = nu_cf(mip)
        Hz_num = nu_val * Bz_num

        dist = sqrt(r_sample**2 + z_val**2)
        if dist < R_sphere:
            Hz_ana = Hz_int_analytical
        else:
            coef = (mu_r - 1) / (mu_r + 2) * R_sphere**3
            Hz_ana = H0 + coef * H0 * (2*z_val**2 - r_sample**2) / dist**5

        err = abs(Hz_num - Hz_ana) / abs(Hz_ana) * 100 if abs(Hz_ana) > 1e-15 else 0
        print(f"{z_val:<10.2f} {Hz_num:<15.6f} {Hz_ana:<15.6f} {err:<10.3f}")
    except Exception as e:
        print(f"{z_val:<10.2f} Error: {e}")

print(f"\nOn radial axis (z = 0.05):")
print(f"{'r (m)':<10} {'Hz_num (A/m)':<15} {'Hz_ana (A/m)':<15} {'Error (%)':<10}")
print("-"*50)

for r_val in [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0]:
    try:
        z_val = 0.05
        mip = mesh(r_val, z_val)
        grad_u = grad(gfu)
        Bz_num = grad_u[0](mip) / r_val
        nu_val = nu_cf(mip)
        Hz_num = nu_val * Bz_num

        dist = sqrt(r_val**2 + z_val**2)
        if dist < R_sphere:
            Hz_ana = Hz_int_analytical
        else:
            coef = (mu_r - 1) / (mu_r + 2) * R_sphere**3
            Hz_ana = H0 + coef * H0 * (2*z_val**2 - r_val**2) / dist**5

        err = abs(Hz_num - Hz_ana) / abs(Hz_ana) * 100 if abs(Hz_ana) > 1e-15 else 0
        print(f"{r_val:<10.2f} {Hz_num:<15.6f} {Hz_ana:<15.6f} {err:<10.3f}")
    except Exception as e:
        print(f"{r_val:<10.2f} Error: {e}")

# ============================================================
# RMS Error
# ============================================================
print("\n" + "="*60)
print("RMS Error Analysis")
print("="*60)

# Interior RMS
err_sum_int = 0
n_int = 0
for z_val in linspace(-0.4, 0.4, 17):
    for r_val in linspace(0.05, 0.4, 8):
        if sqrt(r_val**2 + z_val**2) < R_sphere:
            try:
                mip = mesh(r_val, z_val)
                grad_u = grad(gfu)
                Bz_num = grad_u[0](mip) / r_val
                nu_val = nu_cf(mip)
                Hz_num = nu_val * Bz_num
                Hz_ana = Hz_int_analytical
                err_sum_int += (Hz_num - Hz_ana)**2
                n_int += 1
            except:
                pass

if n_int > 0:
    rms_int = sqrt(err_sum_int / n_int)
    print(f"Interior RMS error: {rms_int:.6e} A/m ({rms_int/H0*100:.2f}% of H0)")
    print(f"  Sampled {n_int} points")

# Exterior RMS
err_sum_ext = 0
n_ext = 0
for z_val in linspace(0.55, 4.5, 20):
    try:
        mip = mesh(0.05, z_val)
        grad_u = grad(gfu)
        Bz_num = grad_u[0](mip) / 0.05
        Hz_num = nu0 * Bz_num

        r_sample = 0.05
        coef = (mu_r - 1) / (mu_r + 2) * R_sphere**3
        dist = sqrt(r_sample**2 + z_val**2)
        Hz_ana = H0 + coef * H0 * (2*z_val**2 - r_sample**2) / dist**5

        err_sum_ext += (Hz_num - Hz_ana)**2
        n_ext += 1
    except:
        pass

if n_ext > 0:
    rms_ext = sqrt(err_sum_ext / n_ext)
    print(f"Exterior RMS error (z-axis): {rms_ext:.6e} A/m ({rms_ext/H0*100:.2f}% of H0)")
    print(f"  Sampled {n_ext} points")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

fig = plt.figure(figsize=(12, 10), dpi=150)

# Compute 2D field
r_grid = linspace(0.02, 2.0, 41)
z_grid = linspace(-2.0, 2.0, 41)
[rr, zz] = meshgrid(r_grid, z_grid)

Br_field = zeros(rr.shape)
Bz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        try:
            mip = mesh(r_grid[ir], z_grid[iz])
            grad_u = grad(gfu)
            Br_field[iz, ir] = -grad_u[1](mip) / r_grid[ir]
            Bz_field[iz, ir] = grad_u[0](mip) / r_grid[ir]
        except:
            Br_field[iz, ir] = nan
            Bz_field[iz, ir] = nan

def draw_sphere(ax):
    theta = linspace(0, pi, 100)
    ax.plot(R_sphere * sin(theta), R_sphere * cos(theta), 'r-', linewidth=2, label='Sphere')

# B field streamlines
ax1 = plt.subplot(2, 2, 1)
Br_plot = Br_field.copy()
Bz_plot = Bz_field.copy()
Br_plot[isnan(Br_plot)] = 0
Bz_plot[isnan(Bz_plot)] = 0

strm = ax1.streamplot(rr, zz, Br_plot, Bz_plot,
                      color='blue', linewidth=0.8, density=1.5,
                      arrowsize=0.6, arrowstyle='->')
draw_sphere(ax1)
ax1.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('$\\mathbf{B}$ Field (A-formulation)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(0, 2)
ax1.set_ylim(-2, 2)
ax1.legend(loc='upper right', fontsize=8)
ax1.grid(True, alpha=0.3)

# |B| contour
ax2 = plt.subplot(2, 2, 2)
B_mag = sqrt(Br_field**2 + Bz_field**2)
B_mag_plot = B_mag.copy()
B_mag_plot[isnan(B_mag_plot)] = 0
levels = linspace(0, B0*3, 20)
contour2 = ax2.contourf(rr, zz, B_mag_plot, levels=levels, cmap='jet', extend='max')
cbar2 = plt.colorbar(contour2, ax=ax2, shrink=0.8)
cbar2.set_label('$|\\mathbf{B}|$ (T)', fontname='Times New Roman', fontsize=10)
draw_sphere(ax2)
ax2.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('$|\\mathbf{B}|$ Contour', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(0, 2)
ax2.set_ylim(-2, 2)
ax2.grid(True, alpha=0.3)

# Hz profile
ax3 = plt.subplot(2, 2, 3)
z_profile = linspace(-4.8, 4.8, 97)
Hz_numerical = zeros(len(z_profile))
Hz_analytical_arr = zeros(len(z_profile))
r_sample = 0.05

for i, z_val in enumerate(z_profile):
    try:
        mip = mesh(r_sample, z_val)
        grad_u = grad(gfu)
        Bz_num = grad_u[0](mip) / r_sample
        nu_val = nu_cf(mip)
        Hz_numerical[i] = nu_val * Bz_num

        dist = sqrt(r_sample**2 + z_val**2)
        if dist < R_sphere:
            Hz_analytical_arr[i] = Hz_int_analytical
        else:
            coef = (mu_r - 1) / (mu_r + 2) * R_sphere**3
            Hz_analytical_arr[i] = H0 + coef * H0 * (2*z_val**2 - r_sample**2) / dist**5
    except:
        Hz_numerical[i] = nan
        Hz_analytical_arr[i] = nan

valid = ~isnan(Hz_numerical)
ax3.plot(z_profile[valid], Hz_numerical[valid], 'k-', linewidth=2, label='NGSolve')
ax3.plot(z_profile, Hz_analytical_arr, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(R_sphere, color='green', linestyle=':', alpha=0.7, label='Sphere boundary')
ax3.axvline(-R_sphere, color='green', linestyle=':', alpha=0.7)
ax3.axhline(Hz_int_analytical, color='orange', linestyle=':', alpha=0.7, label=f'Interior Hz = {Hz_int_analytical:.4f}')
ax3.axhline(H0, color='blue', linestyle=':', alpha=0.7, label=f'H0 = {H0}')
ax3.set_xlabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$H_z$ (A/m)', fontname='Times New Roman', fontsize=10)
ax3.set_title('$H_z$ Profile on $z$-axis', fontname='Times New Roman', fontsize=11)
ax3.legend(loc='best', fontsize=9)
ax3.grid(True, alpha=0.3)

# Error profile
ax4 = plt.subplot(2, 2, 4)
error_profile = zeros(len(z_profile))
for i in range(len(z_profile)):
    if not isnan(Hz_numerical[i]) and not isnan(Hz_analytical_arr[i]) and abs(Hz_analytical_arr[i]) > 1e-15:
        error_profile[i] = (Hz_numerical[i] - Hz_analytical_arr[i]) / Hz_analytical_arr[i] * 100
    else:
        error_profile[i] = nan

valid_err = ~isnan(error_profile)
ax4.plot(z_profile[valid_err], error_profile[valid_err], 'b-', linewidth=2)
ax4.axhline(0, color='gray', linestyle='-', alpha=0.5)
ax4.axvline(R_sphere, color='green', linestyle=':', alpha=0.7)
ax4.axvline(-R_sphere, color='green', linestyle=':', alpha=0.7)
ax4.set_xlabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('Relative Error (%)', fontname='Times New Roman', fontsize=10)
ax4.set_title('$H_z$ Error on $z$-axis', fontname='Times New Roman', fontsize=11)
ax4.grid(True, alpha=0.3)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
