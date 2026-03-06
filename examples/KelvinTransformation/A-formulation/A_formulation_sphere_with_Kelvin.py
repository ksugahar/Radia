"""
Axisymmetric A-formulation for magnetic sphere with Z-OFFSET Kelvin transformation

A-formulation for axisymmetric problems:
  Variable: u = r * A_theta (modified vector potential)
  Equation: -div(nu/r * grad(u)) = 0 (no current source)

Relation to B field:
  Br = -(1/r) * du/dz
  Bz = (1/r) * du/dr

Problem: Magnetic sphere (mu_r) in uniform external field H0 (z-direction)

Z-offset Kelvin transformation:
  - Interior domain: half-circle centered at (0, 0), radius a
  - Exterior domain: half-circle centered at (0, z_offset), radius a
  - r = x is the SAME in both domains (key simplification!)
  - Kelvin factor for nu: (rho'/a)^2

Analytical solution for sphere:
  Interior: Hz = 3/(mu_r + 2) * H0  (uniform)
  Exterior: perturbation field decays as 1/r^3
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("="*60)
print("Axisymmetric A-formulation with Z-OFFSET Kelvin")
print("Magnetic Sphere in Uniform Field")
print("="*60)

# ============================================================
# Parameters
# ============================================================
R_sphere = 0.5          # Sphere radius [m]
a = 1.0                 # Kelvin boundary radius [m]
z_offset = 2.5          # Z-offset for exterior domain [m]
maxh = 0.03             # Mesh size [m]

mu0 = 4*pi*1e-7         # Vacuum permeability [H/m]
nu0 = 1/mu0             # Vacuum reluctivity [m/H]
mu_r = 100              # Relative permeability
H0 = 1.0                # Applied field [A/m]
B0 = mu0 * H0           # Applied flux density [T]

print(f"\nParameters:")
print(f"  Sphere radius: {R_sphere} m")
print(f"  Kelvin boundary radius: a = {a} m")
print(f"  Z-offset: {z_offset} m")
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
# Geometry with Z-offset Kelvin (using full circle -> cutter pattern)
# ============================================================
print("\nCreating geometry with Z-offset Kelvin...")

# Interior domain: FULL circle first, then cut to half
wp1 = WorkPlane()
inner_full = wp1.Circle(a).Face()
inner_full.name = "air_inner"

# Magnetic sphere: FULL circle, then cut
wp2 = WorkPlane()
sphere_full = wp2.Circle(R_sphere).Face()
sphere_full.name = "magnetic"
sphere_full.maxh = maxh/2

# Exterior domain: FULL circle at (0, z_offset), then cut
wp3 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
outer_full = wp3.Circle(a).Face()
outer_full.name = "air_outer"

# GND point at center of exterior domain (represents infinity)
gnd_point = Vertex(Pnt(0, z_offset, 0))
gnd_point.name = "GND"

# Cut to get half-circles (keep r >= 0, i.e., x >= 0)
cutter_left = MoveTo(-a-0.1, -a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()
cutter_left_ext = MoveTo(-a-0.1, z_offset - a - 0.1).Rectangle(a+0.1, 2*a+0.2).Face()

inner_half = inner_full - cutter_left
sphere_half = sphere_full - cutter_left
outer_half = outer_full - cutter_left_ext

# Inner air = inner half - sphere half
inner_air = inner_half - sphere_half
inner_air.name = "air_inner"

# ============================================================
# Name boundaries and identify periodic edges
# ============================================================
print("  Naming boundaries...")

# Inner air edges - find Kelvin arcs using vertex distance
kelvin_inner_edges = []
for edge in inner_air.edges:
    cx = edge.center.x
    cy = edge.center.y
    # Check if both vertices are at distance a from origin
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + v0.p.y**2)
        d1 = sqrt(v1.p.x**2 + v1.p.y**2)
        is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
    except:
        is_kelvin_arc = False

    if cx < 0.01:  # On axis
        edge.name = "axis"
    elif is_kelvin_arc:
        edge.name = "kelvin_int"
        kelvin_inner_edges.append(edge)
    else:
        edge.name = "interface"

# Sphere edges
for edge in sphere_half.edges:
    cx = edge.center.x
    if cx < 0.01:
        edge.name = "axis"
    else:
        edge.name = "sphere_bnd"

# Outer half edges - find Kelvin arcs using vertex distance from (0, z_offset)
kelvin_outer_edges = []
for edge in outer_half.edges:
    cx = edge.center.x
    cy = edge.center.y - z_offset
    # Check if both vertices are at distance a from (0, z_offset)
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
        d1 = sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
        is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
    except:
        is_kelvin_arc = False

    if cx < 0.01:  # On axis
        edge.name = "axis_ext"
    elif is_kelvin_arc:
        edge.name = "kelvin_ext"
        kelvin_outer_edges.append(edge)
    else:
        edge.name = "default"

print(f"  Found {len(kelvin_inner_edges)} interior Kelvin edges")
print(f"  Found {len(kelvin_outer_edges)} exterior Kelvin edges")

# Identify periodic boundaries by matching y-coordinate sign
# (both arcs have the same shape, just offset in y)
matched_pairs = 0
for int_edge in kelvin_inner_edges:
    int_y = int_edge.center.y
    for ext_edge in kelvin_outer_edges:
        ext_y = ext_edge.center.y - z_offset
        # Match by sign of y-coordinate (upper arc with upper arc, lower with lower)
        if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
            matched_pairs += 1
            break

print(f"  Matched {matched_pairs} periodic edge pairs")

# Combine and mesh
shape = Glue([inner_air, sphere_half, outer_half, gnd_point])
geo = OCCGeometry(shape, dim=2)

print(f"  Interior: half-circle centered at (0, 0)")
print(f"  Exterior: half-circle centered at (0, {z_offset})")

# ============================================================
# Mesh
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.5))

print(f"  Number of elements: {mesh.ne}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Finite Element Space
# ============================================================
print("\nSetting up A-formulation...")

# Dirichlet BC: u = 0 on axis (r=0) and at GND (infinity)
fes_before = H1(mesh, order=3, dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
fes = Periodic(fes_before)

# Check if Periodic BC is working
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")

if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

print(f"  Number of DOFs: {fes.ndof}")

u = fes.TrialFunction()
v = fes.TestFunction()

# ============================================================
# Material Properties with Z-offset Kelvin
# ============================================================
# KEY: r = x is the SAME in both interior and exterior domains!
r_weight = IfPos(x - 1e-10, x, 1e-10)

# For exterior domain: rho' = sqrt(x^2 + (y - z_offset)^2)
y_local = y - z_offset
rho_prime = sqrt(x**2 + y_local**2)
rho_prime_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)

# Kelvin transformation factor for reluctivity: nu' = nu0 * (rho'/a)^2
kelvin_factor = (rho_prime_safe / a)**2

# Build reluctivity by region
materials = mesh.GetMaterials()
print(f"\n  Materials: {materials}")

nu_list = []
for mat in materials:
    if "outer" in mat.lower():
        nu_list.append(nu0 * kelvin_factor)
    elif "magnetic" in mat.lower():
        nu_list.append(nu0 / mu_r)
    else:
        nu_list.append(nu0)

nu_cf = CoefficientFunction(nu_list)

print(f"  Key: r = x in both domains (no anisotropy!)")
print(f"  Kelvin factor: (rho'/a)^2")

# ============================================================
# Weak Form (Total potential approach)
# ============================================================
# For A-formulation with uniform field B0 in z-direction:
# Total potential: u = u_s + u_r
# where u_s = B0 * r^2 / 2 gives uniform Bz = B0 in air
#
# The source potential u_s satisfies: -div(nu0/r * grad(u_s)) = 0 in air
# We solve for the total potential u directly
#
# In Kelvin exterior domain:
#   nu' = nu0 * (rho'/a)^2
#   The source potential transforms but we solve total potential directly
#
# Boundary conditions:
#   axis (r=0): u = 0 (natural for u = r*A_theta)
#   GND (infinity): u = 0 (perturbation goes to zero at infinity)
#   Kelvin boundary: periodic BC connects interior and exterior
#
# With periodic BC, the exterior domain sees the source field through
# the matching at the Kelvin boundary

# Source potential for uniform field
us_inner = B0 * x**2 / 2  # gives Bz = B0

# Bilinear form: integral( nu/r * grad(u) . grad(v) dx )
a_form = BilinearForm(fes)
a_form += nu_cf / r_weight * grad(u) * grad(v) * dx

# Linear form: zero for homogeneous problem
f = LinearForm(fes)

print("\nAssembling system...")
a_form.Assemble()
f.Assemble()
print("  System assembled")

# ============================================================
# Solve with Source Potential BC
# ============================================================
print("\nSolving with source potential boundary condition...")

gfu = GridFunction(fes)

# Set boundary condition: u = u_s on Kelvin boundary (interior side)
# This imposes the background field at the computational boundary
# The axis BC u=0 is automatically satisfied

# For the sphere problem, we need to apply u_s at the outer boundary
# With Kelvin transformation, the "outer" is at the Kelvin arc

# Since kelvin_int is not a Dirichlet boundary, we use a different approach:
# Solve the reduced potential problem u_r where u = u_s + u_r
# The equation for u_r in the magnetic region is:
# -div(nu/r * grad(u_r)) = div((nu_0 - nu)/r * grad(u_s))

# For magnetic material:
# RHS = div((nu_0 - nu_mag)/r * grad(u_s)) = (nu_0 - nu_0/mu_r) * div(1/r * grad(u_s))
# Since u_s = B0*x^2/2, grad(u_s) = (B0*x, 0) = (B0*r, 0)
# 1/r * grad(u_s) = (B0, 0)
# div(B0, 0) = 0
# So RHS = 0, which means we need source from boundary integral!

# Actually, the correct approach for reduced potential is:
# a(u_r, v) = -a_s(u_s, v) where a_s uses (nu_0 - nu) instead of nu

# Build source term for reduced potential
# Source coefficient: (nu_0 - nu) in magnetic region only
nu_diff_list = []
for mat in materials:
    if "magnetic" in mat.lower():
        nu_diff_list.append(nu0 - nu0/mu_r)  # nu_0 - nu_mag > 0
    else:
        nu_diff_list.append(0.0)  # zero in air regions

nu_diff_cf = CoefficientFunction(nu_diff_list)

# grad(u_s) = (B0*x, 0)
grad_us = CoefficientFunction((B0 * x, 0))

# Source term: integral( (nu_0 - nu)/r * grad(u_s) . grad(v) dx )
f_source = LinearForm(fes)
f_source += nu_diff_cf / r_weight * InnerProduct(grad_us, grad(v)) * dx("magnetic")
f_source.Assemble()

print(f"  Source term norm: {Norm(f_source.vec):.6e}")

# Solve using residual method
res = gfu.vec.CreateVector()
res.data = f_source.vec - a_form.mat * gfu.vec
gfu.vec.data += a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
print("  Solution converged")

# The result gfu is the reduced potential u_r
# Total field: u_total = u_s + u_r
# Total: Bz = (1/r) * d(u_s + u_r)/dr = B0 + (1/r) * du_r/dr

# ============================================================
# Validation
# ============================================================
print("\n" + "="*60)
print("VALIDATION: Field comparison with analytical solution")
print("="*60)

# Helper function to compute Hz from reduced potential
def compute_Hz(r_val, z_val):
    """Compute Hz = nu * Bz from reduced potential u_r"""
    mip = mesh(r_val, z_val)
    grad_u_r = grad(gfu)
    du_r_dr = grad_u_r[0](mip)

    # Source potential gradient: d(B0*x^2/2)/dx = B0*x
    dus_dr = B0 * r_val

    # Total: Bz = (1/r) * (du_s/dr + du_r/dr)
    Bz_num = (dus_dr + du_r_dr) / r_val

    # Hz = nu * Bz
    nu_val = nu_cf(mip)
    Hz_num = nu_val * Bz_num
    return Hz_num

# In interior domain only (x, y) with y in [-a, a]
r_sample = 0.05
print(f"\nOn z-axis (r = {r_sample} m), Interior domain:")
print(f"{'z (m)':<10} {'Hz_num (A/m)':<15} {'Hz_ana (A/m)':<15} {'Error (%)':<10}")
print("-"*50)

for z_val in [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.95]:
    try:
        Hz_num = compute_Hz(r_sample, z_val)

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

print(f"\nOn radial axis (z = 0.05), Interior domain:")
print(f"{'r (m)':<10} {'Hz_num (A/m)':<15} {'Hz_ana (A/m)':<15} {'Error (%)':<10}")
print("-"*50)

for r_val in [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.95]:
    try:
        z_val = 0.05
        Hz_num = compute_Hz(r_val, z_val)

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
print("RMS Error Analysis (Interior Domain)")
print("="*60)

# Interior RMS (sphere)
err_sum_int = 0
n_int = 0
for z_val in linspace(-0.4, 0.4, 17):
    for r_val in linspace(0.05, 0.4, 8):
        if sqrt(r_val**2 + z_val**2) < R_sphere:
            try:
                Hz_num = compute_Hz(r_val, z_val)
                Hz_ana = Hz_int_analytical
                err_sum_int += (Hz_num - Hz_ana)**2
                n_int += 1
            except:
                pass

if n_int > 0:
    rms_int = sqrt(err_sum_int / n_int)
    print(f"Interior (sphere) RMS error: {rms_int:.6e} A/m ({rms_int/H0*100:.2f}% of H0)")
    print(f"  Sampled {n_int} points")

# Exterior RMS (inner domain, outside sphere)
err_sum_ext = 0
n_ext = 0
r_sample = 0.05
for z_val in linspace(0.55, 0.95, 10):
    try:
        Hz_num = compute_Hz(r_sample, z_val)

        coef = (mu_r - 1) / (mu_r + 2) * R_sphere**3
        dist = sqrt(r_sample**2 + z_val**2)
        Hz_ana = H0 + coef * H0 * (2*z_val**2 - r_sample**2) / dist**5

        err_sum_ext += (Hz_num - Hz_ana)**2
        n_ext += 1
    except:
        pass

if n_ext > 0:
    rms_ext = sqrt(err_sum_ext / n_ext)
    print(f"Exterior (inner domain) RMS error: {rms_ext:.6e} A/m ({rms_ext/H0*100:.2f}% of H0)")
    print(f"  Sampled {n_ext} points")

# ============================================================
# Visualization (3x2 layout matching H-formulation)
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

# Create figure with 3x2 subplots (matching H-formulation)
# Row 1: Interior H (Analytical vs NGSolve)
# Row 2: R-axis and Z-axis profile comparisons
# Row 3: Exterior B and H
fig = plt.figure(figsize=(12, 15), dpi=150)

plot_range = a  # Plot range matches Kelvin radius

# Compute 2D field (interior domain only)
r_grid = linspace(0.02, 0.95, 30)
z_grid = linspace(-0.95, 0.95, 30)
[rr, zz] = meshgrid(r_grid, z_grid)

Hr_field = zeros(rr.shape)
Hz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        try:
            r_val = r_grid[ir]
            z_val = z_grid[iz]
            mip = mesh(r_val, z_val)
            grad_u_r = grad(gfu)
            du_r_dr = grad_u_r[0](mip)
            du_r_dz = grad_u_r[1](mip)
            # Source potential: u_s = B0*x^2/2
            dus_dr = B0 * r_val
            dus_dz = 0
            # Total field: Bz = (1/r)*(du_s/dr + du_r/dr), Br = -(1/r)*(du_s/dz + du_r/dz)
            Br = -(dus_dz + du_r_dz) / r_val
            Bz = (dus_dr + du_r_dr) / r_val
            # H = nu * B (perturbation field = H - Hs)
            nu_val = nu_cf(mip)
            Hr_field[iz, ir] = nu_val * Br - 0  # Hs_r = 0
            Hz_field[iz, ir] = nu_val * Bz - H0  # Hs_z = H0
        except:
            Hr_field[iz, ir] = nan
            Hz_field[iz, ir] = nan

# Analytical perturbation field for interior domain
Hr_analytical = zeros(rr.shape)
Hz_analytical_2d = zeros(rr.shape)
coef_anal = (mu_r - 1) / (mu_r + 2) * R_sphere**3

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        r_val = r_grid[ir]
        z_val = z_grid[iz]
        dist = sqrt(r_val**2 + z_val**2)
        if dist < R_sphere:
            # Interior: uniform perturbation field
            Hr_analytical[iz, ir] = 0
            Hz_analytical_2d[iz, ir] = Hz_int_analytical - H0
        else:
            # Exterior: dipole perturbation field
            Hz_analytical_2d[iz, ir] = coef_anal * H0 * (2*z_val**2 - r_val**2) / dist**5
            Hr_analytical[iz, ir] = coef_anal * H0 * 3*r_val*z_val / dist**5

# Draw magnetic sphere boundary
theta_circle = linspace(0, pi, 100)
r_sphere_plot = R_sphere * sin(theta_circle)
z_sphere_plot = R_sphere * cos(theta_circle)

# Draw Kelvin boundary
r_kelvin_plot = a * sin(theta_circle)
z_kelvin_plot = a * cos(theta_circle)

# Row 1, Col 1: Interior H field (Analytical)
ax1 = plt.subplot(3, 2, 1)
Hr_anal_plot = Hr_analytical.copy()
Hz_anal_plot = Hz_analytical_2d.copy()
Hr_anal_plot[isnan(Hr_anal_plot)] = 0
Hz_anal_plot[isnan(Hz_anal_plot)] = 0
strm1 = ax1.streamplot(rr, zz, Hr_anal_plot, Hz_anal_plot,
                       color='red', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
ax1.fill_betweenx(z_sphere_plot, 0, r_sphere_plot, alpha=0.3, color='lightblue')
ax1.plot(r_sphere_plot, z_sphere_plot, 'r-', linewidth=2, label='Magnetic material')
ax1.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
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
ax2 = plt.subplot(3, 2, 2)
Hr_plot = Hr_field.copy()
Hz_plot = Hz_field.copy()
Hr_plot[isnan(Hr_plot)] = 0
Hz_plot[isnan(Hz_plot)] = 0
strm2 = ax2.streamplot(rr, zz, Hr_plot, Hz_plot,
                       color='black', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
ax2.fill_betweenx(z_sphere_plot, 0, r_sphere_plot, alpha=0.3, color='lightblue')
ax2.plot(r_sphere_plot, z_sphere_plot, 'r-', linewidth=2, label='Magnetic material')
ax2.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
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

# Row 2, Col 1: R-axis profile comparison (perturbation field)
ax3 = plt.subplot(3, 2, 3)
r_profile = linspace(0.05, 0.95, 30)
Hz_numerical_r = zeros(len(r_profile))
Hz_analytical_r = zeros(len(r_profile))
z_sample = 0.05

for i, r_val in enumerate(r_profile):
    try:
        Hz_numerical_r[i] = compute_Hz(r_val, z_sample) - H0  # perturbation
        dist = sqrt(r_val**2 + z_sample**2)
        if dist < R_sphere:
            Hz_analytical_r[i] = Hz_int_analytical - H0
        else:
            Hz_analytical_r[i] = coef_anal * H0 * (2*z_sample**2 - r_val**2) / dist**5
    except:
        Hz_numerical_r[i] = nan
        Hz_analytical_r[i] = nan

ax3.plot(r_profile, Hz_numerical_r, 'k-', linewidth=2, label='NGSolve')
ax3.plot(r_profile, Hz_analytical_r, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(R_sphere, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax3.set_title('R-axis Profile (Perturbation Field)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: Z-axis profile comparison (perturbation field)
ax4 = plt.subplot(3, 2, 4)
z_profile = linspace(-0.95, 0.95, 50)
Hz_pert_numerical_z = zeros(len(z_profile))
Hz_pert_analytical_z = zeros(len(z_profile))
r_sample = 0.05

for i, z_val in enumerate(z_profile):
    try:
        Hz_pert_numerical_z[i] = compute_Hz(r_sample, z_val) - H0  # perturbation
        dist = sqrt(r_sample**2 + z_val**2)
        if dist < R_sphere:
            Hz_pert_analytical_z[i] = Hz_int_analytical - H0
        else:
            Hz_pert_analytical_z[i] = coef_anal * H0 * (2*z_val**2 - r_sample**2) / dist**5
    except:
        Hz_pert_numerical_z[i] = nan
        Hz_pert_analytical_z[i] = nan

valid = ~isnan(Hz_pert_numerical_z)
ax4.plot(z_profile[valid], Hz_pert_numerical_z[valid], 'k-', linewidth=2, label='NGSolve')
ax4.plot(z_profile, Hz_pert_analytical_z, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-R_sphere, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax4.axvline(R_sphere, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax4.set_title('Z-axis Profile (Perturbation Field)', fontname='Times New Roman', fontsize=11)
ax4.minorticks_on()
ax4.tick_params(which='major', direction="in", top=True, right=True)
ax4.tick_params(which='minor', direction="in", top=True, right=True)
ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax4.legend(loc='best', fontsize=9, frameon=False)

# Row 3: Exterior domain fields
# Compute exterior domain field
r_ext = linspace(0.02, 0.95, 25)
z_ext = linspace(-0.95, 0.95, 25)
[rr_ext, zz_ext] = meshgrid(r_ext, z_ext)

Hr_ext = zeros(rr_ext.shape)
Hz_ext = zeros(rr_ext.shape)
Br_ext = zeros(rr_ext.shape)
Bz_ext = zeros(rr_ext.shape)

for nz in range(len(z_ext)):
    for nr in range(len(r_ext)):
        r_from_center = sqrt(r_ext[nr]**2 + z_ext[nz]**2)
        if r_from_center < a - 0.05:  # Inside exterior domain (mapped)
            try:
                # Map to exterior domain coordinates
                mip = mesh(r_ext[nr], z_offset + z_ext[nz])
                grad_u_r = grad(gfu)
                du_r_dr = grad_u_r[0](mip)
                du_r_dz = grad_u_r[1](mip)
                # In exterior domain, source potential needs Kelvin transform
                # But since we're looking at perturbation, just use reduced potential
                r_val = r_ext[nr]
                Br = -du_r_dz / r_val
                Bz = du_r_dr / r_val
                # Kelvin factor for mu
                rho_prime = sqrt(r_ext[nr]**2 + z_ext[nz]**2)
                kelvin_fac = (rho_prime / a)**2
                nu_ext = nu0 * kelvin_fac
                Hr_ext[nz, nr] = nu_ext * Br
                Hz_ext[nz, nr] = nu_ext * Bz
                # B field in exterior (with Kelvin scaling)
                mu_ext = mu0 / kelvin_fac
                Br_ext[nz, nr] = mu_ext * Hr_ext[nz, nr]
                Bz_ext[nz, nr] = mu_ext * Hz_ext[nz, nr]
            except:
                Hr_ext[nz, nr] = nan
                Hz_ext[nz, nr] = nan
                Br_ext[nz, nr] = nan
                Bz_ext[nz, nr] = nan
        else:
            Hr_ext[nz, nr] = nan
            Hz_ext[nz, nr] = nan
            Br_ext[nz, nr] = nan
            Bz_ext[nz, nr] = nan

# Row 3, Col 1: Exterior B field
ax5 = plt.subplot(3, 2, 5)
Br_ext_plot = Br_ext.copy()
Bz_ext_plot = Bz_ext.copy()
Br_ext_plot[isnan(Br_ext_plot)] = 0
Bz_ext_plot[isnan(Bz_ext_plot)] = 0
strm5 = ax5.streamplot(rr_ext, zz_ext, Br_ext_plot, Bz_ext_plot,
                       color='darkblue', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
ax5.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
ax5.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_xlabel("${\\it r'}$ (m)", fontname='Times New Roman', fontsize=10)
plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax5.set_title('Exterior: Flux Density $\\mathbf{B}$', fontname='Times New Roman', fontsize=11)
ax5.set_aspect('equal')
ax5.set_xlim(0, plot_range)
ax5.set_ylim(-plot_range, plot_range)
ax5.minorticks_on()
ax5.tick_params(which='major', direction="in", top=True, right=True)
ax5.tick_params(which='minor', direction="in", top=True, right=True)
ax5.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 3, Col 2: Exterior H field
ax6 = plt.subplot(3, 2, 6)
Hr_ext_plot = Hr_ext.copy()
Hz_ext_plot = Hz_ext.copy()
Hr_ext_plot[isnan(Hr_ext_plot)] = 0
Hz_ext_plot[isnan(Hz_ext_plot)] = 0
strm6 = ax6.streamplot(rr_ext, zz_ext, Hr_ext_plot, Hz_ext_plot,
                       color='darkgreen', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
ax6.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
ax6.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_xlabel("${\\it r'}$ (m)", fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Exterior: Magnetic Field $\\mathbf{H}$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(0, plot_range)
ax6.set_ylim(-plot_range, plot_range)
ax6.minorticks_on()
ax6.tick_params(which='major', direction="in", top=True, right=True)
ax6.tick_params(which='minor', direction="in", top=True, right=True)
ax6.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
