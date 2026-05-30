"""
Axisymmetric A-formulation with Z-OFFSET Kelvin transformation
Single coil magnetostatics problem

Key idea: Offset exterior domain in Z-direction
This keeps r-coordinate (x in 2D) the SAME in both domains!

A-formulation for axisymmetric problems:
  Variable: u = r * A_theta (modified vector potential)
  Equation: -div(nu * grad(u) / r) = J_theta

Kelvin transformation with Z-offset:
  - Interior domain: half-circle centered at (0, 0), radius a
  - Exterior domain: half-circle centered at (0, z_offset), radius a
  - r' = r (unchanged!) - KEY SIMPLIFICATION
  - Kelvin factor for nu: (rho'/a)^2

Problem: Single circular coil at (r_coil, 0)
"""
import os
from numpy import *
from scipy.special import ellipk, ellipe
from ngsolve import *
from netgen.occ import *

print("="*60)
print("Axisymmetric A-formulation with Z-OFFSET Kelvin")
print("Single Coil Problem")
print("="*60)

# ============================================================
# Parameters
# ============================================================
coil_radius = 0.6      # Coil center radius [m]
coil_width = 0.02      # Coil cross-section width [m]
coil_height = 0.02     # Coil cross-section height [m]
a = 1.0                # Kelvin boundary radius [m]
z_offset = 2.5         # Z-offset for exterior domain [m]
maxh = 0.02            # Mesh size [m]

mu0 = 4*pi*1e-7        # Vacuum permeability [H/m]
nu0 = 1/mu0            # Vacuum reluctivity [m/H]
I_total = 1.0          # Total current [A]
J0 = I_total / (coil_width * coil_height)  # Current density [A/m^2]

print(f"\nParameters:")
print(f"  Coil radius: {coil_radius} m")
print(f"  Coil cross-section: {coil_width} x {coil_height} m")
print(f"  Kelvin boundary: a = {a} m")
print(f"  Z-offset: {z_offset} m")
print(f"  Total current: {I_total} A")
print(f"  Current density: {J0:.2e} A/m^2")

# ============================================================
# Geometry with Z-offset
# ============================================================
print("\nCreating geometry with Z-offset...")

# Interior domain: FULL circle first
wp1 = WorkPlane()
inner_full = wp1.Circle(a).Face()
inner_full.name = "air_inner"

# Coil region (full)
coil_full = MoveTo(coil_radius - coil_width/2, -coil_height/2).Rectangle(coil_width, coil_height).Face()
coil_full.name = "coil"
coil_full.maxh = maxh/2

# Exterior domain: FULL circle at (0, z_offset)
wp2 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
outer_full = wp2.Circle(a).Face()
outer_full.name = "air_outer"

# GND point at center of exterior domain
gnd_point = Vertex(Pnt(0, z_offset, 0))
gnd_point.name = "GND"

# Cut to get half-circles (keep r >= 0, i.e., x >= 0)
cutter_left = MoveTo(-a-0.1, -a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()
cutter_left_ext = MoveTo(-a-0.1, z_offset - a - 0.1).Rectangle(a+0.1, 2*a+0.2).Face()

inner_half = inner_full - cutter_left
coil_half = coil_full - cutter_left
outer_half = outer_full - cutter_left_ext

# Subtract coil from inner air
inner_air = inner_half - coil_half
inner_air.name = "air_inner"

# Name boundaries and find Kelvin edges
print("  Naming boundaries...")

# Find Kelvin boundary arcs in inner_air (radius = a, centered at origin)
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

    if cx < 0.01:
        edge.name = "axis"
    elif is_kelvin_arc:
        edge.name = "kelvin_int"
        kelvin_inner_edges.append(edge)
    else:
        edge.name = "default"

for edge in coil_half.edges:
    cx = edge.center.x
    if cx < 0.01:
        edge.name = "axis"
    else:
        edge.name = "coil_bnd"

# Find Kelvin boundary arcs in outer_half (radius = a, centered at (0, z_offset))
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

    if cx < 0.01:
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
shape = Glue([inner_air, coil_half, outer_half, gnd_point])
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

# Kelvin modulation via centralized helper (Nagamine CEFC 2026 canonical):
#   nu_ext = nu_0 * (rho'/a)^2,  rho' = sqrt(x^2 + (y - z_offset)^2)
# See examples/kelvin_transformation/CONVENTION.md.
from radia.kelvin_source import kelvin_nu_factor_axisym_cf, build_material_cf

nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=z_offset, R=a)
nu_cf = build_material_cf(mesh, nu0, nu_kelvin_factor)

print(f"\n  Materials: {mesh.GetMaterials()}")

print(f"  Key: r = x in both domains (no anisotropy!)")
print(f"  Kelvin factor: (rho'/a)^2")

# ============================================================
# Weak Form
# ============================================================
# Bilinear form: integral( nu/r * grad(u) . grad(v) dx )
a_form = BilinearForm(fes)
a_form += nu_cf / r_weight * grad(u) * grad(v) * dx

# Linear form: f(v) = integral( J_theta * v dx )
f = LinearForm(fes)
f += J0 * v * dx("coil")

print("\nAssembling system...")
with TaskManager():
    a_form.Assemble()
    f.Assemble()
    print("  System assembled")

    # ============================================================
    # Solve
    # ============================================================
    print("\nSolving system...")
    gfu = GridFunction(fes)
    gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# ============================================================
# Analytical solution using elliptic integrals
# ============================================================
def CircleEH(OPr, OPz, r_coil):
    """Compute magnetic field for a circular current loop."""
    R = abs(OPr)
    Z = OPz
    r = r_coil

    k2 = abs(4 * r * R) / ((r + R)**2 + Z**2)

    if k2 > 0 and k2 < 1:
        K = ellipk(k2)
        E = ellipe(k2)
        k = sqrt(k2)

        Hr = sign(OPr) * (1/(2*pi)) * (Z/R) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 + R**2 + Z**2) / ((r-R)**2 + Z**2) * E - K)
        Hz = (1/(2*pi)) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 - R**2 - Z**2) / ((r-R)**2 + Z**2) * E + K)
        At = sign(OPr) * sqrt(r/R) * ((2/k - k) * K - 2/k * E) / (2*pi)
    elif k2 == 0:
        K = ellipk(k2)
        E = ellipe(k2)
        Hr = 0.0
        Hz = (1/(2*pi)) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 - R**2 - Z**2) / ((r-R)**2 + Z**2) * E + K)
        At = 0.0
    else:
        Hr = 0.0
        Hz = 0.0
        At = 0.0

    return Hr, Hz, At

def B_analytical(r_obs, z_obs):
    """Compute analytical B field at (r_obs, z_obs) for the coil."""
    if r_obs < 1e-12:
        R = coil_radius
        Hz = I_total * R**2 / (2 * (R**2 + z_obs**2)**1.5)
        return 0.0, mu0 * Hz

    Hr, Hz, At = CircleEH(r_obs, z_obs, coil_radius)
    Br = mu0 * I_total * Hr
    Bz = mu0 * I_total * Hz
    return Br, Bz

def Bz_analytical_on_axis(z_val):
    """Analytical Bz on z-axis."""
    _, Bz = B_analytical(0.001, z_val)
    return Bz

# Sample field on axis
print("\n" + "="*60)
print("VALIDATION: Field on z-axis (r=0.01)")
print("="*60)

z_samples = [0.0, 0.2, 0.4, 0.6, 0.8]
print(f"{'z (m)':<10} {'Bz (NGSolve)':<15} {'Bz (Analytical)':<15} {'Error (%)':<10}")
print("-"*50)

for z_val in z_samples:
    try:
        mip = mesh(0.01, z_val)
        Bz_num = grad(gfu)[0](mip) / 0.01
        Bz_ana = Bz_analytical_on_axis(z_val)
        err = abs(Bz_num - Bz_ana) / abs(Bz_ana) * 100 if abs(Bz_ana) > 1e-15 else 0
        print(f"{z_val:<10.2f} {Bz_num:<15.6e} {Bz_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{z_val:<10.2f} Error: {e}")

# Field at coil center
Bz_center_analytical = mu0 * I_total / (2 * coil_radius)
print(f"\nField at center (r=0.01, z=0):")
try:
    mip = mesh(0.01, 0)
    Bz_center = grad(gfu)[0](mip) / 0.01
    print(f"  NGSolve:    Bz = {Bz_center:.6e} T")
    print(f"  Analytical: Bz = {Bz_center_analytical:.6e} T")
    print(f"  Error: {abs(Bz_center - Bz_center_analytical)/abs(Bz_center_analytical)*100:.3f}%")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# Compute profiles
# ============================================================
print("\nComputing profiles...")

# Z-axis profile
z_profile = linspace(-0.9, 0.9, 91)
Bz_numerical = zeros(len(z_profile))
Bz_analytical_profile = zeros(len(z_profile))

for i, z_val in enumerate(z_profile):
    Bz_analytical_profile[i] = Bz_analytical_on_axis(z_val)
    try:
        mip = mesh(0.01, z_val)
        Bz_numerical[i] = grad(gfu)[0](mip) / 0.01
    except:
        Bz_numerical[i] = nan

# Error statistics
valid_idx = ~isnan(Bz_numerical)
if sum(valid_idx) > 0:
    error = Bz_numerical[valid_idx] - Bz_analytical_profile[valid_idx]
    rms_error = sqrt(mean(error**2))
    max_error = abs(error).max()
    print(f"\n  Profile error statistics:")
    print(f"    RMS error: {rms_error:.6e} T")
    print(f"    Max error: {max_error:.6e} T")
    print(f"    Relative RMS error: {rms_error/abs(Bz_analytical_profile[valid_idx]).mean()*100:.3f}%")

# ============================================================
# 2D field for visualization
# ============================================================
print("\nComputing 2D field...")

r_grid = linspace(0.02, a - 0.05, 41)
z_grid = linspace(-a + 0.05, a - 0.05, 41)
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

# ============================================================
# Save Results
# ============================================================
print("\nSaving results...")

from scipy.io import savemat

mat_data = {
    'rr': rr,
    'zz': zz,
    'Br': Br_field,
    'Bz': Bz_field,
    'z_profile': z_profile,
    'Bz_numerical': Bz_numerical,
    'Bz_analytical': Bz_analytical_profile,
    'coil_radius': coil_radius,
    'a': a,
    'I_total': I_total,
    'mu0': mu0
}

mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, mat_data)
print(f"  MAT file saved to: {mat_file}")

# ============================================================
# Additional profiles
# ============================================================
print("\nComputing additional profiles...")

# R-axis profile (z=0.1m)
z_rprofile = 0.1
r_profile = linspace(0.02, a - 0.05, 51)
Bz_r_numerical = zeros(len(r_profile))
Bz_r_analytical = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    _, Bz_r_analytical[i] = B_analytical(r_val, z_rprofile)
    try:
        mip = mesh(r_val, z_rprofile)
        Bz_r_numerical[i] = grad(gfu)[0](mip) / r_val
    except:
        Bz_r_numerical[i] = nan

# Compute |B| for contour plot
B_mag = sqrt(Br_field**2 + Bz_field**2)

# Error profile
error_profile = zeros(len(z_profile))
for i in range(len(z_profile)):
    if not isnan(Bz_numerical[i]) and abs(Bz_analytical_profile[i]) > 1e-15:
        error_profile[i] = (Bz_numerical[i] - Bz_analytical_profile[i]) / Bz_analytical_profile[i] * 100
    else:
        error_profile[i] = nan

# ============================================================
# Visualization - 3x2 Grid
# ============================================================
print("\nGenerating 3x2 plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

fig = plt.figure(figsize=(12, 15), dpi=150)

def draw_coil(ax):
    coil_x = [coil_radius - coil_width/2, coil_radius + coil_width/2,
              coil_radius + coil_width/2, coil_radius - coil_width/2,
              coil_radius - coil_width/2]
    coil_y = [-coil_height/2, -coil_height/2, coil_height/2, coil_height/2, -coil_height/2]
    ax.fill(coil_x, coil_y, 'red', alpha=0.5)
    ax.plot(coil_x, coil_y, 'r-', linewidth=2)

def draw_kelvin_boundary(ax):
    theta_arc = linspace(-pi/2, pi/2, 100)
    ax.plot(a * cos(theta_arc), a * sin(theta_arc), 'g--', linewidth=2, label='Kelvin boundary')

# Row 1: B field streamlines and |B| contour
ax1 = plt.subplot(3, 2, 1)
strm = ax1.streamplot(rr, zz, Br_field, Bz_field,
                      color='blue', linewidth=0.8, density=1.5,
                      arrowsize=0.6, arrowstyle='->')
draw_coil(ax1)
draw_kelvin_boundary(ax1)
ax1.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('$\\mathbf{B}$ Field (Z-offset Kelvin)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(0, a)
ax1.set_ylim(-a, a)
ax1.legend(loc='upper right', fontsize=8)
ax1.grid(True, alpha=0.3)

ax2 = plt.subplot(3, 2, 2)
B_mag_plot = B_mag.copy()
B_mag_plot[isnan(B_mag_plot)] = 0
levels = linspace(0, nanmax(B_mag)*0.8, 20)
contour2 = ax2.contourf(rr, zz, B_mag_plot*1e6, levels=levels*1e6, cmap='jet', extend='max')
cbar2 = plt.colorbar(contour2, ax=ax2, shrink=0.8)
cbar2.set_label('$|\\mathbf{B}|$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
draw_coil(ax2)
draw_kelvin_boundary(ax2)
ax2.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('$|\\mathbf{B}|$ Contour', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(0, a)
ax2.set_ylim(-a, a)
ax2.grid(True, alpha=0.3)

# Row 2: Z-axis and R-axis profiles
ax3 = plt.subplot(3, 2, 3)
ax3.plot(z_profile, Bz_numerical * 1e6, 'k-', linewidth=2, label='NGSolve')
ax3.plot(z_profile, Bz_analytical_profile * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax3.axhline(Bz_center_analytical * 1e6, color='green', linestyle=':', alpha=0.5, label='Center value')
ax3.set_xlabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$B_z$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
ax3.set_title('$B_z$ Profile on $z$-axis ($r \\approx 0$)', fontname='Times New Roman', fontsize=11)
ax3.legend(loc='best', fontsize=9)
ax3.grid(True, alpha=0.3)

ax4 = plt.subplot(3, 2, 4)
valid_r = ~isnan(Bz_r_numerical)
ax4.plot(r_profile[valid_r], Bz_r_numerical[valid_r] * 1e6, 'k-', linewidth=2, label='NGSolve')
ax4.plot(r_profile, Bz_r_analytical * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(coil_radius, color='green', linestyle=':', alpha=0.7, label=f'Coil ($r = {coil_radius}$ m)')
ax4.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$B_z$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
ax4.set_title(f'$B_z$ Profile on $r$-axis ($z = {z_rprofile}$ m)', fontname='Times New Roman', fontsize=11)
ax4.legend(loc='best', fontsize=9)
ax4.grid(True, alpha=0.3)

# Row 3: Error plot and vector potential
ax5 = plt.subplot(3, 2, 5)
valid_err = ~isnan(error_profile)
ax5.plot(z_profile[valid_err], error_profile[valid_err], 'b-', linewidth=2)
ax5.axhline(0, color='gray', linestyle='-', alpha=0.5)
ax5.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax5.set_xlabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('Relative Error (%)', fontname='Times New Roman', fontsize=10)
ax5.set_title('$B_z$ Error on $z$-axis (vs Analytical)', fontname='Times New Roman', fontsize=11)
ax5.grid(True, alpha=0.3)
valid_idx = ~isnan(Bz_numerical)
if sum(valid_idx) > 0:
    error = Bz_numerical[valid_idx] - Bz_analytical_profile[valid_idx]
    rms_err = sqrt(mean(error**2))
    rel_rms = rms_err / abs(Bz_analytical_profile[valid_idx]).mean() * 100
    ax5.text(0.05, 0.95, f'RMS: {rel_rms:.2f}%', transform=ax5.transAxes,
             fontsize=10, verticalalignment='top', fontname='Times New Roman',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

ax6 = plt.subplot(3, 2, 6)
u_field = zeros(rr.shape)
for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        try:
            mip = mesh(r_grid[ir], z_grid[iz])
            u_field[iz, ir] = gfu(mip)
        except:
            u_field[iz, ir] = nan
u_field_plot = u_field.copy()
u_field_plot[isnan(u_field_plot)] = 0
levels_u = linspace(nanmin(u_field), nanmax(u_field), 20)
contour6 = ax6.contourf(rr, zz, u_field_plot, levels=levels_u, cmap='RdBu_r')
cbar6 = plt.colorbar(contour6, ax=ax6, shrink=0.8)
cbar6.set_label('$u = r A_\\theta$ (Wb/m)', fontname='Times New Roman', fontsize=10)
draw_coil(ax6)
draw_kelvin_boundary(ax6)
ax6.set_xlabel('$r$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('$z$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Vector Potential $u = r A_\\theta$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(0, a)
ax6.set_ylim(-a, a)
ax6.grid(True, alpha=0.3)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
