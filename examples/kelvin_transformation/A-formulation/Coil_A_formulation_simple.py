"""
Axisymmetric A-formulation (WITHOUT Kelvin transformation)
Single coil magnetostatics problem - Baseline validation

A-formulation for axisymmetric problems:
  Variable: u = r * A_theta (modified vector potential)
  Equation: -div(nu * grad(u) / r) = J_theta

  where nu = 1/mu (reluctivity)

Problem: Single circular coil at (r_coil, 0)
  - Coil radius: r_coil
  - Current density: J_theta (concentrated in small area)
  - Analytical: Biot-Savart law on axis
"""
import os
from numpy import *
from scipy.special import ellipk, ellipe
from ngsolve import *
from netgen.occ import *

print("="*60)
print("Axisymmetric A-formulation (No Kelvin)")
print("Single Coil Problem - Baseline Validation")
print("="*60)

# ============================================================
# Parameters
# ============================================================
coil_radius = 0.3      # Coil center radius [m]
coil_width = 0.02      # Coil cross-section width [m]
coil_height = 0.02     # Coil cross-section height [m]
outer_radius = 1.0     # Outer boundary [m]
plot_range = 0.8       # Plot range [m]
maxh = 0.02            # Mesh size [m]

mu0 = 4*pi*1e-7        # Vacuum permeability [H/m]
nu0 = 1/mu0            # Vacuum reluctivity [m/H]
I_total = 1.0          # Total current [A]
J0 = I_total / (coil_width * coil_height)  # Current density [A/m^2]

print(f"\nParameters:")
print(f"  Coil radius: {coil_radius} m")
print(f"  Coil cross-section: {coil_width} x {coil_height} m")
print(f"  Outer boundary: {outer_radius} m")
print(f"  Total current: {I_total} A")
print(f"  Current density: {J0:.2e} A/m^2")

# ============================================================
# Geometry (Axisymmetric: r >= 0, z)
# ============================================================
print("\nCreating geometry...")

# Full domain: half-circle r in [0, outer_radius]
wp = WorkPlane()
outer_circle = wp.Circle(outer_radius).Face()

# Coil region (small rectangle at r=coil_radius, z=0)
coil_rect = MoveTo(coil_radius - coil_width/2, -coil_height/2).Rectangle(coil_width, coil_height).Face()
coil_rect.name = "coil"
coil_rect.maxh = maxh/2

# Cut coil from domain
air = outer_circle - coil_rect
air.name = "air"

# Cut left half (keep r >= 0)
cutter_left = MoveTo(-outer_radius-0.1, -outer_radius-0.1).Rectangle(outer_radius+0.1, 2*outer_radius+0.2).Face()
air = air - cutter_left
coil_rect = coil_rect - cutter_left

# Name boundaries
for edge in air.edges:
    cx = edge.center.x
    cy = edge.center.y
    r_edge = sqrt(cx**2 + cy**2)
    if cx < 0.01:
        edge.name = "axis"
    elif r_edge > outer_radius - 0.01:
        edge.name = "outer"
    else:
        edge.name = "coil_boundary"

for edge in coil_rect.edges:
    if edge.center.x < 0.01:
        edge.name = "axis"
    else:
        edge.name = "coil_boundary"

shape = Glue([air, coil_rect])
geo = OCCGeometry(shape, dim=2)

print(f"  Domain: half-circle r in [0, {outer_radius}]")

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

# Dirichlet BC: u = 0 on axis (r=0) and outer boundary
fes = H1(mesh, order=3, dirichlet="axis|outer")

print(f"  Number of DOFs: {fes.ndof}")
print(f"  FreeDofs: {sum([1 for d in fes.FreeDofs() if d])}")

u = fes.TrialFunction()
v = fes.TestFunction()

# ============================================================
# Material Properties
# ============================================================
nu = nu0  # Uniform reluctivity (vacuum)

# r weight for axisymmetric
r_weight = IfPos(x - 1e-10, x, 1e-10)

# ============================================================
# Weak Form
# ============================================================
# Bilinear form: a(u,v) = integral( nu/r * grad(u) . grad(v) dx )
a_form = BilinearForm(fes)
a_form += nu / r_weight * grad(u) * grad(v) * dx

# Linear form: f(v) = integral( J_theta * v dx )
f = LinearForm(fes)
f += J0 * v * dx("coil")

print("\nAssembling system...")
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
# Analytical solution using elliptic integrals (from CircleEH.m)
# ============================================================
def CircleEH(OPr, OPz, r_coil):
    """
    Compute magnetic field (Hr, Hz) and vector potential (At) for a circular current loop.
    Based on CircleEH.m using complete elliptic integrals.

    Parameters:
        OPr: observation point r-coordinate
        OPz: observation point z-coordinate
        r_coil: coil radius

    Returns:
        Hr: r-component of H field (for unit current I=1)
        Hz: z-component of H field (for unit current I=1)
        At: theta-component of A vector potential (for unit current I=1)
    """
    R = abs(OPr)  # Observation point radial distance
    Z = OPz       # Observation point z coordinate
    r = r_coil    # Coil radius

    # Elliptic integral parameter k^2
    k2 = abs(4 * r * R) / ((r + R)**2 + Z**2)

    if k2 > 0 and k2 < 1:
        K = ellipk(k2)  # Complete elliptic integral of the first kind
        E = ellipe(k2)  # Complete elliptic integral of the second kind
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
    """
    Compute analytical B field at (r_obs, z_obs) for the coil.
    Returns Br, Bz in Tesla.
    """
    if r_obs < 1e-12:
        # On axis: use simplified formula
        R = coil_radius
        Hz = I_total * R**2 / (2 * (R**2 + z_obs**2)**1.5)
        return 0.0, mu0 * Hz

    Hr, Hz, At = CircleEH(r_obs, z_obs, coil_radius)
    # Scale by current
    Br = mu0 * I_total * Hr
    Bz = mu0 * I_total * Hz
    return Br, Bz

def Bz_analytical_on_axis(z_val):
    """Analytical Bz on z-axis using elliptic integral formula."""
    _, Bz = B_analytical(0.001, z_val)  # Use small r to avoid singularity
    return Bz

# Sample field on axis
print("\n" + "="*60)
print("VALIDATION: Field on z-axis (r=0.01)")
print("="*60)

z_samples = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
print(f"{'z (m)':<10} {'Bz (NGSolve)':<15} {'Bz (Analytical)':<15} {'Error (%)':<10}")
print("-"*50)

for z_val in z_samples:
    try:
        mip = mesh(0.01, z_val)
        Bz_num = grad(gfu)[0](mip) / 0.01  # du/dr / r at r=0.01
        Bz_ana = Bz_analytical_on_axis(z_val)
        err = abs(Bz_num - Bz_ana) / abs(Bz_ana) * 100 if abs(Bz_ana) > 1e-15 else 0
        print(f"{z_val:<10.2f} {Bz_num:<15.6e} {Bz_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{z_val:<10.2f} Error: {e}")

# Field at coil center (r=0, z=0)
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
z_profile = linspace(-plot_range, plot_range, 161)
Bz_z_numerical = zeros(len(z_profile))
Bz_z_analytical = zeros(len(z_profile))

for i, z_val in enumerate(z_profile):
    Bz_z_analytical[i] = Bz_analytical_on_axis(z_val)
    try:
        mip = mesh(0.01, z_val)
        Bz_z_numerical[i] = grad(gfu)[0](mip) / 0.01
    except:
        Bz_z_numerical[i] = nan

# R-axis profile (z=0)
r_profile = linspace(0.02, plot_range, 81)
Bz_r_numerical = zeros(len(r_profile))
Bz_r_analytical = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    # Analytical
    _, Bz_r_analytical[i] = B_analytical(r_val, 0.0)
    # Numerical
    try:
        mip = mesh(r_val, 0.0)
        Bz_r_numerical[i] = grad(gfu)[0](mip) / r_val
    except:
        Bz_r_numerical[i] = nan

# Error statistics
valid_idx = ~isnan(Bz_z_numerical)
if sum(valid_idx) > 0:
    error = Bz_z_numerical[valid_idx] - Bz_z_analytical[valid_idx]
    rms_error = sqrt(mean(error**2))
    max_error = abs(error).max()
    print(f"\n  Z-axis profile error statistics:")
    print(f"    RMS error: {rms_error:.6e} T")
    print(f"    Max error: {max_error:.6e} T")
    print(f"    Relative RMS error: {rms_error/abs(Bz_z_analytical[valid_idx]).mean()*100:.3f}%")

# ============================================================
# 2D field for visualization
# ============================================================
print("\nComputing 2D field...")

r_grid = linspace(0.02, plot_range, 51)
z_grid = linspace(-plot_range, plot_range, 51)
[rr, zz] = meshgrid(r_grid, z_grid)

Br_field = zeros(rr.shape)
Bz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        rho_val = sqrt(r_grid[ir]**2 + z_grid[iz]**2)
        if rho_val < outer_radius - 0.05:
            try:
                mip = mesh(r_grid[ir], z_grid[iz])
                grad_u = grad(gfu)
                Br_field[iz, ir] = -grad_u[1](mip) / r_grid[ir]
                Bz_field[iz, ir] = grad_u[0](mip) / r_grid[ir]
            except:
                Br_field[iz, ir] = nan
                Bz_field[iz, ir] = nan
        else:
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
    'Bz_z_numerical': Bz_z_numerical,
    'Bz_z_analytical': Bz_z_analytical,
    'r_profile': r_profile,
    'Bz_r_numerical': Bz_r_numerical,
    'coil_radius': coil_radius,
    'outer_radius': outer_radius,
    'I_total': I_total,
    'mu0': mu0
}

mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, mat_data)
print(f"  MAT file saved to: {mat_file}")

# ============================================================
# Visualization: 3x2 Grid
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

fig = plt.figure(figsize=(12, 15), dpi=150)

# Helper: Draw coil
def draw_coil(ax, color='red', alpha=0.5):
    coil_x = [coil_radius - coil_width/2, coil_radius + coil_width/2,
              coil_radius + coil_width/2, coil_radius - coil_width/2,
              coil_radius - coil_width/2]
    coil_y = [-coil_height/2, -coil_height/2, coil_height/2, coil_height/2, -coil_height/2]
    ax.fill(coil_x, coil_y, color, alpha=alpha)
    ax.plot(coil_x, coil_y, color=color, linewidth=2)

# Row 1, Col 1: B field streamlines
ax1 = plt.subplot(3, 2, 1)
strm1 = ax1.streamplot(rr, zz, Br_field, Bz_field,
                       color='blue', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
draw_coil(ax1, 'red', 0.5)
ax1.plot([coil_radius], [0], 'ro', markersize=8, label='Coil')
plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('$\\mathbf{B}$ Field (Streamlines)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(0, plot_range)
ax1.set_ylim(-plot_range, plot_range)
ax1.legend(loc='upper right', fontsize=8, frameon=False)
ax1.minorticks_on()
ax1.tick_params(which='major', direction="in", top=True, right=True)
ax1.tick_params(which='minor', direction="in", top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 1, Col 2: |B| contour
ax2 = plt.subplot(3, 2, 2)
B_mag = sqrt(Br_field**2 + Bz_field**2) * 1e6  # Convert to µT
levels = linspace(0, nanmax(B_mag)*0.8, 20)
contour2 = ax2.contourf(rr, zz, B_mag, levels=levels, cmap='jet', extend='max')
draw_coil(ax2, 'white', 0.3)
cbar2 = plt.colorbar(contour2, ax=ax2)
cbar2.set_label('$|\\mathbf{B}|$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('$|\\mathbf{B}|$ Contour', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(0, plot_range)
ax2.set_ylim(-plot_range, plot_range)

# Row 2, Col 1: Z-axis profile comparison
ax3 = plt.subplot(3, 2, 3)
ax3.plot(z_profile, Bz_z_numerical * 1e6, 'k-', linewidth=2, label='NGSolve')
ax3.plot(z_profile, Bz_z_analytical * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$B_{z}$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
ax3.set_title('Z-axis Profile ($r \\approx 0$)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: R-axis profile (z=0)
ax4 = plt.subplot(3, 2, 4)
ax4.plot(r_profile, Bz_r_numerical * 1e6, 'k-', linewidth=2, label='NGSolve')
ax4.plot(r_profile, Bz_r_analytical * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(coil_radius, color='green', linestyle=':', linewidth=1.5, alpha=0.7, label='Coil position')
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$B_{z}$ ($\\mu$T)', fontname='Times New Roman', fontsize=10)
ax4.set_title('R-axis Profile ($z=0$)', fontname='Times New Roman', fontsize=11)
ax4.minorticks_on()
ax4.tick_params(which='major', direction="in", top=True, right=True)
ax4.tick_params(which='minor', direction="in", top=True, right=True)
ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax4.legend(loc='best', fontsize=9, frameon=False)

# Row 3, Col 1: Z-axis error
ax5 = plt.subplot(3, 2, 5)
error_percent = (Bz_z_numerical - Bz_z_analytical) / Bz_z_analytical * 100
valid = ~isnan(error_percent) & (abs(Bz_z_analytical) > 1e-12)
ax5.plot(z_profile[valid], error_percent[valid], 'b-', linewidth=1.5)
ax5.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('Relative Error (%)', fontname='Times New Roman', fontsize=10)
ax5.set_title('Z-axis Error vs Analytical', fontname='Times New Roman', fontsize=11)
ax5.minorticks_on()
ax5.tick_params(which='major', direction="in", top=True, right=True)
ax5.tick_params(which='minor', direction="in", top=True, right=True)
ax5.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax5.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)

# Row 3, Col 2: Vector potential (u/r = A_theta)
ax6 = plt.subplot(3, 2, 6)
A_theta = zeros(rr.shape)
for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        try:
            mip = mesh(r_grid[ir], z_grid[iz])
            A_theta[iz, ir] = gfu(mip) / r_grid[ir]
        except:
            A_theta[iz, ir] = nan

levels6 = linspace(nanmin(A_theta), nanmax(A_theta), 20)
contour6 = ax6.contourf(rr, zz, A_theta * 1e6, 20, cmap='coolwarm')
draw_coil(ax6, 'black', 0.3)
cbar6 = plt.colorbar(contour6, ax=ax6)
cbar6.set_label('$A_{\\theta}$ ($\\mu$Wb/m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Vector Potential $A_{\\theta}$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(0, plot_range)
ax6.set_ylim(-plot_range, plot_range)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

try:
    os.startfile(png_file)
except:
    pass

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
