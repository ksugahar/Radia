"""
H-formulation for 2D magnetostatics with perturbation potential
Geometry created internally using OCC (2D circle)
Updated: 2025-11-22
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 2D - OCC Geometry")
print("="*60)

# ============================================================
# Geometry Definition (OCC 2D)
# ============================================================
print("\nCreating geometry...")

# Parameters
circle_radius = 0.5      # Magnetic circle radius [m]
air_inner_radius = 1.0   # Inner air domain radius [m] (outer boundary)
maxh_fine = 0.03         # Fine mesh size [m] (for magnetic circle and inner air)

# Create magnetic circle (finest mesh)
wp = WorkPlane()
mag_circle_shape = wp.Circle(circle_radius).Face()
mag_circle_shape.maxh = maxh_fine

# Create inner air circle and name its boundary (fine mesh)
air_inner_circle_shape = wp.Circle(air_inner_radius).Face()
for edge in air_inner_circle_shape.edges:
    edge.name = "outer"  # Name outer boundary before boolean operation
air_inner_circle_shape.maxh = maxh_fine

# Boolean operations to create two regions (outer air removed)
# Inner air = inner circle - magnetic circle
air_inner_shape = air_inner_circle_shape - mag_circle_shape
air_inner_shape.name = "air_inner"

# Set magnetic material
mag_circle_shape.name = "magnetic"

# Combine into single geometry (only inner air and magnetic circle)
geo = Glue([air_inner_shape, mag_circle_shape])

print(f"Geometry created with two regions:")
print(f"  Magnetic circle radius: {circle_radius} m")
print(f"  Outer boundary radius: {air_inner_radius} m")
print(f"  Fine mesh size: {maxh_fine} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup
# ============================================================
print("\nSetting up H-formulation...")

n = specialcf.normal(mesh.dim)
fes = H1(mesh, order=3)
print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability
mu_d = {"default": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field: H_s = [0, 1] A/m (y-direction in 2D)
Hs = CoefficientFunction((0, 1))
Hsb = BoundaryFromVolumeCF(Hs)

print(f"  Background field: H_s = [0, 1] A/m (y-direction)")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = ∫(∇v)·(mu∇u)dOmega
# Note: 境界項 -∫v(n·mu∇u)dΓ は自然境界条件として省略される
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (PERTURBATION FORMULATION):
# f(v) = ∫(∇v)·(muH_s)dOmega - ∫v(n·muH_s)dΓ
# 注意: Kelvin変換なし（有限領域）では外部境界での境界項が必要
#       外部境界を通る磁束を考慮するため
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                    # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer"))  # 境界項（必須）

a.Assemble()
f.Assemble()

print("  System assembled")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfu = GridFunction(fes)
c = Preconditioner(a, type="local")

solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
           tol=1e-8, printrates=True, maxsteps=10000)

print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# Create evaluation grid (adjusted for outer boundary at r=1.0)
plot_range = 1.1  # Plot range [m]
x = linspace(-plot_range, plot_range, 221)
y = linspace(-plot_range, plot_range, 221)
[xx, yy] = meshgrid(x, y)

# Compute perturbation field: H_pert = -grad(phi)
H = -grad(gfu)
Hx = zeros((shape(xx)))
Hy = zeros((shape(xx)))
phi = zeros((shape(xx)))

for ny in range(len(y)):
    for nx in range(len(x)):
        r = sqrt(x[nx]**2 + y[ny]**2)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(x[nx], y[ny])
                phi[ny, nx] = gfu(mip)
                Hx[ny, nx] = H[0](mip)
                Hy[ny, nx] = H[1](mip)
            except:
                phi[ny, nx] = nan
                Hx[ny, nx] = nan
                Hy[ny, nx] = nan
        else:
            phi[ny, nx] = nan
            Hx[ny, nx] = nan
            Hy[ny, nx] = nan

# Center potential at origin
center_value = gfu(mesh(0, 0))
phi = phi - center_value

print(f"  Potential at origin: {center_value:.6e}")
print(f"  Field at origin: Hy = {H[1](mesh(0,0)):.6f} A/m")

# Expected analytical value (perturbation field interior, y-component)
# For background field [0, 1] in y-direction
Hy_analytical = -1.0 + 2.0/(mu_r + 1)  # 2D formula: -H0 + 2*H0/(mur+1)
print(f"  Analytical (interior): Hy = {Hy_analytical:.6f} A/m")
print(f"  Relative error: {abs(H[1](mesh(0,0)) - Hy_analytical)/abs(Hy_analytical)*100:.3f}%")

# ============================================================
# Profile Comparisons with Analytical Solution
# ============================================================
print("\nComputing axis profiles (perturbation field Hy)...")

profile_range = linspace(-plot_range, plot_range, 221)

# X-axis profile (evaluating Hy along x-axis with background in y-direction)
x_profile = profile_range
Hy_pert_numerical_x = zeros(len(x_profile))
Hy_pert_analytical_x = zeros(len(x_profile))

for i, xval in enumerate(x_profile):
    r = abs(xval)
    if r < air_inner_radius - 0.01:  # Inside mesh domain
        try:
            mip = mesh(xval, 0)
            Hy_pert_numerical_x[i] = H[1](mip)  # Hy component
        except:
            Hy_pert_numerical_x[i] = nan
    else:
        Hy_pert_numerical_x[i] = nan

    if r < circle_radius:
        Hy_pert_analytical_x[i] = -1.0 + 2.0/(mu_r + 1)
    else:
        # 2D exterior with background in y-direction:
        # On x-axis (theta=0 or pi from +x): Hy_pert = -(mur-1)/(mur+1) * (a/r)^2
        Hy_pert_analytical_x[i] = -(mu_r - 1)/(mu_r + 1) * (circle_radius/r)**2

# Y-axis profile (evaluating Hy along y-axis with background in y-direction)
y_profile = profile_range
Hy_pert_numerical_y = zeros(len(y_profile))
Hy_pert_analytical_y = zeros(len(y_profile))

for i, yval in enumerate(y_profile):
    r = abs(yval)
    if r < air_inner_radius - 0.01:  # Inside mesh domain
        try:
            mip = mesh(0, yval)
            Hy_pert_numerical_y[i] = H[1](mip)  # Hy component
        except:
            Hy_pert_numerical_y[i] = nan
    else:
        Hy_pert_numerical_y[i] = nan

    if r < circle_radius:
        Hy_pert_analytical_y[i] = -1.0 + 2.0/(mu_r + 1)
    else:
        # On y-axis (theta=pi/2 from +x): Hy_pert = (mur-1)/(mur+1) * (a/r)^2
        Hy_pert_analytical_y[i] = (mu_r - 1)/(mu_r + 1) * (circle_radius/r)**2

# Error statistics
valid_idx_x = ~isnan(Hy_pert_numerical_x)
interior_idx_x = valid_idx_x & (abs(x_profile) < circle_radius)

print(f"\n  Validation results (X-axis, perturbation field Hy):")
print(f"  -" * 30)

if sum(interior_idx_x) > 0:
    interior_error = Hy_pert_numerical_x[interior_idx_x] - Hy_pert_analytical_x[interior_idx_x]
    max_err_int = max(abs(interior_error))
    rms_err_int = sqrt(mean(interior_error**2))
    rel_err_int = rms_err_int / abs(Hy_pert_analytical_x[interior_idx_x][0]) * 100
    print(f"  Interior (|x| < {circle_radius} m):")
    print(f"    Max error: {max_err_int:.6e} A/m")
    print(f"    RMS error: {rms_err_int:.6e} A/m ({rel_err_int:.3f}%)")

# ============================================================
# Analytical Flux Lines (2D) - Background in y-direction
# ============================================================
print("\nComputing analytical flux lines...")

Hx_analytical = zeros((shape(xx)))
Hy_analytical = zeros((shape(xx)))

for ny in range(len(y)):
    for nx in range(len(x)):
        r = sqrt(x[nx]**2 + y[ny]**2)
        if r < 0.01:
            r = 0.01

        if r < circle_radius:
            # Inside circle: uniform perturbation field
            Hx_analytical[ny, nx] = 0.0
            Hy_analytical[ny, nx] = -1.0 + 2.0/(mu_r + 1)
        else:
            # Outside circle: 2D dipole with background in y-direction [0,1]
            # The 2D dipole solution for background in y-direction:
            # H_pert = C * (a/r)^2 * [sin(2theta), -cos(2theta)] where theta from +x axis
            theta = arctan2(y[ny], x[nx])  # Angle from +x axis
            C = (mu_r - 1)/(mu_r + 1) * (circle_radius/r)**2

            # Correct formula for y-direction background:
            Hx_analytical[ny, nx] = C * sin(2*theta)
            Hy_analytical[ny, nx] = -C * cos(2*theta)  # Note: negative sign

# ============================================================
# Save Results
# ============================================================
print("\nSaving results...")

# Save to .mat file
from scipy.io import savemat
mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, {
    'xx': xx,
    'yy': yy,
    'Hx_analytical': Hx_analytical,
    'Hy_analytical': Hy_analytical,
    'Hx': Hx,
    'Hy': Hy
})
print(f"  MAT file saved to: {mat_file}")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

# Create figure with 2x2 subplots (no exterior domain without Kelvin)
# Row 1: Analytical H and NGSolve H (streamlines)
# Row 2: X-axis and Y-axis profile comparisons
fig = plt.figure(figsize=(12, 10), dpi=150)

# Row 1, Col 1: Interior H field (Analytical)
ax1 = plt.subplot(2, 2, 1)
strm1 = ax1.streamplot(xx, yy, Hx_analytical, Hy_analytical,
                       color='red', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle1 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax1.add_patch(circle1)
ax1.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (Analytical)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(-plot_range, plot_range)
ax1.set_ylim(-plot_range, plot_range)
ax1.minorticks_on()
ax1.tick_params(which='major', direction="in", top=True, right=True)
ax1.tick_params(which='minor', direction="in", top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 1, Col 2: Interior H field (NGSolve)
ax2 = plt.subplot(2, 2, 2)
strm2 = ax2.streamplot(xx, yy, Hx, Hy,
                       color='black', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle2 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax2.add_patch(circle2)
ax2.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (NGSolve)', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(-plot_range, plot_range)
ax2.set_ylim(-plot_range, plot_range)
ax2.minorticks_on()
ax2.tick_params(which='major', direction="in", top=True, right=True)
ax2.tick_params(which='minor', direction="in", top=True, right=True)
ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Apply sign(x) correction for profile plots (flip sign for x<0)
Hy_pert_numerical_x_corrected = Hy_pert_numerical_x * sign(x_profile)
Hy_pert_analytical_x_corrected = Hy_pert_analytical_x * sign(x_profile)
# Handle x=0 case
Hy_pert_numerical_x_corrected[x_profile == 0] = Hy_pert_numerical_x[x_profile == 0]
Hy_pert_analytical_x_corrected[x_profile == 0] = Hy_pert_analytical_x[x_profile == 0]

Hy_pert_numerical_y_corrected = Hy_pert_numerical_y * sign(y_profile)
Hy_pert_analytical_y_corrected = Hy_pert_analytical_y * sign(y_profile)
# Handle y=0 case
Hy_pert_numerical_y_corrected[y_profile == 0] = Hy_pert_numerical_y[y_profile == 0]
Hy_pert_analytical_y_corrected[y_profile == 0] = Hy_pert_analytical_y[y_profile == 0]

# Row 2, Col 1: X-axis profile comparison
ax3 = plt.subplot(2, 2, 3)
ax3.plot(x_profile, Hy_pert_numerical_x_corrected, 'k-', linewidth=2, label='NGSolve')
ax3.plot(x_profile, Hy_pert_analytical_x_corrected, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax3.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$\\mathrm{sign}(x) \\cdot H_{y,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax3.set_title('X-axis Profile ($H_y$ component)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: Y-axis profile comparison
ax4 = plt.subplot(2, 2, 4)
ax4.plot(y_profile, Hy_pert_numerical_y_corrected, 'k-', linewidth=2, label='NGSolve')
ax4.plot(y_profile, Hy_pert_analytical_y_corrected, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax4.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$\\mathrm{sign}(y) \\cdot H_{y,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax4.set_title('Y-axis Profile ($H_y$ component)', fontname='Times New Roman', fontsize=11)
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
