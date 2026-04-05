"""
H-formulation for 2D magnetostatics with perturbation potential and Kelvin transformation
QUADRUPOLE background field: H_s = (x, -y)
Using PERIODIC boundary conditions to connect interior and exterior domains
Based on: 2D_dipole_with_Kelvin.py
Updated: 2025-11-27
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 2D QUADRUPOLE with Kelvin Transformation (Periodic BC)")
print("="*60)

# ============================================================
# Geometry Definition (OCC 2D) - Two separate circular domains
# ============================================================
print("\nCreating geometry with periodic boundary conditions...")

# Parameters
circle_radius = 0.5	  # Magnetic circle radius [m]
kelvin_radius = 1.0	  # Kelvin transformation radius [m] (outer boundary)
maxh_fine = 0.03		# Fine mesh size [m]
plot_range = 1.1      # Plot range [m]

# Offset for exterior domain (placed separately)
offset_x = 3.0  # Offset to place exterior domain away from interior

print(f"Using Kelvin transformation with periodic BC:")
print(f"  - Inner domain: 0 < r < R = {kelvin_radius} m at (0, 0)")
print(f"  - Outer domain: 0 < r' < R = {kelvin_radius} m at ({offset_x}, 0)")
print(f"  - Transformation radius R = {kelvin_radius} m")

# Create WorkPlane
wp = WorkPlane()

# ===== INTERIOR DOMAIN (center at origin) =====
# Magnetic circle
mag_circle = wp.Circle(circle_radius).Face()
mag_circle.maxh = maxh_fine
mag_circle.name = "magnetic"

# Inner air domain (circle_radius < r < kelvin_radius)
inner_circle = wp.Circle(kelvin_radius).Face()
inner_circle.maxh = maxh_fine
# Name Kelvin boundary edge before subtraction (survives Boolean ops)
for edge in inner_circle.edges:
    edge.name = "kelvin_int"
inner_air = inner_circle - mag_circle
inner_air.name = "air_inner"

# ===== EXTERIOR DOMAIN (center at offset position) =====
# Move to offset position
wp = WorkPlane(Axes(Pnt(offset_x, 0, 0), n=Z, h=X))
outer_circle = wp.Circle(kelvin_radius).Face()
outer_circle.maxh = maxh_fine
# Name Kelvin boundary edge before Glue (survives Boolean ops)
for edge in outer_circle.edges:
    edge.name = "kelvin_ext"
outer_circle.name = "air_outer"

# ===== GND VERTEX (center of exterior domain) =====
# This represents r' = 0, which maps to r = infinity in the original space
vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

# Glue all domains (no hole in outer domain)
shape = Glue([inner_air, mag_circle, outer_circle, vertex])

# ===== NAME THE FACES AND EDGES =====
print("\nNaming faces and edges...")
print(f"  Number of faces: {len(shape.faces)}")

# Name the faces
shape.faces[0].name = "air_inner"
shape.faces[1].name = "magnetic"
shape.faces[2].name = "air_outer"

print("\nIdentifying periodic boundaries...")
# Find Kelvin boundary edges by name (named during geometry construction)
kelvin_int_edges = [e for e in shape.edges if e.name == "kelvin_int"]
kelvin_ext_edges = [e for e in shape.edges if e.name == "kelvin_ext"]

print(f"  Found {len(kelvin_int_edges)} interior Kelvin edges")
print(f"  Found {len(kelvin_ext_edges)} exterior Kelvin edges")

if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    for int_edge, ext_edge in zip(kelvin_int_edges, kelvin_ext_edges):
        int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
    print(f"  Periodic identification applied ({len(kelvin_int_edges)} edge pairs)")
else:
    raise RuntimeError(f"Could not find Kelvin boundary edges! "
                       f"(kelvin_int: {len(kelvin_int_edges)}, kelvin_ext: {len(kelvin_ext_edges)})")

# Create geometry
geo = OCCGeometry(shape, dim=2)

print(f"\nGeometry created:")
print(f"  Magnetic circle radius: {circle_radius} m at (0, 0)")
print(f"  Inner air domain: {circle_radius} m < r < {kelvin_radius} m")
print(f"  Outer air domain: 0 < r' < {kelvin_radius} m at ({offset_x}, 0)")
print(f"  Mesh size: {maxh_fine} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
ngmesh = geo.GenerateMesh(maxh=maxh_fine, grading=0.7)
mesh = Mesh(ngmesh)

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup with Periodic BC
# ============================================================
print("\nSetting up H-formulation with Periodic BC...")

# Create finite element space with Periodic BC and Dirichlet BC at GND
fes_before = H1(mesh, order=3, dirichlet_bbnd="GND")

# Check FreeDofs BEFORE Periodic
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])

fes = Periodic(fes_before)  # Apply periodic boundary conditions

# Check FreeDofs AFTER Periodic
freedof_after = sum([1 for d in fes.FreeDofs() if d])

print(f"  Number of DOFs: {fes.ndof}")
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")
if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability
mu_d = {"air_inner": 1*mu0, "air_outer": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field: H_s = (x, -y) A/m (QUADRUPOLE field)
# This satisfies div(H_s) = 0, making it physically valid
# In polar: H_r = r cos(2theta), H_theta = -r sin(2theta)
# Corresponds to potential φ_s = -(1/2)r^2 cos(2theta)
#
# Kelvin transformation for quadrupole field:
# For H_s = (x, -y) in physical space (r > R):
#   In polar coordinates: H_r = r cos(2theta), H_theta = r sin(2theta)
#
# After Kelvin transformation to computational space (r' < R):
#   The field components transform with sign reversal for 2D in-plane components
#   and the radial coordinate transforms as r' = R^2/r
#
# For quadrupole, the transformed field is: H' = (-x', y')
# where (x', y') are coordinates in the exterior (Kelvin) domain

# Interior background field - quadrupole in original coordinates
# For interior domain centered at (0,0)
Hs_x_inner = x  # quadrupole: x component
Hs_y_inner = -y  # quadrupole: -y component

# Detect which domain we're in
is_exterior = IfPos(x - offset_x/2, 1.0, 0.0)  # if x > offset_x/2, we're in exterior

# Exterior background field (Kelvin transformed)
# Coordinates relative to exterior domain center
x_ext = x - offset_x
y_ext = y

# For quadrupole, Kelvin transformation gives: H' = (-x', y')
Hs_x_outer = -x_ext
Hs_y_outer = y_ext

# Background field with domain switching
Hs_x = (1.0 - is_exterior) * Hs_x_inner + is_exterior * Hs_x_outer
Hs_y = (1.0 - is_exterior) * Hs_y_inner + is_exterior * Hs_y_outer

Hs = CoefficientFunction((Hs_x, Hs_y))

print(f"  Background field: H_s = (x, -y) (quadrupole) with Kelvin transformation")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = ∫(∇v)·(mu∇u)dOmega
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (PERTURBATION FORMULATION):
# 注意: Kelvin + 周期BCでは境界項不要
#       周期BCが自動的に連続性を保証するため
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx  # 体積積分のみ

a.Assemble()
f.Assemble()

print("  System assembled")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfu = GridFunction(fes)
c = Preconditioner(a, type="local")

solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat, tol=1e-8, printrates=True, maxsteps=10000)

print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# Compute perturbation field: H_pert = -grad(phi)
H = -grad(gfu)

# Analytical coefficients for quadrupole
# Interior: φ_pert = A r^2 cos2theta, where A = (mu_r-1)/(2(mu_r+1))
# Exterior: φ = B(a⁴/r^2)cos2theta
A_coeff = (mu_r - 1.0)/(2.0*(mu_r + 1.0))
B_coeff = A_coeff  # Same coefficient for continuity

# Evaluate at interior domain on x-axis
try:
    x_test = 0.3
    Hx_origin = H[0](mesh(x_test, 0))
    print(f"  Field at ({x_test}, 0) (interior): Hx = {Hx_origin:.6f} A/m")

    # Expected analytical value (perturbation field interior, x-component at y=0)
    # H_r = -2Ar cos2theta, at theta=0: H_x = H_r = -2*A*r*1
    Hx_analytical = -2.0 * A_coeff * x_test * 1.0  # cos(2*0) = 1
    print(f"  Analytical (interior): Hx = {Hx_analytical:.6f} A/m")
    print(f"  Relative error: {abs(Hx_origin - Hx_analytical)/abs(Hx_analytical)*100:.3f}%")
except:
    print("  Could not evaluate test point")

# ============================================================
# Field Evaluation on Grid (Interior Domain Only)
# ============================================================
print("\nEvaluating field on grid...")

# Grid for plotting (interior domain)
x = linspace(-plot_range, plot_range, 221)
y = linspace(-plot_range, plot_range, 221)
xx, yy = meshgrid(x, y)

# Initialize arrays
phi = zeros((len(y), len(x)))
Hx = zeros((len(y), len(x)))
Hy = zeros((len(y), len(x)))

for ny in range(len(y)):
	for nx in range(len(x)):
		r = sqrt(x[nx]**2 + y[ny]**2)
		# Only evaluate in interior domain (x < offset_x/2) and inside mesh
		if x[nx] < offset_x/2 and r < kelvin_radius - 0.01:
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

# ============================================================
# Profile Evaluation (Interior Domain) - Quadrupole
# ============================================================
print("\nComputing axis profiles (perturbation field)...")

profile_range = linspace(-plot_range, plot_range, 221)

# X-axis profile: Hx vs x at y=0 (quadrupole field)
x_profile = profile_range
Hx_pert_numerical_x = zeros(len(x_profile))
Hx_pert_analytical_x = zeros(len(x_profile))

for i, xval in enumerate(x_profile):
	r = abs(xval)
	if r < kelvin_radius - 0.01:  # Inside mesh domain
		try:
			mip = mesh(xval, 0)
			Hx_pert_numerical_x[i] = H[0](mip)
		except:
			Hx_pert_numerical_x[i] = nan
	else:
		Hx_pert_numerical_x[i] = nan

	# Analytical solution for H_s = (x, -y) at (x, 0): theta=0 or pi
	r = abs(xval)
	theta = arctan2(0, xval)  # theta=0 or pi
	if r < circle_radius:
		# Interior: H_r = -2Ar cos2theta, at y=0: H_x = H_r
		Hx_pert_analytical_x[i] = -2.0 * A_coeff * r * cos(2*theta)
	else:
		# Exterior: H_r = 2B(a⁴/r^3)cos2theta, at y=0: H_x = H_r
		Hx_pert_analytical_x[i] = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)

# Y-axis profile: Hy vs y at x=0 (quadrupole field)
y_profile = profile_range
Hy_pert_numerical_y = zeros(len(y_profile))
Hy_pert_analytical_y = zeros(len(y_profile))

for i, yval in enumerate(y_profile):
	r = abs(yval)
	if r < kelvin_radius - 0.01:  # Inside mesh domain
		try:
			mip = mesh(0, yval)
			Hy_pert_numerical_y[i] = H[1](mip)
		except:
			Hy_pert_numerical_y[i] = nan
	else:
		Hy_pert_numerical_y[i] = nan

	# Analytical solution for H_s = (x, -y) at (0, y): theta=pi/2 or -pi/2
	r = abs(yval)
	theta = arctan2(yval, 0)  # theta=±pi/2
	if r < circle_radius:
		# Interior: H_r = -2Ar cos2theta, H_theta = 2Ar sin2theta
		# At x=0: H_y = H_r sintheta + H_theta costheta
		Hr = -2.0 * A_coeff * r * cos(2*theta)
		Htheta = 2.0 * A_coeff * r * sin(2*theta)
		Hy_pert_analytical_y[i] = Hr * sin(theta) + Htheta * cos(theta)
	else:
		# Exterior: H_r = 2B(a⁴/r^3)cos2theta, H_theta = -2B(a⁴/r^3)sin2theta
		Hr = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)
		Htheta = -2.0 * B_coeff * (circle_radius**4 / r**3) * sin(2*theta)
		Hy_pert_analytical_y[i] = Hr * sin(theta) + Htheta * cos(theta)

# Error statistics
valid_idx_x = ~isnan(Hx_pert_numerical_x)
interior_idx_x = valid_idx_x & (abs(x_profile) < circle_radius)

print(f"\n  Validation results (X-axis, perturbation field Hx):")
print(f"  -" * 30)

if sum(interior_idx_x) > 0:
	interior_error = Hx_pert_numerical_x[interior_idx_x] - Hx_pert_analytical_x[interior_idx_x]
	max_err_int = max(abs(interior_error))
	rms_err_int = sqrt(mean(interior_error**2))
	# Avoid division by zero for interior where analytical might be small
	interior_max = max(abs(Hx_pert_analytical_x[interior_idx_x]))
	if interior_max > 0:
		rel_err_int = rms_err_int / interior_max * 100
	else:
		rel_err_int = 0
	print(f"  Interior (|x| < {circle_radius} m):")
	print(f"	Max error: {max_err_int:.6e} A/m")
	print(f"	RMS error: {rms_err_int:.6e} A/m ({rel_err_int:.3f}%)")

# ============================================================
# Analytical Flux Lines (2D) - Quadrupole
# ============================================================
print("\nComputing analytical flux lines...")

Hx_analytical = zeros(xx.shape)
Hy_analytical = zeros(xx.shape)

for ny in range(len(y)):
	for nx in range(len(x)):
		r = sqrt(x[nx]**2 + y[ny]**2)
		if r < 0.01:
			r = 0.01
		theta = arctan2(y[ny], x[nx])

		if r < circle_radius:
			# Interior: φ_pert = Ar^2cos2theta
			# H_r = -∂φ/∂r = -2Ar cos2theta
			# H_theta = -(1/r)∂φ/∂theta = 2Ar sin2theta
			Hr = -2.0 * A_coeff * r * cos(2*theta)
			Htheta = 2.0 * A_coeff * r * sin(2*theta)
			# Convert to Cartesian
			Hx_analytical[ny, nx] = Hr * cos(theta) - Htheta * sin(theta)
			Hy_analytical[ny, nx] = Hr * sin(theta) + Htheta * cos(theta)
		else:
			# Exterior: φ = B(a⁴/r^2)cos2theta
			# H_r = -∂φ/∂r = 2B(a⁴/r^3)cos2theta
			# H_theta = -(1/r)∂φ/∂theta = 2B(a⁴/r^3)sin2theta
			Hr = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)
			Htheta = 2.0 * B_coeff * (circle_radius**4 / r**3) * sin(2*theta)
			# Convert to Cartesian
			Hx_analytical[ny, nx] = Hr * cos(theta) - Htheta * sin(theta)
			Hy_analytical[ny, nx] = Hr * sin(theta) + Htheta * cos(theta)

# ============================================================
# Save Results
# ============================================================
print("\nSaving results...")

# Save to .mat file (interior domain data saved here, exterior will be added after evaluation)
from scipy.io import savemat
mat_data = {
    'xx': xx,
    'yy': yy,
    'Hx_analytical': Hx_analytical,
    'Hy_analytical': Hy_analytical,
    'Hx': Hx,
    'Hy': Hy
}

# ============================================================
# Evaluate Exterior Domain Field
# ============================================================
print("\nEvaluating exterior domain field...")

# Grid for exterior domain centered at (offset_x, 0)
x_ext = linspace(offset_x - plot_range, offset_x + plot_range, 221)
y_ext = linspace(-plot_range, plot_range, 221)
xx_ext, yy_ext = meshgrid(x_ext, y_ext)

# Evaluate H field in exterior domain
Hx_ext = zeros(xx_ext.shape)
Hy_ext = zeros(xx_ext.shape)

for ny in range(len(y_ext)):
	for nx in range(len(x_ext)):
		# Check if point is inside exterior domain (r' < R from offset center)
		r_from_offset = sqrt((x_ext[nx] - offset_x)**2 + y_ext[ny]**2)
		if r_from_offset < kelvin_radius - 0.05:  # Small margin to avoid boundary
			try:
				mip = mesh(x_ext[nx], y_ext[ny])
				Hx_ext[ny, nx] = H[0](mip)
				Hy_ext[ny, nx] = H[1](mip)
			except:
				Hx_ext[ny, nx] = nan
				Hy_ext[ny, nx] = nan
		else:
			Hx_ext[ny, nx] = nan
			Hy_ext[ny, nx] = nan

# Add exterior domain data to mat_data and save
mat_data['xx_ext'] = xx_ext
mat_data['yy_ext'] = yy_ext
mat_data['Hx_ext'] = Hx_ext
mat_data['Hy_ext'] = Hy_ext
mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, mat_data)
print(f"  MAT file saved to: {mat_file}")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
							  'bf':'serif:bold', 'fontset':'cm'})

# Create figure with 3x2 subplots
# Row 1: Analytical H streamline vs NGSolve H streamline
# Row 2: X-axis and Y-axis profile comparisons
# Row 3: Exterior B and H (NGSolve)
fig = plt.figure(figsize=(12, 15), dpi=150)

# Compute B field for exterior domain (2D: no spatial modulation)
Bx_ext = mu0 * Hx_ext
By_ext = mu0 * Hy_ext

# Row 1, Col 1: Analytical H field streamline
ax1 = plt.subplot(3, 2, 1)
strm1 = ax1.streamplot(xx, yy, Hx_analytical, Hy_analytical,
					   color='red', linewidth=1.0, density=1.5,
					   arrowsize=0.8, arrowstyle='->')
circle1 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
					 alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax1.add_patch(circle1)
kelvin_boundary1 = plt.Circle((0, 0), kelvin_radius, fill=False,
							  edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax1.add_patch(kelvin_boundary1)
ax1.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('Analytical $\\mathbf{H}_{\\mathrm{pert}}$ Streamline', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(-plot_range, plot_range)
ax1.set_ylim(-plot_range, plot_range)
ax1.minorticks_on()
ax1.tick_params(which='major', direction="in", top=True, right=True)
ax1.tick_params(which='minor', direction="in", top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 1, Col 2: NGSolve H field streamline
ax2 = plt.subplot(3, 2, 2)
strm2 = ax2.streamplot(xx, yy, Hx, Hy,
					   color='black', linewidth=1.0, density=1.5,
					   arrowsize=0.8, arrowstyle='->')
circle2 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
					 alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax2.add_patch(circle2)
kelvin_boundary2 = plt.Circle((0, 0), kelvin_radius, fill=False,
							  edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax2.add_patch(kelvin_boundary2)
ax2.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('NGSolve $\\mathbf{H}_{\\mathrm{pert}}$ Streamline', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(-plot_range, plot_range)
ax2.set_ylim(-plot_range, plot_range)
ax2.minorticks_on()
ax2.tick_params(which='major', direction="in", top=True, right=True)
ax2.tick_params(which='minor', direction="in", top=True, right=True)
ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Apply sign(x) correction for X-axis profile analytical solution only
Hx_pert_analytical_x_corrected = Hx_pert_analytical_x * sign(x_profile)
# Handle x=0 case
Hx_pert_analytical_x_corrected[x_profile == 0] = Hx_pert_analytical_x[x_profile == 0]

# Row 2, Col 1: X-axis profile comparison (Hx)
ax3 = plt.subplot(3, 2, 3)
ax3.plot(x_profile, Hx_pert_numerical_x, 'k-', linewidth=2, label='NGSolve')
ax3.plot(x_profile, Hx_pert_analytical_x_corrected, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax3.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$H_{x,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax3.set_title('X-axis Profile ($H_x$ component)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: Y-axis profile comparison (Hy)
ax4 = plt.subplot(3, 2, 4)
ax4.plot(y_profile, Hy_pert_numerical_y, 'k-', linewidth=2, label='NGSolve')
ax4.plot(y_profile, Hy_pert_analytical_y, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax4.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$H_{y,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
ax4.set_title('Y-axis Profile ($H_y$ component)', fontname='Times New Roman', fontsize=11)
ax4.minorticks_on()
ax4.tick_params(which='major', direction="in", top=True, right=True)
ax4.tick_params(which='minor', direction="in", top=True, right=True)
ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax4.legend(loc='best', fontsize=9, frameon=False)

# Row 3, Col 1: Exterior B field
ax5 = plt.subplot(3, 2, 5)
strm5 = ax5.streamplot(xx_ext, yy_ext, Bx_ext, By_ext,
					   color='darkblue', linewidth=1.0, density=1.5,
					   arrowsize=0.8, arrowstyle='->')
circle5 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
					 edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax5.add_patch(circle5)
ax5.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax5.set_title('Exterior: Flux Density $\\mathbf{B}$', fontname='Times New Roman', fontsize=11)
ax5.set_aspect('equal')
ax5.set_xlim(offset_x - plot_range, offset_x + plot_range)
ax5.set_ylim(-plot_range, plot_range)
ax5.minorticks_on()
ax5.tick_params(which='major', direction="in", top=True, right=True)
ax5.tick_params(which='minor', direction="in", top=True, right=True)
ax5.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 3, Col 2: Exterior H field
ax6 = plt.subplot(3, 2, 6)
strm6 = ax6.streamplot(xx_ext, yy_ext, Hx_ext, Hy_ext,
					   color='darkgreen', linewidth=1.0, density=1.5,
					   arrowsize=0.8, arrowstyle='->')
circle6 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
					 edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax6.add_patch(circle6)
ax6.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Exterior: Magnetic Field $\\mathbf{H}$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(offset_x - plot_range, offset_x + plot_range)
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
