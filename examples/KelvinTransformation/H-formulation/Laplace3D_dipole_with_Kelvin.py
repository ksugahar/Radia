"""
3D Laplace equation with Kelvin transformation for dipole field
Conventional formulation (not perturbation)

Boundary conditions:
- Magnetic sphere surface: φ = z (dipole potential)
- Periodic BC at r = R connecting interior and exterior domains

Goal: Verify that φ != 0 at domain boundaries (Kelvin transformation working)
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("3D Laplace Equation with Kelvin Transformation")
print("="*60)

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry...")

# Parameters
sphere_radius = 0.5  # Magnetic sphere radius [m]
kelvin_radius = 1.0  # Kelvin transformation radius [m]
offset_x = 3.0       # Offset for exterior domain
maxh_fine = 0.1      # Mesh size [m]
plot_range = 1.1     # Plot range [m]
box_size = 10.0      # Box size for intersection

# ===== INTERIOR DOMAIN (center at origin) =====
# Magnetic sphere - surface only (no interior mesh)
mag_sphere_surface = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere_surface.bc("dipole_bc")  # Set BC before subtraction

# Inner air domain (sphere_radius < r < kelvin_radius)
inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere.maxh = maxh_fine
# Name Kelvin boundary face before subtraction (survives Boolean ops)
for face in inner_sphere.faces:
    face.name = "kelvin_int"
inner_air = inner_sphere - mag_sphere_surface
inner_air.mat("air_inner")

# ===== EXTERIOR DOMAIN (center at offset position) =====
# Exterior domain: 0 < r' < kelvin_radius
outer_sphere = Sphere(Pnt(offset_x, 0, 0), kelvin_radius)
outer_sphere.maxh = maxh_fine
# Name Kelvin boundary face before Glue (survives Boolean ops)
for face in outer_sphere.faces:
    face.name = "kelvin_ext"
outer_sphere.mat("air_outer")

# ===== GND VERTEX (center of exterior domain) =====
# This represents r' = 0, which maps to r = infinity in the original space
vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

# Glue all domains
geo = Glue([inner_air, outer_sphere, vertex])

# ===== PERIODIC BOUNDARY IDENTIFICATION =====
print("\nIdentifying periodic boundaries...")
print(f"  Number of solids: {len(geo.solids)}")

for i, solid in enumerate(geo.solids):
	print(f"  Solid[{i}] ({solid.name}): {len(solid.faces)} faces")

# Find Kelvin boundary faces by name (named during geometry construction)
kelvin_int_face = None
kelvin_ext_face = None
for solid in geo.solids:
    for face in solid.faces:
        if face.name == "kelvin_int":
            kelvin_int_face = face
            print(f"  Found kelvin_int face in solid '{solid.name}'")
        elif face.name == "kelvin_ext":
            kelvin_ext_face = face
            print(f"  Found kelvin_ext face in solid '{solid.name}'")

if kelvin_int_face is not None and kelvin_ext_face is not None:
    kelvin_int_face.Identify(kelvin_ext_face, "periodic", IdentificationType.PERIODIC)
    print("  Periodic identification applied between kelvin_int and kelvin_ext")
else:
    raise RuntimeError(f"Could not find Kelvin boundary faces!")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# Check if dipole_bc exists
dipole_bc_exists = "dipole_bc" in mesh.GetBoundaries()
print(f"  'dipole_bc' boundary exists: {dipole_bc_exists}")

# ============================================================
# Problem Setup
# ============================================================
print("\nSetting up Laplace equation...")

# Finite element space with Dirichlet BC on dipole surface and GND vertex, plus Periodic BC
fes = H1(mesh, order=2, dirichlet="dipole_bc|GND")
fes = Periodic(fes)

print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties with Kelvin transformation
# Interior: mu = mu_0
# Exterior: mu' = (R^2/r'^2)·mu_0
x_local = x - offset_x
y_local = y
z_local = z
r_prime_sq = x_local**2 + y_local**2 + z_local**2
mu_outer = kelvin_radius**2 / (r_prime_sq + 1e-20) * mu0

mu_d = {"air_inner": 1*mu0, "air_outer": mu_outer}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

print(f"  Permeability:")
print(f"	Interior: mu = mu0")
print(f"	Exterior: mu' = (R^2/r'^2)*mu0")

# ============================================================
# Weak Form (Laplace Equation: div(mu∇φ) = 0)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = ∫mu(∇u)·(∇v) dOmega
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

# Linear form: f(v) = 0 (Laplace equation, no source)
f = LinearForm(fes)

# ============================================================
# Boundary Condition: φ = z/r^2 on magnetic sphere surface
# ============================================================
print("\nApplying boundary condition...")

gfu = GridFunction(fes)

# Set φ = z/r^2 on dipole_bc boundary
r_sq = x**2 + y**2 + z**2
phi_bc = z / r_sq
gfu.Set(phi_bc, definedon=mesh.Boundaries("dipole_bc"))

print(f"  Dipole BC: φ = z/r^2 on sphere surface (r = {sphere_radius} m)")

# Assemble with boundary conditions applied
a.Assemble()
f.Assemble()

# Modify rhs to account for non-homogeneous Dirichlet BC
r = gfu.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec

print("  System assembled with boundary conditions")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

c = Preconditioner(a, type="local")
gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r

print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# Evaluate at key points
print(f"\n  Potential values:")

# On sphere surface (should be φ = z/r^2)
test_points_surface = [
	(0.5, 0, 0, "Sphere (0.5,0,0)"),
	(0, 0, 0.5, "Sphere (0,0,0.5)"),
	(0, 0, -0.5, "Sphere (0,0,-0.5)")
]

print(f"\n  On sphere surface (r={sphere_radius}, should match z/r^2):")
for px, py, pz, label in test_points_surface:
	try:
		phi_val = gfu(mesh(px, py, pz))
		r_val = sqrt(px**2 + py**2 + pz**2)
		phi_analytical = pz / (r_val**2)
		print(f"	{label}: phi = {phi_val:.6f}, z/r^2 = {phi_analytical:.6f}, match = {abs(phi_val-phi_analytical)<0.01}")
	except Exception as e:
		print(f"	{label}: Error - {e}")

# In interior domain
test_points_interior = [
	(0.7, 0, 0),
	(1.0, 0, 0),
	(1.5, 0, 0),
	(0, 0, 1.0),
	(0, 0, 1.5)
]

print(f"\n  In interior domain (air_inner):")
for px, py, pz in test_points_interior:
	try:
		phi_val = gfu(mesh(px, py, pz))
		print(f"	({px:.1f}, {py:.1f}, {pz:.1f}): phi = {phi_val:.6f}")
	except Exception as e:
		print(f"	({px:.1f}, {py:.1f}, {pz:.1f}): Error")

# At periodic boundary
try:
	phi_boundary = gfu(mesh(2.0, 0, 0))
	print(f"\n  At boundary (2.0, 0, 0): phi = {phi_boundary:.6f}")
except:
	print(f"\n  At boundary (2.0, 0, 0): Cannot evaluate")

# Check if φ != 0 at domain boundaries
print(f"\n  VERIFICATION: phi at boundaries should be NON-ZERO for Kelvin transformation")

# ============================================================
# Visualization (4 panels like 2D version)
# ============================================================
print("\nGenerating visualization...")

import matplotlib.pyplot as plt

# Compute field H = -grad(phi)
H = CoefficientFunction((-grad(gfu)[0], -grad(gfu)[1], -grad(gfu)[2]))

# X-Z plane evaluation grid
x_vals = linspace(-plot_range, plot_range, 111)
z_vals = linspace(-plot_range, plot_range, 111)

# X-axis profile (z=0)
x_profile = x_vals
phi_numerical_x = []
phi_analytical_x = []

for xi in x_profile:
    r = abs(xi)
    if r > sphere_radius + 0.01 and r < kelvin_radius - 0.01:
        try:
            phi_num = gfu(mesh(xi, 0, 0))
            phi_numerical_x.append(phi_num)
        except:
            phi_numerical_x.append(nan)
    else:
        phi_numerical_x.append(nan)

    # Analytical: phi = z/r^2 = 0 for z=0
    phi_analytical_x.append(0.0)

phi_numerical_x = array(phi_numerical_x)
phi_analytical_x = array(phi_analytical_x)

# Z-axis profile (x=0)
z_profile = z_vals
phi_numerical_z = []
phi_analytical_z = []

for zi in z_profile:
    r = abs(zi)
    if r > sphere_radius + 0.01 and r < kelvin_radius - 0.01:
        try:
            phi_num = gfu(mesh(0, 0, zi))
            phi_numerical_z.append(phi_num)
        except:
            phi_numerical_z.append(nan)
    else:
        phi_numerical_z.append(nan)

    # Analytical: phi = z/(2*r^3) for dipole field
    if r > sphere_radius:
        phi_analytical_z.append(zi / (2.0 * r**3))
    else:
        phi_analytical_z.append(nan)

phi_numerical_z = array(phi_numerical_z)
phi_analytical_z = array(phi_analytical_z)

# X-Z plane data for flux lines
xx, zz = meshgrid(x_vals, z_vals)
Hx_xz = zeros(xx.shape)
Hz_xz = zeros(xx.shape)

for i in range(len(z_vals)):
    for j in range(len(x_vals)):
        r = sqrt(x_vals[j]**2 + z_vals[i]**2)
        if r > sphere_radius + 0.01 and r < kelvin_radius - 0.01:
            try:
                mip = mesh(x_vals[j], 0, z_vals[i])
                Hx_xz[i, j] = H[0](mip)
                Hz_xz[i, j] = H[2](mip)
            except:
                Hx_xz[i, j] = nan
                Hz_xz[i, j] = nan
        else:
            Hx_xz[i, j] = nan
            Hz_xz[i, j] = nan

# Analytical field (dipole at origin)
# phi = z/(2*r^3), H = -grad(phi)
# Hx = -d/dx[z/(2*r^3)] = 3*x*z/(2*r^5)
# Hz = -d/dz[z/(2*r^3)] = -1/(2*r^3) + 3*z^2/(2*r^5) = -(r^2 - 3*z^2)/(2*r^5)
Hx_analytical = zeros(xx.shape)
Hz_analytical = zeros(xx.shape)
for i in range(len(z_vals)):
    for j in range(len(x_vals)):
        r = sqrt(x_vals[j]**2 + z_vals[i]**2)
        if r > sphere_radius:
            # Dipole field: H = -grad(z/(2*r^3))
            r5 = r**5
            Hx_analytical[i, j] = 3.0 * x_vals[j] * z_vals[i] / (2.0 * r5)
            Hz_analytical[i, j] = -(r**2 - 3.0 * z_vals[i]**2) / (2.0 * r5)
        else:
            Hx_analytical[i, j] = nan
            Hz_analytical[i, j] = nan

# Create figure with 3x2 subplots
# Row 1: Interior B (flux density) and H (magnetic field)
# Row 2: X-axis and Z-axis profile comparisons
# Row 3: Exterior B and H
fig = plt.figure(figsize=(12, 15), dpi=150)

import matplotlib
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

# Compute B field (B = mu0 * H) for interior domain (vacuum)
Bx_xz = mu0 * Hx_xz
Bz_xz = mu0 * Hz_xz

# Create grid for exterior domain centered at (offset_x, 0, 0)
x_ext = linspace(offset_x - plot_range, offset_x + plot_range, 111)
z_ext = linspace(-plot_range, plot_range, 111)
xx_ext, zz_ext = meshgrid(x_ext, z_ext)

# Evaluate H field in exterior domain
Hx_ext = zeros(shape(xx_ext))
Hz_ext = zeros(shape(xx_ext))

for i in range(len(z_ext)):
    for j in range(len(x_ext)):
        # Check if point is inside exterior domain (r' < R from offset center)
        r_from_offset = sqrt((x_ext[j] - offset_x)**2 + z_ext[i]**2)
        if r_from_offset < kelvin_radius - 0.05:  # Small margin to avoid boundary
            try:
                mip = mesh(x_ext[j], 0, z_ext[i])
                Hx_ext[i, j] = H[0](mip)
                Hz_ext[i, j] = H[2](mip)
            except:
                Hx_ext[i, j] = nan
                Hz_ext[i, j] = nan
        else:
            Hx_ext[i, j] = nan
            Hz_ext[i, j] = nan

# Compute B field for exterior domain
Bx_ext = mu0 * Hx_ext
Bz_ext = mu0 * Hz_ext

# Row 1, Col 1: Interior H field (Analytical)
ax1 = plt.subplot(3, 2, 1)
strm1 = ax1.streamplot(xx, zz, Hx_analytical, Hz_analytical,
                       color='red', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle1 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Dipole source')
ax1.add_patch(circle1)
kelvin_boundary1 = plt.Circle((0, 0), kelvin_radius, fill=False,
                              edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax1.add_patch(kelvin_boundary1)
ax1.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax1.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax1.set_title('Interior: $\\mathbf{H}$ (Analytical)', fontname='Times New Roman', fontsize=11)
ax1.set_aspect('equal')
ax1.set_xlim(-plot_range, plot_range)
ax1.set_ylim(-plot_range, plot_range)
ax1.minorticks_on()
ax1.tick_params(which='major', direction="in", top=True, right=True)
ax1.tick_params(which='minor', direction="in", top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 1, Col 2: Interior H field (NGSolve)
ax2 = plt.subplot(3, 2, 2)
strm2 = ax2.streamplot(xx, zz, Hx_xz, Hz_xz,
                       color='black', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle2 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Dipole source')
ax2.add_patch(circle2)
kelvin_boundary2 = plt.Circle((0, 0), kelvin_radius, fill=False,
                              edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax2.add_patch(kelvin_boundary2)
ax2.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax2.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax2.set_title('Interior: $\\mathbf{H}$ (NGSolve)', fontname='Times New Roman', fontsize=11)
ax2.set_aspect('equal')
ax2.set_xlim(-plot_range, plot_range)
ax2.set_ylim(-plot_range, plot_range)
ax2.minorticks_on()
ax2.tick_params(which='major', direction="in", top=True, right=True)
ax2.tick_params(which='minor', direction="in", top=True, right=True)
ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 2, Col 1: X-axis profile
ax3 = plt.subplot(3, 2, 3)
ax3.plot(x_profile, phi_numerical_x, 'k-', linewidth=2, label='NGSolve')
ax3.plot(x_profile, phi_analytical_x, 'r--', linewidth=1.5, label='Analytical')
ax3.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax3.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax3.set_ylabel('$\\phi$', fontname='Times New Roman', fontsize=10)
ax3.set_title('X-axis Profile ($\\phi$ vs $x$ at $z=0$)', fontname='Times New Roman', fontsize=11)
ax3.minorticks_on()
ax3.tick_params(which='major', direction="in", top=True, right=True)
ax3.tick_params(which='minor', direction="in", top=True, right=True)
ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax3.legend(loc='best', fontsize=9, frameon=False)

# Row 2, Col 2: Z-axis profile
ax4 = plt.subplot(3, 2, 4)
ax4.plot(z_profile, phi_numerical_z, 'k-', linewidth=2, label='NGSolve')
ax4.plot(z_profile, phi_analytical_z, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
ax4.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax4.set_ylabel('$\\phi$', fontname='Times New Roman', fontsize=10)
ax4.set_title('Z-axis Profile ($\\phi$ vs $z$ at $x=0$)', fontname='Times New Roman', fontsize=11)
ax4.minorticks_on()
ax4.tick_params(which='major', direction="in", top=True, right=True)
ax4.tick_params(which='minor', direction="in", top=True, right=True)
ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
ax4.legend(loc='best', fontsize=9, frameon=False)

# Row 3, Col 1: Exterior B field
ax5 = plt.subplot(3, 2, 5)
strm5 = ax5.streamplot(xx_ext, zz_ext, Bx_ext, Bz_ext,
                       color='darkblue', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle5 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
                     edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax5.add_patch(circle5)
ax5.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
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
strm6 = ax6.streamplot(xx_ext, zz_ext, Hx_ext, Hz_ext,
                       color='darkgreen', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle6 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
                     edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax6.add_patch(circle6)
ax6.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Exterior: Magnetic Field $\\mathbf{H}$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(offset_x - plot_range, offset_x + plot_range)
ax6.set_ylim(-plot_range, plot_range)
ax6.minorticks_on()
ax6.tick_params(which='major', direction="in", top=True, right=True)
ax6.tick_params(which='minor', direction="in", top=True, right=True)
ax6.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

plt.tight_layout()
plot_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(plot_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {plot_file}")

os.startfile(plot_file)

print("\n" + "="*60)
print("Computation completed successfully")
print("="*60)
