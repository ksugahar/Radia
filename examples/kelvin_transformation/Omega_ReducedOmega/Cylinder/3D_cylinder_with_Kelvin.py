"""
Omega-Reduced Omega Method for 3D Magnetostatics with Kelvin Transformation
Problem: Magnetic cylinder (mu_r=100) in uniform z-directed background field

The cylinder axis is aligned with the z-axis (external field direction).
No analytical solution exists for this geometry.

Based on Omega_ReducedOmega.py implementation:
- Sign convention: B = mu * grad(Omega), H = grad(Omega)
- Source potential: Omega_s = H0 * z (so grad(Omega_s) = H_s)
- Total region (magnetic cylinder): No source term in weak form
- Reduced region (air): Source term from Omega_s

Kelvin transformation:
- Maps infinite exterior domain to finite sphere
- Permeability transformation: mu'(r') = (R/r')^2 * mu0
- Periodic BC couples interior (r=R) with exterior (r'=R)
- Dirichlet BC at exterior center (r'=0 -> r=infinity): Omega = 0

IMPORTANT: 3D Normal Direction
- In NGSolve 3D, specialcf.normal(mesh.dim) on the cylinder boundary points INWARD
  (from air into magnetic cylinder), but Omega-Reduced Omega requires OUTWARD normal.
- Therefore, we NEGATE the normal: normal = -specialcf.normal(mesh.dim)

Author: Claude Code
Date: 2025-01-04
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Omega-Reduced Omega Method - 3D Magnetic Cylinder with Kelvin Transform")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
cylinder_radius = 0.3    # Magnetic cylinder radius [m]
cylinder_height = 1.0    # Magnetic cylinder height [m] (total, centered at z=0)
kelvin_radius = 1.5      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1.0  # [A/m]

# Mesh parameters
maxh_cylinder = 0.06     # Finer mesh for magnetic cylinder
maxh_air = 0.12          # Moderate mesh for air region
fe_order = 2             # Finite element order

# Offset for exterior domain (z-direction)
offset_z = 4.0

print(f"\nProblem parameters:")
print(f"  Cylinder radius: {cylinder_radius} m")
print(f"  Cylinder height: {cylinder_height} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry using OCC...")

half_height = cylinder_height / 2

# Magnetic cylinder at origin, axis aligned with z
mag_cylinder = Cylinder(Pnt(0, 0, -half_height), Z, r=cylinder_radius, h=cylinder_height)
mag_cylinder.mat("magnetic")
mag_cylinder.maxh = maxh_cylinder
# Name the cylinder surfaces
for face in mag_cylinder.faces:
    face.name = "cylinder"

# Inner air domain (sphere minus cylinder)
inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere.maxh = maxh_air
for face in inner_sphere.faces:
    face.name = "kelvin_int"

inner_air = inner_sphere - mag_cylinder
inner_air.mat("air_inner")

# Exterior domain (Kelvin-transformed, centered at z-offset position)
outer_sphere = Sphere(Pnt(0, 0, offset_z), kelvin_radius)
outer_sphere.maxh = maxh_air
outer_sphere.mat("air_outer")
for face in outer_sphere.faces:
    face.name = "kelvin_ext"

# GND vertex at center of exterior domain (represents r'=0 -> r=infinity)
vertex = Vertex(Pnt(0, 0, offset_z))
vertex.name = "GND"

# Glue all domains
geo = Glue([inner_air, mag_cylinder, outer_sphere, vertex])

# Name the solids
geo.solids[0].name = "air_inner"
geo.solids[1].name = "magnetic"
geo.solids[2].name = "air_outer"

# ===== IDENTIFY PERIODIC FACES =====
print("\nIdentifying periodic boundaries...")

print(f"  Number of solids: {len(geo.solids)}")
for i, solid in enumerate(geo.solids):
    print(f"  Solid[{i}] ({solid.name}): {len(solid.faces)} faces")

# Find kelvin_int and kelvin_ext faces by name
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
    print(f"  WARNING: Could not find periodic faces!")

# Generate mesh
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5))

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
fes_before = H1(mesh, order=fe_order, dirichlet="GND")
fes = Periodic(fes_before)  # Apply periodic BC

# Check if Periodic BC is working
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

# ============================================================
# Material Properties
# ============================================================
print("\nSetting up material properties...")

# Distance from exterior domain center (z-offset)
r_prime_sq = x**2 + y**2 + (z - offset_z)**2
r_prime = sqrt(r_prime_sq + 1e-20)

# Kelvin-transformed permeability
mu_kelvin = (kelvin_radius / r_prime)**2 * mu0

mu_dict = {
    "air_inner": mu0,
    "air_outer": mu_kelvin,
    "magnetic": mu_r * mu0
}
Mu = CoefficientFunction([mu_dict[mat] for mat in mesh.GetMaterials()])

print(f"  air_inner: mu = mu0")
print(f"  air_outer: mu = (R/r')^2 * mu0 (Kelvin)")
print(f"  magnetic: mu = {mu_r} * mu0")

# ============================================================
# Weak Form Setup
# ============================================================
print("\nSetting up weak form...")

# Source potential Omega_s = H0 * z
Omega_s = H0 * z

# Source magnetic field and B field
Hs = CoefficientFunction((0.0, 0.0, H0))
Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

# For Kelvin-transformed exterior domain (z-offset):
r_exterior = sqrt(x**2 + y**2 + (z - offset_z)**2 + 1e-20)
Hz_exterior = -(r_exterior / kelvin_radius)**2 * H0
Hs_exterior = CoefficientFunction((0.0, 0.0, Hz_exterior))
Bs_exterior = CoefficientFunction((0.0, 0.0, mu0 * Hz_exterior))

# Bilinear form
a = BilinearForm(fes)
a += Mu * grad(Omega) * grad(psi) * dx("magnetic")
a += Mu * grad(Omega) * grad(psi) * dx("air_inner")
a += Mu * grad(Omega) * grad(psi) * dx("air_outer")
a.Assemble()

# Create GridFunction for solution
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s, BND, mesh.Boundaries("cylinder"))

# Linear form
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
f.Assemble()

# Normal for Neumann BC (NEGATED for correct direction in 3D)
normal = -specialcf.normal(mesh.dim)
f += (normal * Bs) * psi * ds("cylinder")
f.Assemble()

print("  Omega-Reduced Omega formulation:")
print("  Total region (magnetic): H = grad(Omega), no source")
print("  Reduced region (air): Source from Omega_s")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
gfu = gfOmega

print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

fesOt = H1(mesh, order=fe_order, definedon="magnetic")
fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")

Ot = GridFunction(fesOt)
Orr = GridFunction(fesOr)
Oxr = GridFunction(fesOr)

Ot.Set(gfu, VOL, definedon="magnetic")
Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr.Set(Omega_s, BND, mesh.Boundaries("cylinder"))

# B field computation
Bt = grad(Ot) * Mu
Br = (grad(Orr) - grad(Oxr)) * mu0

# Source field Bs is only in Reduced region
Bs_dict = {
    "air_inner": Bs,
    "air_outer": Bs_exterior,
    "magnetic": CoefficientFunction((0.0, 0.0, 0.0))
}
Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

# Total B field
BField = Bt + Br + Bs_cf

# grad(Omega) for validation
grad_Omega = grad(gfu)

# ============================================================
# Perturbation Field Energy Calculation
# ============================================================
print("\nCalculating perturbation field energy...")

# Total region (magnetic cylinder)
H_pert_total = grad(gfu) - Hs
energy_total = Integrate(0.5 * (mu_r * mu0) * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)

# Reduced region (air_inner)
H_pert_reduced = grad(Orr) - grad(Oxr)
energy_reduced = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)

# Kelvin region (air_outer)
H_pert_kelvin = grad(Orr)
energy_kelvin = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)

# Total perturbation energy
energy_total_pert = energy_total + energy_reduced + energy_kelvin

print(f"\nPerturbation field energy:")
print(f"  Cylinder region (magnetic): W_cyl = {energy_total:.6e} J")
print(f"  Air inner region:           W_air = {energy_reduced:.6e} J")
print(f"  Kelvin region (air_outer):  W_kel = {energy_kelvin:.6e} J")
print(f"  Total perturbation energy:  W_tot = {energy_total_pert:.6e} J")
print(f"\nNote: No analytical solution exists for cylinder geometry")

# ============================================================
# Results
# ============================================================
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)

# Sample points along z-axis (inside cylinder)
print(f"\nH_z along z-axis (inside cylinder, r=0):")
for z_val in linspace(-half_height + 0.05, half_height - 0.05, 5):
    try:
        mip = mesh(0.01, 0, z_val)
        Hz = grad_Omega[2](mip)
        print(f"  z={z_val:+.2f}: Hz = {Hz:.6f} A/m")
    except:
        pass

# Sample points along x-axis (z=0)
print(f"\nH_z along x-axis (z=0):")
for x_val in linspace(0.05, kelvin_radius - 0.1, 8):
    try:
        mip = mesh(x_val, 0, 0)
        if x_val < cylinder_radius:
            Hz = grad_Omega[2](mip)
            region = "cylinder"
        else:
            B_val = BField(mip)
            Hz = B_val[2] / mu0
            region = "air"
        print(f"  x={x_val:.2f} ({region}): Hz = {Hz:.6f} A/m")
    except:
        pass

# ============================================================
# VTK Output
# ============================================================
print("\nSaving VTK output...")

vtk_file = os.path.splitext(__file__)[0] + "_result"
VTKOutput(ma=mesh, coefs=[Mu, gfu, grad(gfu)],
          names=["mu", "Omega", "grad_Omega"],
          filename=vtk_file).Do()
print(f"  VTK file saved: {vtk_file}.vtu")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                             'bf': 'serif:bold', 'fontset': 'cm'})

fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=150)

# Plot 1: Z-axis profile (inside cylinder)
ax1 = axes[0, 0]
z_vals = linspace(-kelvin_radius + 0.1, kelvin_radius - 0.1, 80)
Hz_z = []
for zv in z_vals:
    try:
        mip = mesh(0.02, 0, zv)
        in_cylinder = abs(zv) < half_height
        if in_cylinder:
            Hz_z.append(grad_Omega[2](mip))
        else:
            B_val = BField(mip)
            Hz_z.append(B_val[2] / mu0)
    except:
        Hz_z.append(nan)
Hz_z = array(Hz_z)

ax1.plot(z_vals, Hz_z, 'b-', linewidth=2)
ax1.axvline(-half_height, color='red', linestyle='--', alpha=0.7)
ax1.axvline(half_height, color='red', linestyle='--', alpha=0.7, label='Cylinder ends')
ax1.axhline(H0, color='gray', linestyle=':', alpha=0.7, label=f'$H_0$ = {H0} A/m')
ax1.set_xlabel('$z$ (m)', fontsize=11)
ax1.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax1.set_title('Z-axis Profile (x=y=0)', fontsize=12)
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: X-axis profile (z=0)
ax2 = axes[0, 1]
x_vals = linspace(0.02, kelvin_radius - 0.1, 80)
Hz_x = []
for xv in x_vals:
    try:
        mip = mesh(xv, 0, 0)
        if xv < cylinder_radius:
            Hz_x.append(grad_Omega[2](mip))
        else:
            B_val = BField(mip)
            Hz_x.append(B_val[2] / mu0)
    except:
        Hz_x.append(nan)
Hz_x = array(Hz_x)

ax2.plot(x_vals, Hz_x, 'b-', linewidth=2)
ax2.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax2.axhline(H0, color='gray', linestyle=':', alpha=0.7, label=f'$H_0$ = {H0} A/m')
ax2.set_xlabel('$x$ (m)', fontsize=11)
ax2.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax2.set_title('X-axis Profile (y=z=0)', fontsize=12)
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

# Plot 3: Energy summary
ax3 = axes[0, 2]
ax3.axis('off')
energy_text = f"""Perturbation Field Energy

Numerical Results:
  W_cylinder = {energy_total:.4e} J
  W_air      = {energy_reduced + energy_kelvin:.4e} J
  W_total    = {energy_total_pert:.4e} J

Parameters:
  Cylinder radius: {cylinder_radius} m
  Cylinder height: {cylinder_height} m
  mu_r = {mu_r}
  H0 = {H0} A/m

Note: No analytical solution
exists for cylinder geometry.
"""
ax3.text(0.1, 0.5, energy_text, transform=ax3.transAxes, fontsize=10,
         verticalalignment='center', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 4: Bz profile (x-axis)
ax4 = axes[1, 0]
Bz_x = []
for xv in x_vals:
    try:
        mip = mesh(xv, 0, 0)
        B_val = BField(mip)
        Bz_x.append(B_val[2])
    except:
        Bz_x.append(nan)
Bz_x = array(Bz_x)

ax4.plot(x_vals, Bz_x / mu0, 'b-', linewidth=2)
ax4.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax4.set_xlabel('$x$ (m)', fontsize=11)
ax4.set_ylabel('$B_z / \\mu_0$ (A/m)', fontsize=11)
ax4.set_title('$B_z$ Profile (y=z=0)', fontsize=12)
ax4.legend(loc='best', fontsize=9)
ax4.grid(True, alpha=0.3)

# Plot 5: Perturbation Hz (x-axis)
ax5 = axes[1, 1]
Hz_pert_x = []
for xv in x_vals:
    try:
        mip = mesh(xv, 0, 0)
        if xv < cylinder_radius:
            Hz_pert_x.append(grad_Omega[2](mip) - H0)
        else:
            Hz_pert_x.append((grad(Orr) - grad(Oxr))[2](mip))
    except:
        Hz_pert_x.append(nan)
Hz_pert_x = array(Hz_pert_x)

ax5.plot(x_vals, Hz_pert_x, 'b-', linewidth=2)
ax5.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax5.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax5.set_xlabel('$x$ (m)', fontsize=11)
ax5.set_ylabel('$H_{z,pert}$ (A/m)', fontsize=11)
ax5.set_title('Perturbation Field $H_{z,pert}$ (y=z=0)', fontsize=12)
ax5.legend(loc='best', fontsize=9)
ax5.grid(True, alpha=0.3)

# Plot 6: Summary info
ax6 = axes[1, 2]
ax6.axis('off')
summary_text = f"""
3D Omega-Reduced Omega Method + Kelvin
======================================
Problem: Magnetic Cylinder in Uniform Z-Field

Mesh:
  Elements: {mesh.ne}
  DOFs: {fes.ndof}
  Order: {fe_order}

Periodic BC:
  FreeDofs before: {freedof_before}
  FreeDofs after: {freedof_after}
  Reduction: {freedof_before - freedof_after}

Interior Hz (z=0, x=0.1):
"""
try:
    mip = mesh(0.1, 0, 0)
    Hz_center = grad_Omega[2](mip)
    summary_text += f"  Hz = {Hz_center:.6f} A/m"
except:
    summary_text += "  (not available)"

ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
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
