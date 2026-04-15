"""
Omega-Reduced Omega Method for 3D Magnetostatics with Kelvin Transformation
Problem: Magnetic cylinder (mu_r=100) in uniform z-directed background field
Using 1/8 symmetry model (x>=0, y>=0, z>=0)

The cylinder axis is aligned with the z-axis (external field direction).
No analytical solution exists for this geometry.

Symmetry conditions:
- x=0 plane: Neumann BC (dOmega/dn = 0) - natural BC, no explicit setting needed
- y=0 plane: Neumann BC (dOmega/dn = 0) - natural BC, no explicit setting needed
- z=0 plane: Dirichlet BC (Omega = 0) - because Omega_s = H0*z = 0 at z=0

Based on Omega_ReducedOmega.py implementation:
- Sign convention: B = mu * grad(Omega), H = grad(Omega)
- Source potential: Omega_s = H0 * z (so grad(Omega_s) = H_s)

Author: Claude Code
Date: 2025-01-04
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Omega-Reduced Omega Method - 3D Magnetic Cylinder (1/8 Model)")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
cylinder_radius = 0.3    # Magnetic cylinder radius [m]
cylinder_height = 1.0    # Magnetic cylinder height [m] (total, using half: z>=0)
kelvin_radius = 1.5      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1.0  # [A/m]

# Mesh parameters
maxh_cylinder = 0.05     # Finer mesh for magnetic cylinder
maxh_air = 0.10          # Moderate mesh for air region
fe_order = 2             # Finite element order

# Offset for exterior domain (z-direction)
offset_z = 4.0

half_height = cylinder_height / 2  # Only model z >= 0

print(f"\nProblem parameters:")
print(f"  Cylinder radius: {cylinder_radius} m")
print(f"  Cylinder half-height: {half_height} m (1/8 model, z>=0)")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")

# ============================================================
# Geometry Definition (OCC) - 1/8 Model
# ============================================================
print("\nCreating 1/8 geometry using OCC...")

# Create cutting boxes for 1/8 symmetry.
# Keep x >= 0, y >= 0, z >= 0.
# IMPORTANT: cut_x / cut_y must span z up past the OUTER sphere top
# (z = offset_z + kelvin_radius). Without this, the outer sphere's
# upper half is NOT cut and remains a hemisphere -- breaking the
# periodic BC because interior (1/8 sphere) and exterior (hemisphere)
# have incompatible face topology. Symptom: FreeDofs diff = 0 after
# Periodic(fes). Fixed 2026-04-16 for the 1/8 Kelvin-cylinder model.
z_hi = offset_z + kelvin_radius * 2
cut_x = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
            Pnt(0, kelvin_radius*2, z_hi))
cut_y = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
            Pnt(kelvin_radius*2, 0, z_hi))
cut_z = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
            Pnt(kelvin_radius*2, kelvin_radius*2, 0))

# Magnetic cylinder (full, then cut to 1/8)
# Cylinder from z=0 to z=half_height (only upper half)
mag_cylinder_full = Cylinder(Pnt(0, 0, 0), Z, r=cylinder_radius, h=half_height)
mag_cylinder = mag_cylinder_full - cut_x - cut_y  # 1/4 in xy plane

mag_cylinder.mat("magnetic")
mag_cylinder.maxh = maxh_cylinder

# Name faces
for face in mag_cylinder.faces:
    # Check face center to determine type
    fc = face.center
    if abs(fc.x) < 1e-6:
        face.name = "sym_x"  # x=0 symmetry plane
    elif abs(fc.y) < 1e-6:
        face.name = "sym_y"  # y=0 symmetry plane
    elif abs(fc.z) < 1e-6:
        face.name = "sym_z"  # z=0 symmetry plane
    elif abs(fc.z - half_height) < 1e-6:
        face.name = "cylinder"  # Top cap
    else:
        face.name = "cylinder"  # Curved surface

# Inner air domain (1/8 sphere minus 1/8 cylinder)
inner_sphere_full = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere = inner_sphere_full - cut_x - cut_y - cut_z  # 1/8 sphere

inner_sphere.maxh = maxh_air
# Name kelvin boundary first (spherical surface)
for face in inner_sphere.faces:
    fc = face.center
    dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
    if abs(dist - kelvin_radius) < kelvin_radius * 0.1:
        face.name = "kelvin_int"
        print(f"  kelvin_int face at ({fc.x:.2f}, {fc.y:.2f}, {fc.z:.2f})")
    elif abs(fc.x) < 1e-6:
        face.name = "sym_x"
    elif abs(fc.y) < 1e-6:
        face.name = "sym_y"
    elif abs(fc.z) < 1e-6:
        face.name = "sym_z"

inner_air = inner_sphere - mag_cylinder
inner_air.mat("air_inner")

# Name faces after boolean operation
for face in inner_air.faces:
    fc = face.center
    dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
    # Kelvin boundary is the spherical surface
    if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
        face.name = "kelvin_int"
    elif abs(fc.x) < 0.05:
        face.name = "sym_x"
    elif abs(fc.y) < 0.05:
        face.name = "sym_y"
    elif abs(fc.z) < 0.05:
        face.name = "sym_z"
    else:
        face.name = "cylinder"

# Exterior domain (1/8 Kelvin-transformed sphere)
outer_sphere_full = Sphere(Pnt(0, 0, offset_z), kelvin_radius)
# Cut for 1/8: x>=0, y>=0, z>=offset_z (upper octant of the offset sphere)
cut_z_ext = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
                Pnt(kelvin_radius*2, kelvin_radius*2, offset_z))
outer_sphere = outer_sphere_full - cut_x - cut_y - cut_z_ext  # 1/8 sphere

outer_sphere.maxh = maxh_air
outer_sphere.mat("air_outer")
# Name faces
for face in outer_sphere.faces:
    fc = face.center
    dist = sqrt(fc.x**2 + fc.y**2 + (fc.z - offset_z)**2)
    if abs(dist - kelvin_radius) < kelvin_radius * 0.3:  # Spherical boundary
        face.name = "kelvin_ext"
    elif abs(fc.x) < 0.1:
        face.name = "sym_x"
    elif abs(fc.y) < 0.1:
        face.name = "sym_y"
    elif abs(fc.z - offset_z) < 0.1:
        face.name = "sym_z_ext"

# GND vertex at center of exterior domain
vertex = Vertex(Pnt(0, 0, offset_z))
vertex.name = "GND"

# Glue all domains
geo = Glue([inner_air, mag_cylinder, outer_sphere, vertex])

# Name the solids
for i, solid in enumerate(geo.solids):
    if i == 0:
        solid.name = "air_inner"
    elif i == 1:
        solid.name = "magnetic"
    elif i == 2:
        solid.name = "air_outer"

# ===== IDENTIFY PERIODIC FACES =====
print("\nIdentifying periodic boundaries...")

print(f"  Number of solids: {len(geo.solids)}")
for i, solid in enumerate(geo.solids):
    print(f"  Solid[{i}] ({solid.name}): {len(solid.faces)} faces")

# Find kelvin_int and kelvin_ext faces
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

# Apply curved elements
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

# H1 space with:
# - Dirichlet BC at GND (infinity)
# - Dirichlet BC at z=0 symmetry plane (Omega = Omega_s = 0)
# Note: sym_x and sym_y are natural Neumann BC (no explicit setting)
fes_before = H1(mesh, order=fe_order, dirichlet="GND|sym_z|sym_z_ext")
fes = Periodic(fes_before)

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

# Kelvin-transformed permeability (Nagamine CEFC 2026 canonical):
#   mu_ext = mu_0 * (R/r')^2  for 3D spherical (conformal) Kelvin
# See examples/kelvin_transformation/CONVENTION.md.
from radia.kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf

mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=(0.0, 0.0, offset_z),
                                           R=kelvin_radius)
Mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)
# Keep mu_kelvin alias for post-processing
mu_kelvin = mu0 * mu_kelvin_factor

print(f"  air_inner: mu = mu0")
print(f"  air_outer: mu = (R/r')^2 * mu0 [Nagamine CEFC 2026]")
print(f"  magnetic: mu = {mu_r} * mu0")

# ============================================================
# Weak Form Setup
# ============================================================
print("\nSetting up weak form...")

# Source potential
Omega_s = H0 * z

# Source fields
Hs = CoefficientFunction((0.0, 0.0, H0))
Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

# Kelvin-transformed source field
r_exterior = sqrt(x**2 + y**2 + (z - offset_z)**2 + 1e-20)
Hz_exterior = -(r_exterior / kelvin_radius)**2 * H0
Bs_exterior = CoefficientFunction((0.0, 0.0, mu0 * Hz_exterior))

# Bilinear form
a = BilinearForm(fes)
a += Mu * grad(Omega) * grad(psi) * dx("magnetic")
a += Mu * grad(Omega) * grad(psi) * dx("air_inner")
a += Mu * grad(Omega) * grad(psi) * dx("air_outer")
a.Assemble()

# Set Dirichlet BC on cylinder boundary
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s, BND, mesh.Boundaries("cylinder"))

# Linear form
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
f.Assemble()

# Neumann BC (NEGATED normal for 3D)
normal = -specialcf.normal(mesh.dim)
f += (normal * Bs) * psi * ds("cylinder")
f.Assemble()

print("  Omega-Reduced Omega formulation with 1/8 symmetry")

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

Bs_dict = {
    "air_inner": Bs,
    "air_outer": CoefficientFunction((0.0, 0.0, mu0 * Hz_exterior)),
    "magnetic": CoefficientFunction((0.0, 0.0, 0.0))
}
Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

BField = Bt + Br + Bs_cf
grad_Omega = grad(gfu)

# ============================================================
# Perturbation Field Energy (multiply by 8 for full model)
# ============================================================
print("\nCalculating perturbation field energy...")

H_pert_total = grad(gfu) - Hs
energy_total_1_8 = Integrate(0.5 * (mu_r * mu0) * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)

H_pert_reduced = grad(Orr) - grad(Oxr)
energy_reduced_1_8 = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)

H_pert_kelvin = grad(Orr)
energy_kelvin_1_8 = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)

# Full model energy = 8 * (1/8 model energy)
energy_total_full = 8 * energy_total_1_8
energy_reduced_full = 8 * energy_reduced_1_8
energy_kelvin_full = 8 * energy_kelvin_1_8
energy_total_pert = energy_total_full + energy_reduced_full + energy_kelvin_full

print(f"\nPerturbation field energy (1/8 model x 8):")
print(f"  W_cylinder: {energy_total_full:.6e} J")
print(f"  W_air:      {energy_reduced_full + energy_kelvin_full:.6e} J")
print(f"  W_total:    {energy_total_pert:.6e} J")

# ============================================================
# Results
# ============================================================
print("\n" + "=" * 60)
print("RESULTS (1/8 Model)")
print("=" * 60)

# Sample points along z-axis (inside cylinder)
print(f"\nH_z along z-axis (inside cylinder):")
for z_val in linspace(0.05, half_height - 0.05, 5):
    try:
        mip = mesh(0.01, 0.01, z_val)
        Hz = grad_Omega[2](mip)
        print(f"  z={z_val:.2f}: Hz = {Hz:.6f} A/m")
    except:
        pass

# Sample points along x-axis (z=0.1)
print(f"\nH_z along x-axis (z=0.1):")
for x_val in linspace(0.05, kelvin_radius - 0.2, 6):
    try:
        mip = mesh(x_val, 0.01, 0.1)
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

fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

# Plot 1: Z-axis profile
ax1 = axes[0, 0]
z_vals = linspace(0.02, kelvin_radius - 0.1, 60)
Hz_z = []
for zv in z_vals:
    try:
        mip = mesh(0.02, 0.02, zv)
        if zv < half_height:
            Hz_z.append(grad_Omega[2](mip))
        else:
            B_val = BField(mip)
            Hz_z.append(B_val[2] / mu0)
    except:
        Hz_z.append(nan)
Hz_z = array(Hz_z)

ax1.plot(z_vals, Hz_z, 'b-', linewidth=2)
ax1.axvline(half_height, color='red', linestyle='--', alpha=0.7, label='Cylinder top')
ax1.axhline(H0, color='gray', linestyle=':', alpha=0.7)
ax1.set_xlabel('$z$ (m)', fontsize=11)
ax1.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax1.set_title('Z-axis Profile (1/8 model)', fontsize=12)
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: X-axis profile (z=0.1)
ax2 = axes[0, 1]
x_vals = linspace(0.02, kelvin_radius - 0.1, 60)
Hz_x = []
for xv in x_vals:
    try:
        mip = mesh(xv, 0.02, 0.1)
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
ax2.axhline(H0, color='gray', linestyle=':', alpha=0.7)
ax2.set_xlabel('$x$ (m)', fontsize=11)
ax2.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax2.set_title('X-axis Profile at z=0.1 (1/8 model)', fontsize=12)
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

# Plot 3: Summary
ax3 = axes[1, 0]
ax3.axis('off')
summary_text = f"""
1/8 Symmetry Model Results
==========================

Parameters:
  Cylinder radius: {cylinder_radius} m
  Cylinder half-height: {half_height} m
  Kelvin radius: {kelvin_radius} m
  mu_r = {mu_r}
  H0 = {H0} A/m

Mesh (1/8 model):
  Elements: {mesh.ne}
  DOFs: {fes.ndof}
  Order: {fe_order}

Perturbation Energy (full model = 8x):
  W_cylinder: {energy_total_full:.4e} J
  W_air:      {energy_reduced_full + energy_kelvin_full:.4e} J
  W_total:    {energy_total_pert:.4e} J

Periodic BC:
  FreeDofs reduction: {freedof_before - freedof_after}
"""
ax3.text(0.05, 0.95, summary_text, transform=ax3.transAxes,
         fontsize=10, fontfamily='monospace', verticalalignment='top')

# Plot 4: Comparison with full model expected values
ax4 = axes[1, 1]
ax4.axis('off')

# Expected from full model (from previous run)
W_full_expected = 1.659e-05  # From full 3D model

comparison_text = f"""
Comparison with Full Model
==========================

Full 3D model (expected):
  W_total ~= 1.659e-05 J

1/8 model x 8:
  W_total = {energy_total_pert:.4e} J

Ratio: {energy_total_pert / W_full_expected:.4f}

Note: 1/8 model uses symmetry:
  - x=0: Natural Neumann BC
  - y=0: Natural Neumann BC
  - z=0: Dirichlet BC (Omega=0)
"""
ax4.text(0.05, 0.95, comparison_text, transform=ax4.transAxes,
         fontsize=10, fontfamily='monospace', verticalalignment='top')

plt.tight_layout()

# Save figure
png_file = os.path.splitext(__file__)[0] + ".png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved: {png_file}")

try:
    os.startfile(png_file)
except:
    pass

print("\n" + "=" * 60)
print("Computation completed")
print("=" * 60)
