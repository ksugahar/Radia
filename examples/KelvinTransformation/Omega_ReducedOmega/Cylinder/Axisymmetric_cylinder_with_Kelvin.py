"""
Axisymmetric Omega-Reduced Omega Method for Magnetostatics with Kelvin Transformation
Problem: Magnetic cylinder (mu_r=100) in uniform z-directed background field

The cylinder axis is aligned with the z-axis (external field direction).
No analytical solution exists for this geometry.

Based on Omega_ReducedOmega.py implementation:
- Sign convention: B = mu * grad(Omega), H = grad(Omega)
- Source potential: Omega_s = H0 * z (so grad(Omega_s) = H_s)
- Total region (magnetic cylinder): No source term in weak form
- Reduced region (air): Source term from Omega_s

Kelvin transformation:
- Maps infinite exterior domain to finite half-circle
- Permeability transformation: mu'(rho') = (R/rho')^2 * mu0
- Periodic BC couples interior (rho=R) with exterior (rho'=R)

Author: Claude Code
Date: 2025-01-04
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Axisymmetric Omega-Reduced Omega Method with Kelvin Transform")
print("Problem: Magnetic Cylinder in Uniform Z-Field")
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
maxh = 0.04              # Mesh size
fe_order = 3             # Finite element order

# Offset for exterior domain (z-direction for axisymmetric)
offset_z = 4.0

print(f"\nProblem parameters:")
print(f"  Cylinder radius: {cylinder_radius} m")
print(f"  Cylinder height: {cylinder_height} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")

# ============================================================
# Geometry Definition using HALF-CIRCLES (axisymmetric, r >= 0)
# x = r (radial), y = z (axial)
# ============================================================
print("\nCreating geometry with periodic boundary conditions...")

print(f"Using Kelvin transformation with periodic BC:")
print(f"  - Inner domain: 0 <= r < R = {kelvin_radius} m at (0, 0)")
print(f"  - Outer domain: 0 <= r' < R = {kelvin_radius} m at (0, {offset_z})")
print(f"  - Transformation radius R = {kelvin_radius} m")

# ===== INTERIOR DOMAIN (half-circle, r >= 0) =====
# Use Circle().Face() and cut with rectangle (same as sphere approach)
outer_half_int = Circle((0, 0), kelvin_radius).Face()
cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_int = outer_half_int - cutter_int  # Keep x >= 0 part

# Magnetic cylinder (rectangle in axisymmetric coordinates)
# Centered at origin, height along z-axis
half_height = cylinder_height / 2
cylinder_rect = MoveTo(0, -half_height).Rectangle(cylinder_radius, cylinder_height).Face()

# Air region (half-circle minus cylinder)
air_inner = outer_half_int - cylinder_rect

# Name boundaries for air region
# Use threshold based on distance from origin vs cylinder dimensions
for edge in air_inner.edges:
    x_center = edge.center.x
    y_center = edge.center.y
    dist = sqrt(x_center**2 + y_center**2)

    if x_center < 1e-6:  # On z-axis (r = 0)
        edge.name = "axis_int"
    elif abs(dist - kelvin_radius) < kelvin_radius * 0.2:  # On Kelvin boundary
        edge.name = "kelvin_int"
    else:
        # Cylinder surface
        edge.name = "cylinder"
air_inner.faces.name = "air_inner"

# Name cylinder edges
for edge in cylinder_rect.edges:
    if edge.center.x < 1e-6:
        edge.name = "axis_int"
    else:
        edge.name = "cylinder"
cylinder_rect.faces.name = "magnetic"

# ===== EXTERIOR DOMAIN (half-circle, r' >= 0, offset in z-direction) =====
outer_half_ext = Circle((0, offset_z), kelvin_radius).Face()
cutter_ext = MoveTo(-kelvin_radius - 1, offset_z - kelvin_radius - 1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_ext = outer_half_ext - cutter_ext  # Keep x' >= 0 part (r' >= 0)

for edge in outer_half_ext.edges:
    r_center = edge.center.x  # r' coordinate
    if r_center < 1e-6:  # On axis (r' = 0)
        edge.name = "axis_ext"
    else:
        edge.name = "kelvin_ext"
outer_half_ext.faces.name = "air_outer"

# ===== GND VERTEX (center of exterior domain - at infinity) =====
vertex = Vertex(Pnt(0, offset_z, 0))
vertex.name = "GND"

# Glue all domains
shape = Glue([air_inner, cylinder_rect, outer_half_ext, vertex])

# ===== IDENTIFY ALL PERIODIC BOUNDARY EDGE PAIRS =====
print("\nIdentifying periodic boundaries...")

# Find all kelvin edges
kelvin_int_edges = []
kelvin_ext_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_int_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_ext_edges.append(edge)

print(f"  Found {len(kelvin_int_edges)} interior kelvin edges")
print(f"  Found {len(kelvin_ext_edges)} exterior kelvin edges")

# Match edges for z-offset Kelvin
if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    matched_pairs = 0
    for int_edge in kelvin_int_edges:
        int_z = int_edge.center.y  # z coordinate
        for ext_edge in kelvin_ext_edges:
            ext_z = ext_edge.center.y - offset_z  # Relative z from center
            if (int_z > 0 and ext_z > 0) or (int_z < 0 and ext_z < 0):
                int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                print(f"  Identified: int(z={int_z:.3f}) <-> ext(z-offset={ext_z:.3f})")
                matched_pairs += 1
                break
    print(f"  Total matched pairs: {matched_pairs}")

# Create geometry
geo = OCCGeometry(shape, dim=2)

print(f"\nGeometry created:")
print(f"  Magnetic cylinder: radius={cylinder_radius} m, height={cylinder_height} m")
print(f"  Inner air domain: half-circle of radius {kelvin_radius} m")
print(f"  Outer air domain: half-circle of radius {kelvin_radius} m at (0, {offset_z})")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.7))

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
fes_before = H1(mesh, order=fe_order, dirichlet_bbnd="GND")
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

# Coordinate functions (x = r, y = z)
r_coord = x  # Radial coordinate
z_coord = y  # Axial coordinate

# ============================================================
# Material Properties
# ============================================================
print("\nSetting up material properties...")

# Distance squared from exterior domain center (offset in z-direction)
rho_prime_sq = r_coord**2 + (z_coord - offset_z)**2

# Transformed permeability for exterior domain
mu_kelvin = kelvin_radius**2 / (rho_prime_sq + 1e-20) * mu0

mu_dict = {
    "air_inner": mu0,
    "air_outer": mu_kelvin,
    "magnetic": mu_r * mu0
}
Mu = CoefficientFunction([mu_dict[mat] for mat in mesh.GetMaterials()])

print(f"  air_inner: mu = mu0")
print(f"  air_outer: mu = (R/rho')^2 * mu0 (Kelvin)")
print(f"  magnetic: mu = {mu_r} * mu0")

# ============================================================
# Weak Form Setup
# ============================================================
print("\nSetting up weak form...")

# Detect which domain we're in
is_exterior = IfPos(z_coord - offset_z/2, 1.0, 0.0)

# r-weight for axisymmetric formulation
r_weight = IfPos(r_coord - 1e-10, r_coord, 1e-10)

# Source potential: Omega_s = H0 * z
Omega_s = H0 * z_coord

# Source fields
Hs = CoefficientFunction((0.0, H0))
Bs = CoefficientFunction((0.0, mu0 * H0))

# For Kelvin-transformed exterior domain:
rho_prime = sqrt(r_coord**2 + (z_coord - offset_z)**2 + 1e-20)
Hz_exterior = -(rho_prime / kelvin_radius)**2 * H0
Hs_exterior = CoefficientFunction((0.0, Hz_exterior))
Bs_exterior = CoefficientFunction((0.0, mu0 * Hz_exterior))

# ===== Bilinear Form =====
a = BilinearForm(fes)
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("magnetic")
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_inner")
a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_outer")
a.Assemble()

# Dirichlet BC on Total/Reduced interface
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s, BND, mesh.Boundaries("cylinder"))

# Linear form
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * r_weight * dx("air_inner")
f.Assemble()

# Neumann boundary condition
normal = specialcf.normal(mesh.dim)
f += (normal * Bs) * psi * r_weight * ds("cylinder")
f.Assemble()

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
    "magnetic": CoefficientFunction((0.0, 0.0))
}
Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

# Total B field
BField = Bt + Br + Bs_cf

# H field: H = B/mu
HField = BField / Mu

# ============================================================
# Perturbation Field Energy Calculation
# ============================================================
print("\nCalculating perturbation field energy...")

# Perturbation field in Total region: H_pert = grad(Omega) - H_s
H_pert_total = grad(gfu) - Hs

# Perturbation field in Reduced region: H_pert = grad(Omega_r)
# Already defined as grad(Orr) - grad(Oxr) = grad(Omega) - grad(Omega_s)

# Energy in magnetic region (Total)
energy_magnetic = Integrate(
    0.5 * mu_r * mu0 * InnerProduct(H_pert_total, H_pert_total) * r_weight * 2 * pi * dx("magnetic"),
    mesh
)

# Energy in air_inner region (Reduced)
H_pert_air = grad(Orr) - grad(Oxr)
energy_air_inner = Integrate(
    0.5 * mu0 * InnerProduct(H_pert_air, H_pert_air) * r_weight * 2 * pi * dx("air_inner"),
    mesh
)

# Energy in air_outer region (Kelvin)
energy_air_outer = Integrate(
    0.5 * mu_kelvin * InnerProduct(H_pert_air, H_pert_air) * r_weight * 2 * pi * dx("air_outer"),
    mesh
)

energy_total = energy_magnetic + energy_air_inner + energy_air_outer

print(f"\n  Perturbation Field Energy:")
print(f"  " + "-" * 40)
print(f"  W_cylinder (interior): {energy_magnetic:.6e} J")
print(f"  W_air_inner:           {energy_air_inner:.6e} J")
print(f"  W_air_outer (Kelvin):  {energy_air_outer:.6e} J")
print(f"  W_total:               {energy_total:.6e} J")
print(f"  " + "-" * 40)
print(f"  Note: No analytical solution for cylinder geometry")

# ============================================================
# Results
# ============================================================
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)

# Sample points along z-axis (inside cylinder)
print(f"\nH_z along z-axis (r=0, inside cylinder):")
for z_val in linspace(-half_height + 0.05, half_height - 0.05, 5):
    try:
        mip = mesh(0.01, z_val)  # Small r offset to avoid axis singularity
        Hz = grad(gfu)[1](mip)
        print(f"  z={z_val:+.2f}: Hz = {Hz:.6f} A/m")
    except:
        pass

# Sample points along r-axis (z=0)
print(f"\nH_z along r-axis (z=0):")
for r_val in linspace(0.05, kelvin_radius - 0.1, 8):
    try:
        mip = mesh(r_val, 0)
        if r_val < cylinder_radius:
            Hz = grad(gfu)[1](mip)
            region = "cylinder"
        else:
            B_val = BField(mip)
            Hz = B_val[1] / mu0
            region = "air"
        print(f"  r={r_val:.2f} ({region}): Hz = {Hz:.6f} A/m")
    except:
        pass

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                             'bf': 'serif:bold', 'fontset': 'cm'})

fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=150)

# Plot 1: R-axis profile (z=0)
ax1 = axes[0, 0]
r_profile = linspace(0.02, kelvin_radius - 0.05, 100)
Hz_profile = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    try:
        mip = mesh(r_val, 0)
        if r_val < cylinder_radius:
            Hz_profile[i] = grad(gfu)[1](mip)
        else:
            B_val = BField(mip)
            Hz_profile[i] = B_val[1] / mu0
    except:
        Hz_profile[i] = nan

ax1.plot(r_profile, Hz_profile, 'b-', linewidth=2)
ax1.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax1.axhline(H0, color='gray', linestyle=':', alpha=0.7, label=f'$H_0$ = {H0} A/m')
ax1.set_xlabel('$r$ (m)', fontsize=11)
ax1.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax1.set_title('R-axis Profile (z=0)', fontsize=12)
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Z-axis profile (r=0.01)
ax2 = axes[0, 1]
z_profile = linspace(-kelvin_radius + 0.1, kelvin_radius - 0.1, 100)
Hz_zaxis = zeros(len(z_profile))

for i, z_val in enumerate(z_profile):
    try:
        mip = mesh(0.02, z_val)
        r_dist = sqrt(0.02**2 + z_val**2)
        if abs(z_val) < half_height and 0.02 < cylinder_radius:
            Hz_zaxis[i] = grad(gfu)[1](mip)
        else:
            B_val = BField(mip)
            Hz_zaxis[i] = B_val[1] / mu0
    except:
        Hz_zaxis[i] = nan

ax2.plot(z_profile, Hz_zaxis, 'b-', linewidth=2)
ax2.axvline(-half_height, color='red', linestyle='--', alpha=0.7)
ax2.axvline(half_height, color='red', linestyle='--', alpha=0.7, label='Cylinder ends')
ax2.axhline(H0, color='gray', linestyle=':', alpha=0.7, label=f'$H_0$ = {H0} A/m')
ax2.set_xlabel('$z$ (m)', fontsize=11)
ax2.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax2.set_title('Z-axis Profile (r~=0)', fontsize=12)
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

# Plot 3: 2D field visualization
ax3 = axes[0, 2]

r_grid = linspace(0.02, kelvin_radius - 0.05, 40)
z_grid = linspace(-kelvin_radius + 0.05, kelvin_radius - 0.05, 40)
[rr, zz] = meshgrid(r_grid, z_grid)

Hr_field = zeros(rr.shape)
Hz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        r_dist = sqrt(r_grid[ir]**2 + z_grid[iz]**2)
        in_cylinder = (r_grid[ir] < cylinder_radius) and (abs(z_grid[iz]) < half_height)
        if r_dist < kelvin_radius - 0.05:
            try:
                mip = mesh(r_grid[ir], z_grid[iz])
                if in_cylinder:
                    Hr_field[iz, ir] = grad(gfu)[0](mip)
                    Hz_field[iz, ir] = grad(gfu)[1](mip)
                else:
                    B_val = BField(mip)
                    Hr_field[iz, ir] = B_val[0] / mu0
                    Hz_field[iz, ir] = B_val[1] / mu0
            except:
                Hr_field[iz, ir] = nan
                Hz_field[iz, ir] = nan
        else:
            Hr_field[iz, ir] = nan
            Hz_field[iz, ir] = nan

strm = ax3.streamplot(rr, zz, Hr_field, Hz_field,
                      color='black', linewidth=0.8, density=1.5,
                      arrowsize=0.7, arrowstyle='->')

# Draw cylinder
rect = plt.Rectangle((0, -half_height), cylinder_radius, cylinder_height,
                      fill=True, facecolor='lightblue', edgecolor='red',
                      linewidth=2, alpha=0.5)
ax3.add_patch(rect)

# Draw Kelvin boundary
theta_circle = linspace(0, pi, 100)
r_kelvin_plot = kelvin_radius * sin(theta_circle)
z_kelvin_plot = kelvin_radius * cos(theta_circle)
ax3.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')

ax3.set_xlabel('$r$ (m)', fontsize=11)
ax3.set_ylabel('$z$ (m)', fontsize=11)
ax3.set_title('Total Field $\\mathbf{H}$', fontsize=12)
ax3.set_aspect('equal')
ax3.set_xlim(0, kelvin_radius + 0.1)
ax3.set_ylim(-kelvin_radius - 0.1, kelvin_radius + 0.1)

# Plot 4: Bz profile (r-axis)
ax4 = axes[1, 0]
Bz_profile = zeros(len(r_profile))

for i, r_val in enumerate(r_profile):
    try:
        mip = mesh(r_val, 0)
        B_val = BField(mip)
        Bz_profile[i] = B_val[1]
    except:
        Bz_profile[i] = nan

ax4.plot(r_profile, Bz_profile / mu0, 'b-', linewidth=2)
ax4.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax4.set_xlabel('$r$ (m)', fontsize=11)
ax4.set_ylabel('$B_z / \\mu_0$ (A/m)', fontsize=11)
ax4.set_title('$B_z$ Profile (z=0)', fontsize=12)
ax4.legend(loc='best', fontsize=9)
ax4.grid(True, alpha=0.3)

# Plot 5: Perturbation field
ax5 = axes[1, 1]

# Perturbation Hz along r-axis
Hz_pert_profile = zeros(len(r_profile))
for i, r_val in enumerate(r_profile):
    try:
        mip = mesh(r_val, 0)
        if r_val < cylinder_radius:
            H_pert = grad(gfu)(mip) - array([0, H0])
            Hz_pert_profile[i] = H_pert[1]
        else:
            Hz_pert_profile[i] = (grad(Orr) - grad(Oxr))[1](mip)
    except:
        Hz_pert_profile[i] = nan

ax5.plot(r_profile, Hz_pert_profile, 'b-', linewidth=2)
ax5.axvline(cylinder_radius, color='red', linestyle='--', alpha=0.7, label='Cylinder boundary')
ax5.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax5.set_xlabel('$r$ (m)', fontsize=11)
ax5.set_ylabel('$H_{z,pert}$ (A/m)', fontsize=11)
ax5.set_title('Perturbation Field $H_{z,pert}$ (z=0)', fontsize=12)
ax5.legend(loc='best', fontsize=9)
ax5.grid(True, alpha=0.3)

# Plot 6: Summary info
ax6 = axes[1, 2]
ax6.axis('off')

summary_text = f"""
Omega-Reduced Omega Method (Axisymmetric) + Kelvin
==================================================
Problem: Magnetic Cylinder in Uniform Z-Field

Parameters:
  Cylinder radius: {cylinder_radius} m
  Cylinder height: {cylinder_height} m
  Kelvin radius: {kelvin_radius} m
  mu_r = {mu_r}
  H0 = {H0} A/m

Mesh:
  Elements: {mesh.ne}
  DOFs: {fes.ndof}
  Order: {fe_order}

Perturbation Field Energy:
  W_cylinder: {energy_magnetic:.4e} J
  W_air:      {energy_air_inner + energy_air_outer:.4e} J
  W_total:    {energy_total:.4e} J

Note: No analytical solution exists
for cylinder geometry.
"""

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
