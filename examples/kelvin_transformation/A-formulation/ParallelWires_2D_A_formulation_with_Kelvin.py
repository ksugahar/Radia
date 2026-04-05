"""
2D A-formulation with Kelvin transformation for parallel wires
平行2線問題（2次元）のA法 + Kelvin変換

Problem:
  Two parallel wires carrying opposite currents (+I and -I)
  Located at (d/2, 0) and (-d/2, 0)

2D A-formulation:
  Variable: A_z (z-component of vector potential)
  Equation: -div(nu * grad(A_z)) = J_z
  B = curl(A) => Bx = dAz/dy, By = -dAz/dx

Kelvin transformation (full circle):
  - Interior domain: circle centered at origin, radius a
  - Exterior domain: circle centered at (x_offset, 0), radius a
  - Mapping: r' = a^2/r (inversion)
  - For nu: nu' = nu0 * (r'/a)^2 = nu0 * (a/r)^2

Analytical solution for parallel wires:
  A_z = (mu0 * I / 2pi) * ln(r2/r1)
  where r1, r2 are distances from each wire
"""
import os
from numpy import pi, sqrt, cos, sin, linspace, zeros, nan, isnan, meshgrid, array, log, mean, arctan2
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("2D A-formulation with Kelvin Transformation")
print("Parallel Wires Problem")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
wire_distance = 1.4     # Distance between wires (d) [m] -> wires at x=+/-0.8m
wire_radius = 0.02      # Wire cross-section radius [m] - larger radius
a = 1.0                 # Kelvin boundary radius [m]
y_offset = 2.5          # Y-offset for exterior domain [m] (like Z-offset in axisym)
maxh = 0.02             # Mesh size [m] (moderate for debugging)
order = 2               # Polynomial order

mu0 = 4 * pi * 1e-7     # Vacuum permeability [H/m]
nu0 = 1 / mu0           # Vacuum reluctivity [m/H]
I_total = 1.0           # Total current [A]
J0 = I_total / (pi * wire_radius**2)  # Current density [A/m^2]

print(f"\nParameters:")
print(f"  Wire distance: {wire_distance} m")
print(f"  Wire radius: {wire_radius} m")
print(f"  Kelvin boundary: a = {a} m")
print(f"  Y-offset: {y_offset} m")
print(f"  Total current: +/- {I_total} A")
print(f"  Current density: {J0:.2e} A/m^2")

# ============================================================
# Geometry with Y-offset (like Z-offset in axisymmetric)
# KEY: x-coordinate is the SAME in both domains!
# ============================================================
print("\nCreating geometry with Y-offset...")

# Interior domain: full circle at origin
wp1 = WorkPlane()
inner_circle = wp1.Circle(a).Face()
inner_circle.name = "air_inner"

# Wire 1 (positive current) at (+d/2, 0)
wire1_full = MoveTo(wire_distance/2, 0).Circle(wire_radius).Face()
wire1_full.name = "wire_minus"
wire1_full.maxh = maxh / 3

# Wire 2 (negative current) at (-d/2, 0)
wire2_full = MoveTo(-wire_distance/2, 0).Circle(wire_radius).Face()
wire2_full.name = "wire_plus"
wire2_full.maxh = maxh / 3

# Exterior domain: full circle at (0, y_offset) - Y direction offset!
wp2 = WorkPlane(Axes((0, y_offset, 0), n=Z, h=X))
outer_circle = wp2.Circle(a).Face()
outer_circle.name = "air_outer"

# GND point at center of exterior domain (represents infinity)
gnd_point = Vertex(Pnt(0, y_offset, 0))
gnd_point.name = "GND"

# Subtract wires from inner air
inner_air = inner_circle - wire1_full - wire2_full
inner_air.name = "air_inner"

# Name boundaries and find Kelvin edges
print("  Naming boundaries...")

# Find Kelvin boundary arc in inner_air
kelvin_inner_edges = []
for edge in inner_air.edges:
    cx = edge.center.x
    cy = edge.center.y
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + v0.p.y**2)
        d1 = sqrt(v1.p.x**2 + v1.p.y**2)
        is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01
    except:
        is_kelvin_arc = False

    if is_kelvin_arc:
        edge.name = "kelvin_int"
        kelvin_inner_edges.append(edge)
    else:
        edge.name = "wire_bnd"

for edge in wire1_full.edges:
    edge.name = "wire_plus_bnd"

for edge in wire2_full.edges:
    edge.name = "wire_minus_bnd"

# Find Kelvin boundary arc in outer_circle (centered at (0, y_offset))
kelvin_outer_edges = []
for edge in outer_circle.edges:
    cx = edge.center.x
    cy = edge.center.y - y_offset
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + (v0.p.y - y_offset)**2)
        d1 = sqrt(v1.p.x**2 + (v1.p.y - y_offset)**2)
        is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01
    except:
        is_kelvin_arc = False

    if is_kelvin_arc:
        edge.name = "kelvin_ext"
        kelvin_outer_edges.append(edge)

print(f"  Found {len(kelvin_inner_edges)} interior Kelvin edges")
print(f"  Found {len(kelvin_outer_edges)} exterior Kelvin edges")

# Identify periodic boundaries
# For full circle, there should be just one edge each
# Match them directly
matched_pairs = 0
for int_edge in kelvin_inner_edges:
    for ext_edge in kelvin_outer_edges:
        # Full circles - just match them
        int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
        matched_pairs += 1
        print(f"    Matched: int_center=({int_edge.center.x:.2f}, {int_edge.center.y:.2f})")
        print(f"             ext_center=({ext_edge.center.x:.2f}, {ext_edge.center.y:.2f})")
        break

print(f"  Matched {matched_pairs} periodic edge pairs")

# Combine and mesh
shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, gnd_point])
geo = OCCGeometry(shape, dim=2)

print(f"  Interior: full circle centered at (0, 0)")
print(f"  Exterior: full circle centered at (0, {y_offset})")

# ============================================================
# Mesh
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.4)).Curve(2)

print(f"  Number of elements: {mesh.ne}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# Check periodic boundary node correspondence
print("\nChecking periodic boundary correspondence...")
kelvin_int_nodes = []
kelvin_ext_nodes = []
for v in mesh.vertices:
    pt = v.point
    x_v, y_v = pt[0], pt[1]
    r_from_origin = sqrt(x_v**2 + y_v**2)
    r_from_offset = sqrt(x_v**2 + (y_v - y_offset)**2)

    if abs(r_from_origin - a) < 0.01:
        kelvin_int_nodes.append((x_v, y_v))
    if abs(r_from_offset - a) < 0.01:
        kelvin_ext_nodes.append((x_v, y_v))

print(f"  Interior Kelvin nodes: {len(kelvin_int_nodes)}")
print(f"  Exterior Kelvin nodes: {len(kelvin_ext_nodes)}")

# Check a few corresponding points
if len(kelvin_int_nodes) > 0 and len(kelvin_ext_nodes) > 0:
    print("\n  Sample node correspondence:")
    int_sorted = sorted(kelvin_int_nodes, key=lambda p: arctan2(p[1], p[0]))
    ext_sorted = sorted(kelvin_ext_nodes, key=lambda p: arctan2(p[1] - y_offset, p[0]))
    for i in range(min(5, len(int_sorted))):
        xi, yi = int_sorted[i]
        xe, ye = ext_sorted[i]
        theta_i = arctan2(yi, xi)
        theta_e = arctan2(ye - y_offset, xe)
        print(f"    int: ({xi:.4f}, {yi:.4f}) theta={theta_i:.4f}")
        print(f"    ext: ({xe:.4f}, {ye:.4f}) theta={theta_e:.4f}")

# ============================================================
# Finite Element Space
# ============================================================
print("\nSetting up A-formulation...")

# Dirichlet BC at GND (infinity)
fes_before = H1(mesh, order=order, dirichlet_bbnd="GND")
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
# Material Properties with Kelvin
# ============================================================
# (unused - kept for reference)
# x_local = x - x_offset
# rho_prime = sqrt(x_local**2 + y**2)
# rho_prime_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)

# 2D Kelvin transformation with Y-offset (like axisym Z-offset):
#   x' = x (unchanged!)
#   rho' = sqrt(x^2 + (y - y_offset)^2)
#
# Like in axisymmetric case, use (rho'/a)^2 factor

# Build reluctivity by region
materials = mesh.GetMaterials()
print(f"\n  Materials: {materials}")

# For exterior domain: rho' = sqrt(x^2 + (y - y_offset)^2)
y_local_cf = y - y_offset
rho_prime_cf = sqrt(x**2 + y_local_cf**2)
rho_prime_safe_cf = IfPos(rho_prime_cf - 1e-10, rho_prime_cf, 1e-10)

# For 2D Laplace equation: invariant under Kelvin (no factor needed)
nu_cf = nu0

print(f"  Using uniform nu (2D Laplace invariant under Kelvin)")

# ============================================================
# Weak Form
# ============================================================
# 2D A-formulation: -div(nu * grad(A_z)) = J_z
# Bilinear form: integral( nu * grad(u) . grad(v) dx )
a_form = BilinearForm(fes)
a_form += nu_cf * grad(u) * grad(v) * dx

# Linear form: f(v) = integral( J_z * v dx )
# wire_plus: +J0, wire_minus: -J0
f = LinearForm(fes)
f += J0 * v * dx("wire_plus")
f += (-J0) * v * dx("wire_minus")

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
# Analytical solution
# ============================================================
def Az_analytical(x_val, y_val):
    """
    Analytical A_z for two parallel wires with opposite currents.
    A_z = (mu0 * I / 2pi) * ln(r2/r1)
    r1 = distance from wire at (+d/2, 0)
    r2 = distance from wire at (-d/2, 0)
    """
    r1 = sqrt((x_val - wire_distance/2)**2 + y_val**2)
    r2 = sqrt((x_val + wire_distance/2)**2 + y_val**2)

    if r1 < wire_radius or r2 < wire_radius:
        return nan

    if r1 < 1e-12 or r2 < 1e-12:
        return nan

    return (mu0 * I_total / (2 * pi)) * log(r2 / r1)


def B_analytical(x_val, y_val):
    """
    Analytical B field for two parallel wires using Ampere's law.

    Wire 1 at (+d/2, 0) with current +I (flowing in +z):
      H_theta1 = I / (2*pi*r1), direction: counter-clockwise
      phi_hat1 = (-y, x - d/2) / r1

    Wire 2 at (-d/2, 0) with current -I (flowing in -z):
      H_theta2 = I / (2*pi*r2), direction: clockwise (due to -I)
      phi_hat2 = (y, -(x + d/2)) / r2  (opposite direction)
    """
    # Distances from each wire
    r1 = sqrt((x_val - wire_distance/2)**2 + y_val**2)
    r2 = sqrt((x_val + wire_distance/2)**2 + y_val**2)

    if r1 < wire_radius or r2 < wire_radius:
        return nan, nan

    if r1 < 1e-12 or r2 < 1e-12:
        return nan, nan

    # H from wire 1 (at +d/2, 0, current +I in +z direction)
    # H_theta = I/(2*pi*r), phi_hat = (-y, x-d/2)/r (counter-clockwise)
    Hx1 =  I_total / (2 * pi) * (y_val                  ) / (r1**2)
    Hy1 = -I_total / (2 * pi) * (x_val - wire_distance/2) / (r1**2)

    # H from wire 2 (at -d/2, 0, current -I in -z direction)
    # For -I, direction is clockwise: phi_hat = (y, -(x+d/2))/r
    Hx2 = -I_total / (2 * pi) * (y_val                  ) / (r2**2)
    Hy2 =  I_total / (2 * pi) * (x_val + wire_distance/2) / (r2**2)

    # B = mu0 * H
    Bx = mu0 * (Hx1 + Hx2)
    By = mu0 * (Hy1 + Hy2)

    return Bx, By

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# ============================================================
# Validation
# ============================================================
print("\n" + "=" * 60)
print("VALIDATION: Field comparison")
print("=" * 60)

# First, check Az values (not just B)
print(f"\n{'(x,y)':<15} {'Az (NGSolve)':<15} {'Az (Analytical)':<15} {'Diff':<15}")
print("-" * 60)
test_points = [(0, 0.1), (0, 0.5), (0.3, 0.3), (-0.3, 0.3)]
for x_val, y_val in test_points:
    try:
        mip = mesh(x_val, y_val)
        Az_num = gfu(mip)
        Az_ana = Az_analytical(x_val, y_val)
        diff = Az_num - Az_ana
        print(f"({x_val:.1f},{y_val:.1f}){'':<7} {Az_num:<15.6e} {Az_ana:<15.6e} {diff:<15.6e}")
    except Exception as e:
        print(f"({x_val:.1f},{y_val:.1f}){'':<7} Error: {e}")

# Check Az on Kelvin boundary
print(f"\nAz on Kelvin boundary (r={a}):")
print(f"{'theta':<10} {'(x,y)':<20} {'Az (NGSolve)':<15} {'Az (Analytical)':<15}")
print("-" * 60)
for theta_deg in [0, 45, 90, 135, 180]:
    theta_rad = theta_deg * pi / 180
    x_bnd = a * cos(theta_rad) * 0.99  # slightly inside
    y_bnd = a * sin(theta_rad) * 0.99
    try:
        mip = mesh(x_bnd, y_bnd)
        Az_num = gfu(mip)
        Az_ana = Az_analytical(x_bnd, y_bnd)
        print(f"{theta_deg:<10} ({x_bnd:.3f},{y_bnd:.3f}){'':<5} {Az_num:<15.6e} {Az_ana:<15.6e}")
    except Exception as e:
        print(f"{theta_deg:<10} Error: {e}")

# Sample points along y-axis (x=0)
y_samples = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
print(f"\n{'y (m)':<10} {'Bx (NGSolve)':<15} {'Bx (Analytical)':<15} {'Error (%)':<10}")
print("-" * 50)

for y_val in y_samples:
    try:
        mip = mesh(0, y_val)
        grad_Az = grad(gfu)
        Bx_num = grad_Az[1](mip)  # dAz/dy
        Bx_ana, _ = B_analytical(0, y_val)
        err = abs(Bx_num - Bx_ana) / abs(Bx_ana) * 100 if abs(Bx_ana) > 1e-15 else 0
        print(f"{y_val:<10.2f} {Bx_num:<15.6e} {Bx_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{y_val:<10.2f} Error: {e}")

# Field between wires (along x-axis)
print(f"\n{'x (m)':<10} {'By (NGSolve)':<15} {'By (Analytical)':<15} {'Error (%)':<10}")
print("-" * 50)

x_samples = [-0.15, -0.1, 0.0, 0.1, 0.15]
for x_val in x_samples:
    try:
        mip = mesh(x_val, 0.05)
        grad_Az = grad(gfu)
        By_num = -grad_Az[0](mip)  # -dAz/dx
        _, By_ana = B_analytical(x_val, 0.05)
        err = abs(By_num - By_ana) / abs(By_ana) * 100 if abs(By_ana) > 1e-15 else 0
        print(f"{x_val:<10.2f} {By_num:<15.6e} {By_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{x_val:<10.2f} Error: {e}")

# ============================================================
# Compute profiles
# ============================================================
print("\nComputing profiles...")

# Y-axis profile (x=0) - By component
y_profile = linspace(-1.0, 1.0, 101)
By_y_numerical = zeros(len(y_profile))
By_y_analytical = zeros(len(y_profile))

for i, y_val in enumerate(y_profile):
    _, By_y_analytical[i] = B_analytical(0, y_val)
    try:
        mip = mesh(0, y_val)
        By_y_numerical[i] = -grad(gfu)[0](mip)
    except:
        By_y_numerical[i] = nan

# X-axis profile (y=0.1)
x_profile = linspace(-1.0, 1.0, 101)
By_numerical = zeros(len(x_profile))
By_analytical_profile = zeros(len(x_profile))

for i, x_val in enumerate(x_profile):
    _, By_analytical_profile[i] = B_analytical(x_val, 0.1)
    try:
        mip = mesh(x_val, 0.1)
        By_numerical[i] = -grad(gfu)[0](mip)
    except:
        By_numerical[i] = nan

# ============================================================
# 2D field for visualization
# ============================================================
print("\nComputing 2D field...")

x_grid = linspace(-a + 0.05, a - 0.05, 41)
y_grid = linspace(-a + 0.05, a - 0.05, 41)
[xx, yy] = meshgrid(x_grid, y_grid)

Bx_field = zeros(xx.shape)
By_field = zeros(xx.shape)
Az_field = zeros(xx.shape)

for iy in range(len(y_grid)):
    for ix in range(len(x_grid)):
        x_val, y_val = x_grid[ix], y_grid[iy]
        # Skip points inside wires
        r1 = sqrt((x_val - wire_distance/2)**2 + y_val**2)
        r2 = sqrt((x_val + wire_distance/2)**2 + y_val**2)
        if r1 < wire_radius * 1.2 or r2 < wire_radius * 1.2:
            Bx_field[iy, ix] = nan
            By_field[iy, ix] = nan
            Az_field[iy, ix] = nan
            continue

        try:
            mip = mesh(x_val, y_val)
            grad_Az = grad(gfu)
            Bx_field[iy, ix] = grad_Az[1](mip)   # dAz/dy
            By_field[iy, ix] = -grad_Az[0](mip)  # -dAz/dx
            Az_field[iy, ix] = gfu(mip)
        except:
            Bx_field[iy, ix] = nan
            By_field[iy, ix] = nan
            Az_field[iy, ix] = nan

# ============================================================
# Error statistics
# ============================================================
valid_idx = ~isnan(By_y_numerical)
if sum(valid_idx) > 0:
    error = By_y_numerical[valid_idx] - By_y_analytical[valid_idx]
    rms_error = sqrt(mean(error**2))
    max_error = abs(error).max()
    print(f"\n  Y-axis profile (By) error statistics:")
    print(f"    RMS error: {rms_error:.6e} T")
    print(f"    Max error: {max_error:.6e} T")
    if abs(By_y_analytical[valid_idx]).mean() > 1e-15:
        print(f"    Relative RMS error: {rms_error/abs(By_y_analytical[valid_idx]).mean()*100:.3f}%")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})

fig = plt.figure(figsize=(14, 12), dpi=150)


def draw_wires(ax):
    """Draw wire cross-sections."""
    theta = linspace(0, 2*pi, 50)
    # Wire 1 (+)
    ax.fill(wire_distance/2 + wire_radius*cos(theta),
            wire_radius*sin(theta), 'red', alpha=0.5)
    ax.plot(wire_distance/2 + wire_radius*cos(theta),
            wire_radius*sin(theta), 'r-', linewidth=1)
    # Wire 2 (-)
    ax.fill(-wire_distance/2 + wire_radius*cos(theta),
            wire_radius*sin(theta), 'blue', alpha=0.5)
    ax.plot(-wire_distance/2 + wire_radius*cos(theta),
            wire_radius*sin(theta), 'b-', linewidth=1)


def draw_kelvin_boundary(ax):
    """Draw Kelvin boundary circle."""
    theta = linspace(0, 2*pi, 100)
    ax.plot(a * cos(theta), a * sin(theta), 'g--', linewidth=2, label='Kelvin boundary')


# Plot 1: B field streamlines
ax1 = plt.subplot(2, 2, 1)
Bx_plot = Bx_field.copy()
By_plot = By_field.copy()
Bx_plot[isnan(Bx_plot)] = 0
By_plot[isnan(By_plot)] = 0
ax1.streamplot(xx, yy, Bx_plot, By_plot, color='blue', linewidth=0.8, density=1.5)
draw_wires(ax1)
draw_kelvin_boundary(ax1)
ax1.set_xlabel('$x$ (m)')
ax1.set_ylabel('$y$ (m)')
ax1.set_title('$\\mathbf{B}$ Field (Kelvin Transformation)')
ax1.set_aspect('equal')
ax1.set_xlim(-a, a)
ax1.set_ylim(-a, a)
ax1.legend(loc='upper right', fontsize=8)

# Plot 2: |B| contour
ax2 = plt.subplot(2, 2, 2)
B_mag = sqrt(Bx_field**2 + By_field**2)
B_mag_plot = B_mag.copy()
B_mag_plot[isnan(B_mag_plot)] = 0
from numpy import nanmax
levels = linspace(0, nanmax(B_mag)*0.8, 20)
contour2 = ax2.contourf(xx, yy, B_mag_plot*1e6, levels=levels*1e6, cmap='jet', extend='max')
cbar2 = plt.colorbar(contour2, ax=ax2, shrink=0.8)
cbar2.set_label('$|\\mathbf{B}|$ ($\\mu$T)')
draw_wires(ax2)
draw_kelvin_boundary(ax2)
ax2.set_xlabel('$x$ (m)')
ax2.set_ylabel('$y$ (m)')
ax2.set_title('$|\\mathbf{B}|$ Contour')
ax2.set_aspect('equal')
ax2.set_xlim(-a, a)
ax2.set_ylim(-a, a)

# Plot 3: By profile along y-axis (x=0)
ax3 = plt.subplot(2, 2, 3)
valid = ~isnan(By_y_numerical)
ax3.plot(y_profile[valid], By_y_numerical[valid] * 1e6, 'k-', linewidth=2, label='NGSolve')
ax3.plot(y_profile, By_y_analytical * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax3.set_xlabel('$y$ (m)')
ax3.set_ylabel('$B_y$ ($\\mu$T)')
ax3.set_title('$B_y$ Profile on $y$-axis ($x = 0$)')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# Plot 4: By profile along x-axis
ax4 = plt.subplot(2, 2, 4)
valid = ~isnan(By_numerical)
ax4.plot(x_profile[valid], By_numerical[valid] * 1e6, 'k-', linewidth=2, label='NGSolve')
ax4.plot(x_profile, By_analytical_profile * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(wire_distance/2, color='red', linestyle=':', alpha=0.5, label='Wire +')
ax4.axvline(-wire_distance/2, color='blue', linestyle=':', alpha=0.5, label='Wire -')
ax4.set_xlabel('$x$ (m)')
ax4.set_ylabel('$B_y$ ($\\mu$T)')
ax4.set_title('$B_y$ Profile ($y = 0.1$ m)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

plt.tight_layout()

png_file = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

# ============================================================
# VTK Output
# ============================================================
print("\nExporting to VTK...")
vtk_file = os.path.splitext(__file__)[0]
vtk = VTKOutput(mesh, coefs=[gfu], names=["Az"],
                filename=vtk_file, subdivision=2)
vtk.Do()
print(f"  VTK file saved to: {vtk_file}.vtu")

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
