"""
H-formulation for 2D magnetostatics with perturbation potential and Kelvin transformation
Using HALF-CIRCLE geometry to test if Periodic BC works with segmented edges

This is a test to verify whether half-circle geometry causes Periodic BC to fail.
Compare results with 2D_dipole_with_Kelvin.py (full circle version).

Based on: 2D_dipole_with_Kelvin.py
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 2D with Kelvin Transformation (HALF-CIRCLE)")
print("="*60)

# ============================================================
# Geometry Definition (OCC 2D) - Two separate HALF-circular domains
# ============================================================
print("\nCreating HALF-CIRCLE geometry with periodic boundary conditions...")

# Parameters
circle_radius = 0.5   # Magnetic circle radius [m]
kelvin_radius = 1.0   # Kelvin transformation radius [m] (outer boundary)
maxh_fine = 0.03      # Fine mesh size [m]
plot_range = 1.1      # Plot range [m]

# Offset for exterior domain (placed separately)
offset_x = 3.0  # Offset to place exterior domain away from interior

print(f"Using Kelvin transformation with periodic BC (HALF-CIRCLE):")
print(f"  - Inner domain: 0 < r < R = {kelvin_radius} m at (0, 0), x >= 0")
print(f"  - Outer domain: 0 < r' < R = {kelvin_radius} m at ({offset_x}, 0), x' >= 0")
print(f"  - Transformation radius R = {kelvin_radius} m")

# ===== INTERIOR DOMAIN (HALF-circle, x >= 0) =====
# Create full circles first, then cut
outer_half_int = Circle((0, 0), kelvin_radius).Face()
cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_int = outer_half_int - cutter_int  # Keep x >= 0 part

inner_half_int = Circle((0, 0), circle_radius).Face()
inner_half_int = inner_half_int - cutter_int  # Keep x >= 0 part

# Air region (half-annulus)
air_inner = outer_half_int - inner_half_int

# Name boundaries
for edge in air_inner.edges:
    x_center = edge.center.x
    dist = sqrt(edge.center.x**2 + edge.center.y**2)
    if x_center < 1e-6:  # On y-axis (x = 0)
        edge.name = "axis_int"
    elif abs(dist - kelvin_radius) < kelvin_radius * 0.2:
        edge.name = "kelvin_int"
    else:
        edge.name = "circle"
air_inner.faces.name = "air_inner"

# Magnetic material (half-circle)
for edge in inner_half_int.edges:
    if edge.center.x < 1e-6:
        edge.name = "axis_int"
    else:
        edge.name = "circle"
inner_half_int.faces.name = "magnetic"

# ===== EXTERIOR DOMAIN (HALF-circle, x' >= 0) =====
outer_half_ext = Circle((offset_x, 0), kelvin_radius).Face()
cutter_ext = MoveTo(offset_x - kelvin_radius - 1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_ext = outer_half_ext - cutter_ext  # Keep x' >= 0 part

for edge in outer_half_ext.edges:
    x_center = edge.center.x - offset_x
    if x_center < 1e-6:  # On axis (x' = 0)
        edge.name = "axis_ext"
    else:
        edge.name = "kelvin_ext"
outer_half_ext.faces.name = "air_outer"

# ===== GND VERTEX (center of exterior domain) =====
vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

# Glue all domains
shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

# ===== IDENTIFY PERIODIC BOUNDARIES =====
print("\nIdentifying periodic boundaries...")

# Find kelvin boundary edges
kelvin_int_edges = []
kelvin_ext_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_int_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_ext_edges.append(edge)

print(f"  Found {len(kelvin_int_edges)} interior kelvin edges")
print(f"  Found {len(kelvin_ext_edges)} exterior kelvin edges")

# Print edge details
for i, edge in enumerate(kelvin_int_edges):
    print(f"    kelvin_int[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")
for i, edge in enumerate(kelvin_ext_edges):
    print(f"    kelvin_ext[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# Match edges by y-coordinate of center (y>0 with y>0, y<0 with y<0)
if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    matched_pairs = 0
    for int_edge in kelvin_int_edges:
        int_y = int_edge.center.y
        # Find matching exterior edge (same sign of y)
        for ext_edge in kelvin_ext_edges:
            ext_y = ext_edge.center.y
            # Match by sign of y-coordinate
            if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
                int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                print(f"  Identified: int(y={int_y:.3f}) <-> ext(y={ext_y:.3f})")
                matched_pairs += 1
                break
    print(f"  Total matched pairs: {matched_pairs}")

# Create geometry
geo = OCCGeometry(shape, dim=2)

print(f"\nGeometry created (HALF-CIRCLE):")
print(f"  Magnetic circle radius: {circle_radius} m at (0, 0), x >= 0")
print(f"  Inner air domain: {circle_radius} m < r < {kelvin_radius} m, x >= 0")
print(f"  Outer air domain: 0 < r' < {kelvin_radius} m at ({offset_x}, 0), x' >= 0")
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
# Axis boundaries have natural Neumann BC (no constraint needed)
fes_before = H1(mesh, order=3, dirichlet_bbnd="GND")
fes = Periodic(fes_before)

# Check if Periodic BC is working (FreeDofs should decrease)
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")

if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
    print("  This confirms that half-circle geometry breaks Periodic BC.")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability
mu_d = {"air_inner": 1*mu0, "air_outer": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field (2D Kelvin: constant, sign-reversed in exterior)
is_exterior = IfPos(x - offset_x/2, 1.0, 0.0)

Hx_inner = 0.0
Hy_inner = 1.0
Hs_x_outer = 0.0
Hs_y_outer = -1.0

Hs_x = (1.0 - is_exterior) * Hx_inner + is_exterior * Hs_x_outer
Hs_y = (1.0 - is_exterior) * Hy_inner + is_exterior * Hs_y_outer

Hs = CoefficientFunction((Hs_x, Hs_y))

print(f"  Background field: H_s with Kelvin transformation")
print(f"  Relative permeability: μ_r = {mu_r}")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (volume integral only with Kelvin + periodic BC)
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx

a.Assemble()
f.Assemble()

print("  System assembled")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfu = GridFunction(fes)

# Use direct solver for stability
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

print("  Solution converged")

# ============================================================
# Post-processing
# ============================================================
print("\nPost-processing...")

# Compute perturbation field: H_pert = -grad(phi)
H = -grad(gfu)

# Analytical solution (2D)
Hy_analytical = -1.0 + 2.0/(mu_r + 1)

print("\n" + "="*60)
print("RESULTS COMPARISON (HALF-CIRCLE vs ANALYTICAL)")
print("="*60)
print(f"  Analytical interior: Hy = {Hy_analytical:.6f} A/m")
print()

# Evaluate at multiple points
for x_val in [0.1, 0.2, 0.3, 0.4]:
    try:
        Hy_val = H[1](mesh(x_val, 0))
        err = abs(Hy_val - Hy_analytical)/abs(Hy_analytical)*100
        print(f"  x={x_val}: Hy = {Hy_val:.6f} A/m, error = {err:.4f}%")
    except Exception as e:
        print(f"  x={x_val}: Error - {e}")

print()
print("="*60)
print("EXTERIOR:")
print("="*60)

for x_val in [0.6, 0.7, 0.8, 0.9]:
    # Analytical exterior (on x-axis, y=0)
    Hy_ext_analytical = -(mu_r - 1)/(mu_r + 1) * (circle_radius/x_val)**2
    try:
        Hy_val = H[1](mesh(x_val, 0))
        err = abs(Hy_val - Hy_ext_analytical)/abs(Hy_ext_analytical)*100
        print(f"  x={x_val}: Hy = {Hy_val:.6f} A/m (analytical: {Hy_ext_analytical:.6f}), error = {err:.4f}%")
    except Exception as e:
        print(f"  x={x_val}: Error - {e}")

# ============================================================
# Compare with full circle result
# ============================================================
print("\n" + "="*60)
print("CONCLUSION")
print("="*60)
if freedof_before == freedof_after:
    print("  Periodic BC is NOT working with half-circle geometry!")
    print("  This confirms that:")
    print("    - Half-circle creates multiple edge segments")
    print("    - Identify() only connects one segment pair")
    print("    - Full circle is REQUIRED for Periodic BC to work")
else:
    print("  Periodic BC IS working with half-circle geometry.")
    print("  (This would be unexpected based on previous findings)")

print("\n" + "="*60)
print("Computation completed")
print("="*60)
