"""
Debug script to compare geometry edge structures
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Comparing geometry edge structures")
print("=" * 60)

# ============================================================
# Method 1: ArcTo (current parallel wires approach)
# ============================================================
print("\n--- Method 1: ArcTo (Parallel Wires) ---")
wp1 = WorkPlane()
circle_arcto = wp1.MoveTo(0, -1).ArcTo(0, 1, (1, 0)).ArcTo(0, -1, (-1, 0)).Face()

print(f"Face has {len(circle_arcto.edges)} edges")
for i, edge in enumerate(circle_arcto.edges):
    print(f"  Edge {i}: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# ============================================================
# Method 2: Circle - Rectangle (Axisymmetric approach)
# ============================================================
print("\n--- Method 2: Circle - Rectangle (Axisymmetric) ---")
circle_full = Circle((0, 0), 1.0).Face()
cutter = MoveTo(-2, -2).Rectangle(2, 4).Face()
circle_half = circle_full - cutter  # Keep x >= 0 part

print(f"Face has {len(circle_half.edges)} edges")
for i, edge in enumerate(circle_half.edges):
    print(f"  Edge {i}: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# ============================================================
# Now let's check the mesh edge pairing after meshing
# ============================================================
print("\n" + "=" * 60)
print("Testing mesh generation and Periodic identification")
print("=" * 60)

# Create two domains with ArcTo method
a = 1.0
y_offset = 2.5

wp1 = WorkPlane()
inner_circle = wp1.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
inner_circle.name = "air_inner"

wp2 = WorkPlane(Axes((0, y_offset, 0), n=Z, h=X))
outer_circle = wp2.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
outer_circle.name = "air_outer"

vertex = Vertex(Pnt(0, 0, 0))
vertex.name = "fix_point"

# Name and identify edges
kelvin_inner_edges = []
for edge in inner_circle.edges:
    edge.name = "kelvin_int"
    kelvin_inner_edges.append(edge)
    print(f"Inner edge: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

kelvin_outer_edges = []
for edge in outer_circle.edges:
    edge.name = "kelvin_ext"
    kelvin_outer_edges.append(edge)
    print(f"Outer edge: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# GLUE first, then identify (like Axisymmetric)
shape = Glue([inner_circle, outer_circle, vertex])

# Find edges from glued shape
kelvin_inner_edges = []
kelvin_outer_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_inner_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_outer_edges.append(edge)

print(f"After Glue: {len(kelvin_inner_edges)} inner, {len(kelvin_outer_edges)} outer kelvin edges")

# Identify by x-coordinate sign
print("\nIdentifying periodic boundaries...")
for int_edge in kelvin_inner_edges:
    int_x = int_edge.center.x
    for ext_edge in kelvin_outer_edges:
        ext_x = ext_edge.center.x
        if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
            print(f"  Matching: inner ({int_x:.3f}, {int_edge.center.y:.3f}) <-> outer ({ext_x:.3f}, {ext_edge.center.y:.3f})")
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
            break

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.3)).Curve(2)

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")

# Check boundary elements
print("\nBoundary info:")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# Check if periodic is working by testing Periodic FES
fes_before = H1(mesh, order=2, dirichlet_bbnd="fix_point")
fes_periodic = Periodic(fes_before)
print(f"\n  FES without Periodic: {fes_before.ndof} DOFs")
print(f"  FES with Periodic: {fes_periodic.ndof} DOFs")
print(f"  DOF reduction: {fes_before.ndof - fes_periodic.ndof}")

if fes_before.ndof == fes_periodic.ndof:
    print("  WARNING: Periodic BC not reducing DOFs!")
else:
    print("  OK: Periodic BC is working")
