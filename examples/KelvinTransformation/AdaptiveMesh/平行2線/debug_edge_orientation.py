"""
Debug script to check edge orientation for periodic matching
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Testing edge orientation for periodic matching")
print("=" * 60)

a = 1.0
y_offset = 2.5

# Create circles with ArcTo (2 arcs each)
wp1 = WorkPlane()
inner_circle = wp1.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
inner_circle.name = "air_inner"

wp2 = WorkPlane(Axes((0, y_offset, 0), n=Z, h=X))
outer_circle = wp2.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
outer_circle.name = "air_outer"

vertex = Vertex(Pnt(0, 0, 0))
vertex.name = "fix_point"

# Name edges before glue
for edge in inner_circle.edges:
    edge.name = "kelvin_int"

for edge in outer_circle.edges:
    edge.name = "kelvin_ext"

# Glue
shape = Glue([inner_circle, outer_circle, vertex])

# Find edges and check vertices
print("\nEdge vertex information:")
print("-" * 60)

kelvin_inner_edges = []
kelvin_outer_edges = []

for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_inner_edges.append(edge)
        v0, v1 = edge.vertices
        print(f"Inner edge: ({v0.p.x:.3f}, {v0.p.y:.3f}) -> ({v1.p.x:.3f}, {v1.p.y:.3f}), center=({edge.center.x:.3f}, {edge.center.y:.3f})")
    elif edge.name == "kelvin_ext":
        kelvin_outer_edges.append(edge)
        v0, v1 = edge.vertices
        print(f"Outer edge: ({v0.p.x:.3f}, {v0.p.y:.3f}) -> ({v1.p.x:.3f}, {v1.p.y:.3f}), center=({edge.center.x:.3f}, {edge.center.y:.3f})")

print(f"\nTotal: {len(kelvin_inner_edges)} inner, {len(kelvin_outer_edges)} outer")

# Check if edges go in the same direction
print("\n" + "=" * 60)
print("Checking edge directions for matching:")
print("=" * 60)

for i, int_edge in enumerate(kelvin_inner_edges):
    int_v0, int_v1 = int_edge.vertices
    int_x = int_edge.center.x

    for ext_edge in kelvin_outer_edges:
        ext_v0, ext_v1 = ext_edge.vertices
        ext_x = ext_edge.center.x

        if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
            print(f"\nPair {i+1} (x {'>' if int_x > 0 else '<'} 0):")
            print(f"  Inner: ({int_v0.p.x:.3f}, {int_v0.p.y:.3f}) -> ({int_v1.p.x:.3f}, {int_v1.p.y:.3f})")
            print(f"  Outer: ({ext_v0.p.x:.3f}, {ext_v0.p.y:.3f}) -> ({ext_v1.p.x:.3f}, {ext_v1.p.y:.3f})")

            # Check if directions match (inner y direction vs outer y direction)
            int_dy = int_v1.p.y - int_v0.p.y
            ext_dy = ext_v1.p.y - ext_v0.p.y
            print(f"  Inner dy: {int_dy:.3f}, Outer dy: {ext_dy:.3f}")

            if int_dy * ext_dy > 0:
                print(f"  -> Directions MATCH")
            else:
                print(f"  -> Directions OPPOSITE - may need orientation reversal!")
            break

# Now test with identification
print("\n" + "=" * 60)
print("Testing Periodic after identification:")
print("=" * 60)

# Identify
for int_edge in kelvin_inner_edges:
    int_x = int_edge.center.x
    for ext_edge in kelvin_outer_edges:
        ext_x = ext_edge.center.x
        if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
            break

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.3)).Curve(2)

fes_before = H1(mesh, order=2, dirichlet_bbnd="fix_point")
fes_periodic = Periodic(fes_before)

print(f"FES without Periodic: {fes_before.ndof} DOFs")
print(f"FES with Periodic: {fes_periodic.ndof} DOFs")
print(f"DOF reduction: {fes_before.ndof - fes_periodic.ndof}")

# Let's also check what identifications are in the mesh
print(f"\nMesh identifications: {mesh.GetIdentifications()}")
