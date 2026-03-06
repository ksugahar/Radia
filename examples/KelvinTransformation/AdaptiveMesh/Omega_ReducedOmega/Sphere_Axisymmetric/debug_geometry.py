"""
Debug script to check Periodic in Axisymmetric case
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Testing Axisymmetric Periodic identification")
print("=" * 60)

kelvin_radius = 1.0
offset_x = 3.0

# INTERIOR DOMAIN (half-circle)
outer_half_int = Circle((0, 0), kelvin_radius).Face()
cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_int = outer_half_int - cutter_int

# EXTERIOR DOMAIN
outer_half_ext = Circle((offset_x, 0), kelvin_radius).Face()
cutter_ext = MoveTo(offset_x - kelvin_radius - 1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_ext = outer_half_ext - cutter_ext

# Name edges
for edge in outer_half_int.edges:
    x_center = edge.center.x
    if x_center < 1e-6:
        edge.name = "axis_int"
    else:
        edge.name = "kelvin_int"
        print(f"Inner kelvin edge: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

for edge in outer_half_ext.edges:
    x_center = edge.center.x - offset_x
    if x_center < 1e-6:
        edge.name = "axis_ext"
    else:
        edge.name = "kelvin_ext"
        print(f"Outer kelvin edge: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

outer_half_int.faces.name = "air_inner"
outer_half_ext.faces.name = "air_outer"

vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

shape = Glue([outer_half_int, outer_half_ext, vertex])

# Identify periodic
kelvin_int_edges = []
kelvin_ext_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_int_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_ext_edges.append(edge)

print(f"\nFound {len(kelvin_int_edges)} interior kelvin edges")
print(f"Found {len(kelvin_ext_edges)} exterior kelvin edges")

if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    for int_edge in kelvin_int_edges:
        int_z = int_edge.center.y
        for ext_edge in kelvin_ext_edges:
            ext_z = ext_edge.center.y
            if (int_z > 0 and ext_z > 0) or (int_z < 0 and ext_z < 0):
                print(f"  Matching: ({int_edge.center.x:.3f}, {int_edge.center.y:.3f}) <-> ({ext_edge.center.x:.3f}, {ext_edge.center.y:.3f})")
                int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                break

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.3)).Curve(2)

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")

# Check Periodic DOF reduction
fes_before = H1(mesh, order=2, dirichlet_bbnd="GND")
fes_periodic = Periodic(fes_before)
print(f"\nFES without Periodic: {fes_before.ndof} DOFs")
print(f"FES with Periodic: {fes_periodic.ndof} DOFs")
print(f"DOF reduction: {fes_before.ndof - fes_periodic.ndof}")

if fes_before.ndof == fes_periodic.ndof:
    print("WARNING: Periodic BC not reducing DOFs!")
else:
    print("OK: Periodic BC is working")
