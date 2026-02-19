"""
Debug script to check FreeDofs reduction with Periodic
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Testing FreeDofs reduction with Periodic")
print("=" * 60)

a = 1.0
y_offset = 2.5
wire_distance = 1.4
wire_radius = 0.02

# Create geometry
wp1 = WorkPlane()
inner_circle = wp1.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
inner_circle.name = "air_inner"

wire1_full = MoveTo(wire_distance/2, 0).Circle(wire_radius).Face()
wire1_full.name = "wire_minus"
wire2_full = MoveTo(-wire_distance/2, 0).Circle(wire_radius).Face()
wire2_full.name = "wire_plus"

wp2 = WorkPlane(Axes((0, y_offset, 0), n=Z, h=X))
outer_circle = wp2.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
outer_circle.name = "air_outer"

fix_point = Vertex(Pnt(wire_distance/2 + wire_radius, 0, 0))
fix_point.name = "fix_point"

inner_air = inner_circle - wire1_full - wire2_full
inner_air.name = "air_inner"

# Name edges
for edge in inner_air.edges:
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + v0.p.y**2)
        d1 = sqrt(v1.p.x**2 + v1.p.y**2)
        is_kelvin = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01
    except:
        is_kelvin = False
    edge.name = "kelvin_int" if is_kelvin else "wire_bnd"

for edge in outer_circle.edges:
    edge.name = "kelvin_ext"

# Glue first
shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, fix_point])

# Find and identify edges
kelvin_int = [e for e in shape.edges if e.name == "kelvin_int"]
kelvin_ext = [e for e in shape.edges if e.name == "kelvin_ext"]

print(f"Found {len(kelvin_int)} interior, {len(kelvin_ext)} exterior kelvin edges")

for int_edge in kelvin_int:
    int_x = int_edge.center.x
    for ext_edge in kelvin_ext:
        ext_x = ext_edge.center.x
        if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
            break

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.5)).Curve(2)

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")

# Test for order=3
order = 3
print(f"\n=== Testing Order = {order} with Refinement ===")

def check_freedofs(mesh, label):
    fes_before = H1(mesh, order=order, dirichlet_bbnd="fix_point")
    fes_after = Periodic(fes_before)
    freedofs_before = sum(1 for i in range(fes_before.ndof) if fes_before.FreeDofs()[i])
    freedofs_after = sum(1 for i in range(fes_after.ndof) if fes_after.FreeDofs()[i])
    reduction = freedofs_before - freedofs_after
    print(f"{label}: ndof={fes_after.ndof}, FreeDofs reduction={reduction}")
    return fes_after, reduction

# Initial mesh
fes, red0 = check_freedofs(mesh, "Initial mesh")

# Refine some elements
for el in mesh.Elements():
    mesh.SetRefinementFlag(el, False)

# Mark first 50 elements
for i, el in enumerate(mesh.Elements()):
    if i < 50:
        mesh.SetRefinementFlag(el, True)

mesh.Refine()
fes, red1 = check_freedofs(mesh, "After Refine 1")

# Refine again
for el in mesh.Elements():
    mesh.SetRefinementFlag(el, False)
for i, el in enumerate(mesh.Elements()):
    if i < 100:
        mesh.SetRefinementFlag(el, True)

mesh.Refine()
fes, red2 = check_freedofs(mesh, "After Refine 2")

# Refine again
for el in mesh.Elements():
    mesh.SetRefinementFlag(el, False)
for i, el in enumerate(mesh.Elements()):
    if i < 200:
        mesh.SetRefinementFlag(el, True)

mesh.Refine()
fes, red3 = check_freedofs(mesh, "After Refine 3")

print(f"\nFreeDofs reductions: {red0} -> {red1} -> {red2} -> {red3}")
if red1 > red0 and red2 > red1:
    print("OK: FreeDofs reduction increases with refinement (more boundary DOFs)")
elif red1 == 0 or red2 == 0:
    print("WARNING: Periodic identification lost after refinement!")
