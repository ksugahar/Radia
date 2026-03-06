"""
Debug script to test periodic boundary behavior after mesh.Refine()
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

# Parameters
wire_distance = 1.4
wire_radius = 0.02
a = 1.0
y_offset = 2.5
maxh_initial = 0.5
order = 3

mu0 = 4 * pi * 1e-7
nu0 = 1 / mu0
I_total = 1.0
J0 = I_total / (pi * wire_radius**2)

def create_geometry():
    """Create the parallel wires geometry with Kelvin transformation.
    Using Circle - Rectangle method (same as Axisymmetric) to ensure proper periodic BCs.
    """
    # Interior domain: circle cut in half using rectangle (same as Axisymmetric)
    inner_circle_full = Circle((0, 0), a).Face()
    # Don't cut - use full circle for 2D problem

    # Actually, let's try: cut at y=0 to create 2 arcs
    # Use rectangle to cut left half (x < 0)
    # NO - in 2D we need full circle. Let's try another approach.

    # The key insight: in Axisymmetric, Circle creates 2 edges naturally (upper and lower half)
    # Let's check if our ArcTo method creates proper periodic matching

    # Use standard ArcTo method
    wp1 = WorkPlane()
    inner_circle = wp1.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
    inner_circle.name = "air_inner"

    # Wires
    wire1_full = MoveTo(wire_distance/2, 0).Circle(wire_radius).Face()
    wire1_full.name = "wire_minus"
    wire2_full = MoveTo(-wire_distance/2, 0).Circle(wire_radius).Face()
    wire2_full.name = "wire_plus"

    # Exterior domain: 2 arcs split at x=0
    wp2 = WorkPlane(Axes((0, y_offset, 0), n=Z, h=X))
    outer_circle = wp2.MoveTo(0, -a).ArcTo(0, a, (a, 0)).ArcTo(0, -a, (-a, 0)).Face()
    outer_circle.name = "air_outer"

    # Fix point
    fix_point = Vertex(Pnt(wire_distance/2 + wire_radius, 0, 0))
    fix_point.name = "fix_point"

    # Subtract wires from inner air
    inner_air = inner_circle - wire1_full - wire2_full
    inner_air.name = "air_inner"

    # Name boundaries
    kelvin_inner_edges = []
    for edge in inner_air.edges:
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

    kelvin_outer_edges = []
    for edge in outer_circle.edges:
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

    # Identify periodic boundaries
    print(f"Found {len(kelvin_inner_edges)} interior kelvin edges")
    print(f"Found {len(kelvin_outer_edges)} exterior kelvin edges")

    matched_pairs = 0
    for int_edge in kelvin_inner_edges:
        int_x = int_edge.center.x
        for ext_edge in kelvin_outer_edges:
            ext_x = ext_edge.center.x
            if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
                int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
                matched_pairs += 1
                break

    print(f"Matched {matched_pairs} periodic edge pairs")

    shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, fix_point])
    return OCCGeometry(shape, dim=2)

def check_periodic_info(mesh, label=""):
    """Check periodic vertex mapping."""
    print(f"\n=== {label} ===")

    # Check periodic identification using NT_VERTEX and NT_EDGE
    try:
        periodic_verts = mesh.GetPeriodicNodePairs(NT_VERTEX)
        print(f"Periodic vertex pairs: {len(periodic_verts)}")
    except:
        periodic_verts = []
        print("Could not get periodic vertex pairs")

    try:
        periodic_edges = mesh.GetPeriodicNodePairs(NT_EDGE)
        print(f"Periodic edge pairs: {len(periodic_edges)}")
    except:
        periodic_edges = []
        print("Could not get periodic edge pairs")

    return len(periodic_verts), len(periodic_edges)


def solve_and_get_error(mesh, order):
    """Solve and return error estimator."""
    fes_before = H1(mesh, order=order, dirichlet_bbnd="fix_point")
    fes = Periodic(fes_before)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    a_form = BilinearForm(fes)
    a_form += nu0 * grad(u) * grad(v) * dx
    a_form.Assemble()

    f = LinearForm(fes)
    f += J0 * v * dx("wire_plus")
    f += (-J0) * v * dx("wire_minus")
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # Error estimator
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    W = Integrate(nu0 * (B * B), mesh)
    if W < 1e-20:
        W = 1.0

    H = nu0 * B
    Hfes = HCurl(mesh, order=order+1, type1=False)
    Hint = GridFunction(Hfes)
    Hint.Set(H)

    err = ((Hint - H) * (Hint - H)) / nu0 / W
    element_errors = Integrate(err, mesh, element_wise=True)
    total_error = sqrt(sum(element_errors))

    return fes.ndof, total_error, element_errors, gfu


print("=" * 60)
print("Creating geometry...")
geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.4)).Curve(2)

print(f"\nInitial mesh: {mesh.ne} elements, {mesh.nv} vertices")

# Check initial periodic info
v0, e0 = check_periodic_info(mesh, "Before Refinement")

# Solve
ndof, error, elem_errors, gfu = solve_and_get_error(mesh, order)
print(f"\nIteration 1: DOFs={ndof}, Error={error:.6e}")

# Refine
theta = 0.5
sorted_indices = sorted(range(len(elem_errors)), key=lambda i: elem_errors[i], reverse=True)
total_error = sum(elem_errors)
marked_error = 0
marked = []
for idx in sorted_indices:
    marked.append(idx)
    marked_error += elem_errors[idx]
    if marked_error >= theta * total_error:
        break

for el in mesh.Elements():
    mesh.SetRefinementFlag(el, False)
for el_nr in marked:
    mesh.SetRefinementFlag(ElementId(VOL, el_nr), True)

mesh.Refine()

# Check periodic info after refinement
v1, e1 = check_periodic_info(mesh, "After Refinement")

# Solve again
ndof2, error2, elem_errors2, gfu2 = solve_and_get_error(mesh, order)
print(f"\nIteration 2: DOFs={ndof2}, Error={error2:.6e}")

# Summary
print("\n" + "=" * 60)
print("Summary:")
print(f"  Before Refine: {v0} vertex pairs, {e0} edge pairs")
print(f"  After Refine:  {v1} vertex pairs, {e1} edge pairs")
print(f"  Error change: {error:.6e} -> {error2:.6e}")
if error2 > error:
    print("  WARNING: Error INCREASED after refinement!")
print("=" * 60)
