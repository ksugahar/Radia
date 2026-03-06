"""
Debug script to test periodic boundary behavior in Axisymmetric case
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

# Parameters
sphere_radius = 0.5
kelvin_radius = 1.0
mu_r = 100
mu0 = 4 * pi * 1e-7
H0 = 1.0
maxh_initial = 0.3
order = 4
offset_x = 3.0

def create_geometry():
    """Create geometry with periodic boundary conditions."""
    print("Creating half-circle geometry...")

    # INTERIOR DOMAIN (half-circle, r >= 0)
    outer_half_int = Circle((0, 0), kelvin_radius).Face()
    cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
    outer_half_int = outer_half_int - cutter_int

    inner_half_int = Circle((0, 0), sphere_radius).Face()
    inner_half_int = inner_half_int - cutter_int

    air_inner = outer_half_int - inner_half_int

    for edge in air_inner.edges:
        x_center = edge.center.x
        dist = sqrt(edge.center.x**2 + edge.center.y**2)
        if x_center < 1e-6:
            edge.name = "axis_int"
        elif dist > (kelvin_radius + sphere_radius) / 2:
            edge.name = "kelvin_int"
        else:
            edge.name = "sphere"
    air_inner.faces.name = "air_inner"

    for edge in inner_half_int.edges:
        if edge.center.x < 1e-6:
            edge.name = "axis_int"
        else:
            edge.name = "sphere"
    inner_half_int.faces.name = "magnetic"

    # EXTERIOR DOMAIN
    outer_half_ext = Circle((offset_x, 0), kelvin_radius).Face()
    cutter_ext = MoveTo(offset_x - kelvin_radius - 1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
    outer_half_ext = outer_half_ext - cutter_ext

    for edge in outer_half_ext.edges:
        x_center = edge.center.x - offset_x
        if x_center < 1e-6:
            edge.name = "axis_ext"
        else:
            edge.name = "kelvin_ext"
    outer_half_ext.faces.name = "air_outer"

    vertex = Vertex(Pnt(offset_x, 0, 0))
    vertex.name = "GND"

    shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

    # Identify periodic boundaries
    kelvin_int_edges = []
    kelvin_ext_edges = []
    for edge in shape.edges:
        if edge.name == "kelvin_int":
            kelvin_int_edges.append(edge)
        elif edge.name == "kelvin_ext":
            kelvin_ext_edges.append(edge)

    print(f"Found {len(kelvin_int_edges)} interior kelvin edges")
    print(f"Found {len(kelvin_ext_edges)} exterior kelvin edges")

    if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
        matched_pairs = 0
        for int_edge in kelvin_int_edges:
            int_z = int_edge.center.y
            for ext_edge in kelvin_ext_edges:
                ext_z = ext_edge.center.y
                if (int_z > 0 and ext_z > 0) or (int_z < 0 and ext_z < 0):
                    int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                    matched_pairs += 1
                    break
        print(f"Matched {matched_pairs} periodic edge pairs")

    return OCCGeometry(shape, dim=2)

def check_periodic_info(mesh, label=""):
    """Check periodic vertex mapping."""
    print(f"\n=== {label} ===")

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
    r_coord = x
    z_coord = y
    rho_prime_sq = (x - offset_x)**2 + y**2

    mu_inner = IfPos(x - sphere_radius, mu0, mu_r * mu0)
    mu_outer = (kelvin_radius**2 / rho_prime_sq) * mu0

    mu_cf = mesh.MaterialCF({"air_inner": mu_inner, "magnetic": mu_inner, "air_outer": mu_outer})

    fes_before = H1(mesh, order=order, dirichlet_bbnd="GND")
    fes = Periodic(fes_before)

    Omega = fes.TrialFunction()
    psi = fes.TestFunction()

    a_form = BilinearForm(fes)
    a_form += 2 * pi * r_coord * mu_cf * grad(Omega) * grad(psi) * dx("air_inner|magnetic")
    a_form += 2 * pi * (x - offset_x) * mu_cf * grad(Omega) * grad(psi) * dx("air_outer")
    a_form.Assemble()

    f = LinearForm(fes)
    f += 2 * pi * r_coord * H0 * mu_cf * grad(psi)[1] * dx("air_inner|magnetic")
    f += 2 * pi * (x - offset_x) * H0 * mu_cf * grad(psi)[1] * dx("air_outer")
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # Error estimator
    H_r = -grad(gfu)[0]
    H_z = -grad(gfu)[1] + H0

    W = Integrate(2 * pi * r_coord * mu_cf * (H_r**2 + H_z**2), mesh, definedon=mesh.Materials("air_inner|magnetic"))
    if W < 1e-20:
        W = 1.0

    B_r = mu_cf * H_r
    B_z = mu_cf * H_z

    Bfes = HDiv(mesh, order=order+1)
    Bint = GridFunction(Bfes)
    Bint.Set((B_r, B_z))

    err = (2 * pi * r_coord / mu_cf) * ((Bint[0] - B_r)**2 + (Bint[1] - B_z)**2) / W
    element_errors = Integrate(err, mesh, element_wise=True)
    total_error = sqrt(sum(element_errors))

    return fes.ndof, total_error, element_errors, gfu


print("=" * 60)
print("Creating geometry...")
geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.4)).Curve(2)

print(f"\nInitial mesh: {mesh.ne} elements, {mesh.nv} vertices")

v0, e0 = check_periodic_info(mesh, "Before Refinement")

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

v1, e1 = check_periodic_info(mesh, "After Refinement")

ndof2, error2, elem_errors2, gfu2 = solve_and_get_error(mesh, order)
print(f"\nIteration 2: DOFs={ndof2}, Error={error2:.6e}")

print("\n" + "=" * 60)
print("Summary:")
print(f"  Before Refine: {v0} vertex pairs, {e0} edge pairs")
print(f"  After Refine:  {v1} vertex pairs, {e1} edge pairs")
print(f"  Error change: {error:.6e} -> {error2:.6e}")
if error2 > error:
    print("  WARNING: Error INCREASED after refinement!")
else:
    print("  OK: Error decreased as expected")
print("=" * 60)
