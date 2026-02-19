"""
Debug script to compare solution quality at different orders
"""
from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Comparing solution quality across orders")
print("=" * 60)

# Parameters
a = 1.0
y_offset = 2.5
wire_distance = 1.4
wire_radius = 0.02
mu0 = 4 * pi * 1e-7
nu0 = 1 / mu0
I_total = 1.0
J0 = I_total / (pi * wire_radius**2)

def B_analytical(x_val, y_val):
    """Analytical B field for two parallel wires."""
    r1 = sqrt((x_val - wire_distance/2)**2 + y_val**2)
    r2 = sqrt((x_val + wire_distance/2)**2 + y_val**2)

    if r1 < wire_radius or r2 < wire_radius:
        return float('nan'), float('nan')
    if r1 < 1e-12 or r2 < 1e-12:
        return float('nan'), float('nan')

    Hx1 =  I_total / (2 * pi) * y_val / (r1**2)
    Hy1 = -I_total / (2 * pi) * (x_val - wire_distance/2) / (r1**2)
    Hx2 = -I_total / (2 * pi) * y_val / (r2**2)
    Hy2 =  I_total / (2 * pi) * (x_val + wire_distance/2) / (r2**2)

    Bx = mu0 * (Hx1 + Hx2)
    By = mu0 * (Hy1 + Hy2)
    return Bx, By

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

shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, fix_point])

kelvin_int = [e for e in shape.edges if e.name == "kelvin_int"]
kelvin_ext = [e for e in shape.edges if e.name == "kelvin_ext"]

for int_edge in kelvin_int:
    int_x = int_edge.center.x
    for ext_edge in kelvin_ext:
        ext_x = ext_edge.center.x
        if (int_x > 0 and ext_x > 0) or (int_x < 0 and ext_x < 0):
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
            break

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.3)).Curve(2)

print(f"Mesh: {mesh.ne} elements")

# Test points
test_points = [(0, 0.1), (0, 0.3), (0, 0.5), (0, 0.7), (0.3, 0.3)]

for order in [2, 3, 4]:
    print(f"\n{'='*60}")
    print(f"Order = {order}")
    print("=" * 60)

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

    print(f"\nField comparison at test points:")
    print(f"{'Point':<15} {'By_num':<15} {'By_ana':<15} {'Error(%)':<10}")
    print("-" * 55)

    for px, py in test_points:
        try:
            mip = mesh(px, py)
            By_num = -grad(gfu)[0](mip)
            _, By_ana = B_analytical(px, py)
            err = abs(By_num - By_ana) / abs(By_ana) * 100
            print(f"({px}, {py})"[:15].ljust(15) + f"{By_num:<15.6e} {By_ana:<15.6e} {err:<10.3f}")
        except Exception as e:
            print(f"({px}, {py}): Error - {e}")

    # Check solution continuity at Kelvin boundary
    print(f"\nKelvin boundary continuity check:")
    theta_vals = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]  # Various angles on circle
    for theta in theta_vals:
        # Point on inner Kelvin boundary (radius a from origin)
        x_int = a * 0.9 * sqrt(1 - (theta/2)**2)  # Inside boundary
        y_int = a * 0.9 * theta / 2

        # Point on outer Kelvin boundary
        x_ext = a * 0.9 * sqrt(1 - (theta/2)**2)
        y_ext = y_offset + a * 0.9 * theta / 2

        try:
            mip_int = mesh(x_int, y_int)
            mip_ext = mesh(x_ext, y_ext)
            Az_int = gfu(mip_int)
            Az_ext = gfu(mip_ext)
            diff = abs(Az_int - Az_ext)
            print(f"  theta={theta:.2f}: Az_int={Az_int:.6e}, Az_ext={Az_ext:.6e}, diff={diff:.6e}")
        except:
            pass
