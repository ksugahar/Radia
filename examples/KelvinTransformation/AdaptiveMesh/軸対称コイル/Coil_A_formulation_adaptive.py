"""
Axisymmetric A-formulation with Z-OFFSET Kelvin transformation
Single coil magnetostatics problem WITH ADAPTIVE MESH REFINEMENT

Adaptive algorithm:
  1. ZZ-type error estimator (flux recovery)
  2. Doerfler marking strategy
  3. Laplacian mesh smoothing after refinement
  4. VTU output at each iteration

A-formulation for axisymmetric problems:
  Variable: u = r * A_theta (modified vector potential)
  Equation: -div(nu * grad(u) / r) = J_theta

Kelvin transformation with Z-offset:
  - Interior domain: half-circle centered at (0, 0), radius a
  - Exterior domain: half-circle centered at (0, z_offset), radius a
  - r' = r (unchanged!) - KEY SIMPLIFICATION
  - Kelvin factor for nu: (rho'/a)^2
"""
import os
from numpy import pi, sqrt, cos, sin, linspace, zeros, nan, isnan, meshgrid, array, log, mean, sign
from scipy.special import ellipk, ellipe
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Axisymmetric A-formulation with Z-OFFSET Kelvin")
print("ADAPTIVE MESH REFINEMENT")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
coil_radius = 0.6      # Coil center radius [m]
coil_width = 0.02      # Coil cross-section width [m]
coil_height = 0.02     # Coil cross-section height [m]
a = 1.0                # Kelvin boundary radius [m]
z_offset = 2.5         # Z-offset for exterior domain [m]
maxh_initial = 0.3     # Initial mesh size [m] (coarser for demonstration)
order = 2              # Polynomial order

mu0 = 4 * pi * 1e-7    # Vacuum permeability [H/m]
nu0 = 1 / mu0          # Vacuum reluctivity [m/H]
I_total = 1.0          # Total current [A]
J0 = I_total / (coil_width * coil_height)  # Current density [A/m^2]

# Adaptive mesh parameters
max_iterations = 6         # Number of refinement iterations
theta = 0.5                # Doerfler marking parameter

# Smoothing parameters
smoothing_iterations = 3   # Number of Laplacian smoothing iterations
smoothing_weight = 0.4     # Relaxation weight (0 < w < 1)

print(f"\nParameters:")
print(f"  Coil radius: {coil_radius} m")
print(f"  Coil cross-section: {coil_width} x {coil_height} m")
print(f"  Kelvin boundary: a = {a} m")
print(f"  Z-offset: {z_offset} m")
print(f"  Total current: {I_total} A")
print(f"  Current density: {J0:.2e} A/m^2")
print(f"\nAdaptive mesh parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Max iterations: {max_iterations}")
print(f"  Doerfler theta: {theta}")
print(f"  Smoothing iterations: {smoothing_iterations}")
print(f"  Smoothing weight: {smoothing_weight}")


# ============================================================
# Geometry with Z-offset (function to regenerate)
# ============================================================
def create_geometry():
    """Create the coil geometry with Kelvin transformation."""
    # Interior domain: FULL circle first
    wp1 = WorkPlane()
    inner_full = wp1.Circle(a).Face()
    inner_full.name = "air_inner"

    # Coil region (full)
    coil_full = MoveTo(coil_radius - coil_width/2, -coil_height/2).Rectangle(coil_width, coil_height).Face()
    coil_full.name = "coil"
    coil_full.maxh = maxh_initial / 2

    # Exterior domain: FULL circle at (0, z_offset)
    wp2 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp2.Circle(a).Face()
    outer_full.name = "air_outer"

    # GND point at center of exterior domain
    gnd_point = Vertex(Pnt(0, z_offset, 0))
    gnd_point.name = "GND"

    # Cut to get half-circles (keep r >= 0, i.e., x >= 0)
    cutter_left = MoveTo(-a - 0.1, -a - 0.1).Rectangle(a + 0.1, 2 * a + 0.2).Face()
    cutter_left_ext = MoveTo(-a - 0.1, z_offset - a - 0.1).Rectangle(a + 0.1, 2 * a + 0.2).Face()

    inner_half = inner_full - cutter_left
    coil_half = coil_full - cutter_left
    outer_half = outer_full - cutter_left_ext

    # Subtract coil from inner air
    inner_air = inner_half - coil_half
    inner_air.name = "air_inner"

    # Name boundaries and find Kelvin edges
    kelvin_inner_edges = []
    for edge in inner_air.edges:
        cx = edge.center.x
        cy = edge.center.y
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = sqrt(v1.p.x**2 + v1.p.y**2)
            is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
        except:
            is_kelvin_arc = False

        if cx < 0.01:
            edge.name = "axis"
        elif is_kelvin_arc:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        else:
            edge.name = "default"

    for edge in coil_half.edges:
        cx = edge.center.x
        if cx < 0.01:
            edge.name = "axis"
        else:
            edge.name = "coil_bnd"

    kelvin_outer_edges = []
    for edge in outer_half.edges:
        cx = edge.center.x
        cy = edge.center.y - z_offset
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
        except:
            is_kelvin_arc = False

        if cx < 0.01:
            edge.name = "axis_ext"
        elif is_kelvin_arc:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)
        else:
            edge.name = "default"

    # Identify periodic boundaries
    for int_edge in kelvin_inner_edges:
        int_y = int_edge.center.y
        for ext_edge in kelvin_outer_edges:
            ext_y = ext_edge.center.y - z_offset
            if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
                int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
                break

    shape = Glue([inner_air, coil_half, outer_half, gnd_point])
    geo = OCCGeometry(shape, dim=2)

    return geo


# ============================================================
# Solve A-formulation
# ============================================================
def solve_A_formulation(mesh, order):
    """Solve A-formulation on given mesh."""
    fes_before = H1(mesh, order=order, dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
    fes = Periodic(fes_before)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # r-weight
    r_weight = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin transformation factor
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_prime_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_factor = (rho_prime_safe / a)**2

    # Build reluctivity by region
    materials = mesh.GetMaterials()
    nu_list = []
    for mat in materials:
        if "outer" in mat.lower():
            nu_list.append(nu0 * kelvin_factor)
        else:
            nu_list.append(nu0)
    nu_cf = CoefficientFunction(nu_list)

    # Bilinear form
    a_form = BilinearForm(fes)
    a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
    a_form.Assemble()

    # Linear form
    f = LinearForm(fes)
    f += J0 * v * dx("coil")
    f.Assemble()

    # Solve
    gfu = GridFunction(fes)
    gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    return fes, gfu, nu_cf, r_weight


# ============================================================
# ZZ Error Estimator
# ============================================================
def compute_error_estimator(mesh, fes, gfu, nu_cf, r_weight):
    """
    Compute ZZ-type error estimator using flux recovery.
    For A-formulation: flux = nu/r * grad(u)
    """
    flux = nu_cf / r_weight * grad(gfu)

    fes_flux = HDiv(mesh, order=fes.globalorder - 1)
    gf_flux = GridFunction(fes_flux)
    gf_flux.Set(flux)

    err = (flux - gf_flux) * (flux - gf_flux)
    element_errors = Integrate(err, mesh, element_wise=True)

    return element_errors


# ============================================================
# Doerfler Marking
# ============================================================
def mark_elements_doerfler(mesh, element_errors, theta):
    """Doerfler marking strategy."""
    sorted_indices = sorted(range(len(element_errors)),
                           key=lambda i: element_errors[i], reverse=True)

    total_error = sum(element_errors)
    marked_error = 0
    marked = []

    for idx in sorted_indices:
        marked.append(idx)
        marked_error += element_errors[idx]
        if marked_error >= theta * total_error:
            break

    return marked


# ============================================================
# Laplacian Smoothing
# ============================================================
def laplacian_smoothing(mesh, n_iter, weight):
    """
    Laplacian mesh smoothing with boundary constraints.
    Interior vertices: move toward neighbor centroid.
    Boundary/corner vertices: keep fixed.
    """
    print(f"    Applying Laplacian smoothing ({n_iter} iterations)...")

    ngmesh = mesh.ngmesh
    pnts = ngmesh.Points()
    n_vertices = len(pnts)

    # Build vertex adjacency
    neighbors = {i: set() for i in range(1, n_vertices + 1)}
    for el in ngmesh.Elements2D():
        verts = [v.nr for v in el.vertices]
        for i, v in enumerate(verts):
            for j, w in enumerate(verts):
                if i != j:
                    neighbors[v].add(w)

    # Identify boundary vertices
    boundary_vertices = set()
    boundary_edges = {}

    for seg in ngmesh.Elements1D():
        verts = [v.nr for v in seg.vertices]
        for v in verts:
            boundary_vertices.add(v)
            if v not in boundary_edges:
                boundary_edges[v] = []
        if len(verts) == 2:
            boundary_edges[verts[0]].append(verts[1])
            boundary_edges[verts[1]].append(verts[0])

    # Corner vertices
    corner_vertices = set()
    for v in boundary_vertices:
        if len(boundary_edges.get(v, [])) != 2:
            corner_vertices.add(v)

    # Smoothing iterations
    for iteration in range(n_iter):
        new_positions = {}

        for vid in range(1, n_vertices + 1):
            px = pnts[vid][0]
            py = pnts[vid][1]

            if vid in corner_vertices:
                new_positions[vid] = (px, py)
            elif vid in boundary_vertices:
                new_positions[vid] = (px, py)
            else:
                if len(neighbors[vid]) > 0:
                    avg_x = sum(pnts[n][0] for n in neighbors[vid]) / len(neighbors[vid])
                    avg_y = sum(pnts[n][1] for n in neighbors[vid]) / len(neighbors[vid])
                    new_x = px + weight * (avg_x - px)
                    new_y = py + weight * (avg_y - py)
                    new_positions[vid] = (new_x, new_y)
                else:
                    new_positions[vid] = (px, py)

        for vid in range(1, n_vertices + 1):
            pnts[vid][0] = new_positions[vid][0]
            pnts[vid][1] = new_positions[vid][1]

    ngmesh.Update()
    print(f"    Smoothing completed")


# ============================================================
# VTK Output
# ============================================================
def output_vtk(mesh, gfu, iteration):
    """Output mesh and solution to VTK file."""
    vtk_filename = f"coil_iter_{iteration:02d}"
    vtk = VTKOutput(mesh, coefs=[gfu], names=["u_rAtheta"],
                    filename=vtk_filename, subdivision=0)
    vtk.Do()
    return vtk_filename + ".vtu"


# ============================================================
# Analytical solution
# ============================================================
def CircleEH(OPr, OPz, r_coil):
    """Compute magnetic field for a circular current loop."""
    R = abs(OPr)
    Z = OPz
    r = r_coil

    k2 = abs(4 * r * R) / ((r + R)**2 + Z**2)

    if k2 > 0 and k2 < 1:
        K = ellipk(k2)
        E = ellipe(k2)
        Hr = sign(OPr) * (1/(2*pi)) * (Z/R) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 + R**2 + Z**2) / ((r-R)**2 + Z**2) * E - K)
        Hz = (1/(2*pi)) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 - R**2 - Z**2) / ((r-R)**2 + Z**2) * E + K)
    elif k2 == 0:
        K = ellipk(k2)
        E = ellipe(k2)
        Hr = 0.0
        Hz = (1/(2*pi)) / sqrt((r+R)**2 + Z**2) * \
             ((r**2 - R**2 - Z**2) / ((r-R)**2 + Z**2) * E + K)
    else:
        Hr = 0.0
        Hz = 0.0

    return Hr, Hz


def Bz_analytical_on_axis(z_val):
    """Analytical Bz on z-axis."""
    if abs(z_val) < 1e-12:
        return mu0 * I_total / (2 * coil_radius)
    Hr, Hz = CircleEH(0.001, z_val, coil_radius)
    return mu0 * I_total * Hz


# ============================================================
# Main: Adaptive Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating geometry...")
print("=" * 60)

geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.5))

print(f"\nInitial mesh:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")

# History tracking
history = {
    'ndof': [],
    'elements': [],
    'error': [],
    'Bz_error_percent': []
}

print("\n" + "=" * 60)
print("Starting Adaptive Mesh Refinement")
print("=" * 60)

for iteration in range(max_iterations):
    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}/{max_iterations}")
    print("=" * 60)

    # Solve
    fes, gfu, nu_cf, r_weight = solve_A_formulation(mesh, order)

    # Compute error estimator
    element_errors = compute_error_estimator(mesh, fes, gfu, nu_cf, r_weight)
    total_error = sqrt(sum(element_errors))

    # Compute Bz error at center
    try:
        mip = mesh(0.01, 0)
        Bz_center = grad(gfu)[0](mip) / 0.01
        Bz_analytical = mu0 * I_total / (2 * coil_radius)
        Bz_error_percent = abs(Bz_center - Bz_analytical) / abs(Bz_analytical) * 100
    except:
        Bz_error_percent = nan

    # Record history
    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)
    history['Bz_error_percent'].append(Bz_error_percent)

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")
    print(f"  Bz error at center: {Bz_error_percent:.3f}%")

    # Output VTK
    vtk_file = output_vtk(mesh, gfu, iteration)
    print(f"  VTK saved: {vtk_file}")

    # Mark and refine (except last iteration)
    if iteration < max_iterations - 1:
        marked = mark_elements_doerfler(mesh, element_errors, theta)
        print(f"  Marked elements: {len(marked)}")

        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, False)
        for el_nr in marked:
            mesh.SetRefinementFlag(ElementId(VOL, el_nr), True)

        mesh.Refine()

        # Apply Laplacian smoothing
        laplacian_smoothing(mesh, smoothing_iterations, smoothing_weight)


# ============================================================
# Final Validation
# ============================================================
print("\n" + "=" * 60)
print("VALIDATION: Field on z-axis (r=0.01)")
print("=" * 60)

z_samples = [0.0, 0.2, 0.4, 0.6, 0.8]
print(f"{'z (m)':<10} {'Bz (NGSolve)':<15} {'Bz (Analytical)':<15} {'Error (%)':<10}")
print("-" * 50)

for z_val in z_samples:
    try:
        mip = mesh(0.01, z_val)
        Bz_num = grad(gfu)[0](mip) / 0.01
        Bz_ana = Bz_analytical_on_axis(z_val)
        err = abs(Bz_num - Bz_ana) / abs(Bz_ana) * 100 if abs(Bz_ana) > 1e-15 else 0
        print(f"{z_val:<10.2f} {Bz_num:<15.6e} {Bz_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{z_val:<10.2f} Error: {e}")


# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Final Statistics")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'Bz Err(%)':<10}")
print("-" * 48)
for i in range(len(history['ndof'])):
    print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
          f"{history['error'][i]:<12.4e} {history['Bz_error_percent'][i]:<10.3f}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
      f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")
print(f"  Bz error: {history['Bz_error_percent'][0]:.3f}% -> {history['Bz_error_percent'][-1]:.3f}%")


# ============================================================
# Compute profiles for visualization
# ============================================================
print("\nComputing profiles...")

# Z-axis profile
z_profile = linspace(-0.9, 0.9, 91)
Bz_numerical = zeros(len(z_profile))
Bz_analytical_profile = zeros(len(z_profile))

for i, z_val in enumerate(z_profile):
    Bz_analytical_profile[i] = Bz_analytical_on_axis(z_val)
    try:
        mip = mesh(0.01, z_val)
        Bz_numerical[i] = grad(gfu)[0](mip) / 0.01
    except:
        Bz_numerical[i] = nan

# 2D field
r_grid = linspace(0.02, a - 0.05, 41)
z_grid = linspace(-a + 0.05, a - 0.05, 41)
[rr, zz] = meshgrid(r_grid, z_grid)

Br_field = zeros(rr.shape)
Bz_field = zeros(rr.shape)

for iz in range(len(z_grid)):
    for ir in range(len(r_grid)):
        try:
            mip = mesh(r_grid[ir], z_grid[iz])
            grad_u = grad(gfu)
            Br_field[iz, ir] = -grad_u[1](mip) / r_grid[ir]
            Bz_field[iz, ir] = grad_u[0](mip) / r_grid[ir]
        except:
            Br_field[iz, ir] = nan
            Bz_field[iz, ir] = nan


# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})

fig = plt.figure(figsize=(14, 12), dpi=150)


def draw_coil(ax):
    coil_x = [coil_radius - coil_width/2, coil_radius + coil_width/2,
              coil_radius + coil_width/2, coil_radius - coil_width/2,
              coil_radius - coil_width/2]
    coil_y = [-coil_height/2, -coil_height/2, coil_height/2, coil_height/2, -coil_height/2]
    ax.fill(coil_x, coil_y, 'red', alpha=0.5)
    ax.plot(coil_x, coil_y, 'r-', linewidth=2)


def draw_kelvin_boundary(ax):
    theta_arc = linspace(-pi/2, pi/2, 100)
    ax.plot(a * cos(theta_arc), a * sin(theta_arc), 'g--', linewidth=2, label='Kelvin boundary')


# Plot 1: Final mesh (inner domain only)
ax1 = plt.subplot(2, 2, 1)
for el in mesh.Elements():
    vertices = [mesh[v].point for v in el.vertices]
    # Only plot inner domain (y < z_offset/2)
    avg_y = sum(v[1] for v in vertices) / len(vertices)
    if avg_y < z_offset / 2:
        xs = [v[0] for v in vertices] + [vertices[0][0]]
        ys = [v[1] for v in vertices] + [vertices[0][1]]
        ax1.plot(xs, ys, 'b-', linewidth=0.3)

draw_coil(ax1)
draw_kelvin_boundary(ax1)
ax1.set_xlabel('$r$ (m)')
ax1.set_ylabel('$z$ (m)')
ax1.set_title(f'Final Mesh ({mesh.ne} elements, with smoothing)')
ax1.set_aspect('equal')
ax1.set_xlim(-0.05, a + 0.05)
ax1.set_ylim(-a - 0.05, a + 0.05)
ax1.legend(loc='upper right', fontsize=8)

# Plot 2: B field streamlines
ax2 = plt.subplot(2, 2, 2)
Br_plot = Br_field.copy()
Bz_plot = Bz_field.copy()
Br_plot[isnan(Br_plot)] = 0
Bz_plot[isnan(Bz_plot)] = 0
ax2.streamplot(rr, zz, Br_plot, Bz_plot, color='blue', linewidth=0.8, density=1.5)
draw_coil(ax2)
draw_kelvin_boundary(ax2)
ax2.set_xlabel('$r$ (m)')
ax2.set_ylabel('$z$ (m)')
ax2.set_title('$\\mathbf{B}$ Field (Z-offset Kelvin)')
ax2.set_aspect('equal')
ax2.set_xlim(0, a)
ax2.set_ylim(-a, a)

# Plot 3: Convergence
ax3 = plt.subplot(2, 2, 3)
ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Adaptive')
ndof_ref = array(history['ndof'])
err_ref = history['error'][0] * (ndof_ref[0] / ndof_ref) ** (order / 2)
ax3.loglog(ndof_ref, err_ref, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/2}})$')
ax3.set_xlabel('DOFs')
ax3.set_ylabel('Error Estimator')
ax3.set_title('Convergence')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xticks([10**3, 10**4])
ax3.set_xticklabels(['$10^3$', '$10^4$'])

# Plot 4: Bz profile comparison
ax4 = plt.subplot(2, 2, 4)
valid = ~isnan(Bz_numerical)
ax4.plot(z_profile[valid], Bz_numerical[valid] * 1e6, 'k-', linewidth=2, label='NGSolve (Adaptive)')
ax4.plot(z_profile, Bz_analytical_profile * 1e6, 'r--', linewidth=1.5, label='Analytical')
ax4.set_xlabel('$z$ (m)')
ax4.set_ylabel('$B_z$ ($\\mu$T)')
ax4.set_title('$B_z$ Profile on $z$-axis ($r \\approx 0$)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

plt.tight_layout()

png_file = "coil_adaptive_convergence.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
