"""
2D A-formulation with Kelvin transformation for parallel wires
平行2線問題（2次元）のA法 + Kelvin変換 WITH METRIC-BASED REMESHING

Key improvement: Element-wise error-based local mesh size control
  - Computes ideal mesh size for each element based on its error
  - Uses the formula: h_ideal = h_current * (eta_target / eta_element)^(1/(p+1))
  - Generates meshsizefile for Netgen with local size constraints
  - Allows mesh coarsening where error is small

Algorithm:
  1. Solve on current mesh
  2. Compute ZZ error estimator per element
  3. For each element, compute ideal mesh size based on its error
  4. Write local size field to meshsizefile and regenerate mesh
  5. Repeat for 8 iterations

2D A-formulation:
  Variable: A_z (z-component of vector potential)
  Equation: -div(nu * grad(A_z)) = J_z
  B = curl(A) => Bx = dAz/dy, By = -dAz/dx
"""
import os
import glob

# Delete existing .png and .mat files in the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
for ext in ['*.png', '*.mat']:
    for f in glob.glob(os.path.join(script_dir, ext)):
        os.remove(f)
        print(f"Deleted: {f}")

from numpy import pi, sqrt, cos, sin, linspace, zeros, nan, isnan, meshgrid, array, log, mean, arctan2
from ngsolve import *
from ngsolve import TaskManager
from netgen.occ import *
from netgen.meshing import MeshingParameters
import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})

print("=" * 60)
print("2D A-formulation with Kelvin Transformation")
print("Parallel Wires - METRIC-BASED REMESHING")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
wire_distance = 1.4     # Distance between wires (d) [m] -> wires at x=+/-0.7m
wire_radius = 0.02      # Wire cross-section radius [m]
a = 1.0                 # Kelvin boundary radius [m]
y_offset = 2.5          # Y-offset for exterior domain [m]
maxh_initial = 0.3      # Initial mesh size [m]
order = 3               # Polynomial order

mu0 = 4 * pi * 1e-7     # Vacuum permeability [H/m]
nu0 = 1 / mu0           # Vacuum reluctivity [m/H]
I_total = 1.0           # Total current [A]
J0 = I_total / (pi * wire_radius**2)  # Current density [A/m^2]

# Adaptive mesh parameters
max_iterations = 10     # Stop after 10 iterations
h_min = 0.0005          # Minimum mesh size [m]
h_max = 0.3             # Maximum mesh size [m]
grading = 0.3           # Mesh grading parameter
eta_target_factor = 0.5 # Target error reduction factor per iteration

print(f"\nParameters:")
print(f"  Wire distance: {wire_distance} m")
print(f"  Wire radius: {wire_radius} m")
print(f"  Kelvin boundary: a = {a} m")
print(f"  Y-offset: {y_offset} m")
print(f"  Total current: +/- {I_total} A")
print(f"\nMetric-based remeshing parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations")
print(f"  h_min: {h_min} m, h_max: {h_max} m")


# ============================================================
# Create base geometry
# ============================================================
def create_base_geometry():
    """Create the parallel wires geometry with Kelvin transformation."""
    inner_circle = Circle((0, 0), a).Face()
    inner_circle.name = "air_inner"

    wire1_full = MoveTo(wire_distance/2, 0).Circle(wire_radius).Face()
    wire1_full.name = "wire_minus"

    wire2_full = MoveTo(-wire_distance/2, 0).Circle(wire_radius).Face()
    wire2_full.name = "wire_plus"

    outer_circle = Circle((0, y_offset), a).Face()
    outer_circle.name = "air_outer"

    fix_point = Vertex(Pnt(wire_distance/2 + wire_radius, 0, 0))
    fix_point.name = "fix_point"

    inner_air = inner_circle - wire1_full - wire2_full
    inner_air.name = "air_inner"

    for edge in inner_air.edges:
        dist = sqrt(edge.center.x**2 + edge.center.y**2)
        if abs(dist - a) < 0.01:
            edge.name = "kelvin_int"
        else:
            edge.name = "wire_bnd"

    for edge in wire1_full.edges:
        edge.name = "wire_minus_bnd"
    for edge in wire2_full.edges:
        edge.name = "wire_plus_bnd"
    for edge in outer_circle.edges:
        edge.name = "kelvin_ext"

    shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, fix_point])

    kelvin_inner_edges = []
    kelvin_outer_edges = []
    for edge in shape.edges:
        dist_from_origin = sqrt(edge.start.x**2 + edge.start.y**2)
        dist_from_outer = sqrt(edge.start.x**2 + (edge.start.y - y_offset)**2)

        if abs(dist_from_origin - a) < 0.01:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        elif abs(dist_from_outer - a) < 0.01:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)

    for int_edge in kelvin_inner_edges:
        for ext_edge in kelvin_outer_edges:
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)

    return shape


# ============================================================
# Generate mesh with local size field
# ============================================================
def generate_mesh_with_local_sizes(shape, local_sizes=None, default_maxh=0.3):
    import tempfile

    geo = OCCGeometry(shape, dim=2)
    mp = MeshingParameters(maxh=default_maxh, grading=grading)

    if local_sizes is not None and len(local_sizes) > 0:
        temp_dir = tempfile.gettempdir()
        meshsize_file = os.path.join(temp_dir, "meshsize_field.txt")
        with open(meshsize_file, 'w') as f:
            f.write(f"{len(local_sizes)}\n")
            for (px, py, h) in local_sizes:
                h_clamped = max(h_min, min(h_max, h))
                f.write(f"{px} {py} 0 {h_clamped}\n")
            f.write("0\n")
        ngmesh = geo.GenerateMesh(mp=mp, meshsizefilename=meshsize_file)
    else:
        ngmesh = geo.GenerateMesh(mp=mp)

    mesh = Mesh(ngmesh).Curve(2)
    return mesh


# ============================================================
# Solve A-formulation
# ============================================================
def solve_A_formulation(mesh, order):
    fes_before = H1(mesh, order=order, dirichlet_bbnd="fix_point")
    fes = Periodic(fes_before)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # 2D in-plane: Kelvin factor = 1 (kelvin_factor_2d_inplane_cf)
    nu_cf = nu0

    a_form = BilinearForm(fes)
    a_form += nu_cf * grad(u) * grad(v) * dx
    with TaskManager():
        a_form.Assemble()

        f = LinearForm(fes)
        f += J0 * v * dx("wire_plus")
        f += (-J0) * v * dx("wire_minus")
        f.Assemble()

        gfu = GridFunction(fes)
        gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        return fes, gfu, nu_cf


# ============================================================
# Compute error and local sizes
# ============================================================
def compute_error_and_local_sizes(mesh, fes, gfu, nu_cf, solution_order, prev_local_sizes=None):
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    W = Integrate(nu_cf * (B * B), mesh)
    if W < 1e-20:
        W = 1.0

    H = nu_cf * B
    recovery_order = solution_order - 1
    Hfes = HCurl(mesh, order=recovery_order, type1=False)
    Hint = GridFunction(Hfes)
    Hint.Set(H)

    err = ((Hint - H) * (Hint - H)) / nu_cf / W
    element_errors = Integrate(err, mesh, element_wise=True)

    nonzero_errors = [e for e in element_errors if e > 1e-20]
    eta_mean = sum(nonzero_errors) / len(nonzero_errors) if nonzero_errors else 1.0
    eta_target = eta_mean * eta_target_factor

    local_sizes = []
    exponent = 1.0 / (solution_order + 1)

    prev_size_dict = {}
    if prev_local_sizes is not None:
        for (px, py, h) in prev_local_sizes:
            key = (round(px, 4), round(py, 4))
            if key in prev_size_dict:
                prev_size_dict[key] = min(prev_size_dict[key], h)
            else:
                prev_size_dict[key] = h

    for el in mesh.Elements():
        verts = el.vertices
        pts = [mesh[v].point for v in verts]
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)

        if len(pts) >= 3:
            x0, y0 = pts[0][0], pts[0][1]
            x1, y1 = pts[1][0], pts[1][1]
            x2, y2 = pts[2][0], pts[2][1]
            area = 0.5 * abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))
            h_current = sqrt(area)
        else:
            h_current = 0.1

        eta = element_errors[el.nr]

        if eta > 1e-20:
            ratio = eta_target / eta
            ratio = max(0.3, min(3.0, ratio))
            h_ideal = h_current * (ratio ** exponent)
        else:
            h_ideal = min(h_current * 1.5, h_max)

        key = (round(cx, 4), round(cy, 4))
        if key in prev_size_dict:
            h_prev = prev_size_dict[key]
            h_ideal = min(h_ideal, h_prev * 1.1)

        h_ideal = max(h_min, min(h_max, h_ideal))
        local_sizes.append((cx, cy, h_ideal))

    n_wire_pts = 20
    for theta in linspace(0, 2*pi, n_wire_pts, endpoint=False):
        wx = wire_distance/2 + wire_radius * 1.2 * cos(theta)
        wy = wire_radius * 1.2 * sin(theta)
        local_sizes.append((wx, wy, h_min * 3))
        wx = -wire_distance/2 + wire_radius * 1.2 * cos(theta)
        wy = wire_radius * 1.2 * sin(theta)
        local_sizes.append((wx, wy, h_min * 3))

    total_error = sqrt(sum(element_errors))

    return element_errors, local_sizes, total_error


# ============================================================
# Analytical solution and numerical energy
# ============================================================
def B_analytical(x_val, y_val):
    r1 = sqrt((x_val - wire_distance/2)**2 + y_val**2)
    r2 = sqrt((x_val + wire_distance/2)**2 + y_val**2)

    if r1 < wire_radius or r2 < wire_radius:
        return nan, nan
    if r1 < 1e-12 or r2 < 1e-12:
        return nan, nan

    Hx1 =  I_total / (2 * pi) * y_val / (r1**2)
    Hy1 = -I_total / (2 * pi) * (x_val - wire_distance/2) / (r1**2)
    Hx2 = -I_total / (2 * pi) * y_val / (r2**2)
    Hy2 =  I_total / (2 * pi) * (x_val + wire_distance/2) / (r2**2)

    Bx = mu0 * (Hx1 + Hx2)
    By = mu0 * (Hy1 + Hy2)
    return Bx, By


def compute_numerical_energy(mesh, gfu, nu_cf):
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    W = 0.5 * Integrate(nu_cf * (B * B), mesh)
    return W


# ============================================================
# Plotting helpers
# ============================================================
theta_w = linspace(0, 2*pi, 50)
theta_k = linspace(0, 2*pi, 100)


def draw_wires(ax):
    ax.fill(wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'red', alpha=0.5)
    ax.plot(wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'r-', linewidth=1)
    ax.fill(-wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'blue', alpha=0.5)
    ax.plot(-wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'b-', linewidth=1)


def draw_kelvin_boundary(ax, y_center=0):
    ax.plot(a * cos(theta_k), y_center + a * sin(theta_k), 'g--', linewidth=2, label='Kelvin boundary')


def generate_iteration_plot(iter_idx, n_elements, history, inner_error_map, outer_error_map, gfu, mesh):
    """Generate and save plot for a single iteration."""
    fig = plt.figure(figsize=(14, 14), dpi=150)

    norm = Normalize(vmin=-10, vmax=-2)
    cmap_obj = plt.colormaps['jet']

    # ===== Top-left: Inner domain ZZ error map + flux lines =====
    ax1 = plt.subplot(2, 2, 1)

    if len(inner_error_map) > 0:
        polygons = [item[0] for item in inner_error_map]
        err_values = array([item[1] for item in inner_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='none', alpha=0.9)
        ax1.add_collection(pc)

        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar1 = plt.colorbar(sm, ax=ax1)
        cbar1.set_label('$\\log_{10}$(ZZ Error)')

    n_grid = 100
    x_grid = linspace(-a * 0.98, a * 0.98, n_grid)
    y_grid = linspace(-a * 0.98, a * 0.98, n_grid)
    X, Y = meshgrid(x_grid, y_grid)
    Az_grid = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            x_pt, y_pt = X[i, j], Y[i, j]
            r = sqrt(x_pt**2 + y_pt**2)
            if r < a - 0.01:
                r1 = sqrt((x_pt - wire_distance/2)**2 + y_pt**2)
                r2 = sqrt((x_pt + wire_distance/2)**2 + y_pt**2)
                if r1 > wire_radius and r2 > wire_radius:
                    try:
                        Az_grid[i, j] = gfu(mesh(x_pt, y_pt))
                    except:
                        Az_grid[i, j] = nan
                else:
                    Az_grid[i, j] = nan
            else:
                Az_grid[i, j] = nan

    Az_min = Az_grid[~isnan(Az_grid)].min() if not all(isnan(Az_grid.flatten())) else 0
    Az_max = Az_grid[~isnan(Az_grid)].max() if not all(isnan(Az_grid.flatten())) else 1
    n_levels = 20
    levels = linspace(Az_min, Az_max, n_levels)
    ax1.contour(X, Y, Az_grid, levels=levels, colors='k', linewidths=0.5, alpha=0.7)

    draw_wires(ax1)
    draw_kelvin_boundary(ax1, 0)
    ax1.set_xlabel('$x$ (m)')
    ax1.set_ylabel('$y$ (m)')
    ax1.set_title(f'Inner Domain: ZZ Error + Flux Lines ({len(inner_error_map)} elements)')
    ax1.set_aspect('equal')
    ax1.set_xlim(-a - 0.05, a + 0.05)
    ax1.set_ylim(-a - 0.05, a + 0.05)

    # ===== Top-right: Outer domain ZZ error map + flux lines =====
    ax2 = plt.subplot(2, 2, 2)

    if len(outer_error_map) > 0:
        polygons = [item[0] for item in outer_error_map]
        err_values = array([item[1] for item in outer_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='none', alpha=0.9)
        ax2.add_collection(pc)

        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar2 = plt.colorbar(sm, ax=ax2)
        cbar2.set_label('$\\log_{10}$(ZZ Error)')

    x_grid_out = linspace(-a * 0.98, a * 0.98, n_grid)
    y_grid_out = linspace(y_offset - a * 0.98, y_offset + a * 0.98, n_grid)
    X_out, Y_out = meshgrid(x_grid_out, y_grid_out)
    Az_grid_out = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            x_pt, y_pt = X_out[i, j], Y_out[i, j]
            r_out = sqrt(x_pt**2 + (y_pt - y_offset)**2)
            if r_out < a - 0.01:
                try:
                    Az_grid_out[i, j] = gfu(mesh(x_pt, y_pt))
                except:
                    Az_grid_out[i, j] = nan
            else:
                Az_grid_out[i, j] = nan

    Az_out_valid = Az_grid_out[~isnan(Az_grid_out)]
    if len(Az_out_valid) > 0:
        Az_min_out = Az_out_valid.min()
        Az_max_out = Az_out_valid.max()
        levels_out = linspace(Az_min_out, Az_max_out, n_levels)
        ax2.contour(X_out, Y_out, Az_grid_out, levels=levels_out, colors='k', linewidths=0.5, alpha=0.7)

    draw_kelvin_boundary(ax2, y_offset)
    ax2.set_xlabel('$x$ (m)')
    ax2.set_ylabel('$y$ (m)')
    ax2.set_title(f'Outer Domain: ZZ Error + Flux Lines ({len(outer_error_map)} elements)')
    ax2.set_aspect('equal')
    ax2.set_xlim(-a - 0.05, a + 0.05)
    ax2.set_ylim(y_offset - a - 0.05, y_offset + a + 0.05)

    # ===== Bottom-left: DOF vs Error =====
    ax3 = plt.subplot(2, 2, 3)

    ax3.loglog(history['ndof'][:iter_idx+1], history['error'][:iter_idx+1], 'ko-', linewidth=2, markersize=6, label='Metric-based')
    ax3.loglog(history['ndof'][iter_idx], history['error'][iter_idx], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    N_ref = 1e5
    err_ref_point = 3e-4
    ndof_line = array([1e3, 1e6])
    err_line = err_ref_point * (N_ref / ndof_line) ** (order / 2)
    ax3.loglog(ndof_line, err_line, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/2}})$')

    ax3.set_xlabel('DOFs')
    ax3.set_ylabel('Error Estimator')
    ax3.set_title('DOF vs Error Estimator')
    ax3.set_xlim(1e3, 1e6)
    ax3.set_ylim(1e-4, 1e-1)
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # ===== Bottom-right: DOF vs Magnetic Energy =====
    ax4 = plt.subplot(2, 2, 4)

    ax4.semilogx(history['ndof'][:iter_idx+1], history['energy'][:iter_idx+1], 'gs-', linewidth=2, markersize=6, label='Magnetic Energy')
    ax4.semilogx(history['ndof'][iter_idx], history['energy'][iter_idx], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    ax4.set_xlabel('DOFs')
    ax4.set_ylabel('Magnetic Energy (J/m)')
    ax4.set_title('DOF vs Magnetic Energy')
    ax4.set_xlim(1e3, 1e6)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    fig.suptitle(f'Solution order={order}, Iteration {iter_idx+1}: {n_elements} elements, DOFs={history["ndof"][iter_idx]}, Error={history["error"][iter_idx]:.2e}',
                 fontsize=14, fontweight='bold', y=1.02)

    png_file = os.path.join(script_dir, f"ParallelWires_2D_metric_iter_{iter_idx+1:02d}.png")
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close()
    return png_file


def save_iteration_mat(iter_idx, history):
    mat_filename = os.path.join(script_dir, f"ParallelWires_2D_metric_iter_{iter_idx+1:02d}.mat")
    mat_data = {
        'ndof': array(history['ndof']),
        'elements': array(history['elements']),
        'error': array(history['error']),
        'By_error_percent': array(history['By_error_percent']),
        'energy': array(history['energy'])
    }
    sio.savemat(mat_filename, mat_data)
    return mat_filename


# ============================================================
# Main: Metric-based Remeshing Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating base geometry...")
print("=" * 60)

base_shape = create_base_geometry()

print("\nGenerating initial mesh...")
mesh = generate_mesh_with_local_sizes(base_shape, local_sizes=None, default_maxh=maxh_initial)

print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")

history = {
    'ndof': [],
    'elements': [],
    'error': [],
    'By_error_percent': [],
    'energy': []
}

print("\n" + "=" * 60)
print("Starting Metric-Based Remeshing")
print("=" * 60)

local_sizes = None

for iteration in range(max_iterations):
    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}")
    print("=" * 60)

    fes, gfu, nu_cf = solve_A_formulation(mesh, order)

    element_errors, new_local_sizes, total_error = compute_error_and_local_sizes(
        mesh, fes, gfu, nu_cf, order, local_sizes
    )

    try:
        mip = mesh(0, 0.1)
        By_num = -grad(gfu)[0](mip)
        _, By_ana = B_analytical(0, 0.1)
        By_error_percent = abs(By_num - By_ana) / abs(By_ana) * 100
    except:
        By_error_percent = nan

    W_energy = compute_numerical_energy(mesh, gfu, nu_cf)

    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)
    history['By_error_percent'].append(By_error_percent)
    history['energy'].append(W_energy)

    # Compute element-wise ZZ error maps
    inner_error_map = []
    outer_error_map = []

    for el_nr, el in enumerate(mesh.Elements()):
        el_verts = [mesh[v].point for v in el.vertices]
        cy = sum(v[1] for v in el_verts) / len(el_verts)

        mat_name = str(el.mat)
        if 'wire' in mat_name:
            continue

        err_val = element_errors[el_nr]
        if err_val > 1e-20:
            log_err = log(err_val) / log(10)
        else:
            log_err = -20

        poly = [(v[0], v[1]) for v in el_verts]

        if cy < y_offset / 2:
            inner_error_map.append((poly, log_err))
        else:
            outer_error_map.append((poly, log_err))

    h_values = [h for (_, _, h) in new_local_sizes]
    h_min_actual = min(h_values)
    h_max_actual = max(h_values)

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")
    print(f"  By error at (0, 0.1): {By_error_percent:.3f}%")

    png_file = generate_iteration_plot(iteration, mesh.ne, history, inner_error_map, outer_error_map, gfu, mesh)
    mat_file = save_iteration_mat(iteration, history)
    print(f"  PNG saved: {png_file}")
    print(f"  MAT saved: {mat_file}")

    # Regenerate mesh with new size field (skip on last iteration)
    if iteration < max_iterations - 1:
        local_sizes = new_local_sizes
        base_shape = create_base_geometry()
        mesh = generate_mesh_with_local_sizes(base_shape, local_sizes, default_maxh=h_max_actual)
        print(f"  New mesh: {mesh.ne} elements, {mesh.nv} vertices")


# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Convergence History")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'Energy (J/m)':<14}")
print("-" * 52)
for i in range(len(history['ndof'])):
    print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
          f"{history['error'][i]:<12.4e} {history['energy'][i]:<14.6e}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
if history['error'][-1] > 0:
    print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
          f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")
print(f"  Energy: {history['energy'][0]:.6e} -> {history['energy'][-1]:.6e} J/m")

# Save final MAT file with full convergence history
final_mat = os.path.join(script_dir, "ParallelWires_2D_metric_remesh_with_Kelvin.mat")
sio.savemat(final_mat, {
    'ndof': array(history['ndof']),
    'elements': array(history['elements']),
    'error': array(history['error']),
    'By_error_percent': array(history['By_error_percent']),
    'energy': array(history['energy'])
})
print(f"\nFinal MAT saved: {final_mat}")

# Open the final iteration plot
final_png = os.path.join(script_dir, f"ParallelWires_2D_metric_iter_{max_iterations:02d}.png")
if os.path.exists(final_png):
    os.startfile(final_png)

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
