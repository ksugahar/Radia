"""
2D A-formulation with Kelvin transformation for parallel wires
平行2線問題（2次元）のA法 + Kelvin変換 WITH ADAPTIVE MESH REFINEMENT

Adaptive algorithm:
  1. ZZ-type error estimator (flux recovery)
  2. Doerfler marking strategy
  3. Red-green refinement
  4. VTU and PNG output at each iteration

2D A-formulation:
  Variable: A_z (z-component of vector potential)
  Equation: -div(nu * grad(A_z)) = J_z
  B = curl(A) => Bx = dAz/dy, By = -dAz/dx

Kelvin transformation (full circle):
  - Interior domain: circle centered at origin, radius a
  - Exterior domain: circle centered at (0, y_offset), radius a
  - For 2D Laplace: invariant under Kelvin (no factor needed)

Stop condition: DOF >= 1e6 only
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
from netgen.occ import *
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
print("Parallel Wires - ADAPTIVE MESH REFINEMENT")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
wire_distance = 1.4     # Distance between wires (d) [m] -> wires at x=+/-0.7m
wire_radius = 0.02      # Wire cross-section radius [m]
a = 1.0                 # Kelvin boundary radius [m]
y_offset = 2.5          # Y-offset for exterior domain [m]
maxh_initial = 0.5      # Initial mesh size [m] (coarser for adaptive refinement)
order = 2               # Polynomial order

mu0 = 4 * pi * 1e-7     # Vacuum permeability [H/m]
nu0 = 1 / mu0           # Vacuum reluctivity [m/H]
I_total = 1.0           # Total current [A]
J0 = I_total / (pi * wire_radius**2)  # Current density [A/m^2]

# Adaptive mesh parameters
max_iterations = 8      # Stop after this many iterations
theta = 0.5             # Doerfler marking parameter

print(f"\nParameters:")
print(f"  Wire distance: {wire_distance} m")
print(f"  Wire radius: {wire_radius} m")
print(f"  Kelvin boundary: a = {a} m")
print(f"  Y-offset: {y_offset} m")
print(f"  Total current: +/- {I_total} A")
print(f"\nAdaptive mesh parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations")
print(f"  Doerfler theta: {theta}")


# ============================================================
# Geometry with fix_point on wire boundary (instead of GND)
# ============================================================
def create_geometry():
    """Create the parallel wires geometry with Kelvin transformation."""
    # Interior domain: full circle at origin
    inner_circle = Circle((0, 0), a).Face()
    inner_circle.name = "air_inner"

    # Wire 1 (negative current) at (+d/2, 0)
    wire1_full = MoveTo(wire_distance/2, 0).Circle(wire_radius).Face()
    wire1_full.name = "wire_minus"
    wire1_full.maxh = maxh_initial / 3

    # Wire 2 (positive current) at (-d/2, 0)
    wire2_full = MoveTo(-wire_distance/2, 0).Circle(wire_radius).Face()
    wire2_full.name = "wire_plus"
    wire2_full.maxh = maxh_initial / 3

    # Exterior domain: full circle at (0, y_offset)
    outer_circle = Circle((0, y_offset), a).Face()
    outer_circle.name = "air_outer"

    # Fix point on wire boundary for uniqueness (instead of GND at infinity)
    fix_point = Vertex(Pnt(wire_distance/2 + wire_radius, 0, 0))
    fix_point.name = "fix_point"

    # Subtract wires from inner air
    inner_air = inner_circle - wire1_full - wire2_full
    inner_air.name = "air_inner"

    # Name boundaries
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

    # Glue all shapes
    shape = Glue([inner_air, wire1_full, wire2_full, outer_circle, fix_point])

    # Identify periodic boundaries (Kelvin BC) by distance from centers
    # After Glue, edge names may not be preserved, so use geometry to identify
    # Note: Use edge.start instead of edge.center because full Circle() edges
    # return their geometric center (0,0), not a point on the edge itself
    kelvin_inner_edges = []
    kelvin_outer_edges = []
    for edge in shape.edges:
        # Check distance from origin for inner kelvin boundary (use edge.start)
        dist_from_origin = sqrt(edge.start.x**2 + edge.start.y**2)
        # Check distance from y_offset for outer kelvin boundary (use edge.start)
        dist_from_outer = sqrt(edge.start.x**2 + (edge.start.y - y_offset)**2)

        if abs(dist_from_origin - a) < 0.01:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        elif abs(dist_from_outer - a) < 0.01:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)

    print(f"  Found {len(kelvin_inner_edges)} interior kelvin edges")
    print(f"  Found {len(kelvin_outer_edges)} exterior kelvin edges")

    # Simple 1-to-1 identification (both circles have 1 edge each)
    for int_edge in kelvin_inner_edges:
        for ext_edge in kelvin_outer_edges:
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)

    print(f"  Matched {len(kelvin_inner_edges)} periodic edge pairs")

    geo = OCCGeometry(shape, dim=2)

    return geo


# ============================================================
# Solve A-formulation
# ============================================================
def solve_A_formulation(mesh, order):
    """Solve 2D A-formulation on given mesh."""
    # Dirichlet at fix_point on wire boundary
    fes_before = H1(mesh, order=order, dirichlet_bbnd="fix_point")
    fes = Periodic(fes_before)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    nu_cf = nu0

    a_form = BilinearForm(fes)
    a_form += nu_cf * grad(u) * grad(v) * dx
    a_form.Assemble()

    f = LinearForm(fes)
    f += J0 * v * dx("wire_plus")
    f += (-J0) * v * dx("wire_minus")
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    return fes, gfu, nu_cf


# ============================================================
# ZZ Error Estimator
# ============================================================
def compute_error_estimator(mesh, fes, gfu, nu_cf, solution_order):
    """Compute ZZ-type error estimator using HCurl flux recovery."""
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
# Analytical solution
# ============================================================
def B_analytical(x_val, y_val):
    """Analytical B field for two parallel wires."""
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
    """Compute numerical magnetic energy from FEM solution.

    W = (1/2) * integral(nu * |B|^2) = (1/2) * integral(nu * |grad(Az)|^2)
    """
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    W = 0.5 * Integrate(nu_cf * (B * B), mesh)
    return W


# ============================================================
# Plotting helpers
# ============================================================
theta_w = linspace(0, 2*pi, 50)
theta_k = linspace(0, 2*pi, 100)


def draw_wires(ax):
    """Draw wire cross-sections."""
    ax.fill(wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'red', alpha=0.5)
    ax.plot(wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'r-', linewidth=1)
    ax.fill(-wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'blue', alpha=0.5)
    ax.plot(-wire_distance/2 + wire_radius*cos(theta_w),
            wire_radius*sin(theta_w), 'b-', linewidth=1)


def draw_kelvin_boundary(ax, y_center=0):
    """Draw Kelvin boundary circle."""
    ax.plot(a * cos(theta_k), y_center + a * sin(theta_k), 'g--', linewidth=2, label='Kelvin boundary')


def generate_iteration_plot(iter_idx, n_elements, history, inner_error_map, outer_error_map, gfu, mesh):
    """Generate and save plot for a single iteration.

    Layout:
    - Top-left: Inner domain ZZ error map (log scale) + flux lines + colorbar
    - Top-right: Outer domain ZZ error map (log scale) + colorbar
    - Bottom-left: DOF vs Error Estimator convergence curve
    - Bottom-right: DOF vs Magnetic Energy
    """
    fig = plt.figure(figsize=(14, 14), dpi=150)

    # Colormap settings for error maps
    norm = Normalize(vmin=-10, vmax=-2)
    cmap_obj = plt.colormaps['jet']

    # ===== Top-left: Inner domain ZZ error map + flux lines =====
    ax1 = plt.subplot(2, 2, 1)

    # Plot error map
    if len(inner_error_map) > 0:
        polygons = [item[0] for item in inner_error_map]
        err_values = array([item[1] for item in inner_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='none', alpha=0.9)
        ax1.add_collection(pc)

        # Add colorbar for inner domain
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar1 = plt.colorbar(sm, ax=ax1)
        cbar1.set_label('$\\log_{10}$(ZZ Error)')

    # Draw flux lines (contour of Az)
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
                # Check if inside wire
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

    # Draw flux lines
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

    # ===== Top-right: Outer domain ZZ error map =====
    ax2 = plt.subplot(2, 2, 2)

    if len(outer_error_map) > 0:
        polygons = [item[0] for item in outer_error_map]
        err_values = array([item[1] for item in outer_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='none', alpha=0.9)
        ax2.add_collection(pc)

        # Add colorbar for outer domain
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar2 = plt.colorbar(sm, ax=ax2)
        cbar2.set_label('$\\log_{10}$(ZZ Error)')

    # Draw flux lines for outer domain
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

    # Draw flux lines for outer domain
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

    # ===== Bottom-left: DOF vs Error convergence curve =====
    ax3 = plt.subplot(2, 2, 3)

    ax3.loglog(history['ndof'][:iter_idx+1], history['error'][:iter_idx+1], 'ko-', linewidth=2, markersize=6, label='Adaptive')
    ax3.loglog(history['ndof'][iter_idx], history['error'][iter_idx], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    # Reference line
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

    png_file = os.path.join(script_dir, f"ParallelWires_2D_adaptive_iter_{iter_idx+1:02d}.png")
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.close()
    return png_file


def save_iteration_mat(iter_idx, history):
    """Save MAT file with history up to current iteration."""
    mat_filename = os.path.join(script_dir, f"ParallelWires_2D_adaptive_iter_{iter_idx+1:02d}.mat")
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
# Main: Adaptive Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating geometry...")
print("=" * 60)

geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.4)).Curve(2)

print(f"\nInitial mesh:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")

# History tracking
history = {
    'ndof': [],
    'elements': [],
    'error': [],
    'By_error_percent': [],
    'energy': []
}

# Store mesh data and solutions for each iteration
# Note: NGSolve mesh is modified in-place by Refine(), so we save vertex/element data
mesh_data_list = []  # List of (vertices, elements, n_elements, n_vertices) tuples

print("\n" + "=" * 60)
print("Starting Adaptive Mesh Refinement")
print("=" * 60)

for iteration in range(max_iterations):
    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}")
    print("=" * 60)

    # Solve
    fes, gfu, nu_cf = solve_A_formulation(mesh, order)

    # Compute error estimator
    element_errors = compute_error_estimator(mesh, fes, gfu, nu_cf, order)
    total_error = sqrt(sum(element_errors))

    # Separate inner and outer domain errors for debugging
    inner_errors = []
    outer_errors = []
    for el_nr, el in enumerate(mesh.Elements()):
        el_verts = [mesh[v].point for v in el.vertices]
        cy = sum(v[1] for v in el_verts) / len(el_verts)
        if cy < y_offset / 2:
            inner_errors.append(element_errors[el_nr])
        else:
            outer_errors.append(element_errors[el_nr])
    inner_total = sqrt(sum(inner_errors)) if inner_errors else 0
    outer_total = sqrt(sum(outer_errors)) if outer_errors else 0
    print(f"  Inner domain error: {inner_total:.6e} ({len(inner_errors)} elements)")
    print(f"  Outer domain error: {outer_total:.6e} ({len(outer_errors)} elements)")

    # Compute By error at center
    try:
        mip = mesh(0, 0.1)
        By_num = -grad(gfu)[0](mip)
        _, By_ana = B_analytical(0, 0.1)
        By_error_percent = abs(By_num - By_ana) / abs(By_ana) * 100
    except:
        By_error_percent = nan

    # Compute magnetic energy
    W_energy = compute_numerical_energy(mesh, gfu, nu_cf)

    # Record history
    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)
    history['By_error_percent'].append(By_error_percent)
    history['energy'].append(W_energy)

    # Store mesh data (vertices and elements) - mesh is modified in-place by Refine()
    vertices = [(mesh[v].point[0], mesh[v].point[1]) for v in mesh.vertices]
    elements = [([mesh[v].point for v in el.vertices], str(el.mat)) for el in mesh.Elements()]
    mesh_data_list.append((vertices, elements, mesh.ne, mesh.nv))

    # Compute element-wise ZZ error maps for inner and outer domains
    inner_error_map = []  # List of (polygon, log_error) tuples
    outer_error_map = []  # List of (polygon, log_error) tuples

    for el_nr, el in enumerate(mesh.Elements()):
        el_verts = [mesh[v].point for v in el.vertices]
        cy = sum(v[1] for v in el_verts) / len(el_verts)

        # Skip wire elements
        mat_name = str(el.mat)
        if 'wire' in mat_name:
            continue

        # Get ZZ error for this element
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

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")
    print(f"  By error at (0, 0.1): {By_error_percent:.3f}%")

    # Generate PNG and save MAT for this iteration
    png_file = generate_iteration_plot(iteration, mesh.ne, history, inner_error_map, outer_error_map, gfu, mesh)
    mat_file = save_iteration_mat(iteration, history)
    print(f"  PNG saved: {png_file}")
    print(f"  MAT saved: {mat_file}")

    # Mark and refine (skip on last iteration)
    if iteration < max_iterations - 1:
        marked = mark_elements_doerfler(mesh, element_errors, theta)
        print(f"  Marked elements: {len(marked)}")

        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, False)
        for el_nr in marked:
            mesh.SetRefinementFlag(ElementId(VOL, el_nr), True)

        mesh.Refine()


# ============================================================
# Final Validation
# ============================================================
print("\n" + "=" * 60)
print("VALIDATION: Final Field comparison")
print("=" * 60)

y_samples = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
print(f"\n{'y (m)':<10} {'Bx (NGSolve)':<15} {'Bx (Analytical)':<15} {'Error (%)':<10}")
print("-" * 50)

for y_val in y_samples:
    try:
        mip = mesh(0, y_val)
        grad_Az = grad(gfu)
        Bx_num = grad_Az[1](mip)
        Bx_ana, _ = B_analytical(0, y_val)
        err = abs(Bx_num - Bx_ana) / abs(Bx_ana) * 100 if abs(Bx_ana) > 1e-15 else 0
        print(f"{y_val:<10.2f} {Bx_num:<15.6e} {Bx_ana:<15.6e} {err:<10.3f}")
    except Exception as e:
        print(f"{y_val:<10.2f} Error: {e}")


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
final_mat = os.path.join(script_dir, "ParallelWires_2D_adaptive_with_Kelvin.mat")
sio.savemat(final_mat, {
    'ndof': array(history['ndof']),
    'elements': array(history['elements']),
    'error': array(history['error']),
    'By_error_percent': array(history['By_error_percent']),
    'energy': array(history['energy'])
})
print(f"\nFinal MAT saved: {final_mat}")

# Open the final iteration plot
final_png = os.path.join(script_dir, f"ParallelWires_2D_adaptive_iter_{max_iterations:02d}.png")
if os.path.exists(final_png):
    os.startfile(final_png)

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
