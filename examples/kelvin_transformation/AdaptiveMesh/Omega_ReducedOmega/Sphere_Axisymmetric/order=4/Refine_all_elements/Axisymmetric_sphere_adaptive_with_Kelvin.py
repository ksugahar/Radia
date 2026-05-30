"""
Axisymmetric Omega Method for Magnetostatics with Kelvin Transformation
WITH UNIFORM MESH REFINEMENT (All Elements Refined)

Problem: Magnetic sphere (mu_r=100) in uniform z-directed background field

Formulation (Total Scalar Potential):
- H = -grad(Omega)
- From div(B) = 0: div(mu * grad(Omega)) = 0

Kelvin transformation:
- Maps infinite exterior domain to finite half-circle
- Permeability transformation: mu'(rho') = (R/rho')^2 * mu0
- Periodic BC couples interior (rho=R) with exterior (rho'=R)

Uniform refinement algorithm:
  1. Solve on current mesh
  2. Compute ZZ-type error estimator for monitoring (not for marking)
  3. Refine ALL elements uniformly
  4. VTU output at each iteration

Stop condition: DOF >= max_dof
"""
import os
import glob

# Delete existing output files
script_dir = os.path.dirname(os.path.abspath(__file__))
print("Deleting existing output files...")
deleted_count = 0
for ext in ['*.png', '*.vtu', '*.mat']:
    for f in glob.glob(os.path.join(script_dir, ext)):
        try:
            os.remove(f)
            print(f"  Deleted: {f}")
            deleted_count += 1
        except Exception as e:
            print(f"  Failed to delete {f}: {e}")
print(f"Deleted {deleted_count} files.")

from numpy import pi, sqrt, linspace, zeros, nan, isnan, meshgrid, array, log, log10, sin, cos
from ngsolve import *
from netgen.occ import *
import scipy.io as sio
from radia.kelvin_source import kelvin_mu_factor_axisym_cf, build_material_cf

# Import matplotlib for plotting
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})

print("=" * 60)
print("Axisymmetric Omega Method with Kelvin Transform")
print("ADAPTIVE MESH REFINEMENT (Local Refinement)")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
sphere_radius = 0.5      # Magnetic sphere radius [m]
kelvin_radius = 1.0      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1.0  # [A/m]

# Mesh parameters
maxh_initial = 0.5       # Initial mesh size (adjusted for ~100 elements)
order = 4                # Finite element order

# Offset for exterior domain (z-direction for axisymmetric)
offset_z = 3.0

# Uniform refinement parameters
max_iterations = 6       # Stop after 6 iterations

print(f"\nProblem parameters:")
print(f"  Sphere radius: {sphere_radius} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")
print(f"\nUniform refinement parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations")
print(f"  Refinement: ALL elements (uniform)")


# ============================================================
# Geometry Definition using HALF-CIRCLES (axisymmetric, r >= 0)
# ============================================================
def create_geometry():
    """Create geometry with periodic boundary conditions."""
    print("\nCreating half-circle geometry with periodic boundary conditions...")

    # ===== INTERIOR DOMAIN (half-circle, r >= 0) =====
    outer_half_int = Circle((0, 0), kelvin_radius).Face()
    cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
    outer_half_int = outer_half_int - cutter_int  # Keep x >= 0 part

    inner_half_int = Circle((0, 0), sphere_radius).Face()
    inner_half_int = inner_half_int - cutter_int  # Keep x >= 0 part

    # Air region (half-annulus)
    air_inner = outer_half_int - inner_half_int

    # Name boundaries
    for edge in air_inner.edges:
        x_center = edge.center.x
        dist = sqrt(edge.center.x**2 + edge.center.y**2)
        if x_center < 1e-6:  # On z-axis (r = 0)
            edge.name = "axis_int"
        elif abs(dist - kelvin_radius) < kelvin_radius * 0.2:
            edge.name = "kelvin_int"
        else:
            edge.name = "sphere"
    air_inner.faces.name = "air_inner"

    # Magnetic material (half-circle)
    for edge in inner_half_int.edges:
        if edge.center.x < 1e-6:
            edge.name = "axis_int"
        else:
            edge.name = "sphere"
    inner_half_int.faces.name = "magnetic"

    # ===== EXTERIOR DOMAIN (half-circle, z-offset for axisymmetric) =====
    outer_half_ext = Circle((0, offset_z), kelvin_radius).Face()
    cutter_ext = MoveTo(-kelvin_radius - 1, offset_z - kelvin_radius - 1).Rectangle(kelvin_radius + 1, 2 * kelvin_radius + 2).Face()
    outer_half_ext = outer_half_ext - cutter_ext  # Keep x >= 0 part (r >= 0)

    for edge in outer_half_ext.edges:
        r_center = edge.center.x  # r coordinate
        if r_center < 1e-6:  # On axis (r = 0)
            edge.name = "axis_ext"
        else:
            edge.name = "kelvin_ext"
    outer_half_ext.faces.name = "air_outer"

    # ===== GND VERTEX (center of exterior domain) =====
    vertex = Vertex(Pnt(0, offset_z, 0))
    vertex.name = "GND"

    # Glue all domains
    shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

    # ===== IDENTIFY ALL PERIODIC BOUNDARY EDGE PAIRS =====
    print("Identifying periodic boundaries...")

    # Find all kelvin edges
    kelvin_int_edges = []
    kelvin_ext_edges = []
    for edge in shape.edges:
        if edge.name == "kelvin_int":
            kelvin_int_edges.append(edge)
        elif edge.name == "kelvin_ext":
            kelvin_ext_edges.append(edge)

    print(f"  Found {len(kelvin_int_edges)} interior kelvin edges")
    print(f"  Found {len(kelvin_ext_edges)} exterior kelvin edges")

    # Match edges by r-coordinate (x) of center
    # Interior: z > 0 or z < 0, Exterior: (z - offset_z) > 0 or < 0
    if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
        matched_pairs = 0
        for int_edge in kelvin_int_edges:
            int_z = int_edge.center.y  # z in interior
            for ext_edge in kelvin_ext_edges:
                ext_z_local = ext_edge.center.y - offset_z  # z' in exterior (local coord)
                # Match: both positive z or both negative z
                if (int_z > 0 and ext_z_local > 0) or (int_z < 0 and ext_z_local < 0):
                    int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                    matched_pairs += 1
                    break
        print(f"  Total matched pairs: {matched_pairs}")

    return OCCGeometry(shape, dim=2)


# ============================================================
# Solve Omega formulation (Omega-Reduced Omega with Dirichlet Lifting)
# ============================================================
def solve_omega_formulation(mesh, fe_order):
    """
    Solve Omega-Reduced Omega formulation on given mesh.
    Following EMPY_Analysis/Omega_ReducedOmega.py convention.

    - Total region (magnetic): H = grad(Omega)
    - Reduced region (air): H = H_s - grad(Omega), Omega is perturbation
    """
    import numpy as np_arr

    # H1 space with Dirichlet BC at GND (infinity)
    fes_before = H1(mesh, order=fe_order, dirichlet_bbnd="GND")
    fes = Periodic(fes_before)  # Apply periodic BC

    # Trial and test functions
    Omega = fes.TrialFunction()
    psi = fes.TestFunction()

    # Coordinate functions (x = r, y = z)
    r_coord = x  # Radial coordinate
    z_coord = y  # Axial coordinate

    # Distance squared from exterior domain center (z-offset for axisymmetric)
    rho_prime_sq = r_coord**2 + (z_coord - offset_z)**2

    # Transformed permeability for exterior domain (Kelvin) via centralized helper
    mu_kelvin_factor = kelvin_mu_factor_axisym_cf(
        z_offset=offset_z, R=kelvin_radius,
        r_coord=r_coord, z_coord=z_coord,
    )
    Mu = build_material_cf(
        mesh, mu0, mu_kelvin_factor,
        outer_keyword="air_outer",
        overrides={"magnetic": mu_r * mu0},
    )
    mu_kelvin = mu0 * mu_kelvin_factor  # alias for post-processing

    # r-weight for axisymmetric formulation (z-offset: r_weight = r for all domains)
    r_weight = r_coord
    r_weight_inner = r_coord  # For boundary integrals

    # ========================================
    # Omega-Reduced Omega Formulation
    # (Following Omega_ReducedOmega.py)
    # ========================================

    # Source potential: Omega_s = H0 * z (so grad(Omega_s) = (0, H0) = H_s)
    Omega_s = H0 * z_coord

    # Source magnetic field H_s = (0, H0) and B_s = mu0 * H_s
    Hs = CoefficientFunction((0.0, H0))
    Bs = CoefficientFunction((0.0, mu0 * H0))

    # For Kelvin-transformed exterior domain:
    rho_prime = sqrt(rho_prime_sq + 1e-20)
    Hz_exterior = -(rho_prime / kelvin_radius)**2 * H0
    Bs_exterior = CoefficientFunction((0.0, mu0 * Hz_exterior))

    # ===== Bilinear Form =====
    a = BilinearForm(fes)
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("magnetic")
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_inner")
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_outer")
    a.Assemble()

    # ===== Dirichlet BC on Total/Reduced interface =====
    gfOmega = GridFunction(fes)
    gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # ===== Linear Form =====
    f = LinearForm(fes)
    # Source term in Reduced region only (air_inner)
    # NOTE: Do NOT include air_outer (Kelvin region) in source term
    f += Mu * grad(gfOmega) * grad(psi) * r_weight * dx("air_inner")
    f.Assemble()

    # Extract FreeDofs part of f
    fcut = np_arr.array(f.vec.FV())[fes.FreeDofs()]
    np_arr.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

    # Neumann boundary condition on sphere
    normal = specialcf.normal(mesh.dim)
    f += (normal * Bs) * psi * r_weight_inner * ds("sphere")
    f.Assemble()

    # Solve
    gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # ===== Compute full B field for error estimation =====
    # Following Omega_ReducedOmega.py post-processing
    fesOt = H1(mesh, order=fe_order, definedon="magnetic")
    fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")

    Ot = GridFunction(fesOt)
    Orr = GridFunction(fesOr)
    Oxr = GridFunction(fesOr)

    Ot.Set(gfOmega, VOL, definedon="magnetic")
    Orr.Set(gfOmega, VOL, definedon="air_inner|air_outer")
    Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # B field: Bt in Total, Br + Bs in Reduced
    Bt = grad(Ot) * Mu  # Total region
    Br = (grad(Orr) - grad(Oxr)) * mu0  # Reduced region perturbation

    # Source field Bs per region
    Bs_dict = {
        "air_inner": Bs,
        "air_outer": Bs_exterior,
        "magnetic": CoefficientFunction((0.0, 0.0))
    }
    Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

    # Total B field (for error estimation)
    BField = Bt + Br + Bs_cf

    return fes, gfOmega, Mu, r_weight, BField


# ============================================================
# ZZ Error Estimator with H(div) Recovery
# ============================================================
def compute_error_estimator(mesh, fes, gfu, mu, r_weight, BField):
    """
    Compute ZZ-type error estimator using H(div) flux recovery.
    Uses full B field (BField) for error estimation.
    """
    flux = BField  # Use full B field for error estimation

    # Total magnetic energy for normalization (using B field norm)
    W = Integrate(InnerProduct(flux, flux) * r_weight, mesh)
    if abs(W) < 1e-20:
        W = 1.0

    # H(div) recovery (order - 1 for ZZ estimator)
    recovery_order = max(1, fes.globalorder - 1)
    fes_flux = HDiv(mesh, order=recovery_order)
    gf_flux = GridFunction(fes_flux)
    gf_flux.Set(flux)

    # Error = difference between direct and recovered flux
    err = (flux - gf_flux) * (flux - gf_flux) * r_weight / abs(W)
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
# VTK Output
# ============================================================
def output_vtk(mesh, gfu, iteration):
    """Output mesh and solution to VTK file."""
    vtk_filename = os.path.join(script_dir, f"axi_omega_iter_{iteration:02d}")
    vtk = VTKOutput(mesh, coefs=[gfu], names=["Omega_pert"],
                    filename=vtk_filename, subdivision=0)
    vtk.Do()
    return vtk_filename + ".vtu"


# ============================================================
# Main: Adaptive Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating geometry...")
print("=" * 60)

geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.5))
with TaskManager():
	mesh.Curve(order)  # Curve the mesh for higher-order elements

print(f"\nInitial mesh:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# History tracking
history = {
    'ndof': [],
    'elements': [],
    'error': [],
    'Hz_error_percent': [],
    'energy': []
}

# Analytical solution for interior (total field)
Hz_total_analytical = 3.0 / (mu_r + 2) * H0

# Analytical perturbation field energy (3D energy, NOT divided by 2*pi)
#
# 3D analytical energy (from 3D_sphere_with_Kelvin.py):
# - Inside sphere: W_in = (1/2) * mu_r * mu0 * |H_pert|^2 * (4*pi/3) * a^3
# - Outside sphere (dipole): W_out = mu0 * m^2 / (12*pi*a^3) where m = 4*pi*a^3*(mu_r-1)/(mu_r+2)*H0
#
# For axisymmetric calculation with 2*pi multiplier, we get 3D energy directly
H_pert_analytical = -(mu_r - 1) / (mu_r + 2) * H0
V_sphere = (4.0/3.0) * pi * sphere_radius**3
W_in_analytical = 0.5 * mu_r * mu0 * H_pert_analytical**2 * V_sphere
m_dipole = 4 * pi * sphere_radius**3 * (mu_r - 1) / (mu_r + 2) * H0
W_out_analytical = mu0 * m_dipole**2 / (12 * pi * sphere_radius**3)
W_analytical = W_in_analytical + W_out_analytical

# Theta for circles
theta_circle = linspace(0, pi, 100)


def draw_sphere_and_kelvin(ax, z_offset=0, show_legend=True):
    """Draw sphere and Kelvin boundary."""
    r_sphere = sphere_radius * sin(theta_circle)
    z_sphere = sphere_radius * cos(theta_circle)
    if z_offset == 0:
        ax.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
        label1 = 'Magnetic sphere' if show_legend else None
        ax.plot(r_sphere, z_sphere, 'r-', linewidth=2, label=label1)

    r_kelvin = kelvin_radius * sin(theta_circle)
    z_kelvin = kelvin_radius * cos(theta_circle) + z_offset
    label2 = 'Kelvin boundary' if show_legend else None
    ax.plot(r_kelvin, z_kelvin, 'g--', linewidth=1.5, label=label2)


def generate_iteration_plot(iter_num, mesh, gfu, element_errors, history, n_elements, n_vertices):
    """Generate 2x2 plot for current iteration and save MAT file with plotting data."""
    fig = plt.figure(figsize=(14, 14), dpi=150)

    # Separate inner and outer domain elements with error information
    inner_error_map = []  # List of (polygon, log_err)
    outer_error_map = []
    inner_elements = []
    outer_elements = []
    inner_el_idx = []
    outer_el_idx = []

    for el_idx, el in enumerate(mesh.Elements()):
        el_verts = [mesh[v].point for v in el.vertices]
        cx = sum(v[0] for v in el_verts) / len(el_verts)
        cy = sum(v[1] for v in el_verts) / len(el_verts)

        # Get element error
        err_val = element_errors[el_idx]
        if err_val > 1e-20:
            log_err = log10(err_val)
        else:
            log_err = -20

        poly = [(v[0], v[1]) for v in el_verts]

        if cy < offset_z / 2:
            inner_elements.append(el_verts)
            inner_el_idx.append(el_idx)
            inner_error_map.append((poly, log_err))
        else:
            outer_elements.append(el_verts)
            outer_el_idx.append(el_idx)
            outer_error_map.append((poly, log_err))

    # ===== Compute B field for streamplot (flux lines) =====
    # In axisymmetric (r, z) coordinates:
    # B = mu * H = mu * (-grad(Omega) + H_s)
    # Use NGSolve's grad(gfu) directly for accuracy

    # Create grid for inner domain
    n_grid = 50  # Coarser grid for streamplot (faster and cleaner)
    r_inner = linspace(0.02, kelvin_radius - 0.02, n_grid)
    z_inner = linspace(-kelvin_radius + 0.02, kelvin_radius - 0.02, n_grid)
    R_grid, Z_grid = meshgrid(r_inner, z_inner)

    # Compute H field using NGSolve's gradient
    H_grad = -grad(gfu)  # H_pert = -grad(Omega)

    # B field arrays
    B_r_inner = zeros((n_grid, n_grid))
    B_z_inner = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_grid[i, j], Z_grid[i, j]
            # Check if inside domain (half-circle)
            if r_pt**2 + z_pt**2 < (kelvin_radius - 0.01)**2 and r_pt > 0.01:
                try:
                    mip = mesh(r_pt, z_pt)
                    # Get H_pert from NGSolve
                    H_pert_r = H_grad[0](mip)
                    H_pert_z = H_grad[1](mip)
                    # Total H = H_pert + H_s, where H_s = (0, H0)
                    H_r = H_pert_r
                    H_z = H_pert_z + H0
                    # Determine mu at this point
                    rho = sqrt(r_pt**2 + z_pt**2)
                    if rho < sphere_radius:
                        mu_val = mu_r * mu0
                    else:
                        mu_val = mu0
                    # B = mu * H
                    B_r_inner[i, j] = mu_val * H_r
                    B_z_inner[i, j] = mu_val * H_z
                except:
                    B_r_inner[i, j] = nan
                    B_z_inner[i, j] = nan
            else:
                B_r_inner[i, j] = nan
                B_z_inner[i, j] = nan

    # Create grid for outer domain (Kelvin transformed, z-offset)
    r_outer = linspace(0.02, kelvin_radius - 0.02, n_grid)
    z_outer = linspace(offset_z - kelvin_radius + 0.02, offset_z + kelvin_radius - 0.02, n_grid)
    R_out_grid, Z_out_grid = meshgrid(r_outer, z_outer)

    B_r_outer = zeros((n_grid, n_grid))
    B_z_outer = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_out_grid[i, j], Z_out_grid[i, j]
            z_local = z_pt - offset_z  # Local z' coordinate
            # Check if inside domain (half-circle centered at z = offset_z)
            if r_pt**2 + z_local**2 < (kelvin_radius - 0.01)**2 and r_pt > 0.01:
                try:
                    mip = mesh(r_pt, z_pt)
                    # Get H_pert from NGSolve
                    H_pert_r = H_grad[0](mip)
                    H_pert_z = H_grad[1](mip)
                    # Kelvin transformed source field: H_s' = -(rho'/R)^2 * H0
                    rho_prime = sqrt(r_pt**2 + z_local**2)
                    H_s_z = -(rho_prime / kelvin_radius)**2 * H0
                    # Total H
                    H_r = H_pert_r
                    H_z = H_pert_z + H_s_z
                    # Kelvin transformed mu: mu' = (R/rho')^2 * mu0
                    mu_val = (kelvin_radius / (rho_prime + 1e-20))**2 * mu0
                    # B = mu * H
                    B_r_outer[i, j] = mu_val * H_r
                    B_z_outer[i, j] = mu_val * H_z
                except:
                    B_r_outer[i, j] = nan
                    B_z_outer[i, j] = nan
            else:
                B_r_outer[i, j] = nan
                B_z_outer[i, j] = nan

    # ===== Save MAT file with all plotting data =====
    # Convert element data to arrays for MAT file
    inner_el_verts_r = []
    inner_el_verts_z = []
    for el_verts in inner_elements:
        inner_el_verts_r.append([v[0] for v in el_verts])
        inner_el_verts_z.append([v[1] for v in el_verts])

    outer_el_verts_r = []
    outer_el_verts_z = []
    for el_verts in outer_elements:
        outer_el_verts_r.append([v[0] for v in el_verts])
        outer_el_verts_z.append([v[1] for v in el_verts])

    inner_errors = array([element_errors[i] for i in inner_el_idx]) if len(inner_el_idx) > 0 else array([])

    mat_iter_data = {
        'iter_num': iter_num,
        'n_elements': n_elements,
        'n_vertices': n_vertices,
        'ndof': history['ndof'][-1],
        'error': history['error'][-1],
        'Hz_error_percent': history['Hz_error_percent'][-1],
        # History up to this iteration
        'history_ndof': array(history['ndof']),
        'history_error': array(history['error']),
        'history_Hz_error_percent': array(history['Hz_error_percent']),
        # Inner domain mesh
        'inner_n_elements': len(inner_elements),
        'inner_el_verts_r': array(inner_el_verts_r, dtype=object),
        'inner_el_verts_z': array(inner_el_verts_z, dtype=object),
        'inner_element_errors': inner_errors,
        # Outer domain mesh
        'outer_n_elements': len(outer_elements),
        'outer_el_verts_r': array(outer_el_verts_r, dtype=object),
        'outer_el_verts_z': array(outer_el_verts_z, dtype=object),
        # B field grids for streamplot
        'R_grid_inner': R_grid,
        'Z_grid_inner': Z_grid,
        'B_r_inner': B_r_inner,
        'B_z_inner': B_z_inner,
        'R_grid_outer': R_out_grid,
        'Z_grid_outer': Z_out_grid,
        'B_r_outer': B_r_outer,
        'B_z_outer': B_z_outer,
        # Parameters
        'sphere_radius': sphere_radius,
        'kelvin_radius': kelvin_radius,
        'offset_z': offset_z,
        'order': order
    }
    mat_iter_file = os.path.join(script_dir, f"axi_omega_iter_{iter_num:02d}.mat")
    sio.savemat(mat_iter_file, mat_iter_data)
    print(f"  MAT saved: {mat_iter_file}")

    # Colormap settings for error visualization
    norm = Normalize(vmin=-12, vmax=-3)
    cmap_obj = plt.colormaps['jet']

    # ===== Top-left: Inner domain ZZ error map (no flux lines, white mesh edges) =====
    ax1 = plt.subplot(2, 2, 1)
    if len(inner_error_map) > 0:
        polygons = [item[0] for item in inner_error_map]
        err_values = array([item[1] for item in inner_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='white', linewidths=0.3, alpha=0.9)
        ax1.add_collection(pc)
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar1 = plt.colorbar(sm, ax=ax1)
        cbar1.set_label('$\\log_{10}$(ZZ Error)')

    draw_sphere_and_kelvin(ax1, 0, show_legend=True)
    ax1.set_xlabel('$r$ (m)')
    ax1.set_ylabel('$z$ (m)')
    ax1.set_title(f'Inner Domain: ZZ Error ({len(inner_elements)} elements)')
    ax1.set_aspect('equal')
    ax1.set_xlim(-0.05, kelvin_radius + 0.05)
    ax1.set_ylim(-kelvin_radius - 0.05, kelvin_radius + 0.05)

    # ===== Top-right: Outer domain ZZ error map (no flux lines, white mesh edges) =====
    ax2 = plt.subplot(2, 2, 2)
    if len(outer_error_map) > 0:
        polygons = [item[0] for item in outer_error_map]
        err_values = array([item[1] for item in outer_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='white', linewidths=0.3, alpha=0.9)
        ax2.add_collection(pc)
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar2 = plt.colorbar(sm, ax=ax2)
        cbar2.set_label('$\\log_{10}$(ZZ Error)')

    draw_sphere_and_kelvin(ax2, offset_z, show_legend=False)
    ax2.set_xlabel("$r'$ (m)")
    ax2.set_ylabel("$z'$ (m)")
    ax2.set_title(f'Outer Domain: ZZ Error ({len(outer_elements)} elements)')
    ax2.set_aspect('equal')
    ax2.set_xlim(-0.05, kelvin_radius + 0.05)
    ax2.set_ylim(offset_z - kelvin_radius - 0.05, offset_z + kelvin_radius + 0.05)

    # ===== Bottom-left: DOF vs Error convergence curve =====
    ax3 = plt.subplot(2, 2, 3)
    ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Uniform (All Elements)')
    ax3.loglog(history['ndof'][-1], history['error'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    ndof_line = array([1e2, 1e6])
    err_ref_point = history['error'][0]
    N_ref = history['ndof'][0]
    err_line = err_ref_point * (N_ref / ndof_line) ** (order / 2)
    ax3.loglog(ndof_line, err_line, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/2}})$')

    ax3.set_xlabel('DOFs')
    ax3.set_ylabel('Error Estimator')
    ax3.set_title(f'DOFs vs Error (Iteration {iter_num})')
    ax3.set_xlim(1e2, 1e6)
    ax3.set_ylim(1e-5, 1e-1)
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # ===== Bottom-right: DOF vs Magnetic Energy =====
    ax4 = plt.subplot(2, 2, 4)

    ax4.semilogx(history['ndof'], history['energy'], 'gs-', linewidth=2, markersize=6, label='Magnetic Energy')
    ax4.axhline(y=W_analytical, color='r', linestyle='--', linewidth=2, label=f'Analytical: {W_analytical:.4e}')
    ax4.semilogx(history['ndof'][-1], history['energy'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    ax4.set_xlabel('DOFs')
    ax4.set_ylabel('Magnetic Energy (J/m)')
    ax4.set_title('DOF vs Magnetic Energy')
    ax4.set_xlim(1e2, 1e6)
    ax4.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    fig.suptitle(f'Iteration {iter_num}: {n_elements} elements, DOFs={history["ndof"][-1]}, Error={history["error"][-1]:.2e}',
                 fontsize=14, fontweight='bold', y=1.02)

    png_file = os.path.join(script_dir, f"axi_omega_iter_{iter_num:02d}.png")
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    print(f"  PNG saved: {png_file}")
    plt.close()

print("\n" + "=" * 60)
print("Starting Adaptive Mesh Refinement")
print("=" * 60)

for iteration in range(max_iterations):
    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}")
    print("=" * 60)

    # Solve
    fes, gfu, mu, r_weight, BField = solve_omega_formulation(mesh, order)

    # Compute error estimator (using full B field)
    element_errors = compute_error_estimator(mesh, fes, gfu, mu, r_weight, BField)
    total_error = sqrt(sum(element_errors))

    # Compute Hz error inside sphere
    # In Total region (magnetic): H = grad(Omega), so Hz = grad(gfu)[1] directly
    try:
        Hz_total_val = grad(gfu)[1](mesh(0.1, 0))
        Hz_error_percent = abs(Hz_total_val - Hz_total_analytical) / abs(Hz_total_analytical) * 100
    except:
        Hz_error_percent = nan

    # Compute magnetic energy of perturbation field
    # Following 3D_sphere_with_Kelvin.py implementation
    #
    # Total region (magnetic): H_pert = grad(Omega) - H_s, B_pert = mu * H_pert
    # Reduced region (air_inner): H_pert = grad(Omega_r) = grad(Orr) - grad(Oxr)
    # Kelvin region (air_outer): H_pert = grad(Orr)
    #
    # Energy: W = (1/2) * integral(B_pert * H_pert) dV

    # Source field H_s = grad(Omega_s) = (0, H0)
    H_s_cf = CoefficientFunction((0.0, H0))

    # Need separate GridFunctions for accurate perturbation field in Reduced region
    # Following 3D_sphere_with_Kelvin.py lines 360-373
    from ngsolve import VOL
    from ngsolve import TaskManager
    Omega_s = H0 * y  # y is z in axisymmetric

    fesOt_energy = H1(mesh, order=fes.globalorder, definedon="magnetic")
    fesOr_energy = H1(mesh, order=fes.globalorder, definedon="air_inner|air_outer")

    Ot_energy = GridFunction(fesOt_energy)
    Orr_energy = GridFunction(fesOr_energy)
    Oxr_energy = GridFunction(fesOr_energy)

    Ot_energy.Set(gfu, VOL, definedon="magnetic")
    Orr_energy.Set(gfu, VOL, definedon="air_inner|air_outer")
    Oxr_energy.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # --- Total region (magnetic sphere) ---
    # H_pert = grad(Omega) - H_s
    # B_pert = mu_r * mu0 * H_pert
    # Energy = 2*pi * (1/2) * mu_r * mu0 * |H_pert|^2 (3D energy)
    H_pert_total = grad(Ot_energy) - H_s_cf
    W_total = 2 * pi * 0.5 * (mu_r * mu0) * Integrate(InnerProduct(H_pert_total, H_pert_total) * r_weight,
                                               mesh, definedon=mesh.Materials("magnetic"))

    # --- Reduced region (air_inner) ---
    # H_pert = grad(Orr) - grad(Oxr) (perturbation from source potential)
    H_pert_reduced = grad(Orr_energy) - grad(Oxr_energy)
    W_reduced = 2 * pi * 0.5 * mu0 * Integrate(InnerProduct(H_pert_reduced, H_pert_reduced) * r_weight,
                                        mesh, definedon=mesh.Materials("air_inner"))

    # --- Kelvin region (air_outer) ---
    # H_pert = grad(Orr)
    # Using mu_kelvin for Kelvin-transformed region (centralized helper)
    mu_kelvin_energy = mu0 * kelvin_mu_factor_axisym_cf(
        z_offset=offset_z, R=kelvin_radius,
    )
    H_pert_kelvin = grad(Orr_energy)
    W_kelvin = 2 * pi * 0.5 * Integrate(mu_kelvin_energy * InnerProduct(H_pert_kelvin, H_pert_kelvin) * r_weight,
                                mesh, definedon=mesh.Materials("air_outer"))

    # Total perturbation energy
    W_energy = W_total + W_reduced + W_kelvin

    # Record history
    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)
    history['Hz_error_percent'].append(Hz_error_percent)
    history['energy'].append(W_energy)

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")
    print(f"  Hz error inside sphere: {Hz_error_percent:.4f}%")
    print(f"  Magnetic energy: {W_energy:.6e} J/m")

    # Output VTK
    vtk_file = output_vtk(mesh, gfu, iteration)
    print(f"  VTK saved: {vtk_file}")

    # Generate PNG for this iteration (with current mesh and solution)
    generate_iteration_plot(iteration + 1, mesh, gfu, element_errors, history, mesh.ne, mesh.nv)

    # Refine ALL elements (uniform refinement) - skip on last iteration
    if iteration < max_iterations - 1:
        print(f"  Refining ALL {mesh.ne} elements")
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, True)
        mesh.Refine()
        mesh.Curve(order)  # Re-curve mesh after refinement


# ============================================================
# Final Validation
# ============================================================
print("\n" + "=" * 60)
print("FINAL VALIDATION")
print("=" * 60)

print(f"\nAnalytical Hz_total inside sphere: {Hz_total_analytical:.6f} A/m")

print("\n" + "=" * 60)
print("INTERIOR (magnetic sphere):")
print("In Total region: H = grad(Omega)")
print("=" * 60)

for r_val in [0.1, 0.2, 0.3, 0.4]:
    try:
        # In Total region: H = grad(Omega), so Hz = grad(gfu)[1] directly
        Hz_total_val = grad(gfu)[1](mesh(r_val, 0))
        err = abs(Hz_total_val - Hz_total_analytical) / abs(Hz_total_analytical) * 100
        print(f"  r={r_val}: Hz_total = {Hz_total_val:.6f}, error = {err:.4f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")

print("\n" + "=" * 60)
print("EXTERIOR (air region):")
print("In Reduced region: H = H_s + grad(Omega_r) - grad(Omega_s)")
print("=" * 60)

# Post-processing for exterior - need to compute H = H_s + grad(Omega_r) - grad(Omega_s)
fesOr_val = H1(mesh, order=order, definedon="air_inner|air_outer")
Orr_val = GridFunction(fesOr_val)
Oxr_val = GridFunction(fesOr_val)
Omega_s_val = H0 * y  # y is z in axisymmetric

Orr_val.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr_val.Set(Omega_s_val, BND, mesh.Boundaries("sphere"))

for r_val in [0.6, 0.7, 0.8, 0.9]:
    Hz_ext_analytical = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r_val)**3)
    try:
        # In Reduced region: Hz = H0 + (grad(Orr) - grad(Oxr))[1]
        Hz_pert = (grad(Orr_val) - grad(Oxr_val))[1](mesh(r_val, 0))
        Hz_total_val = H0 + Hz_pert
        err = abs(Hz_total_val - Hz_ext_analytical) / abs(Hz_ext_analytical) * 100
        print(f"  r={r_val}: Hz_total = {Hz_total_val:.6f} (analytical: {Hz_ext_analytical:.6f}), error = {err:.4f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")


# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Convergence History")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'Hz Err(%)':<10}")
print("-" * 48)
for i in range(len(history['ndof'])):
    print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
          f"{history['error'][i]:<12.4e} {history['Hz_error_percent'][i]:<10.4f}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
if history['error'][-1] > 0:
    print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
          f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")
print(f"  Hz error: {history['Hz_error_percent'][0]:.4f}% -> {history['Hz_error_percent'][-1]:.4f}%")

# Save history to .mat file
mat_filename = os.path.join(script_dir, os.path.splitext(os.path.basename(__file__))[0] + ".mat")
mat_data = {
    'ndof': array(history['ndof']),
    'elements': array(history['elements']),
    'error': array(history['error']),
    'Hz_error_percent': array(history['Hz_error_percent']),
    'energy': array(history['energy']),
    'order': order
}
sio.savemat(mat_filename, mat_data)
print(f"\nData saved to: {mat_filename}")

# Open the final iteration plot
final_png = os.path.join(script_dir, f"axi_omega_iter_{max_iterations:02d}.png")
if os.path.exists(final_png):
    try:
        os.startfile(final_png)
    except:
        pass

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
